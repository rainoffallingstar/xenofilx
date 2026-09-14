#!/usr/bin/env python3
"""Run one Note 4 BS-seq gradient cell and evaluate the three manuscript pipelines.

The cell contract (fixture labels, sizes, SHA-256) is read from
``benchmark/note4/mixture-cells.json`` and the fixtures are pulled from the published
``fallingstar10/otter-data`` dataset. The three evaluated pipelines are:

* ``modern_xenofilx``                - strand-aware bisulfite NM
* ``legacy_picard_xenofilter``       - Picard ``SetNmMdAndUqTags`` with ``IS_BISULFITE_SEQUENCE=true``
                                       followed by R ``XenofilteR``
* ``conventional_picard_xenofilter`` - the unpatched control (``IS_BISULFITE_SEQUENCE`` omitted)

All tools run through ``enva`` except the statically linked ``xenofilx`` binary.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import subprocess
import sys
import time
import urllib.request

DOWNLOAD_ATTEMPTS = 5
SCRIPT_DIRECTORY = pathlib.Path(__file__).resolve().parent


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", required=True, help="Path to benchmark/note4/mixture-cells.json")
    parser.add_argument("--cell", required=True, help="Gradient cell label, e.g. human0p5-replicate1")
    parser.add_argument("--fixture-base-url", required=True, help="Hugging Face dataset revision URL")
    parser.add_argument("--reference-directory", required=True, help="Directory holding hg19.fa and mm10.fa")
    parser.add_argument("--xenofilx", required=True, help="Path to the compiled xenofilx binary")
    parser.add_argument("--enva", required=True, help="Path to the compiled enva binary")
    parser.add_argument("--environment", required=True, help="enva environment name for the analysis toolchain")
    parser.add_argument("--work-directory", required=True, help="Directory holding per-cell evidence")
    parser.add_argument(
        "--fixture-cache",
        default="ci-artifacts/fixture-cache",
        help="Shared content-addressed cache so repeated runs of one cell download it once",
    )
    parser.add_argument("--gate6-scripts", default=str(SCRIPT_DIRECTORY), help="Directory holding the gate6 helpers")
    parser.add_argument("--threads", type=int, default=4, help="Thread count for xenofilx and XenofilteR")
    parser.add_argument("--sort-memory", default="4G", help="samtools sort memory for xenofilx")
    parser.add_argument("--mm-threshold", type=int, default=None, help="Override the contract mismatch threshold")
    parser.add_argument("--unmapped-penalty", type=int, default=None, help="Override the contract unmapped penalty")
    parser.add_argument("--skip-conventional-control", action="store_true", help="Skip the unpatched Picard control arm")
    return parser.parse_args()


def load_cell(contract_path: pathlib.Path, label: str) -> tuple[dict, dict]:
    contract = json.loads(contract_path.read_text())
    for cell in contract["cells"]:
        if cell["label"] == label:
            return contract, cell
    raise SystemExit(f"cell {label!r} is not present in {contract_path}")


def sha256_of(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_verified(url: str, destination: pathlib.Path, expected_sha256: str, expected_size: int) -> None:
    if destination.exists() and destination.stat().st_size == expected_size:
        if sha256_of(destination) == expected_sha256:
            print(f"reuse {destination} ({expected_size} bytes)")
            return
    for attempt in range(1, DOWNLOAD_ATTEMPTS + 1):
        try:
            with urllib.request.urlopen(url, timeout=120) as response, destination.open("wb") as handle:
                while True:
                    chunk = response.read(1024 * 1024)
                    if not chunk:
                        break
                    handle.write(chunk)
            if destination.stat().st_size != expected_size:
                raise ValueError(f"size mismatch: {destination.stat().st_size} != {expected_size}")
            computed = sha256_of(destination)
            if computed != expected_sha256:
                raise ValueError(f"sha256 mismatch: {computed} != {expected_sha256}")
            print(f"downloaded {destination} ({expected_size} bytes)")
            return
        except Exception as error:  # noqa: BLE001 - the retry loop reports every failure mode
            print(f"download attempt {attempt} failed for {destination}: {error}", file=sys.stderr)
            time.sleep(3 * attempt)
    raise SystemExit(f"could not download {url}")


def run(command: list[str], log_path: pathlib.Path) -> None:
    print("+ " + " ".join(command), flush=True)
    with log_path.open("w") as log_handle:
        completed = subprocess.run(command, stdout=log_handle, stderr=subprocess.STDOUT, check=False)
    if completed.returncode != 0:
        print(log_path.read_text()[-4000:], file=sys.stderr)
        raise SystemExit(f"command failed with exit code {completed.returncode}: {' '.join(command)}")


def enva_command(arguments: argparse.Namespace, command: str) -> list[str]:
    return [arguments.enva, "run", "--name", arguments.environment, "--command", command]


def main() -> None:
    arguments = parse_arguments()
    contract_path = pathlib.Path(arguments.contract)
    contract, cell = load_cell(contract_path, arguments.cell)
    reference_directory = pathlib.Path(arguments.reference_directory)
    parameters = dict(contract["parameters"])
    if arguments.mm_threshold is not None:
        parameters["mm_threshold"] = arguments.mm_threshold
    if arguments.unmapped_penalty is not None:
        parameters["unmapped_penalty"] = arguments.unmapped_penalty
    cell_directory_name = arguments.cell
    if arguments.mm_threshold is not None or arguments.unmapped_penalty is not None:
        cell_directory_name = (
            f"{arguments.cell}-t{parameters['mm_threshold']}-p{parameters['unmapped_penalty']}"
        )
    work_directory = pathlib.Path(arguments.work_directory) / cell_directory_name
    log_directory = work_directory / "logs"
    log_directory.mkdir(parents=True, exist_ok=True)

    # Fixtures are content addressed in a shared cache so a sweep that re-runs the same cell
    # (for example the 12 parameter combinations of one replicate) downloads it only once.
    fixture_cache = pathlib.Path(arguments.fixture_cache)
    fixture_cache.mkdir(parents=True, exist_ok=True)

    dataset_prefix = f"{arguments.fixture_base_url}/xenofilx/mixture/bs-pdx/{arguments.cell}"
    fixtures = {
        "graft": (f"{dataset_prefix}/graft.hg19.bam", cell["graft_hg19_bam"]),
        "host": (f"{dataset_prefix}/host.mm10.bam", cell["host_mm10_bam"]),
        "truth": (f"{dataset_prefix}/truth.tsv", cell["truth_tsv"]),
    }
    local_paths = {}
    for role, (url, metadata) in fixtures.items():
        suffix = ".tsv" if role == "truth" else ".bam"
        local_path = fixture_cache / f"{arguments.cell}.{role}{suffix}"
        download_verified(url, local_path, metadata["sha256"], metadata["size_bytes"])
        local_paths[role] = local_path

    graft_reference = reference_directory / f"{contract['references']['graft']['label']}.fa"
    host_reference = reference_directory / f"{contract['references']['host']['label']}.fa"
    legacy_script = pathlib.Path(arguments.gate6_scripts) / "gate6-xenofilter-legacy.R"
    evaluator = pathlib.Path(arguments.gate6_scripts) / "gate6-mixture-evaluate.py"

    modern_output_directory = work_directory / "modern-xenofilx"
    modern_output_directory.mkdir(exist_ok=True)
    modern_bam = modern_output_directory / "mixture_graft_Filtered.bam"
    run(
        [
            arguments.xenofilx, "run",
            "--graft", str(local_paths["graft"]),
            "--host", str(local_paths["host"]),
            "--output", str(modern_output_directory),
            "--output-names", modern_bam.name,
            "--graft-ref", str(graft_reference),
            "--host-ref", str(host_reference),
            "--mm-threshold", str(parameters["mm_threshold"]),
            "--unmapped-penalty", str(parameters["unmapped_penalty"]),
            "--threads", str(arguments.threads),
            "--sort-memory", arguments.sort_memory,
            "--recalculate-nm",
            "--bisulfite",
        ],
        log_directory / "modern-xenofilx.log",
    )

    arms = {
        "legacy_picard_xenofilter": "true",
    }
    if not arguments.skip_conventional_control:
        arms["conventional_picard_xenofilter"] = "false"

    tool_outputs = [f"modern_xenofilx={modern_bam}"]
    for arm_name, bisulfite_flag in arms.items():
        arm_directory = work_directory / arm_name
        arm_directory.mkdir(exist_ok=True)
        picard_graft = arm_directory / "graft.picard-nm.bam"
        picard_host = arm_directory / "host.picard-nm.bam"
        extra = " IS_BISULFITE_SEQUENCE=true" if bisulfite_flag == "true" else ""
        run(
            enva_command(
                arguments,
                f"picard SetNmMdAndUqTags I={local_paths['graft']} O={picard_graft} "
                f"R={graft_reference}{extra}",
            ),
            log_directory / f"picard-{arm_name}-graft.log",
        )
        run(
            enva_command(
                arguments,
                f"picard SetNmMdAndUqTags I={local_paths['host']} O={picard_host} "
                f"R={host_reference}{extra}",
            ),
            log_directory / f"picard-{arm_name}-host.log",
        )
        destination = arm_directory / "destination"
        destination.mkdir(exist_ok=True)
        run(
            enva_command(
                arguments,
                f"Rscript {legacy_script} {picard_graft} {picard_host} {destination} "
                f"{arguments.threads} {parameters['mm_threshold']} {parameters['unmapped_penalty']} "
                f"NM {arm_directory / 'xenofilter-run.json'}",
            ),
            log_directory / f"xenofilter-{arm_name}.log",
        )
        filtered = sorted((destination / "Filtered_bams").glob("*_Filtered.bam"))
        if len(filtered) != 1:
            raise SystemExit(f"{arm_name}: expected one filtered BAM, found {len(filtered)}")
        tool_outputs.append(f"{arm_name}={filtered[0]}")

    report_path = work_directory / "source-performance.json"
    tool_output_arguments = " ".join(f"--tool-output {entry}" for entry in tool_outputs)
    run(
        enva_command(
            arguments,
            f"python3 {evaluator} --truth {local_paths['truth']} --samtools samtools "
            f"{tool_output_arguments} --require-known-selection "
            f"--graft-source human --report {report_path}",
        ),
        log_directory / "evaluation.log",
    )
    print(f"cell {arguments.cell} complete: {report_path}")


if __name__ == "__main__":
    main()
