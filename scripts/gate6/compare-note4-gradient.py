#!/usr/bin/env python3
"""Compare Note 4 BS-seq gradient metrics against the published Supplementary Table S6.

Reads every ``source-performance.json`` produced by
``run-note4-gradient-cell.py``, aggregates the modern ``xenofilx`` and prototype
Picard + ``XenofilteR`` arms per human fraction, and compares the means with
``benchmark/note4/table-s6-gradient.tsv``.

The comparison fails only for fractions that were actually evaluated, so the
representative push/PR slice validates the cells it covers while the full dispatch
sweep validates all 43.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import re
import statistics
import sys

MODERN_TOOL = "modern_xenofilx"
LEGACY_TOOL = "legacy_picard_xenofilter"

METRIC_NAMES = {
    "modern_recall": (MODERN_TOOL, "sensitivity_graft_recall"),
    "legacy_recall": (LEGACY_TOOL, "sensitivity_graft_recall"),
    "modern_spec": (MODERN_TOOL, "specificity_host_recall"),
    "legacy_spec": (LEGACY_TOOL, "specificity_host_recall"),
    "modern_ppv": (MODERN_TOOL, "precision_graft_ppv"),
    "legacy_ppv": (LEGACY_TOOL, "precision_graft_ppv"),
}


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metrics-directory", required=True, help="Directory holding per-cell evidence")
    parser.add_argument("--expected", required=True, help="Path to benchmark/note4/table-s6-gradient.tsv")
    parser.add_argument("--report", required=True, help="Path for the JSON comparison report")
    parser.add_argument("--tolerance", type=float, default=0.05, help="Allowed deviation in percentage points")
    return parser.parse_args()


def fraction_from_label(label: str) -> float:
    if label.startswith("human0-replicate"):
        return 0.0
    match = re.match(r"human0p(\d+)-replicate\d+", label)
    if not match:
        raise ValueError(f"unrecognised cell label: {label}")
    return float("0." + match.group(1))


def load_expected(path: pathlib.Path) -> dict[float, dict[str, float | None]]:
    expected: dict[float, dict[str, float | None]] = {}
    lines = path.read_text().splitlines()
    header = lines[0].split("\t")
    for line in lines[1:]:
        if not line.strip():
            continue
        values = line.split("\t")
        row = dict(zip(header, values))
        fraction = float(row["fraction"])
        expected[fraction] = {
            name: (None if row[name] == "N/A" else float(row[name])) for name in METRIC_NAMES
        }
    return expected


def collect_metrics(metrics_directory: pathlib.Path) -> dict[float, list[dict]]:
    grouped: dict[float, list[dict]] = {}
    for report_path in sorted(metrics_directory.glob("*/source-performance.json")):
        payload = json.loads(report_path.read_text())
        grouped.setdefault(fraction_from_label(report_path.parent.name), []).append(payload)
    return grouped


def mean_metric(payloads: list[dict], tool: str, key: str) -> float | None:
    values = [
        payload["tools"][tool][key]
        for payload in payloads
        if tool in payload.get("tools", {}) and payload["tools"][tool].get(key) is not None
    ]
    return statistics.fmean(values) if values else None


def main() -> None:
    arguments = parse_arguments()
    expected = load_expected(pathlib.Path(arguments.expected))
    observed = collect_metrics(pathlib.Path(arguments.metrics_directory))
    if not observed:
        raise SystemExit(f"no per-cell source-performance.json found under {arguments.metrics_directory}")

    comparisons = []
    failures = []
    for fraction in sorted(observed):
        payloads = observed[fraction]
        row: dict[str, object] = {"fraction": fraction, "evaluated_cells": len(payloads), "metrics": {}}
        expected_row = expected.get(fraction)
        for name, (tool, key) in METRIC_NAMES.items():
            value = mean_metric(payloads, tool, key)
            record: dict[str, object] = {"observed": None if value is None else round(value * 100, 4)}
            if value is not None and expected_row and expected_row[name] is not None:
                reference = expected_row[name]
                deviation = round((value * 100) - reference, 4)
                record["expected"] = reference
                record["deviation_percentage_points"] = deviation
                if abs(deviation) > arguments.tolerance:
                    failures.append(f"fraction {fraction:g} {name}: {deviation:+.4f} pp")
            row["metrics"][name] = record
        comparisons.append(row)

    report = {
        "schema_version": "otter.xenofilx-note4-gradient-comparison/v1",
        "tolerance_percentage_points": arguments.tolerance,
        "evaluated_fraction_count": len(observed),
        "comparisons": comparisons,
        "failures": failures,
    }
    pathlib.Path(arguments.report).write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"evaluated_fractions": len(observed), "failures": failures}, indent=2))
    if failures:
        print("gradient comparison failed beyond tolerance", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
