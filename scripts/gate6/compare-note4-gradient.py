#!/usr/bin/env python3
"""Compare Note 4 BS-seq gradient metrics against the recorded per-cell expectations.

``xenofilx`` selection is deterministic, so every evaluated cell must reproduce the metric
values recorded by the original evaluation run. The per-cell expectations live in
``benchmark/note4/per-cell-expected-s6.json`` and cover all 43 Supplementary Table S6 cells;
``benchmark/note4/table-s6-gradient.tsv`` holds the aggregated table for reference.

The comparison validates whatever subset of cells the workflow actually ran, so the
representative push/PR slice and the full dispatch sweep use the same code path.
"""

from __future__ import annotations

import argparse
import json
import pathlib
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
    parser.add_argument("--expected", required=True, help="Path to per-cell-expected-s6.json")
    parser.add_argument("--report", required=True, help="Path for the JSON comparison report")
    parser.add_argument(
        "--tolerance",
        type=float,
        default=0.01,
        help="Allowed deviation in percentage points; the pipeline is deterministic",
    )
    return parser.parse_args()


def load_expected(path: pathlib.Path) -> dict[str, dict[str, str]]:
    return json.loads(path.read_text())


def load_observed(metrics_directory: pathlib.Path) -> dict[str, dict]:
    observed = {}
    for report_path in sorted(metrics_directory.glob("*/source-performance.json")):
        label = report_path.parent.name
        observed[label] = json.loads(report_path.read_text())
    return observed


def observed_value(payload: dict, metric: str) -> float | None:
    tool, key = METRIC_NAMES[metric]
    tools = payload.get("tools", {})
    if tool not in tools:
        return None
    value = tools[tool].get(key)
    return None if value is None else round(value * 100, 4)


def main() -> None:
    arguments = parse_arguments()
    expected = load_expected(pathlib.Path(arguments.expected))
    observed = load_observed(pathlib.Path(arguments.metrics_directory))
    if not observed:
        raise SystemExit(f"no per-cell source-performance.json found under {arguments.metrics_directory}")

    unknown = sorted(label for label in observed if label not in expected)
    if unknown:
        raise SystemExit(f"evaluated cells are not part of the recorded contract: {unknown}")

    comparisons = []
    failures = []
    for label in sorted(observed):
        payload = observed[label]
        row: dict[str, object] = {"cell": label, "metrics": {}}
        for metric in METRIC_NAMES:
            value = observed_value(payload, metric)
            record: dict[str, object] = {"observed": value}
            reference_text = expected[label].get(metric)
            if value is not None and reference_text not in (None, "N/A"):
                reference = float(reference_text)
                deviation = round(value - reference, 4)
                record["expected"] = reference
                record["deviation_percentage_points"] = deviation
                if abs(deviation) > arguments.tolerance:
                    failures.append(f"{label} {metric}: {deviation:+.4f} pp")
            row["metrics"][metric] = record
        comparisons.append(row)

    report = {
        "schema_version": "otter.xenofilx-note4-per-cell-comparison/v1",
        "tolerance_percentage_points": arguments.tolerance,
        "evaluated_cell_count": len(observed),
        "comparisons": comparisons,
        "failures": failures,
    }
    pathlib.Path(arguments.report).write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"evaluated_cells": sorted(observed), "failures": failures}, indent=2))
    if failures:
        print("per-cell comparison failed beyond tolerance", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
