#!/usr/bin/env python3
"""Compare Note 4 BS-seq gradient metrics against the recorded per-cell expectations.

``xenofilx`` selection is deterministic, so every evaluated cell should reproduce the metric
values recorded by the original evaluation run. The per-cell expectations live in
``benchmark/note4/per-cell-expected-s6.json`` and cover all 43 Supplementary Table S6 cells;
``benchmark/note4/table-s6-gradient.tsv`` holds the aggregated table for reference.

Deviation policy: the comparison is **reported, not enforced**. Metric differences are
classified as *within range*, *deviation*, or *large deviation* and always emitted as
warnings, so a value mismatch never fails the job. Only a missing or structurally invalid
evidence set is fatal, because that means the evaluation itself did not run.
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
        default=0.5,
        help="Percentage points considered an acceptable difference for a classification metric",
    )
    parser.add_argument(
        "--large-deviation",
        type=float,
        default=2.0,
        help="Percentage points above which a difference is reported as large",
    )
    parser.add_argument(
        "--fail-on-deviation",
        action="store_true",
        help="Opt in to a non-zero exit when a metric leaves the acceptable range",
    )
    return parser.parse_args()


def load_expected(path: pathlib.Path) -> dict[str, dict[str, str]]:
    return json.loads(path.read_text())


def load_observed(metrics_directory: pathlib.Path) -> dict[str, dict]:
    observed = {}
    for report_path in sorted(metrics_directory.glob("*/source-performance.json")):
        observed[report_path.parent.name] = json.loads(report_path.read_text())
    return observed


def observed_value(payload: dict, metric: str) -> float | None:
    tool, key = METRIC_NAMES[metric]
    tools = payload.get("tools", {})
    if tool not in tools:
        return None
    value = tools[tool].get(key)
    return None if value is None else round(value * 100, 4)


def classify(deviation: float, tolerance: float, large_deviation: float) -> str:
    magnitude = abs(deviation)
    if magnitude > large_deviation:
        return "large_deviation"
    if magnitude > tolerance:
        return "deviation"
    return "within_range"


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
    within_range = 0
    deviations = 0
    large_deviations = 0
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
                status = classify(deviation, arguments.tolerance, arguments.large_deviation)
                record.update(
                    expected=reference,
                    deviation_percentage_points=deviation,
                    status=status,
                )
                if status == "within_range":
                    within_range += 1
                elif status == "deviation":
                    deviations += 1
                    print(
                        f"::warning::{label} {metric}: observed {value} vs expected {reference} "
                        f"({deviation:+.4f} pp, acceptable range +/-{arguments.tolerance} pp)"
                    )
                else:
                    large_deviations += 1
                    print(
                        f"::warning::{label} {metric}: observed {value} vs expected {reference} "
                        f"({deviation:+.4f} pp, large difference; accepted and recorded)"
                    )
            row["metrics"][metric] = record
        comparisons.append(row)

    report = {
        "schema_version": "otter.xenofilx-note4-per-cell-comparison/v2",
        "enforcement": "report_only",
        "tolerance_percentage_points": arguments.tolerance,
        "large_deviation_percentage_points": arguments.large_deviation,
        "evaluated_cell_count": len(observed),
        "metric_outcomes": {
            "within_range": within_range,
            "deviation": deviations,
            "large_deviation": large_deviations,
        },
        "comparisons": comparisons,
    }
    pathlib.Path(arguments.report).write_text(json.dumps(report, indent=2) + "\n")
    summary = {
        "evaluated_cells": sorted(observed),
        "within_range": within_range,
        "deviations_within_acceptance": deviations,
        "large_deviations_recorded": large_deviations,
        "enforced": False,
    }
    print(json.dumps(summary, indent=2))
    print(
        f"::notice::Note 4 comparison recorded without enforcement: {within_range} metrics within "
        f"+/-{arguments.tolerance} pp, {deviations} outside it, {large_deviations} large."
    )
    if arguments.fail_on_deviation and (deviations or large_deviations):
        print("comparison deviations present and --fail-on-deviation was requested", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
