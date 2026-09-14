#!/usr/bin/env python3
"""Aggregate a Xenofilx parameter grid into per-assay metrics and a recommendation."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


METRIC_KEYS = (
    "accuracy",
    "sensitivity_graft_recall",
    "specificity_host_recall",
    "precision_graft_ppv",
    "f1_graft",
)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid-directory", required=True, type=Path)
    parser.add_argument("--assay", required=True)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument(
        "--selection-metric",
        default="f1_graft",
        choices=METRIC_KEYS,
        help="Metric used to select the recommended grid cell",
    )
    return parser.parse_args()


def load_cell(path: Path) -> dict[str, object]:
    document = json.loads(path.read_text(encoding="utf-8"))
    tool = document["tools"]["modern_xenofilx"]
    return {key: tool[key] for key in METRIC_KEYS}


def main() -> int:
    arguments = parse_arguments()
    grid_paths = sorted(arguments.grid_directory.glob("*.metrics.json"))
    if not grid_paths:
        raise SystemExit(f"no grid metrics found in {arguments.grid_directory}")

    cells = []
    for grid_path in grid_paths:
        stem = grid_path.name.removesuffix(".metrics.json")
        # Expected stem: t<threshold>-p<penalty>
        threshold_part, penalty_part = stem.split("-")
        cell = {
            "label": stem,
            "mm_threshold": int(threshold_part.removeprefix("t")),
            "unmapped_penalty": int(penalty_part.removeprefix("p")),
            "metrics": load_cell(grid_path),
        }
        cells.append(cell)

    best = max(cells, key=lambda cell: float(cell["metrics"][arguments.selection_metric]))
    report = {
        "schema_version": "otter.xenofilx-parameter-grid/v1",
        "assay": arguments.assay,
        "selection_metric": arguments.selection_metric,
        "grid_size": len(cells),
        "cells": cells,
        "recommendation": {
            "mm_threshold": best["mm_threshold"],
            "unmapped_penalty": best["unmapped_penalty"],
            "selection_metric": arguments.selection_metric,
            "selection_value": best["metrics"][arguments.selection_metric],
            "scope": (
                f"Valid only for the {arguments.assay} fixture and assay; it is not a "
                "cross-assay recommendation."
            ),
        },
    }

    arguments.report.parent.mkdir(parents=True, exist_ok=True)
    arguments.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report["recommendation"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
