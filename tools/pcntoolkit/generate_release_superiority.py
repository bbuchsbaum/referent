# /// script
# requires-python = "==3.12.11"
# dependencies = [
#   "numpy==2.4.6",
#   "pandas==3.0.5",
#   "pcntoolkit==1.3.0",
#   "scipy==1.18.1",
# ]
# ///
"""Generate independent PCNtoolkit predictions for the release benchmark."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from generate_fixtures import COMPARATOR_VERSION, fit_scenario, scenario_frame, stable_csv


def md5(path: Path) -> str:
    return hashlib.md5(path.read_bytes(), usedforsecurity=False).hexdigest()


def write_release_superiority(
    output: Path, replicates: int, seed_start: int
) -> None:
    if replicates < 5:
        raise ValueError("release superiority requires at least five replicates")
    output.mkdir(parents=True, exist_ok=True)
    input_parts: list[pd.DataFrame] = []
    prediction_parts: list[pd.DataFrame] = []
    warning_receipts: list[dict[str, object]] = []
    seeds = list(range(seed_start, seed_start + replicates))

    for replicate, seed in enumerate(seeds, start=1):
        frame = scenario_frame(
            "skew_heavy",
            n_train=500,
            n_validation=100,
            n_test=500,
            seed=seed,
        )
        frame["row_id"] = frame["row_id"].map(
            lambda value: f"rep-{replicate:03d}-{value}"
        )
        frame.insert(0, "seed", seed)
        frame.insert(0, "replicate", replicate)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            prediction = fit_scenario(frame, "skew_heavy")
        numeric_evidence = prediction[
            ["observed", "median", "predictive_sd", "centile", "z", "log_density",
             "q05", "q25", "q50", "q75", "q95"]
        ].to_numpy()
        if not np.isfinite(numeric_evidence).all():
            raise RuntimeError(f"non-finite PCNtoolkit evidence in replicate {replicate}")
        warning_receipts.append(
            {
                "replicate": replicate,
                "count": len(caught),
                "messages": [str(item.message) for item in caught],
            }
        )
        prediction.insert(0, "seed", seed)
        prediction.insert(0, "replicate", replicate)
        input_parts.append(frame)
        prediction_parts.append(prediction)

    inputs = pd.concat(input_parts, ignore_index=True)
    predictions = pd.concat(prediction_parts, ignore_index=True)
    input_path = output / "release_inputs.csv"
    prediction_path = output / "release_pcntoolkit_predictions.csv"
    stable_csv(inputs, input_path)
    stable_csv(predictions, prediction_path)

    versions = {
        package: importlib.metadata.version(package)
        for package in ("pcntoolkit", "numpy", "pandas", "scipy")
    }
    if versions["pcntoolkit"] != COMPARATOR_VERSION:
        raise RuntimeError(
            f"expected PCNtoolkit {COMPARATOR_VERSION}, got {versions['pcntoolkit']}"
        )
    receipt = {
        "schema_version": "1.0.0",
        "comparator": {"package": "pcntoolkit", "version": COMPARATOR_VERSION},
        "scenario": "skew_heavy",
        "replicates": replicates,
        "seeds": seeds,
        "split_sizes": {"train": 500, "validation": 100, "test": 500},
        "selection_rule": "prespecified SHASH after family validation",
        "fit_warnings": warning_receipts,
        "dependencies": versions,
        "files": {
            input_path.name: {
                "md5": md5(input_path),
                "rows": len(inputs),
                "columns": list(inputs.columns),
            },
            prediction_path.name: {
                "md5": md5(prediction_path),
                "rows": len(predictions),
                "columns": list(predictions.columns),
            },
        },
    }
    (output / "release_generation_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--replicates", type=int, default=20)
    parser.add_argument("--seed-start", type=int, default=20260824)
    args = parser.parse_args()
    write_release_superiority(args.output.resolve(), args.replicates, args.seed_start)


if __name__ == "__main__":
    main()
