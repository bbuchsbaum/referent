# /// script
# requires-python = "==3.12.11"
# dependencies = [
#   "numpy==2.4.6",
#   "pandas==3.0.5",
#   "pcntoolkit==1.3.0",
#   "scipy==1.18.1",
# ]
# ///
"""Regenerate deterministic PCNtoolkit 1.3.0 comparison fixtures.

Run from the repository root with:

    uv run tools/pcntoolkit/generate_fixtures.py

The ordinary R test suite reads the checked-in results and never invokes this
script. The manifest intentionally has no wall-clock timestamp so two clean
runs in the pinned environment are byte-for-byte identical.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from pcntoolkit import NormativeModel
from pcntoolkit.dataio.norm_data import NormData
from pcntoolkit.math_functions.basis_function import BsplineBasisFunction
from pcntoolkit.math_functions.shash import S, S_inv, m1m2
from pcntoolkit.regression_model.blr import BLR


COMPARATOR_VERSION = "1.3.0"
SCHEMA_VERSION = "1.1.0"
SKEW_HEAVY_LBFGSB_EPSILON = 0.01
SEEDS = {
    "null_linear": 1101,
    "linear_gaussian": 1102,
    "multiple_encoded": 1103,
    "nonlinear_heteroskedastic": 1104,
    "log_linear": 1105,
    "balanced_site": 1106,
    "skew_heavy": 1107,
    "unequal_site": 1108,
    "covariate_shift": 1109,
}
ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = ROOT / "tests" / "fixtures" / "pcntoolkit" / "v1.3.0"
CENTILES = np.array([0.05, 0.25, 0.5, 0.75, 0.95], dtype=float)


def stable_csv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False, float_format="%.17g", lineterminator="\n")


def md5(path: Path) -> str:
    return hashlib.md5(path.read_bytes(), usedforsecurity=False).hexdigest()


def gaussian_semantics() -> pd.DataFrame:
    rows: list[dict[str, float | str]] = []
    cases = [
        ("standard", 0.0, 1.0, -0.35, 0.73),
        ("shifted", 3.25, 0.7, 2.1, 0.01),
        ("wide", -4.0, 5.5, 8.0, 0.999),
        ("lower_tail", 1.0, 2.0, -15.0, 1e-12),
        ("upper_tail", -2.0, 0.4, 2.5, 1 - 1e-12),
    ]
    for case, mu, sigma, y, p in cases:
        rows.append(
            {
                "case": case,
                "mu": mu,
                "sigma": sigma,
                "y": y,
                "p": p,
                "density": stats.norm.pdf(y, loc=mu, scale=sigma),
                "log_density": stats.norm.logpdf(y, loc=mu, scale=sigma),
                "cdf": stats.norm.cdf(y, loc=mu, scale=sigma),
                "log_cdf": stats.norm.logcdf(y, loc=mu, scale=sigma),
                "upper_tail": stats.norm.sf(y, loc=mu, scale=sigma),
                "log_upper_tail": stats.norm.logsf(y, loc=mu, scale=sigma),
                "quantile": stats.norm.ppf(p, loc=mu, scale=sigma),
                "z": (y - mu) / sigma,
            }
        )
    return pd.DataFrame(rows)


def shashb_semantics() -> pd.DataFrame:
    rows: list[dict[str, float | str]] = []
    cases = [
        ("identity", 0.0, 1.0, 0.0, 1.0, -0.4, 0.1),
        ("skew", 2.0, 1.3, 0.6, 1.0, 3.2, 0.5),
        ("heavy", -1.0, 0.8, -0.35, 0.7, -2.4, 0.95),
        ("light", 4.0, 2.2, 0.25, 1.5, 5.1, 0.99),
        ("tail", -3.0, 0.45, 0.8, 0.65, 0.5, 1 - 1e-8),
    ]
    for case, mu_b, sigma_b, epsilon, delta, y, p in cases:
        mean_s, raw_second_s = m1m2(epsilon, delta)
        var_s = raw_second_s - mean_s**2
        remapped = ((y - mu_b) / sigma_b) * np.sqrt(var_s) + mean_s
        gaussian_score = S(remapped, epsilon, delta)
        log_density = (
            stats.norm.logpdf(gaussian_score)
            + np.log(delta)
            + 0.5 * np.log1p(gaussian_score**2)
            - 0.5 * np.log1p(remapped**2)
            + 0.5 * np.log(var_s)
            - np.log(sigma_b)
        )
        raw_quantile = S_inv(stats.norm.ppf(p), epsilon, delta)
        quantile = mu_b + sigma_b * (raw_quantile - mean_s) / np.sqrt(var_s)
        sigma_r = sigma_b / (delta * np.sqrt(var_s))
        mu_r = mu_b - sigma_b * mean_s / np.sqrt(var_s)
        rows.append(
            {
                "case": case,
                "mu_b": mu_b,
                "sigma_b": sigma_b,
                "epsilon": epsilon,
                "delta": delta,
                "y": y,
                "p": p,
                "moment_1": mean_s,
                "moment_2": raw_second_s,
                "mu_r": mu_r,
                "sigma_r": sigma_r,
                "density": np.exp(log_density),
                "log_density": log_density,
                "cdf": stats.norm.cdf(gaussian_score),
                "log_cdf": stats.norm.logcdf(gaussian_score),
                "upper_tail": stats.norm.sf(gaussian_score),
                "log_upper_tail": stats.norm.logsf(gaussian_score),
                "quantile": quantile,
                "z": gaussian_score,
            }
        )
    return pd.DataFrame(rows)


def metric_semantics() -> pd.DataFrame:
    observed = np.array([-2.0, -0.2, 0.1, 0.7, 1.4, 2.8, 5.0])
    predicted = np.array([-1.5, -0.4, 0.3, 0.5, 1.8, 2.1, 4.2])
    z = np.array([-1.8, -0.9, -0.2, 0.1, 0.5, 1.4, 2.7])
    pit = stats.norm.cdf(z)
    residual = observed - predicted
    grid = CENTILES
    mace = np.mean([abs(q - np.mean(pit <= q)) for q in grid])
    pit_a = np.array([0.97, 0.99])
    pit_b = np.array([0.02, 0.08, 0.22, 0.38, 0.55, 0.71, 0.86, 0.96])
    mace_a = np.mean([abs(q - np.mean(pit_a <= q)) for q in grid])
    mace_b = np.mean([abs(q - np.mean(pit_b <= q)) for q in grid])
    mace_pooled = np.mean(
        [abs(q - np.mean(np.concatenate([pit_a, pit_b]) <= q)) for q in grid]
    )
    z_sample_standardized = (z - np.mean(z)) / np.std(z, ddof=1)
    values = {
        "rmse": np.sqrt(np.mean(residual**2)),
        "smse_population": np.mean(residual**2) / np.var(observed, ddof=0),
        "ev": 1 - np.var(residual, ddof=1) / np.var(observed, ddof=1),
        "rho_spearman": stats.spearmanr(observed, predicted).statistic,
        "pearson": stats.pearsonr(observed, predicted).statistic,
        "mace_no_batch": mace,
        "skew_bias_corrected": stats.skew(z, bias=False),
        "skew_unadjusted": stats.skew(z, bias=True),
        "skew_referent_sample_sd": np.mean(z_sample_standardized**3),
        "excess_kurtosis_bias_corrected": stats.kurtosis(z, fisher=True, bias=False),
        "excess_kurtosis_unadjusted": stats.kurtosis(z, fisher=True, bias=True),
        "excess_kurtosis_referent_sample_sd": np.mean(z_sample_standardized**4) - 3,
        "mace_equal_batch_average": np.mean([mace_a, mace_b]),
        "mace_pooled_rows": mace_pooled,
    }
    return pd.DataFrame({"metric": list(values), "value": list(values.values())})


def scenario_frame(
    name: str,
    n_train: int = 500,
    n_validation: int = 100,
    n_test: int = 501,
    seed: int | None = None,
) -> pd.DataFrame:
    rng = np.random.default_rng(SEEDS[name] if seed is None else seed)
    n = n_train + n_validation + n_test
    split_counts = (("train", n_train), ("validation", n_validation), ("test", n_test))
    split = np.concatenate([np.repeat(part, count) for part, count in split_counts])
    row_id = [
        f"{name}-{part}-{i:04d}"
        for part, count in split_counts
        for i in range(count)
    ]
    if name == "covariate_shift":
        x1 = np.concatenate(
            [
                rng.uniform(-2.0, 0.75, n_train),
                rng.uniform(0.25, 2.0, n_validation + n_test),
            ]
        )
    else:
        x1 = rng.uniform(-2.0, 2.0, n)
    x2 = rng.normal(size=n)
    sex_m = rng.integers(0, 2, n)
    if name == "unequal_site":
        site = rng.choice(
            np.array(["site-1", "site-2", "site-3", "site-4"]),
            size=n,
            p=np.array([0.65, 0.2, 0.1, 0.05]),
        )
    else:
        site = np.array([f"site-{i % 4 + 1}" for i in rng.permutation(n)])
    noise = rng.normal(size=n)

    if name == "null_linear":
        y = 2.0 + 0.8 * noise
    elif name == "linear_gaussian":
        y = 0.7 + 1.4 * x1 + 0.75 * noise
    elif name == "multiple_encoded":
        y = -0.3 + 0.9 * x1 - 0.55 * x2 + 0.6 * sex_m + 0.7 * noise
    elif name == "nonlinear_heteroskedastic":
        sigma = np.exp(-0.25 + 0.28 * x1)
        y = 0.4 + np.sin(1.4 * x1) + 0.2 * x1**2 + sigma * noise
    elif name == "log_linear":
        y = np.exp(1.2 + 0.35 * x1 + 0.25 * noise)
    elif name == "balanced_site":
        shifts = {"site-1": -0.8, "site-2": -0.2, "site-3": 0.35, "site-4": 0.9}
        y = 1.0 + 0.65 * x1 + np.array([shifts[s] for s in site]) + 0.65 * noise
    elif name == "skew_heavy":
        mu = 0.4 + np.sin(1.35 * x1) + 0.15 * x1**2
        sigma = np.exp(-0.35 + 0.22 * x1)
        epsilon, delta = 0.9, 0.7
        y = mu + sigma * delta * np.sinh((np.arcsinh(noise) + epsilon) / delta)
    elif name == "unequal_site":
        shifts = {"site-1": -0.7, "site-2": -0.1, "site-3": 0.5, "site-4": 1.1}
        scales = {"site-1": 0.6, "site-2": 0.8, "site-3": 1.0, "site-4": 1.25}
        y = 0.8 + 0.7 * x1 + np.array([shifts[s] for s in site]) + np.array(
            [scales[s] for s in site]
        ) * noise
    elif name == "covariate_shift":
        y = -0.2 + 1.1 * x1 + 0.7 * noise
    else:
        raise ValueError(f"unknown scenario: {name}")

    return pd.DataFrame(
        {
            "scenario": name,
            "row_id": row_id,
            "split": split,
            "x1": x1,
            "x2": x2,
            "sex_M": sex_m,
            "site": site,
            "y": y,
        }
    )


def fit_scenario(frame: pd.DataFrame, name: str) -> pd.DataFrame:
    train = frame.loc[frame["split"] == "train"].copy()
    test = frame.loc[frame["split"] == "test"].copy()
    covariates = {
        "null_linear": ["x1"],
        "linear_gaussian": ["x1"],
        "multiple_encoded": ["x1", "x2", "sex_M"],
        "nonlinear_heteroskedastic": ["x1"],
        "log_linear": ["x1"],
        "balanced_site": ["x1"],
        "skew_heavy": ["x1"],
        "unequal_site": ["x1"],
        "covariate_shift": ["x1"],
    }[name]
    batch_effects = ["site"] if name in ("balanced_site", "unequal_site") else None
    response = "y_model"
    if name == "log_linear":
        train[response] = np.log(train["y"])
        test[response] = np.log(test["y"])
    else:
        train[response] = train["y"]
        test[response] = test["y"]

    kwargs: dict[str, object] = {"n_iter": 300, "tol": 1e-6}
    if name in ("nonlinear_heteroskedastic", "skew_heavy"):
        kwargs.update(
            heteroskedastic=True,
            basis_function_mean=BsplineBasisFunction(basis_column=0, degree=3, nknots=5),
            basis_function_var=BsplineBasisFunction(basis_column=0, degree=3, nknots=4),
        )
    if name in ("balanced_site", "unequal_site"):
        kwargs.update(fixed_effect=True, fixed_effect_var=True)
    if name == "skew_heavy":
        # PCNtoolkit 1.3.0 defaults to a very large 0.1 finite-difference
        # step. On the registered skew-heavy matrix that makes otherwise
        # valid L-BFGS-B fits probe ill-conditioned posterior matrices. A
        # 0.01 step retains the pinned model and optimizer while avoiding
        # those invalid probes across the locked fresh-seed replicate set.
        kwargs["l_bfgs_b_epsilon"] = SKEW_HEAVY_LBFGSB_EPSILON

    fit_data = NormData.from_dataframe(
        f"{name}-train",
        train,
        covariates=covariates,
        batch_effects=batch_effects,
        response_vars=[response],
        subject_ids="row_id",
    )
    test_data = NormData.from_dataframe(
        f"{name}-test",
        test,
        covariates=covariates,
        batch_effects=batch_effects,
        response_vars=[response],
        subject_ids="row_id",
    )
    model = NormativeModel(
        BLR(**kwargs),
        savemodel=False,
        evaluate_model=False,
        saveresults=False,
        saveplots=False,
        inscaler="standardize",
        outscaler="standardize",
        name=f"fixture-{name}",
    )
    model.fit(fit_data)
    predicted = model.predict(test_data)
    yhat = predicted["Yhat"].sel(response_vars=response).values
    z = predicted["Z"].sel(response_vars=response).values
    logp = predicted["logp"].sel(response_vars=response).values
    # NormativeModel maps Yhat/centiles back through outscaler but leaves logp
    # on the standardised outcome scale. A density on the response scale needs
    # the inverse-standardisation Jacobian. numpy's ddof=0 matches the scaler.
    logp = logp - np.log(np.std(train[response].to_numpy(), ddof=0))
    quantiles = {
        q: predicted["centiles"].sel(response_vars=response, centile=q).values
        for q in CENTILES
    }

    if name == "log_linear":
        yhat = np.exp(yhat)
        quantiles = {q: np.exp(values) for q, values in quantiles.items()}
        logp = logp - np.log(test["y"].to_numpy())

    q05, q95 = quantiles[0.05], quantiles[0.95]
    pred_sd = (q95 - q05) / (2 * stats.norm.ppf(0.95))
    out = pd.DataFrame(
        {
            "scenario": name,
            "row_id": test["row_id"].to_numpy(),
            "observed": test["y"].to_numpy(),
            "median": yhat,
            "predictive_sd": pred_sd,
            "centile": stats.norm.cdf(z),
            "z": z,
            "log_density": logp,
            "q05": q05,
            "q25": quantiles[0.25],
            "q50": quantiles[0.5],
            "q75": quantiles[0.75],
            "q95": q95,
        }
    )
    return out


def write_fixture(output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    frames = [scenario_frame(name) for name in SEEDS]
    inputs = pd.concat(frames, ignore_index=True)
    predictions = pd.concat(
        [fit_scenario(frame, name) for name, frame in zip(SEEDS, frames)],
        ignore_index=True,
    )

    products = {
        "semantic_gaussian.csv": gaussian_semantics(),
        "semantic_shashb.csv": shashb_semantics(),
        "semantic_metrics.csv": metric_semantics(),
        "fitted_inputs.csv": inputs,
        "fitted_predictions.csv": predictions,
    }
    for filename, frame in products.items():
        stable_csv(frame, output / filename)

    file_receipts = {
        filename: {
            "md5": md5(output / filename),
            "rows": int(len(frame)),
            "columns": list(frame.columns),
        }
        for filename, frame in products.items()
    }
    versions = {
        package: importlib.metadata.version(package)
        for package in ("pcntoolkit", "numpy", "pandas", "scipy")
    }
    if versions["pcntoolkit"] != COMPARATOR_VERSION:
        raise RuntimeError(f"expected PCNtoolkit {COMPARATOR_VERSION}, got {versions['pcntoolkit']}")
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "comparator": {"package": "pcntoolkit", "version": COMPARATOR_VERSION},
        "generator": "tools/pcntoolkit/generate_fixtures.py",
        "python": platform.python_version(),
        "dependencies": versions,
        "seeds": SEEDS,
        "centiles": CENTILES.tolist(),
        "files": file_receipts,
        "scenario_classes": {
            "null_linear": "matched_estimator",
            "linear_gaussian": "matched_estimator",
            "multiple_encoded": "matched_estimator",
            "nonlinear_heteroskedastic": "same_estimand",
            "log_linear": "matched_estimator",
            "balanced_site": "same_estimand",
            "skew_heavy": "same_estimand",
            "unequal_site": "same_estimand",
            "covariate_shift": "matched_estimator",
        },
        "optimizer_controls": {
            "skew_heavy": {
                "optimizer": "l-bfgs-b",
                "l_bfgs_b_epsilon": SKEW_HEAVY_LBFGSB_EPSILON,
            }
        },
        "normalisations": {
            "log_density": "PCNtoolkit logp minus log(population SD of the fitted outcome); log-transformed scenarios also minus log(y)",
            "predictive_sd": "(q95 - q05) / (2 * qnorm(0.95)) on the reported outcome scale",
        },
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    write_fixture(args.output.resolve())


if __name__ == "__main__":
    main()
