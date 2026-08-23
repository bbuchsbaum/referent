# /// script
# requires-python = "==3.12.11"
# dependencies = [
#   "arviz==0.23.4",
#   "numpy==2.4.6",
#   "nutpie==0.16.8",
#   "pandas==3.0.5",
#   "pcntoolkit==1.3.0",
#   "pymc==5.28.5",
#   "scipy==1.18.1",
#   "xarray==2026.7.0",
# ]
# ///
"""Run stochastic HBR and new-site transfer evidence.

This is a release/manual runner, not an ordinary package test. It emits both
PCNtoolkit's draw-averaged outputs and proper posterior-predictive mixture
quantities, convergence receipts, and site-stratified metrics.
"""

from __future__ import annotations

import argparse
import importlib.metadata
import json
from pathlib import Path

import arviz as az
import numpy as np
import pandas as pd
import pymc as pm
from scipy import optimize, special, stats

from pcntoolkit import NormativeModel
from pcntoolkit.dataio.norm_data import NormData
from pcntoolkit.regression_model.hbr import HBR


def simulate(seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows: list[dict[str, object]] = []
    shifts = {"site-1": -0.7, "site-2": -0.2, "site-3": 0.25, "site-4": 0.7, "site-5": 1.1}
    scales = {"site-1": 0.65, "site-2": 0.75, "site-3": 0.85, "site-4": 0.95, "site-5": 1.1}
    for site in shifts:
        n = 80 if site != "site-5" else 100
        x = rng.uniform(-2, 2, n)
        y = 0.5 + 0.8 * x + shifts[site] + scales[site] * rng.normal(size=n)
        if site == "site-5":
            split = np.repeat(["adapt", "transport_test"], [25, 75])
        else:
            split = np.repeat(["reference_train", "observed_site_test"], [50, 30])
        for i in range(n):
            rows.append(
                {
                    "row_id": f"{site}-{split[i]}-{i:04d}",
                    "site": site,
                    "split": split[i],
                    "x": x[i],
                    "y": y[i],
                }
            )
    return pd.DataFrame(rows)


def norm_data(name: str, frame: pd.DataFrame) -> NormData:
    return NormData.from_dataframe(
        name,
        frame.copy(),
        covariates=["x"],
        batch_effects=["site"],
        response_vars=["y"],
        subject_ids="row_id",
    )


def convergence(model: NormativeModel, stage: str) -> tuple[pd.DataFrame, dict[str, object]]:
    idata = model["y"].idata
    # The fitted idata does not contain the temporary per-subject variables
    # created during prediction, so the full posterior is the convergence set.
    summary = az.summary(idata, fmt="wide")
    summary.insert(0, "parameter", summary.index.astype(str))
    summary.insert(0, "stage", stage)
    divergences = int(idata.sample_stats["diverging"].sum().item())
    max_rhat = float(np.nanmax(summary["r_hat"].to_numpy()))
    min_bulk = float(np.nanmin(summary["ess_bulk"].to_numpy()))
    min_tail = float(np.nanmin(summary["ess_tail"].to_numpy()))
    receipt = {
        "stage": stage,
        "divergences": divergences,
        "max_rhat": max_rhat,
        "min_ess_bulk": min_bulk,
        "min_ess_tail": min_tail,
        "pass": bool(divergences == 0 and max_rhat <= 1.01 and min_bulk >= 400 and min_tail >= 400),
    }
    return summary.reset_index(drop=True), receipt


def posterior_normal_draws(model: NormativeModel, frame: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    data = norm_data("posterior-mixture", frame)
    model.preprocess(data)
    response_data = data.sel({"response_vars": "y"})
    X, be, _, Y, _ = model.extract_data(response_data)
    hbr = model["y"]
    pymc_model = hbr.likelihood.create_model_with_data(X, be, hbr.be_maps, Y)
    params = hbr.likelihood.compile_params(pymc_model, X, be, hbr.be_maps, Y)
    var_names = [f"{name}_per_subject" for name in params]
    with pymc_model:
        for name, (value, dims) in params.items():
            pm.Deterministic(f"{name}_per_subject", value, dims=dims)
        posterior = pm.sample_posterior_predictive(
            hbr.idata, extend_inferencedata=False, var_names=var_names, progressbar=False
        )
    extracted = az.extract(posterior, "posterior_predictive", var_names=var_names)
    arrays = {
        name: hbr.extract_and_reshape(extracted, len(frame), f"{name}_per_subject").values
        for name in params
    }
    return arrays["mu"], arrays["sigma"]


def mixture_quantile(mu: np.ndarray, sigma: np.ndarray, p: float) -> np.ndarray:
    out = np.empty(mu.shape[0])
    for i in range(mu.shape[0]):
        lo = np.min(mu[i] - 12 * sigma[i])
        hi = np.max(mu[i] + 12 * sigma[i])
        out[i] = optimize.brentq(
            lambda q: np.mean(stats.norm.cdf(q, mu[i], sigma[i])) - p,
            lo,
            hi,
        )
    return out


def predict_with_mixture(
    model: NormativeModel, frame: pd.DataFrame, lane: str, outcome_mean: float, outcome_sd: float
) -> pd.DataFrame:
    predicted = model.predict(norm_data(f"{lane}-test", frame))
    reported_z = predicted["Z"].sel(response_vars="y").values
    reported_q05 = predicted["centiles"].sel(response_vars="y", centile=0.05).values
    reported_q50 = predicted["centiles"].sel(response_vars="y", centile=0.5).values
    reported_q95 = predicted["centiles"].sel(response_vars="y", centile=0.95).values
    reported_logp = predicted["logp"].sel(response_vars="y").values - np.log(outcome_sd)

    mu, sigma = posterior_normal_draws(model, frame)
    y_scaled = (frame["y"].to_numpy() - outcome_mean) / outcome_sd
    component_z = (y_scaled[:, None] - mu) / sigma
    mixture_cdf = np.mean(stats.norm.cdf(component_z), axis=1)
    mixture_logp = special.logsumexp(stats.norm.logpdf(component_z) - np.log(sigma), axis=1) - np.log(
        sigma.shape[1]
    ) - np.log(outcome_sd)
    mixture_q05 = outcome_mean + outcome_sd * mixture_quantile(mu, sigma, 0.05)
    mixture_q50 = outcome_mean + outcome_sd * mixture_quantile(mu, sigma, 0.5)
    mixture_q95 = outcome_mean + outcome_sd * mixture_quantile(mu, sigma, 0.95)
    return pd.DataFrame(
        {
            "lane": lane,
            "row_id": frame["row_id"].to_numpy(),
            "site": frame["site"].to_numpy(),
            "observed": frame["y"].to_numpy(),
            "reported_mean_draw_z": reported_z,
            "reported_mean_draw_log_density": reported_logp,
            "reported_mean_draw_q05": reported_q05,
            "reported_mean_draw_q50": reported_q50,
            "reported_mean_draw_q95": reported_q95,
            "mixture_centile": mixture_cdf,
            "mixture_z": stats.norm.ppf(mixture_cdf),
            "mixture_log_density": mixture_logp,
            "mixture_q05": mixture_q05,
            "mixture_q50": mixture_q50,
            "mixture_q95": mixture_q95,
        }
    )


def site_summary(predictions: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (lane, site), group in predictions.groupby(["lane", "site"], sort=True):
        u = group["mixture_centile"].to_numpy()
        grid = np.array([0.05, 0.25, 0.5, 0.75, 0.95])
        rows.append(
            {
                "method": "pcntoolkit_hbr_mixture",
                "lane": lane,
                "site": site,
                "n": len(group),
                "mean_log_score": group["mixture_log_density"].mean(),
                "median_rmse": np.sqrt(np.mean((group["observed"] - group["mixture_q50"]) ** 2)),
                "coverage90": np.mean(
                    (group["observed"] >= group["mixture_q05"])
                    & (group["observed"] <= group["mixture_q95"])
                ),
                "mace": np.mean([abs(q - np.mean(u <= q)) for q in grid]),
                "mean_z": group["mixture_z"].mean(),
                "var_z": group["mixture_z"].var(ddof=1),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--draws", type=int, default=1500)
    parser.add_argument("--tune", type=int, default=500)
    parser.add_argument("--chains", type=int, default=4)
    parser.add_argument("--cores", type=int, default=4)
    parser.add_argument("--nuts-sampler", default="nutpie")
    parser.add_argument("--seed", type=int, default=20260823)
    parser.add_argument("--progressbar", action="store_true")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    data = simulate(args.seed)
    reference_train = data[data["split"] == "reference_train"].copy()
    observed_test = data[data["split"] == "observed_site_test"].copy()
    adaptation = data[data["split"] == "adapt"].copy()
    transport_test = data[data["split"] == "transport_test"].copy()
    outcome_mean = float(reference_train["y"].mean())
    outcome_sd = float(reference_train["y"].std(ddof=0))

    hbr = HBR(
        draws=args.draws,
        tune=args.tune,
        chains=args.chains,
        cores=args.cores,
        nuts_sampler=args.nuts_sampler,
        progressbar=args.progressbar,
    )
    base = NormativeModel(
        hbr, savemodel=False, evaluate_model=False, saveresults=False, saveplots=False,
        inscaler="standardize", outscaler="standardize", name="hbr-site-evidence"
    )
    base.fit(norm_data("reference-train", reference_train))
    base_convergence, base_receipt = convergence(base, "base")
    observed_predictions = predict_with_mixture(
        base, observed_test, "observed_site", outcome_mean, outcome_sd
    )

    transferred = base.transfer(
        norm_data("site-5-adaptation", adaptation),
        save_dir=str(args.output / "transfer_model"),
        draws=args.draws,
        tune=args.tune,
        chains=args.chains,
        cores=args.cores,
        nuts_sampler=args.nuts_sampler,
        progressbar=args.progressbar,
    )
    transfer_convergence, transfer_receipt = convergence(transferred, "transfer")
    transport_predictions = predict_with_mixture(
        transferred, transport_test, "transferred_site", outcome_mean, outcome_sd
    )
    predictions = pd.concat([observed_predictions, transport_predictions], ignore_index=True)
    convergence_table = pd.concat([base_convergence, transfer_convergence], ignore_index=True)

    data.to_csv(args.output / "site_data.csv", index=False, float_format="%.17g")
    predictions.to_csv(args.output / "pcntoolkit_predictions.csv", index=False, float_format="%.17g")
    site_summary(predictions).to_csv(args.output / "pcntoolkit_site_summary.csv", index=False, float_format="%.17g")
    convergence_table.to_csv(args.output / "hbr_convergence.csv", index=False, float_format="%.17g")
    receipt = {
        "comparator": "pcntoolkit",
        "version": importlib.metadata.version("pcntoolkit"),
        "draws": args.draws,
        "tune": args.tune,
        "chains": args.chains,
        "cores": args.cores,
        "nuts_sampler": args.nuts_sampler,
        "requested_seed": args.seed,
        "hbr_public_random_seed_argument": False,
        "aggregation": {
            "reported_z": "mean of draw-specific z",
            "reported_quantile": "mean of draw-specific quantiles",
            "reported_logp": "mean of draw-specific log densities",
            "mixture_outputs": "CDF and density averaged over posterior draws, then transformed",
        },
        "convergence_thresholds": {
            "max_rhat": 1.01,
            "min_ess_bulk": 400,
            "min_ess_tail": 400,
            "divergences": 0,
        },
        "stages": [base_receipt, transfer_receipt],
        "all_stages_pass": bool(base_receipt["pass"] and transfer_receipt["pass"]),
    }
    (args.output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
