#!/usr/bin/env python
"""Compute a compact sensitivity/data-priority screen for PFAS-FLUX.

The screen is intentionally lightweight. It recomputes the selected receptor
load ratio using the same triangular/log-triangular parameter ranges used in
the benchmark, then ranks interpretable aggregate drivers by Spearman rank
correlation. It is a screening metric, not a substitute for a calibrated
variance-decomposition analysis.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def triangular(rng: np.random.Generator, low: float, mode: float, high: float, n: int) -> np.ndarray:
    low, mode, high = float(low), float(mode), float(high)
    low = max(low, 0.0)
    high = max(high, low)
    mode = min(max(mode, low), high)
    if high == low:
        return np.full(n, low)
    return rng.triangular(low, mode, high, n)


def log_triangular(rng: np.random.Generator, low: float, mode: float, high: float, n: int) -> np.ndarray:
    low, mode, high = float(low), float(mode), float(high)
    if low <= 0 or mode <= 0 or high <= 0:
        return triangular(rng, low, mode, high, n)
    return 10 ** rng.triangular(np.log10(low), np.log10(mode), np.log10(high), n)


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute PFAS-FLUX sensitivity/data-priority screen.")
    parser.add_argument("--input-dir", default="data/synthetic_case_archive")
    parser.add_argument("--output-file", default="outputs/tables/sensitivity_rank_summary.csv")
    parser.add_argument("--selected-plane-id", default="CP3")
    parser.add_argument("--selected-receptor-id", default="R3")
    parser.add_argument("--selected-group-id", default="G3")
    parser.add_argument("--n-realizations", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20260703)
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    n = args.n_realizations
    rng = np.random.default_rng(args.seed)

    cp = pd.read_csv(input_dir / "control_plane_cells.csv")
    kz = pd.read_csv(input_dir / "assumed_K_zones.csv")
    rec = pd.read_csv(input_dir / "receptor_flows.csv")
    groups = pd.read_csv(input_dir / "analyte_groups.csv")

    group = groups.loc[groups["group_id"] == args.selected_group_id].iloc[0]
    included = [x.strip().upper().replace(" ", "") for x in str(group["included_analytes"]).replace(",", ";").split(";")]
    # The benchmark sensitivity currently supports PFOS + PFHxS and groups that contain these analytes.
    analytes = []
    if "PFOS" in included:
        analytes.append("PFOS_ug_L")
    if "PFHXS" in included:
        analytes.append("PFHxS_ug_L")
    if not analytes:
        raise ValueError("Sensitivity script currently requires a selected group containing PFOS and/or PFHxS.")

    cp_sel = cp.loc[cp["plane_id"] == args.selected_plane_id].merge(kz, on="k_zone_id", suffixes=("", "_kzone"))
    if cp_sel.empty:
        raise ValueError(f"No control-plane cells found for {args.selected_plane_id}")

    measured_md = np.zeros(n)
    qA_terms = []
    c_weight_terms = []
    K_terms = []
    grad_terms = []
    thick_terms = []

    for _, cell in cp_sel.iterrows():
        K = log_triangular(rng, cell["K_low_m_d"], cell["K_central_m_d"], cell["K_high_m_d"], n)
        grad = triangular(rng, cell["gradient_low"], cell["gradient_normal_central"], cell["gradient_high"], n)
        thick = triangular(rng, cell["effective_thickness_low_m"], cell["effective_thickness_central_m"], cell["effective_thickness_high_m"], n)
        width = float(cell["y_max_m"] - cell["y_min_m"])
        area = width * thick
        qA = K * grad * area
        csum = np.zeros(n)
        for analyte in analytes:
            c0 = float(cell[analyte])
            csum += triangular(rng, 0.5 * c0, c0, 2.0 * c0, n)
        measured_md += csum * qA * 1.0e-3
        qA_terms.append(qA)
        c_weight_terms.append(csum * qA)
        K_terms.append(K)
        grad_terms.append(grad)
        thick_terms.append(thick)

    fcomp = triangular(rng, group["F_comp_low"], group["F_comp_central"], group["F_comp_high"], n)
    adjusted_md = measured_md * fcomp

    receptor = rec.loc[rec["receptor_id"] == args.selected_receptor_id].iloc[0]
    Q_mix = triangular(rng, receptor["Q_low_m3_d"], receptor["Q_central_m3_d"], receptor["Q_high_m3_d"], n)
    f_mix = triangular(rng, receptor["mixing_factor_low"], receptor["mixing_factor_central"], receptor["mixing_factor_high"], n)
    criterion_ug_L = float(receptor["criterion_placeholder_PFOS_plus_PFHxS_ng_L"]) / 1000.0
    allowable = criterion_ug_L * Q_mix * f_mix * 1.0e-3
    load_ratio = adjusted_md / allowable

    hydraulic_capacity = np.sum(qA_terms, axis=0)
    flux_weighted_concentration = np.sum(c_weight_terms, axis=0) / hydraulic_capacity
    drivers = {
        "receptor flow": Q_mix,
        "hydraulic capacity": hydraulic_capacity,
        "PFOS+PFHxS concentration": flux_weighted_concentration,
        "normal gradient": np.mean(grad_terms, axis=0),
        "F_comp factor": fcomp,
        "mixing factor": f_mix,
        "saturated thickness": np.mean(thick_terms, axis=0),
    }

    rows = []
    for name, values in drivers.items():
        rho, _ = spearmanr(values, load_ratio)
        rows.append({
            "output": f"{args.selected_receptor_id} load ratio",
            "driver": name,
            "spearman_r": float(rho),
            "abs_spearman_r": float(abs(rho)),
            "n_realizations": n,
            "note": "screening data-priority metric; not a formal variance decomposition",
        })

    out = pd.DataFrame(rows).sort_values("abs_spearman_r", ascending=False)
    output_file = Path(args.output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output_file, index=False)
    print(f"Wrote {output_file}")


if __name__ == "__main__":
    main()
