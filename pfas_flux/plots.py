"""Plotting utilities for the PFAS technical note workflow."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def _save_figure(fig: plt.Figure, out_dir: str | Path, stem: str) -> List[Path]:
    """Save a figure as PNG and SVG."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    paths = []
    for ext in ("png", "svg"):
        path = out / f"{stem}.{ext}"
        fig.savefig(path, bbox_inches="tight", dpi=300)
        paths.append(path)
    plt.close(fig)
    return paths


def plot_control_plane_group_bars(
    group_summary: pd.DataFrame,
    out_dir: str | Path,
    selected_group_id: str = "G3",
) -> List[Path]:
    """Bar plot of measured and F_comp-adjusted group mass discharge by plane."""
    df = group_summary[group_summary["group_id"] == selected_group_id].copy()
    if df.empty:
        return []
    df = df.sort_values("plane_id")
    x = np.arange(len(df))
    width = 0.38
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.bar(
        x - width / 2,
        df["measured_group_mass_discharge_g_d"],
        width,
        label="Measured target-suite load",
    )
    ax.bar(
        x + width / 2,
        df["F_comp_adjusted_mass_discharge_g_d"],
        width,
        label="F_comp scenario-adjusted load",
    )
    ax.set_xticks(x)
    ax.set_xticklabels(df["plane_id"].tolist())
    ax.set_ylabel("Mass discharge (g/d)")
    ax.set_xlabel("Control plane")
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    return _save_figure(fig, out_dir, f"figure_group_{selected_group_id}_mass_discharge_by_plane")


def plot_analyte_bars(plane_analyte: pd.DataFrame, out_dir: str | Path) -> List[Path]:
    """Grouped bar plot of analyte-specific mass discharge by plane."""
    if plane_analyte.empty:
        return []
    pivot = plane_analyte.pivot_table(
        index="plane_id", columns="analyte", values="mass_discharge_g_d", aggfunc="sum"
    ).fillna(0.0)
    pivot = pivot.sort_index()
    fig, ax = plt.subplots(figsize=(8, 4.8))
    bottom = np.zeros(len(pivot))
    x = np.arange(len(pivot))
    for analyte in pivot.columns:
        values = pivot[analyte].values
        ax.bar(x, values, bottom=bottom, label=analyte)
        bottom += values
    ax.set_xticks(x)
    ax.set_xticklabels(pivot.index.tolist())
    ax.set_ylabel("Mass discharge (g/d)")
    ax.set_xlabel("Control plane")
    ax.legend(ncol=3, fontsize=8)
    ax.grid(axis="y", alpha=0.3)
    return _save_figure(fig, out_dir, "figure_stacked_analyte_mass_discharge_by_plane")


def plot_uncertainty_group_errorbars(
    uncertainty_group: pd.DataFrame,
    out_dir: str | Path,
    selected_group_id: str = "G3",
) -> List[Path]:
    """P05-P95 uncertainty plot for selected group by control plane."""
    df = uncertainty_group[uncertainty_group["group_id"] == selected_group_id].copy()
    if df.empty:
        return []
    df = df.sort_values("plane_id")
    x = np.arange(len(df))
    median = df["F_comp_adjusted_mass_discharge_g_d_p50"].values
    low = median - df["F_comp_adjusted_mass_discharge_g_d_p05"].values
    high = df["F_comp_adjusted_mass_discharge_g_d_p95"].values - median
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.errorbar(x, median, yerr=[low, high], fmt="o", capsize=5)
    ax.set_xticks(x)
    ax.set_xticklabels(df["plane_id"].tolist())
    ax.set_ylabel("F_comp-adjusted mass discharge (g/d)")
    ax.set_xlabel("Control plane")
    ax.grid(axis="y", alpha=0.3)
    return _save_figure(fig, out_dir, f"figure_uncertainty_{selected_group_id}_by_plane")


def plot_receptor_probability(
    receptor_uncertainty: pd.DataFrame,
    out_dir: str | Path,
) -> List[Path]:
    """Plot receptor exceedance probability using synthetic placeholder allowable loads."""
    if receptor_uncertainty.empty:
        return []
    df = receptor_uncertainty.sort_values("receptor_id")
    x = np.arange(len(df))
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.bar(x, df["P_exceed_placeholder_allowable"].values)
    ax.set_xticks(x)
    ax.set_xticklabels(df["receptor_id"].tolist())
    ax.set_ylim(0, 1)
    ax.set_ylabel("P(load > placeholder allowable load)")
    ax.set_xlabel("Synthetic receptor")
    ax.grid(axis="y", alpha=0.3)
    return _save_figure(fig, out_dir, "figure_receptor_exceedance_probability")


def plot_all(
    deterministic: dict,
    uncertainty: dict,
    out_dir: str | Path,
    selected_group_id: str = "G3",
) -> List[Path]:
    """Generate the default technical note figure set."""
    paths: List[Path] = []
    paths.extend(
        plot_control_plane_group_bars(
            deterministic["control_plane_group_summary"], out_dir, selected_group_id
        )
    )
    paths.extend(plot_analyte_bars(deterministic["control_plane_analyte_summary"], out_dir))
    paths.extend(
        plot_uncertainty_group_errorbars(
            uncertainty["uncertainty_control_plane_group_summary"], out_dir, selected_group_id
        )
    )
    paths.extend(plot_receptor_probability(uncertainty["uncertainty_receptor_summary"], out_dir))
    return paths
