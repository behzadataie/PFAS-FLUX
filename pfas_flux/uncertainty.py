"""Monte Carlo uncertainty calculations for the PFAS technical note workflow."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd

from .calculations import parse_analyte_group
from .schema import ANALYTE_COLUMNS, COLUMN_ANALYTE_MAP, UNITS


@dataclass(frozen=True)
class MonteCarloOptions:
    """Monte Carlo settings for the technical note workflow."""

    n: int = 10000
    seed: int = 20260703
    concentration_low_multiplier: float = 0.5
    concentration_high_multiplier: float = 2.0
    selected_group_id: str = "G3"
    use_fcomp_adjusted_for_receptors: bool = True


def _triangular(
    rng: np.random.Generator, low: float, mode: float, high: float, n: int
) -> np.ndarray:
    """Sample a triangular distribution with defensive bounds."""
    low = float(low)
    mode = float(mode)
    high = float(high)
    if not np.isfinite(low) or not np.isfinite(mode) or not np.isfinite(high):
        raise ValueError(f"Non-finite triangular inputs: low={low}, mode={mode}, high={high}")
    low = max(low, 0.0)
    high = max(high, low)
    mode = min(max(mode, low), high)
    if high == low:
        return np.full(n, low)
    return rng.triangular(low, mode, high, size=n)


def _log_triangular(
    rng: np.random.Generator, low: float, mode: float, high: float, n: int
) -> np.ndarray:
    """Sample a log-space triangular distribution for positive hydraulic K."""
    low = float(low)
    mode = float(mode)
    high = float(high)
    if low <= 0 or mode <= 0 or high <= 0:
        return _triangular(rng, low, mode, high, n)
    log_low = np.log10(low)
    log_mode = np.log10(mode)
    log_high = np.log10(high)
    log_mode = min(max(log_mode, log_low), log_high)
    return np.power(10.0, rng.triangular(log_low, log_mode, log_high, size=n))


def _quantiles(values: np.ndarray) -> Dict[str, float]:
    """Return summary statistics for a numeric vector."""
    values = np.asarray(values, dtype=float)
    return {
        "mean": float(np.nanmean(values)),
        "p05": float(np.nanpercentile(values, 5)),
        "p50": float(np.nanpercentile(values, 50)),
        "p95": float(np.nanpercentile(values, 95)),
        "cv": float(np.nanstd(values) / np.nanmean(values)) if np.nanmean(values) else np.nan,
    }


def _sample_concentration(
    rng: np.random.Generator,
    central: float,
    n: int,
    low_multiplier: float,
    high_multiplier: float,
) -> np.ndarray:
    """Sample a concentration uncertainty distribution around the central value."""
    central = float(central)
    if not np.isfinite(central) or central <= 0:
        return np.zeros(n)
    low = central * low_multiplier
    high = central * high_multiplier
    return _triangular(rng, low, central, high, n)


def _sample_fcomp(rng: np.random.Generator, group: pd.Series, n: int) -> np.ndarray:
    return _triangular(
        rng,
        float(group.get("F_comp_low", 1.0)),
        float(group.get("F_comp_central", 1.0)),
        float(group.get("F_comp_high", group.get("F_comp_central", 1.0))),
        n,
    )


def _prepare_control_plane_cells(
    control_plane_cells: pd.DataFrame, assumed_k_zones: pd.DataFrame
) -> pd.DataFrame:
    cp = control_plane_cells.copy()
    kz = assumed_k_zones.copy()
    for df in (cp, kz):
        for col in df.columns:
            if col.endswith("_m_d") or col.endswith("_m") or col.endswith("_m2") or "gradient" in col or col in ANALYTE_COLUMNS:
                df[col] = pd.to_numeric(df[col], errors="coerce")
    merged = cp.merge(kz, on="k_zone_id", how="left", suffixes=("", "_kzone"))
    merged["cell_width_m"] = merged["y_max_m"] - merged["y_min_m"]
    merged["gradient_mode"] = merged["gradient_normal_central"].fillna(merged["gradient_central"])
    return merged


def control_plane_uncertainty(
    control_plane_cells: pd.DataFrame,
    assumed_k_zones: pd.DataFrame,
    analyte_groups: pd.DataFrame,
    options: MonteCarloOptions,
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[Tuple[str, str, str], np.ndarray]]:
    """Run Monte Carlo calculations for control-plane loads.

    Returns
    -------
    analyte_summary, group_summary, samples
        ``samples`` stores arrays by (kind, plane_id, label), where kind is
        "analyte", "group_measured" or "group_adjusted".
    """
    rng = np.random.default_rng(options.seed)
    n = options.n
    cp = _prepare_control_plane_cells(control_plane_cells, assumed_k_zones)

    samples: Dict[Tuple[str, str, str], np.ndarray] = {}

    for _, cell in cp.iterrows():
        k = _log_triangular(
            rng,
            cell["K_low_m_d"],
            cell["K_central_m_d"],
            cell["K_high_m_d"],
            n,
        )
        grad = _triangular(
            rng,
            cell["gradient_low"],
            cell["gradient_mode"],
            cell["gradient_high"],
            n,
        )
        thick = _triangular(
            rng,
            cell["effective_thickness_low_m"],
            cell["effective_thickness_central_m"],
            cell["effective_thickness_high_m"],
            n,
        )
        area = cell["cell_width_m"] * thick
        q = k * grad
        for col in ANALYTE_COLUMNS:
            if col not in cp.columns:
                continue
            conc = _sample_concentration(
                rng,
                cell.get(col, 0.0),
                n,
                options.concentration_low_multiplier,
                options.concentration_high_multiplier,
            )
            md = conc * q * area * UNITS.conversion_factor_to_g_d
            label = COLUMN_ANALYTE_MAP[col]
            key = ("analyte", str(cell["plane_id"]), label)
            samples[key] = samples.get(key, np.zeros(n)) + md

    analyte_rows = []
    for (kind, plane_id, analyte), values in samples.items():
        if kind != "analyte":
            continue
        row = {"plane_id": plane_id, "analyte": analyte}
        row.update({f"mass_discharge_g_d_{k}": v for k, v in _quantiles(values).items()})
        row["n_realizations"] = n
        analyte_rows.append(row)
    analyte_summary = pd.DataFrame(analyte_rows).sort_values(["plane_id", "analyte"])

    group_rows = []
    for _, group in analyte_groups.iterrows():
        group_id = group["group_id"]
        included_cols = parse_analyte_group(group["included_analytes"])
        included_labels = [COLUMN_ANALYTE_MAP[col] for col in included_cols]
        fcomp = _sample_fcomp(rng, group, n)
        for plane_id in sorted(cp["plane_id"].astype(str).unique()):
            measured = np.zeros(n)
            for label in included_labels:
                measured += samples.get(("analyte", plane_id, label), np.zeros(n))
            adjusted = measured * fcomp
            samples[("group_measured", plane_id, group_id)] = measured
            samples[("group_adjusted", plane_id, group_id)] = adjusted
            row = {
                "plane_id": plane_id,
                "group_id": group_id,
                "included_analytes": ";".join(included_labels),
                "excluded_information": group.get("excluded_information", ""),
                "n_realizations": n,
            }
            for prefix, values in (
                ("measured_group_mass_discharge_g_d", measured),
                ("F_comp_adjusted_mass_discharge_g_d", adjusted),
            ):
                row.update({f"{prefix}_{k}": v for k, v in _quantiles(values).items()})
            group_rows.append(row)
    group_summary = pd.DataFrame(group_rows).sort_values(["plane_id", "group_id"])
    return analyte_summary, group_summary, samples


def surface_water_uncertainty(
    surface_water_nodes: pd.DataFrame,
    analyte_groups: pd.DataFrame,
    options: MonteCarloOptions,
) -> Tuple[pd.DataFrame, Dict[Tuple[str, str, str], np.ndarray]]:
    """Run uncertainty calculations for surface-water/stormwater node loads."""
    rng = np.random.default_rng(options.seed + 17)
    n = options.n
    sw = surface_water_nodes.copy()
    for col in sw.columns:
        if col.startswith("flow_") or col in ANALYTE_COLUMNS:
            sw[col] = pd.to_numeric(sw[col], errors="coerce")

    analyte_samples: Dict[Tuple[str, str, str], np.ndarray] = {}
    for _, node in sw.iterrows():
        flow = _triangular(
            rng, node["flow_low_m3_d"], node["flow_central_m3_d"], node["flow_high_m3_d"], n
        )
        for col in ANALYTE_COLUMNS:
            c = _sample_concentration(
                rng,
                node.get(col, 0.0),
                n,
                options.concentration_low_multiplier,
                options.concentration_high_multiplier,
            )
            label = COLUMN_ANALYTE_MAP[col]
            load = c * flow * UNITS.conversion_factor_to_g_d
            analyte_samples[("surface_analyte", str(node["node_id"]), label)] = load

    rows = []
    samples: Dict[Tuple[str, str, str], np.ndarray] = {}
    samples.update(analyte_samples)
    for _, group in analyte_groups.iterrows():
        group_id = group["group_id"]
        included_cols = parse_analyte_group(group["included_analytes"])
        labels = [COLUMN_ANALYTE_MAP[col] for col in included_cols]
        fcomp = _sample_fcomp(rng, group, n)
        for node_id in sorted(sw["node_id"].astype(str).unique()):
            measured = np.zeros(n)
            for label in labels:
                measured += analyte_samples.get(("surface_analyte", node_id, label), np.zeros(n))
            adjusted = measured * fcomp
            samples[("surface_group_measured", node_id, group_id)] = measured
            samples[("surface_group_adjusted", node_id, group_id)] = adjusted
            node_type = sw.loc[sw["node_id"].astype(str) == node_id, "node_type"].iloc[0]
            row = {
                "node_id": node_id,
                "node_type": node_type,
                "group_id": group_id,
                "included_analytes": ";".join(labels),
                "n_realizations": n,
                "unit_note": "g/d if flow is m3/d; g/event if flow is an event volume",
            }
            for prefix, values in (
                ("measured_group_load", measured),
                ("F_comp_adjusted_group_load", adjusted),
            ):
                row.update({f"{prefix}_{k}": v for k, v in _quantiles(values).items()})
            rows.append(row)
    return pd.DataFrame(rows).sort_values(["node_id", "group_id"]), samples


def well_uncertainty(
    wells: pd.DataFrame,
    receptor_flows: pd.DataFrame,
    analyte_groups: pd.DataFrame,
    options: MonteCarloOptions,
) -> Dict[Tuple[str, str, str], np.ndarray]:
    """Run uncertainty calculations for receptor links to wells."""
    rng = np.random.default_rng(options.seed + 23)
    n = options.n
    wells_num = wells.copy()
    for col in ANALYTE_COLUMNS:
        if col in wells_num.columns:
            wells_num[col] = pd.to_numeric(wells_num[col], errors="coerce")
    well_map = wells_num.set_index("well_id")
    samples: Dict[Tuple[str, str, str], np.ndarray] = {}
    for _, rec in receptor_flows.iterrows():
        link = str(rec.get("linked_plane_or_node", ""))
        if link not in well_map.index:
            continue
        q = _triangular(
            rng,
            float(rec["Q_low_m3_d"]),
            float(rec["Q_central_m3_d"]),
            float(rec["Q_high_m3_d"]),
            n,
        )
        well = well_map.loc[link]
        for _, group in analyte_groups.iterrows():
            group_id = group["group_id"]
            cols = parse_analyte_group(group["included_analytes"])
            c_sum = np.zeros(n)
            for col in cols:
                c_sum += _sample_concentration(
                    rng,
                    float(well.get(col, 0.0)),
                    n,
                    options.concentration_low_multiplier,
                    options.concentration_high_multiplier,
                )
            measured = c_sum * q * UNITS.conversion_factor_to_g_d
            fcomp = _sample_fcomp(rng, group, n)
            adjusted = measured * fcomp
            samples[("well_group_measured", link, group_id)] = measured
            samples[("well_group_adjusted", link, group_id)] = adjusted
    return samples


def receptor_uncertainty(
    receptor_flows: pd.DataFrame,
    control_samples: Dict[Tuple[str, str, str], np.ndarray],
    surface_samples: Dict[Tuple[str, str, str], np.ndarray],
    well_samples: Dict[Tuple[str, str, str], np.ndarray],
    options: MonteCarloOptions,
) -> pd.DataFrame:
    """Estimate P(current load > allowable load) by receptor."""
    rng = np.random.default_rng(options.seed + 31)
    n = options.n
    rows = []
    group_id = options.selected_group_id
    for _, rec in receptor_flows.iterrows():
        link = str(rec["linked_plane_or_node"])
        q = _triangular(rng, rec["Q_low_m3_d"], rec["Q_central_m3_d"], rec["Q_high_m3_d"], n)
        mix = _triangular(
            rng,
            rec["mixing_factor_low"],
            rec["mixing_factor_central"],
            rec["mixing_factor_high"],
            n,
        )
        criterion_ug_L = float(rec["criterion_placeholder_PFOS_plus_PFHxS_ng_L"]) / 1000.0
        allowable = criterion_ug_L * q * mix * UNITS.conversion_factor_to_g_d
        if link.startswith("CP"):
            kind = "group_adjusted" if options.use_fcomp_adjusted_for_receptors else "group_measured"
            load = control_samples.get((kind, link, group_id), np.full(n, np.nan))
            source = "control-plane mass discharge"
        elif link.startswith("SW"):
            kind = (
                "surface_group_adjusted"
                if options.use_fcomp_adjusted_for_receptors
                else "surface_group_measured"
            )
            load = surface_samples.get((kind, link, group_id), np.full(n, np.nan))
            source = "surface-water or stormwater load"
        else:
            kind = "well_group_adjusted" if options.use_fcomp_adjusted_for_receptors else "well_group_measured"
            load = well_samples.get((kind, link, group_id), np.full(n, np.nan))
            source = "extracted-water load"
        ratio = load / allowable
        row = {
            "receptor_id": rec["receptor_id"],
            "receptor_type": rec["receptor_type"],
            "linked_plane_or_node": link,
            "group_id": group_id,
            "load_source": source,
            "load_type": "F_comp adjusted scenario"
            if options.use_fcomp_adjusted_for_receptors
            else "measured target suite",
            "P_exceed_placeholder_allowable": float(np.nanmean(load > allowable)),
            "n_realizations": n,
            "warning": "allowable load uses placeholder synthetic criterion; not a regulatory decision",
        }
        for prefix, values in (("current_load_g_d", load), ("allowable_load_g_d", allowable), ("load_ratio", ratio)):
            row.update({f"{prefix}_{k}": v for k, v in _quantiles(values).items()})
        rows.append(row)
    return pd.DataFrame(rows).sort_values("receptor_id")


def uncertainty_driver_proxy(
    control_plane_cells: pd.DataFrame,
    assumed_k_zones: pd.DataFrame,
    analyte_groups: pd.DataFrame,
    options: MonteCarloOptions,
) -> pd.DataFrame:
    """Create a simple value-of-information proxy for uncertainty drivers.

    This is not a formal variance decomposition. It ranks the input ranges that
    are most likely to drive uncertainty in a sparse-data technical note example.
    A formal Sobol or regression-based sensitivity analysis can be added later.
    """
    cp = _prepare_control_plane_cells(control_plane_cells, assumed_k_zones)
    rows = []
    for _, group in analyte_groups.iterrows():
        group_id = group["group_id"]
        f_low = float(group.get("F_comp_low", 1.0))
        f_high = float(group.get("F_comp_high", 1.0))
        for plane_id, sub in cp.groupby("plane_id"):
            scores = []
            for _, cell in sub.iterrows():
                scores.append(("hydraulic conductivity K", np.log10(cell["K_high_m_d"] / cell["K_low_m_d"])))
                scores.append(("hydraulic gradient", np.log10(cell["gradient_high"] / cell["gradient_low"])))
                scores.append(
                    (
                        "control-plane thickness/area",
                        np.log10(cell["effective_thickness_high_m"] / cell["effective_thickness_low_m"]),
                    )
                )
            scores.append(
                (
                    "concentration interpolation/analytical variability",
                    np.log10(options.concentration_high_multiplier / options.concentration_low_multiplier),
                )
            )
            if f_low > 0:
                scores.append(("analyte completeness factor F_comp", np.log10(f_high / f_low)))
            scores_df = pd.DataFrame(scores, columns=["driver", "log10_range_score"])
            ranked = scores_df.groupby("driver", as_index=False)["log10_range_score"].max()
            ranked = ranked.sort_values("log10_range_score", ascending=False)
            rows.append(
                {
                    "plane_id": plane_id,
                    "group_id": group_id,
                    "dominant_driver_proxy": ranked.iloc[0]["driver"],
                    "dominant_log10_range_score": float(ranked.iloc[0]["log10_range_score"]),
                    "second_driver_proxy": ranked.iloc[1]["driver"] if len(ranked) > 1 else "",
                    "second_log10_range_score": float(ranked.iloc[1]["log10_range_score"]) if len(ranked) > 1 else np.nan,
                    "note": "screening proxy only; use formal sensitivity analysis for real decisions",
                }
            )
    return pd.DataFrame(rows).sort_values(["plane_id", "group_id"])


def monte_carlo_workflow(tables: Dict[str, pd.DataFrame], options: MonteCarloOptions) -> Dict[str, pd.DataFrame]:
    """Run all technical note Monte Carlo calculations."""
    analyte_summary, group_summary, control_samples = control_plane_uncertainty(
        tables["control_plane_cells"], tables["assumed_K_zones"], tables["analyte_groups"], options
    )
    sw_group_summary, surface_samples = surface_water_uncertainty(
        tables["surface_water_nodes"], tables["analyte_groups"], options
    )
    well_samples = well_uncertainty(
        tables["wells"], tables["receptor_flows"], tables["analyte_groups"], options
    )
    receptor_summary = receptor_uncertainty(
        tables["receptor_flows"], control_samples, surface_samples, well_samples, options
    )
    drivers = uncertainty_driver_proxy(
        tables["control_plane_cells"], tables["assumed_K_zones"], tables["analyte_groups"], options
    )
    return {
        "uncertainty_control_plane_analyte_summary": analyte_summary,
        "uncertainty_control_plane_group_summary": group_summary,
        "uncertainty_surface_water_group_summary": sw_group_summary,
        "uncertainty_receptor_summary": receptor_summary,
        "uncertainty_driver_proxy": drivers,
    }
