"""Deterministic calculations for the PFAS technical note workflow."""

from __future__ import annotations

import re
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd

from .schema import ANALYTE_COLUMNS, ANALYTE_COLUMN_MAP, COLUMN_ANALYTE_MAP, UNITS


def _to_number(series: pd.Series) -> pd.Series:
    """Convert a pandas Series to numeric values with NaN for bad entries."""
    return pd.to_numeric(series, errors="coerce")


def normalise_analyte_token(token: str) -> str | None:
    """Convert an analyte token from analyte_groups.csv to a concentration column.

    The CSV deliberately stores human-readable text such as "6:2 FTS". This
    function maps that text to the machine-readable concentration column.
    Unknown text is ignored rather than treated as an error because some group
    definitions may include descriptive words such as "scenario terminal
    products".
    """
    clean = token.strip()
    if not clean:
        return None
    clean = clean.replace("and", ";") if clean.lower() == "and" else clean
    key = clean.upper().replace(" ", "")
    lookup = {k.upper().replace(" ", ""): v for k, v in ANALYTE_COLUMN_MAP.items()}
    return lookup.get(key)


def parse_analyte_group(included_analytes: str) -> List[str]:
    """Parse a human-readable analyte group into concentration columns."""
    # Split on semicolons or commas, and also capture the common phrase
    # "6:2 FTS and scenario terminal products" by keeping known analytes only.
    pieces = re.split(r"[;,]", str(included_analytes))
    columns: List[str] = []
    for piece in pieces:
        # First try whole piece.
        col = normalise_analyte_token(piece)
        if col and col not in columns:
            columns.append(col)
            continue
        # Then search for known analyte names inside longer descriptive text.
        for label, candidate_col in ANALYTE_COLUMN_MAP.items():
            if label.upper() in piece.upper() and candidate_col not in columns:
                columns.append(candidate_col)
    return columns


def control_plane_cell_results(
    control_plane_cells: pd.DataFrame, assumed_k_zones: pd.DataFrame
) -> pd.DataFrame:
    """Calculate deterministic central-case mass flux and mass discharge by cell.

    Formula implemented:
        q_i = K_i * I_i
        Md_a,i = C_a,i * q_i * A_i * 1e-3

    where C is in ug/L, q is in m/d, A is in m2, and Md is in g/d.
    """
    k = assumed_k_zones.copy()
    cp = control_plane_cells.copy()

    numeric_cols = [
        "K_central_m_d",
        "gradient_central",
        "area_m2",
        "gradient_normal_central",
        "y_min_m",
        "y_max_m",
    ] + ANALYTE_COLUMNS
    for df in (k, cp):
        for col in numeric_cols:
            if col in df.columns:
                df[col] = _to_number(df[col])

    merged = cp.merge(k, on="k_zone_id", how="left", suffixes=("", "_kzone"))
    if merged["K_central_m_d"].isna().any():
        bad = merged.loc[merged["K_central_m_d"].isna(), "k_zone_id"].unique()
        raise ValueError(f"K zones missing K_central_m_d: {bad}")

    merged["gradient_used"] = merged["gradient_normal_central"].fillna(
        merged["gradient_central"]
    )
    merged["q_m_d"] = merged["K_central_m_d"] * merged["gradient_used"]
    merged["Q_m3_d"] = merged["q_m_d"] * merged["area_m2"]

    rows = []
    id_cols = [
        "plane_id",
        "cell_id",
        "plane_type",
        "k_zone_id",
        "x_m",
        "y_min_m",
        "y_max_m",
        "area_m2",
        "K_central_m_d",
        "gradient_used",
        "q_m_d",
        "Q_m3_d",
    ]
    for _, row in merged.iterrows():
        for col in ANALYTE_COLUMNS:
            if col not in merged.columns:
                continue
            concentration = row[col]
            if pd.isna(concentration):
                continue
            md_g_d = concentration * row["q_m_d"] * row["area_m2"] * UNITS.conversion_factor_to_g_d
            j_g_m2_d = concentration * row["q_m_d"] * UNITS.conversion_factor_to_g_d
            out = {name: row.get(name, np.nan) for name in id_cols}
            out.update(
                {
                    "analyte": COLUMN_ANALYTE_MAP[col],
                    "concentration_ug_L": concentration,
                    "mass_flux_g_m2_d": j_g_m2_d,
                    "mass_discharge_g_d": md_g_d,
                    "mass_discharge_mg_d": md_g_d * 1000.0,
                    "calculation_basis": "deterministic central case",
                    "formula": "Md_g_d = C_ug_L * K_m_d * gradient * area_m2 * 1e-3",
                }
            )
            rows.append(out)
    return pd.DataFrame(rows)


def plane_analyte_summary(cell_results: pd.DataFrame) -> pd.DataFrame:
    """Aggregate deterministic mass discharge by control plane and analyte."""
    group_cols = ["plane_id", "plane_type", "analyte"]
    out = (
        cell_results.groupby(group_cols, as_index=False)
        .agg(
            mass_discharge_g_d=("mass_discharge_g_d", "sum"),
            mass_discharge_mg_d=("mass_discharge_mg_d", "sum"),
            contributing_cells=("cell_id", "nunique"),
            mean_q_m_d=("q_m_d", "mean"),
            total_Q_m3_d=("Q_m3_d", "sum"),
        )
        .sort_values(["plane_id", "analyte"])
    )
    return out


def group_mass_discharge(
    plane_analyte: pd.DataFrame, analyte_groups: pd.DataFrame
) -> pd.DataFrame:
    """Calculate measured and completeness-adjusted group mass discharge.

    The adjusted value is a scenario calculation. It is not a measurement.
    """
    rows = []
    for _, group in analyte_groups.iterrows():
        group_id = group["group_id"]
        included_cols = parse_analyte_group(group["included_analytes"])
        included_labels = [COLUMN_ANALYTE_MAP[col] for col in included_cols]
        f_low = float(group.get("F_comp_low", 1.0))
        f_central = float(group.get("F_comp_central", 1.0))
        f_high = float(group.get("F_comp_high", f_central))
        subset = plane_analyte[plane_analyte["analyte"].isin(included_labels)].copy()
        for (plane_id, plane_type), sdf in subset.groupby(["plane_id", "plane_type"]):
            measured_g_d = float(sdf["mass_discharge_g_d"].sum())
            rows.append(
                {
                    "plane_id": plane_id,
                    "plane_type": plane_type,
                    "group_id": group_id,
                    "included_analytes": ";".join(included_labels),
                    "excluded_information": group.get("excluded_information", ""),
                    "F_comp_low": f_low,
                    "F_comp_central": f_central,
                    "F_comp_high": f_high,
                    "measured_group_mass_discharge_g_d": measured_g_d,
                    "F_comp_adjusted_mass_discharge_g_d": measured_g_d * f_central,
                    "adjusted_is_measurement": False,
                    "note": "F_comp-adjusted load is a scenario value, not a measured target-suite load.",
                }
            )
    return pd.DataFrame(rows).sort_values(["plane_id", "group_id"])


def surface_water_loads(surface_water_nodes: pd.DataFrame) -> pd.DataFrame:
    """Calculate surface-water or stormwater loads for each node and analyte.

    If the flow column is a daily flow, the result is g/d. If the user has
    placed an event volume in the flow field, the result is g/event. The code
    preserves the archive terminology and flags the unit interpretation for the
    user to confirm in the manuscript.
    """
    sw = surface_water_nodes.copy()
    for col in ["flow_low_m3_d", "flow_central_m3_d", "flow_high_m3_d"] + ANALYTE_COLUMNS:
        if col in sw.columns:
            sw[col] = _to_number(sw[col])

    rows = []
    for _, row in sw.iterrows():
        for col in ANALYTE_COLUMNS:
            if col not in sw.columns:
                continue
            concentration = row[col]
            if pd.isna(concentration):
                continue
            for scenario in ("low", "central", "high"):
                flow = row.get(f"flow_{scenario}_m3_d", np.nan)
                load = concentration * flow * UNITS.conversion_factor_to_g_d
                rows.append(
                    {
                        "node_id": row["node_id"],
                        "node_type": row.get("node_type", ""),
                        "scenario": scenario,
                        "analyte": COLUMN_ANALYTE_MAP[col],
                        "concentration_ug_L": concentration,
                        "flow_m3_per_supplied_time_unit": flow,
                        "load_g_per_supplied_time_unit": load,
                        "unit_note": "g/d if flow is m3/d; g/event if flow is an event volume",
                    }
                )
    return pd.DataFrame(rows)


def surface_water_group_loads(
    surface_loads: pd.DataFrame, analyte_groups: pd.DataFrame
) -> pd.DataFrame:
    """Aggregate surface-water/stormwater loads by node, scenario and PFAS group."""
    rows = []
    for _, group in analyte_groups.iterrows():
        included_cols = parse_analyte_group(group["included_analytes"])
        included_labels = [COLUMN_ANALYTE_MAP[col] for col in included_cols]
        f_central = float(group.get("F_comp_central", 1.0))
        subset = surface_loads[surface_loads["analyte"].isin(included_labels)].copy()
        for (node_id, node_type, scenario), sdf in subset.groupby(
            ["node_id", "node_type", "scenario"]
        ):
            measured = float(sdf["load_g_per_supplied_time_unit"].sum())
            rows.append(
                {
                    "node_id": node_id,
                    "node_type": node_type,
                    "scenario": scenario,
                    "group_id": group["group_id"],
                    "included_analytes": ";".join(included_labels),
                    "measured_group_load_g_per_supplied_time_unit": measured,
                    "F_comp_central": f_central,
                    "F_comp_adjusted_group_load_g_per_supplied_time_unit": measured * f_central,
                    "unit_note": "g/d if flow is m3/d; g/event if flow is an event volume",
                }
            )
    return pd.DataFrame(rows).sort_values(["node_id", "scenario", "group_id"])


def well_group_loads(
    wells: pd.DataFrame, receptor_flows: pd.DataFrame, analyte_groups: pd.DataFrame
) -> pd.DataFrame:
    """Estimate extraction/load at well-type receptor links.

    This is used only for receptor comparisons where a receptor is linked to a
    water-use bore rather than to a control plane. It is not a plume mass
    discharge; it is an extracted-water load under the assumed receptor flow.
    """
    rows = []
    well_map = wells.set_index("well_id")
    for _, rec in receptor_flows.iterrows():
        link = str(rec.get("linked_plane_or_node", ""))
        if link not in well_map.index:
            continue
        well = well_map.loc[link]
        q = float(rec.get("Q_central_m3_d", np.nan))
        for _, group in analyte_groups.iterrows():
            included_cols = parse_analyte_group(group["included_analytes"])
            included_labels = [COLUMN_ANALYTE_MAP[col] for col in included_cols]
            c_sum = 0.0
            for col in included_cols:
                c_sum += float(well.get(col, 0.0))
            measured = c_sum * q * UNITS.conversion_factor_to_g_d
            rows.append(
                {
                    "well_id": link,
                    "receptor_id": rec["receptor_id"],
                    "group_id": group["group_id"],
                    "included_analytes": ";".join(included_labels),
                    "assumed_extraction_m3_d": q,
                    "measured_group_concentration_ug_L": c_sum,
                    "measured_extracted_load_g_d": measured,
                    "F_comp_central": float(group.get("F_comp_central", 1.0)),
                    "F_comp_adjusted_extracted_load_g_d": measured
                    * float(group.get("F_comp_central", 1.0)),
                    "note": "extracted-water load, not a control-plane mass discharge",
                }
            )
    return pd.DataFrame(rows)


def receptor_allowable_loads(receptor_flows: pd.DataFrame) -> pd.DataFrame:
    """Convert receptor concentration criteria and flow assumptions to loads.

    The criterion column in the synthetic archive is a placeholder. Users must
    replace it with a project-appropriate criterion before using the output for
    any real decision.
    """
    rows = []
    for _, row in receptor_flows.iterrows():
        criterion_ng_L = float(row["criterion_placeholder_PFOS_plus_PFHxS_ng_L"])
        criterion_ug_L = criterion_ng_L / 1000.0
        for scenario in ("low", "central", "high"):
            q = float(row[f"Q_{scenario}_m3_d"])
            mf = float(row[f"mixing_factor_{scenario}"])
            allowed = criterion_ug_L * q * mf * UNITS.conversion_factor_to_g_d
            rows.append(
                {
                    "receptor_id": row["receptor_id"],
                    "receptor_type": row["receptor_type"],
                    "linked_plane_or_node": row["linked_plane_or_node"],
                    "scenario": scenario,
                    "criterion_ng_L": criterion_ng_L,
                    "criterion_ug_L": criterion_ug_L,
                    "Q_m3_d": q,
                    "mixing_factor": mf,
                    "allowable_load_g_d": allowed,
                    "warning": "criterion is a placeholder in the synthetic archive; replace before real use",
                }
            )
    return pd.DataFrame(rows)


def receptor_comparison(
    group_plane: pd.DataFrame,
    surface_group: pd.DataFrame,
    well_group: pd.DataFrame,
    receptor_allowable: pd.DataFrame,
    selected_group_id: str = "G3",
    use_adjusted: bool = False,
) -> pd.DataFrame:
    """Compare current load with central receptor allowable load."""
    allowed = receptor_allowable[receptor_allowable["scenario"] == "central"].copy()
    rows = []
    for _, rec in allowed.iterrows():
        link = str(rec["linked_plane_or_node"])
        load = np.nan
        source = "not found"
        if link.startswith("CP"):
            sub = group_plane[
                (group_plane["plane_id"] == link) & (group_plane["group_id"] == selected_group_id)
            ]
            if not sub.empty:
                col = (
                    "F_comp_adjusted_mass_discharge_g_d"
                    if use_adjusted
                    else "measured_group_mass_discharge_g_d"
                )
                load = float(sub.iloc[0][col])
                source = "control-plane mass discharge"
        elif link.startswith("SW"):
            sub = surface_group[
                (surface_group["node_id"] == link)
                & (surface_group["scenario"] == "central")
                & (surface_group["group_id"] == selected_group_id)
            ]
            if not sub.empty:
                col = (
                    "F_comp_adjusted_group_load_g_per_supplied_time_unit"
                    if use_adjusted
                    else "measured_group_load_g_per_supplied_time_unit"
                )
                load = float(sub.iloc[0][col])
                source = "surface-water or stormwater load"
        else:
            sub = well_group[
                (well_group["well_id"] == link) & (well_group["group_id"] == selected_group_id)
            ]
            if not sub.empty:
                col = (
                    "F_comp_adjusted_extracted_load_g_d"
                    if use_adjusted
                    else "measured_extracted_load_g_d"
                )
                load = float(sub.iloc[0][col])
                source = "extracted-water load"
        allowed_load = float(rec["allowable_load_g_d"])
        rows.append(
            {
                "receptor_id": rec["receptor_id"],
                "receptor_type": rec["receptor_type"],
                "linked_plane_or_node": link,
                "group_id": selected_group_id,
                "load_source": source,
                "load_type": "F_comp adjusted scenario" if use_adjusted else "measured target suite",
                "current_load_g_d": load,
                "allowable_load_g_d": allowed_load,
                "load_ratio": load / allowed_load if allowed_load > 0 and not pd.isna(load) else np.nan,
                "exceeds_placeholder_allowable": bool(load > allowed_load)
                if not pd.isna(load)
                else False,
                "warning": "allowable load uses placeholder synthetic criterion; not a regulatory decision",
            }
        )
    return pd.DataFrame(rows)


def basis_audit(tables: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Summarise value-basis labels in the archive."""
    rows = []
    for table_name, df in tables.items():
        basis_cols = [c for c in df.columns if c.endswith("_basis") or c.endswith("basis")]
        for col in basis_cols:
            counts = df[col].fillna("missing").value_counts(dropna=False)
            for value, count in counts.items():
                rows.append(
                    {
                        "table": table_name,
                        "basis_column": col,
                        "basis_label": value,
                        "count": int(count),
                    }
                )
    return pd.DataFrame(rows).sort_values(["table", "basis_column", "basis_label"])


def deterministic_workflow(
    tables: Dict[str, pd.DataFrame], selected_group_id: str = "G3"
) -> Dict[str, pd.DataFrame]:
    """Run the deterministic central-case technical note calculations."""
    cell = control_plane_cell_results(tables["control_plane_cells"], tables["assumed_K_zones"])
    plane = plane_analyte_summary(cell)
    groups = group_mass_discharge(plane, tables["analyte_groups"])
    sw = surface_water_loads(tables["surface_water_nodes"])
    sw_groups = surface_water_group_loads(sw, tables["analyte_groups"])
    well_groups = well_group_loads(
        tables["wells"], tables["receptor_flows"], tables["analyte_groups"]
    )
    allowable = receptor_allowable_loads(tables["receptor_flows"])
    receptor_measured = receptor_comparison(
        groups, sw_groups, well_groups, allowable, selected_group_id, use_adjusted=False
    )
    receptor_adjusted = receptor_comparison(
        groups, sw_groups, well_groups, allowable, selected_group_id, use_adjusted=True
    )
    audit = basis_audit(tables)
    return {
        "control_plane_cell_results": cell,
        "control_plane_analyte_summary": plane,
        "control_plane_group_summary": groups,
        "surface_water_analyte_loads": sw,
        "surface_water_group_loads": sw_groups,
        "well_group_loads": well_groups,
        "receptor_allowable_loads": allowable,
        "receptor_comparison_measured": receptor_measured,
        "receptor_comparison_fcomp_adjusted": receptor_adjusted,
        "archive_basis_audit": audit,
    }
