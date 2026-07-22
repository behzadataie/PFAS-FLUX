"""Input schema and analyte mapping for the PFAS technical note code."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List


ANALYTE_COLUMN_MAP: Dict[str, str] = {
    "PFOS": "PFOS_ug_L",
    "PFHXS": "PFHxS_ug_L",
    "PFHxS": "PFHxS_ug_L",
    "PFOA": "PFOA_ug_L",
    "PFBA": "PFBA_ug_L",
    "PFBS": "PFBS_ug_L",
    "6:2 FTS": "FTS_6_2_ug_L",
    "6:2FTS": "FTS_6_2_ug_L",
    "FTS_6_2": "FTS_6_2_ug_L",
    "FTS 6:2": "FTS_6_2_ug_L",
}

# Canonical analyte labels used in output tables.
COLUMN_ANALYTE_MAP: Dict[str, str] = {
    "PFOS_ug_L": "PFOS",
    "PFHxS_ug_L": "PFHxS",
    "PFOA_ug_L": "PFOA",
    "PFBA_ug_L": "PFBA",
    "PFBS_ug_L": "PFBS",
    "FTS_6_2_ug_L": "6:2 FTS",
}

ANALYTE_COLUMNS: List[str] = list(COLUMN_ANALYTE_MAP.keys())


REQUIRED_FILES = [
    "wells.csv",
    "surface_water_nodes.csv",
    "source_areas.csv",
    "assumed_K_zones.csv",
    "control_plane_cells.csv",
    "receptor_flows.csv",
    "analyte_groups.csv",
    "archive_manifest.csv",
]


REQUIRED_COLUMNS = {
    "assumed_K_zones.csv": [
        "k_zone_id",
        "K_low_m_d",
        "K_central_m_d",
        "K_high_m_d",
        "gradient_low",
        "gradient_central",
        "gradient_high",
        "effective_thickness_low_m",
        "effective_thickness_central_m",
        "effective_thickness_high_m",
    ],
    "control_plane_cells.csv": [
        "plane_id",
        "cell_id",
        "plane_type",
        "y_min_m",
        "y_max_m",
        "area_m2",
        "k_zone_id",
        "gradient_normal_central",
        "PFOS_ug_L",
        "PFHxS_ug_L",
    ],
    "surface_water_nodes.csv": [
        "node_id",
        "node_type",
        "flow_low_m3_d",
        "flow_central_m3_d",
        "flow_high_m3_d",
        "PFOS_ug_L",
        "PFHxS_ug_L",
    ],
    "receptor_flows.csv": [
        "receptor_id",
        "receptor_type",
        "linked_plane_or_node",
        "Q_low_m3_d",
        "Q_central_m3_d",
        "Q_high_m3_d",
        "mixing_factor_low",
        "mixing_factor_central",
        "mixing_factor_high",
        "criterion_placeholder_PFOS_plus_PFHxS_ng_L",
    ],
    "analyte_groups.csv": [
        "group_id",
        "included_analytes",
        "excluded_information",
        "F_comp_low",
        "F_comp_central",
        "F_comp_high",
    ],
    "wells.csv": [
        "well_id",
        "role",
        "x_m",
        "y_m",
        "k_zone_id",
        "PFOS_ug_L",
        "PFHxS_ug_L",
    ],
}


@dataclass(frozen=True)
class UnitConvention:
    """Mass-discharge unit convention used by this code.

    Inputs are expected to use:
    - concentration: micrograms per litre (ug/L)
    - Darcy flux: metres per day (m/d)
    - area: square metres (m2)

    Then C_ug_L * q_m_d * A_m2 * 1e-3 = grams per day.
    """

    concentration: str = "ug/L"
    darcy_flux: str = "m/d"
    area: str = "m2"
    mass_discharge: str = "g/d"
    conversion_factor_to_g_d: float = 1.0e-3


UNITS = UnitConvention()
