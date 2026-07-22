"""Minimal regression tests for the PFAS technical note code."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd

from pfas_flux.calculations import control_plane_cell_results, deterministic_workflow
from pfas_flux.io import read_archive


def test_unit_conversion_single_cell():
    cp = pd.DataFrame(
        {
            "plane_id": ["CPX"],
            "cell_id": ["C01"],
            "plane_type": ["test"],
            "x_m": [0.0],
            "y_min_m": [0.0],
            "y_max_m": [10.0],
            "area_m2": [100.0],
            "k_zone_id": ["KZ"],
            "gradient_normal_central": [0.01],
            "PFOS_ug_L": [10.0],
            "PFHxS_ug_L": [0.0],
            "PFOA_ug_L": [0.0],
            "PFBA_ug_L": [0.0],
            "PFBS_ug_L": [0.0],
            "FTS_6_2_ug_L": [0.0],
        }
    )
    kz = pd.DataFrame(
        {
            "k_zone_id": ["KZ"],
            "K_central_m_d": [1.0],
            "gradient_central": [0.01],
            "K_low_m_d": [1.0],
            "K_high_m_d": [1.0],
            "gradient_low": [0.01],
            "gradient_high": [0.01],
            "effective_thickness_low_m": [10.0],
            "effective_thickness_central_m": [10.0],
            "effective_thickness_high_m": [10.0],
        }
    )
    out = control_plane_cell_results(cp, kz)
    pfos = out[out["analyte"] == "PFOS"].iloc[0]
    # 10 ug/L * 1 m/d * 0.01 * 100 m2 * 1e-3 = 0.01 g/d
    assert abs(pfos["mass_discharge_g_d"] - 0.01) < 1e-12


def test_archive_runs():
    root = Path(__file__).resolve().parents[1]
    tables = read_archive(root / "data" / "synthetic_case_archive")
    out = deterministic_workflow(tables, selected_group_id="G3")
    assert not out["control_plane_cell_results"].empty
    assert not out["control_plane_group_summary"].empty
    assert "receptor_comparison_measured" in out
