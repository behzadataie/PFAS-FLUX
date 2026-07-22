#!/usr/bin/env python3
"""Run the PFAS technical note control-plane flux workflow.

This script is the executable entry point for the technical note code package.
It reads the realistic synthetic CSV archive, calculates deterministic and
Monte Carlo mass-discharge outputs, and writes manuscript-ready CSV tables and
figures.

Example
-------
python run_pfas_flux.py \
    --input-dir data/synthetic_case_archive \
    --output-dir outputs \
    --n-realizations 10000 \
    --selected-group-id G3
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict

import pandas as pd

from pfas_flux.calculations import deterministic_workflow
from pfas_flux.io import ensure_output_dirs, read_archive, write_table
from pfas_flux.plots import plot_all
from pfas_flux.uncertainty import MonteCarloOptions, monte_carlo_workflow


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="PFAS technical note control-plane mass-discharge workflow"
    )
    parser.add_argument(
        "--input-dir",
        default="data/synthetic_case_archive",
        help="Directory containing the technical note CSV archive.",
    )
    parser.add_argument(
        "--output-dir",
        default="outputs",
        help="Directory for output tables, figures and logs.",
    )
    parser.add_argument(
        "--selected-group-id",
        default="G3",
        help="Analyte group used for receptor decision summaries. G3 is PFOS+PFHxS in the synthetic archive.",
    )
    parser.add_argument(
        "--n-realizations",
        type=int,
        default=10000,
        help="Number of Monte Carlo realizations.",
    )
    parser.add_argument("--seed", type=int, default=20260703, help="Random seed.")
    parser.add_argument(
        "--concentration-low-multiplier",
        type=float,
        default=0.5,
        help="Low multiplier for concentration uncertainty distributions.",
    )
    parser.add_argument(
        "--concentration-high-multiplier",
        type=float,
        default=2.0,
        help="High multiplier for concentration uncertainty distributions.",
    )
    parser.add_argument(
        "--measured-only-receptors",
        action="store_true",
        help="Use measured target-suite load instead of F_comp-adjusted load for receptor P_exceed.",
    )
    parser.add_argument(
        "--skip-figures", action="store_true", help="Do not create PNG/SVG figures."
    )
    return parser.parse_args()


def write_json(path: Path, data: Dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)


def main() -> None:
    args = parse_args()
    paths = ensure_output_dirs(args.output_dir)
    tables = read_archive(args.input_dir)

    deterministic = deterministic_workflow(tables, selected_group_id=args.selected_group_id)
    for name, df in deterministic.items():
        write_table(df, paths["tables"], f"{name}.csv")

    options = MonteCarloOptions(
        n=args.n_realizations,
        seed=args.seed,
        concentration_low_multiplier=args.concentration_low_multiplier,
        concentration_high_multiplier=args.concentration_high_multiplier,
        selected_group_id=args.selected_group_id,
        use_fcomp_adjusted_for_receptors=not args.measured_only_receptors,
    )
    uncertainty = monte_carlo_workflow(tables, options)
    for name, df in uncertainty.items():
        write_table(df, paths["tables"], f"{name}.csv")

    figure_paths = []
    if not args.skip_figures:
        figure_paths = plot_all(
            deterministic,
            uncertainty,
            paths["figures"],
            selected_group_id=args.selected_group_id,
        )

    run_metadata = {
        "input_dir": str(Path(args.input_dir).resolve()),
        "output_dir": str(Path(args.output_dir).resolve()),
        "selected_group_id": args.selected_group_id,
        "n_realizations": args.n_realizations,
        "seed": args.seed,
        "concentration_low_multiplier": args.concentration_low_multiplier,
        "concentration_high_multiplier": args.concentration_high_multiplier,
        "receptor_load_type": "measured target suite"
        if args.measured_only_receptors
        else "F_comp adjusted scenario",
        "tables_written": sorted([p.name for p in paths["tables"].glob("*.csv")]),
        "figures_written": sorted([p.name for p in figure_paths]),
        "unit_convention": "C_ug_L * K_m_d * gradient * area_m2 * 1e-3 = Md_g_d",
        "important_warning": "The receptor criterion in the synthetic archive is a placeholder. Replace it before any real decision.",
    }
    write_json(paths["logs"] / "run_metadata.json", run_metadata)

    print("PFAS technical note workflow complete.")
    print(f"Tables:  {paths['tables']}")
    print(f"Figures: {paths['figures']}")
    print(f"Log:     {paths['logs'] / 'run_metadata.json'}")


if __name__ == "__main__":
    main()
