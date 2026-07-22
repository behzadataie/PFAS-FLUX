# PFAS-FLUX

PFAS-FLUX is a reproducible control-plane workflow for converting sparse PFAS groundwater monitoring data into compound-specific mass discharge, analyte-group mass discharge, receptor-load exceedance probability, and hydrogeologic data-priority metrics.

This repository accompanies the technical note:

> **PFAS-FLUX: An open control-plane workflow for PFAS mass discharge and receptor-load uncertainty in groundwater**

The workflow does **not** introduce a new definition of mass flux. It implements established groundwater mass-discharge calculations and adds a PFAS-specific reporting layer for:

- compound-explicit and group mass discharge;
- target-suite scope and scenario completeness factors;
- receptor allowable-load comparisons;
- Monte Carlo uncertainty propagation; and
- screening-level data-priority metrics.

The included benchmark is a realistic synthetic example. It contains no real coordinates, site names, bore identifiers, property identifiers, or proprietary field data.

## Repository layout

```text
PFAS-FLUX-HydrogeologyJournal/
├── pfas_flux/                  # Python package
│   ├── calculations.py          # deterministic mass-flux/mass-discharge calculations
│   ├── io.py                    # CSV reading, output writing, and archive validation
│   ├── plots.py                 # manuscript and supplementary plotting functions
│   ├── schema.py                # expected columns, analyte names, and unit conventions
│   └── uncertainty.py           # Monte Carlo sampling and uncertainty summaries
├── data/synthetic_case_archive/
│   ├── wells.csv
│   ├── surface_water_nodes.csv
│   ├── source_areas.csv
│   ├── assumed_K_zones.csv
│   ├── control_plane_cells.csv
│   ├── receptor_flows.csv
│   ├── analyte_groups.csv
│   └── archive_manifest.csv
├── scripts/
│   ├── run_example_bash.sh
│   ├── run_example_windows.ps1
│   └── compute_sensitivity_summary.py
├── outputs/                    # example reproducible outputs
├── figures/                    # manuscript figures as SVG and PNG
├── docs/                       # file dictionary, Monte Carlo design, release notes
├── manuscript/                 # manuscript and supplementary DOCX files
├── tests/                      # unit tests
├── run_pfas_flux.py   # command-line runner retained for compatibility
└── pyproject.toml
```

## Quick start

Create a virtual environment if desired, then install dependencies:

```bash
python -m pip install -r requirements.txt
```

Run the benchmark:

```bash
python run_pfas_flux.py \
  --input-dir data/synthetic_case_archive \
  --output-dir outputs \
  --selected-group-id G3 \
  --n-realizations 10000 \
  --seed 20260703
```

On Windows PowerShell, use one line:

```powershell
python .\run_pfas_flux.py --input-dir .\data\synthetic_case_archive --output-dir .\outputs --selected-group-id G3 --n-realizations 10000 --seed 20260703
```

Run the compact sensitivity/data-priority screen:

```bash
python scripts/compute_sensitivity_summary.py \
  --input-dir data/synthetic_case_archive \
  --output-file outputs/tables/sensitivity_rank_summary.csv \
  --selected-plane-id CP3 \
  --selected-receptor-id R3 \
  --selected-group-id G3 \
  --n-realizations 10000 \
  --seed 20260703
```

Run tests:

```bash
pytest -q
```

## Outputs

The workflow writes CSV tables and figures under `outputs/`:

- `control_plane_cell_results.csv`
- `control_plane_analyte_summary.csv`
- `control_plane_group_summary.csv`
- `receptor_allowable_loads.csv`
- `receptor_comparison_measured.csv`
- `receptor_comparison_fcomp_adjusted.csv`
- `uncertainty_control_plane_group_summary.csv`
- `uncertainty_receptor_summary.csv`
- `uncertainty_driver_proxy.csv`
- `sensitivity_rank_summary.csv`

The placeholder receptor criteria in the synthetic archive are for method demonstration only. Replace them before any real site use.

## Monte Carlo design

The default uncertainty implementation samples:

- hydraulic conductivity using a log-triangular distribution from low, central, and high values;
- hydraulic gradient, effective saturated thickness, concentration, receptor flow, and mixing factor using triangular distributions;
- `F_comp` as a triangular scenario factor; and
- concentrations as triangular values around the central concentration using default 0.5x and 2.0x multipliers.

Input variables are sampled independently in the benchmark. Site applications may require correlated hydraulic parameters, calibrated concentration interpolation, or formal variance decomposition.

## Citation

Use the citation in `CITATION.cff`. If the repository is archived on Zenodo, replace the placeholder DOI in the manuscript and citation file with the Zenodo DOI.

## License

MIT License. See `LICENSE`.
