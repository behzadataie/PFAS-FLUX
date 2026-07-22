# GitHub and Zenodo release checklist

Before making the repository public:

1. Replace `repository-code` in `CITATION.cff` with the real GitHub URL.
2. Confirm the selected software license with all co-authors.
3. Run `pytest -q`.
4. Run `python run_pfas_flux.py --input-dir data/synthetic_case_archive --output-dir outputs --selected-group-id G3 --n-realizations 10000 --seed 20260703`.
5. Run `python scripts/compute_sensitivity_summary.py --input-dir data/synthetic_case_archive --output-file outputs/tables/sensitivity_rank_summary.csv --selected-plane-id CP3 --selected-receptor-id R3 --selected-group-id G3`.
6. Confirm that outputs in `outputs/tables` and `outputs/figures` regenerate cleanly.
7. Replace placeholder receptor criteria before using the workflow for any real case.
8. Archive a tagged release on Zenodo and update the manuscript Code and data availability section with the final DOI.
