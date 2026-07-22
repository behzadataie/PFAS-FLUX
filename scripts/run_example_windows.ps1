$ErrorActionPreference = "Stop"

python -m pip install --upgrade pip
python -m pip install -r .\requirements.txt

python .\run_pfas_flux.py `
  --input-dir .\data\synthetic_case_archive `
  --output-dir .\outputs `
  --selected-group-id G3 `
  --n-realizations 10000 `
  --seed 20260703

python .\scripts\compute_sensitivity_summary.py `
  --input-dir .\data\synthetic_case_archive `
  --output-file .\outputs\tables\sensitivity_rank_summary.csv `
  --selected-plane-id CP3 `
  --selected-receptor-id R3 `
  --selected-group-id G3 `
  --n-realizations 10000 `
  --seed 20260703

Write-Host "Done. Tables are in .\outputs\tables and figures are in .\outputs\figures"
