# GitHub and Zenodo release instructions

1. Create a clean public GitHub repository named `PFAS-FLUX` or similar.
2. Upload the repository contents exactly as provided.
3. Confirm the example run succeeds on a clean machine:

   ```bash
   python -m pip install -r requirements.txt
   pytest -q
   python run_pfas_flux.py --input-dir data/synthetic_case_archive --output-dir outputs --selected-group-id G3 --n-realizations 10000 --seed 20260703
   ```

4. Enable the GitHub-Zenodo integration or manually upload a release ZIP to Zenodo.
5. Create a tagged release, for example `v0.3.0`.
6. Copy the Zenodo DOI into:
   - `CITATION.cff`
   - manuscript Code and data availability section
   - supplementary information
7. Freeze the final repository used for review. Do not overwrite it after submission unless a revised manuscript requires a new tagged release.
