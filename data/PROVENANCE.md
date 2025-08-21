# Data Provenance (CSV → Figures / Results)

This document records how CSVs were produced and which figures/results they support.

## Environment
- Python 3.12 (Anaconda); interpreter `/home/user/anaconda3/bin/python3`
- NumPy 1.26.4 (`/home/user/anaconda3/lib/python3.12/site-packages`)
- Exact package versions are pinned in `requirements.txt` and `environment.yml`.

## Reproduction commands
Run from the repository root with `PYTHONPATH=.`:

- Main run (creates `data/run.csv`)
PYTHONPATH=. python -m scripts.triad_full --n 100000 --csv data/run.csv

css
Copy
Edit

- Summary sweep (creates `data/summary.csv`)
PYTHONPATH=. python -m scripts.regen_full_csv --n 200000 --out data/summary.csv

bash
Copy
Edit

## CSV files and their use
- `data/triad_convergence_data.csv` → convergence plots of \(F(n)\) (linear/log–log).
- `data/triad_convergence_data_exactpi.csv` → **convergence_loglog_exactpi.png**.
- `data/triad_pi_sensitivity.csv` → **kappa_alpha_beta_heatmap.png**.
- `data/triad_balance_table_to_1e6.csv` → **kappa_levelsets_alpha_beta_n1e6.png**.
- `data/triad_balance_to_1e6.csv` → auxiliary balance tables/curves (optional).
- `data/label_ref_report.csv` → auxiliary label/reference report (optional).
- `data/summary.csv` → **fig2_scaled_remainder.png** (scaled remainder diagnostic).
- `data/run.csv` → sample run used in quickstart/sanity checks.

> Note: the manuscript embeds only PNG figures; it does not read CSVs directly.
> Keeping the CSVs here ensures transparency and exact reproducibility.

## Integrity checks
Create checksums for the current CSVs:
```bash
for f in data/*.csv; do sha256sum "$f"; done > data/SHA256SUMS.txt
Verify later with:

bash
Copy
Edit
sha256sum -c data/SHA256SUMS.txt
Repository & release
Repo: https://github.com/andronacheion/An-Experimental-Study-of-a-Combinatorial-Arithmetic-Entropy-Triad-and-Its-Limit-Constant

Release used in the manuscript: v1.0.0 (commit 17b77ce546cba40adc786ff6248ba219aa36b633)

SHA256(release-v1.0.0.zip): 399d99176ee4ec8b279a2dc184b406f1880a347f28f4716b4939c0c752c530e0

Zenodo DOI: to be added at the proof stage
