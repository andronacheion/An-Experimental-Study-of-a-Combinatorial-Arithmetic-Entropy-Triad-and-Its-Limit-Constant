# An Experimental Study of a Combinatorial–Arithmetic–Entropy Triad and Its Limit Constant

Reproducible source for a study of a combinatorial–arithmetic–entropy triad.  
We give closed forms for two components, unconditional tail bounds for the third, and high-precision evidence for a limit constant \(\kappa\) with rate \(O((\ln n)/n)\).

## Quickstart
```bash
# Python 3.12 (Anaconda), NumPy 1.26.4 recommended
pip install -r requirements.txt
# or: conda env create -f environment.yml

# Reproduce the main run (creates data/run.csv)
PYTHONPATH=. python -m scripts.triad_full --n 100000 --csv data/run.csv

# Regenerate the summary sweep (creates data/summary.csv)
PYTHONPATH=. python -m scripts.regen_full_csv --n 200000 --out data/summary.csv
Sanity: for N=100000, F(N) ≈ -0.5266594632.

Project layout
bash
Copy
Edit
src/        # library modules (triad, tails)
scripts/    # runnable scripts (sanity checks, CSV generators)
paper/      # LaTeX manuscript (final_56.tex)
figures/    # figures included in the paper
data/       # small outputs (CSV) + PROVENANCE.md
Reproducibility
Interpreter: /home/user/anaconda3/bin/python3

NumPy: 1.26.4 (/home/user/anaconda3/lib/python3.12/site-packages)

Exact versions are pinned in requirements.txt and environment.yml.

Data & Code
GitHub: https://github.com/andronacheion/An-Experimental-Study-of-a-Combinatorial-Arithmetic-Entropy-Triad-and-Its-Limit-Constant

Release: v1.0.0 (commit 17b77ce546cba40adc786ff6248ba219aa36b633)

SHA256(release-v1.0.0.zip): 399d99176ee4ec8b279a2dc184b406f1880a347f28f4716b4939c0c752c530e0

Zenodo DOI: to be added at the proof stage

Figures
See figures/README.md for the list used in the manuscript and license info (CC BY 4.0 for figures/data).

Data provenance
See data/PROVENANCE.md for CSV→figure mapping, generation commands, and checksums.

Citation
If you use this code, please cite the paper and this repository (see CITATION.cff).

License
Code: MIT (see LICENSE)

Figures & data: CC BY 4.0

Contact
Corresponding author: andronache.c.ion@proton.me
