#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, time
import numpy as np
from pathlib import Path

# Local imports
from src.triad import sieve_totients, sieve_prime_counts, Xi_terms, E_terms, Omega_terms

def kahan_prefix_sum(arr: np.ndarray) -> np.ndarray:
    s = 0.0
    c = 0.0
    out = np.zeros_like(arr, dtype=np.float64)
    for i in range(len(arr)):
        y = float(arr[i]) - c
        t = s + y
        c = (t - s) - y
        s = t
        out[i] = s
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=1_000_000, help="Upper limit N")
    ap.add_argument("--out", type=str, default="triad_balance_full_to_1e6.csv", help="Output CSV path")
    args = ap.parse_args()

    N = int(args.n)
    t0 = time.time()
    phi = sieve_totients(N)
    pi_vals = sieve_prime_counts(N)
    Tx = Xi_terms(N)
    Te = E_terms(N)
    Tom = Omega_terms(N, phi, pi_vals)

    Xi_pref = kahan_prefix_sum(Tx)
    E_pref  = kahan_prefix_sum(Te)
    Om_pref = kahan_prefix_sum(Tom)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["n","Xi","Omega","E","F"])
        for n in range(1, N+1):
            xi = float(Xi_pref[n])
            om = float(Om_pref[n])
            ee = float(E_pref[n])
            w.writerow([n, xi, om, ee, xi + om - ee])

    print(f"[done] Wrote {out_path.resolve()} in {time.time()-t0:.2f}s")

if __name__ == "__main__":
    main()
