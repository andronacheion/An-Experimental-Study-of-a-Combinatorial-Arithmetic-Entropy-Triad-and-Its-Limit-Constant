# scripts/make_dataset.py
# Produce CSV (optionally .gz) with rows: n, Xi, Omega, E, F
from __future__ import annotations
import argparse, csv, gzip, sys
from pathlib import Path
from time import perf_counter
from src.triad import sieve_phi, sieve_pi, stream_prefix_values

def main():
    ap = argparse.ArgumentParser(description="Generate triad balance dataset up to N.")
    ap.add_argument("--max-n", type=int, default=1_000_000, help="Upper limit N (default: 1e6).")
    ap.add_argument("--stride", type=int, default=1, help="Write every 'stride'-th n (default: 1).")
    ap.add_argument("--out", type=Path, default=Path("data/triad_balance_to_1e6.csv"),
                    help="Output CSV path ('.gz' -> gzip).")
    ap.add_argument("--progress", action="store_true", help="Print progress every ~1%%.")
    args = ap.parse_args()

    N = int(args.max_n)
    stride = max(int(args.stride), 1)
    out_path: Path = args.out
    out_path.parent.mkdir(parents=True, exist_ok=True)

    t0 = perf_counter()
    if args.progress:
        print(f"[info] Building phi, pi up to N={N} ...")
    phi = sieve_phi(N)
    pi  = sieve_pi(N)
    t1 = perf_counter()
    if args.progress:
        print(f"[info] Sieves done in {t1 - t0:.2f}s")

    # Choose writer (plain CSV or gz)
    if out_path.suffix.endswith(".gz"):
        f = gzip.open(out_path, "wt", newline="")
    else:
        f = open(out_path, "w", newline="")

    with f:
        w = csv.writer(f)
        w.writerow(["n", "Xi", "Omega", "E", "F"])
        next_tick = 0.01
        for (n, xi, om, e_val, f_val) in stream_prefix_values(N, phi, pi):
            if (n % stride) == 0:
                w.writerow([n, f"{xi:.15f}", f"{om:.15f}", f"{e_val:.15f}", f"{f_val:.15f}"])
            if args.progress and n >= int(N * next_tick):
                print(f"[progress] n={n}/{N} ({100.0*n/N:.0f}%)")
                next_tick += 0.01

    t2 = perf_counter()
    if args.progress:
        print(f"[done] Wrote {out_path} in {t2 - t1:.2f}s (total {t2 - t0:.2f}s)")

if __name__ == "__main__":
    # Allow `python -m scripts.make_dataset --progress` from repo root
    # or `python scripts/make_dataset.py --progress`
    if __package__ is None and __name__ == "__main__":
        # Make triad importable when running as a script
        sys.path.append(str(Path(__file__).resolve().parents[1] / "src"))
        from src.triad import sieve_phi, sieve_pi, stream_prefix_values
    main()
