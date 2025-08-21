# scripts/sanity_check.py
from __future__ import annotations
import argparse, sys
from pathlib import Path
from math import log
from typing import Dict, Tuple

# Make src importable both with `python -m scripts.sanity_check` and `python scripts/sanity_check.py`
if __package__ is None and __name__ == "__main__":
    sys.path.append(str(Path(__file__).resolve().parents[1] / "src"))

from triad import sieve_phi, sieve_pi, stream_prefix_values
from tails import triad_tail_bound

# Valorile de referință din draft (opțional, doar pentru coloană ΔF vs REF)
REF_F: Dict[int, float] = {
    10_000:   -0.5255069390,
    20_000:   -0.5261163974,
    50_000:   -0.5265148581,
    100_000:  -0.5266594632,
    200_000:  -0.5267366742,
    1_000_000:-0.5268040299,
}

def compute_values(points, N_ref: int, D_omega: float, progress: bool):
    """
    Returnează:
      vals: dict n -> (Xi, Omega, E, F)
      R:    dict n -> (R_Xi, R_E, R_Omega, R_total)
      kappa_interval: (kappa_min, kappa_max), construit din F(N_ref) ± R_total(N_ref)
    """
    N = max(max(points), N_ref)
    if progress:
        print(f"[info] Construiesc phi, pi până la N={N} ...")
    phi = sieve_phi(N)
    pi  = sieve_pi(N)
    if progress:
        print("[info] Sieves gata.")

    need = set(points) | {N_ref}
    vals: Dict[int, Tuple[float,float,float,float]] = {}
    R: Dict[int, Tuple[float,float,float,float]] = {}

    # Agregăm valorile printr-un singur stream
    next_tick = 0.01
    for (n, xi, om, e_val, f_val) in stream_prefix_values(N, phi, pi):
        if n in need:
            vals[n] = (xi, om, e_val, f_val)
            R[n] = triad_tail_bound(n, D_omega=D_omega)
        if progress and n >= int(N * next_tick):
            print(f"[progress] n={n}/{N} ({100.0*n/N:.0f}%)")
            next_tick += 0.01
        if len(vals) == len(need):
            break

    if N_ref not in vals:
        raise RuntimeError("Nu am capturat F(N_ref) din stream; verifică N_ref <= N.")
    _, _, _, F_ref = vals[N_ref]
    _, _, _, Rtot_ref = R[N_ref]
    kappa_interval = (F_ref - Rtot_ref, F_ref + Rtot_ref)
    return vals, R, kappa_interval

def main():
    ap = argparse.ArgumentParser(
        description="Sanity check + tail certificates pentru Triad Balance la praguri fixe."
    )
    ap.add_argument("--points", type=str,
                    default="10000,20000,50000,100000,200000,1000000",
                    help="Lista de n separate prin virgulă (implicit: 10k,20k,50k,100k,200k,1e6).")
    ap.add_argument("--nref", type=int, default=200_000,
                    help="N_ref pentru certificarea lui kappa din F(N_ref) ± R_total(N_ref) (implicit: 200000).")
    ap.add_argument("--D-omega", type=float, default=10.0,
                    help="Constanta D_omega în limita pentru R_Omega (implicit: 10.0).")
    ap.add_argument("--progress", action="store_true", help="Afișează progres la 1%.")
    args = ap.parse_args()

    points = sorted({int(x) for x in args.points.split(",") if x.strip()})
    vals, R, (kappa_min, kappa_max) = compute_values(points, args.nref, args.__dict__["D_omega"], args.progress)
    kappa_mid = 0.5 * (kappa_min + kappa_max)
    one_over_ln2 = 1.0 / log(2.0)

    print("\n[cert] Interval pentru kappa din N_ref={}:".format(args.nref))
    print("       kappa ∈ [{:.12f}, {:.12f}]  (latime ≈ {:.3e})".format(
        kappa_min, kappa_max, (kappa_max - kappa_min)))

    header = (
        f"{'n':>10}  {'F(n)':>14}  {'R_tot(n)':>11}  {'|F-k_mid|':>11}  "
        f"{'OK?':>4}  {'scaled rem.':>12}  {'ΔF vs REF':>12}"
    )
    print("\n" + header)
    print("-" * len(header))

    for n in points:
        xi, om, e_val, f_val = vals[n]
        rx, re, ro, rtot = R[n]
        # Verificare certificat: |F(n) - kappa| <= R_tot(n)  și  |kappa - kappa_mid| <= R_tot(N_ref)
        # Deci garantat: |F(n) - kappa_mid| <= R_tot(n) + R_tot(N_ref)
        ok = abs(f_val - kappa_mid) <= (rtot + (kappa_max - kappa_min) / 2.0)
        scaled = (n * (f_val - kappa_mid) / log(n)) if n > 1 else float('nan')
        ref = REF_F.get(n, float('nan'))
        delta = f_val - ref if ref == ref else float('nan')
        print(f"{n:10d}  {f_val:14.10f}  {rtot:11.3e}  {abs(f_val-kappa_mid):11.3e}  "
              f"{'OK' if ok else 'FAIL':>4}  {scaled:12.6f}  "
              f"{(delta if ref == ref else float('nan')):12.2e}")

    print("\nNote:")
    print("• OK înseamnă că |F(n) - kappa_mid| ≤ R_tot(n) + R_tot(N_ref) (compatibilitate cu certificatul).")
    print("• scaled remainder = n (F(n) - kappa_mid) / ln n; te aștepți ≈ 1/ln 2 ≈ {:.6f} pe n mare.".format(one_over_ln2))
    print("• ΔF vs REF folosește valorile din draft (dacă există pentru n-ul respectiv).")

if __name__ == "__main__":
    main()
