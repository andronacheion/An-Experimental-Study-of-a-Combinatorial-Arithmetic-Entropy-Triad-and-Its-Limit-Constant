#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Triad reproducibility script (see paper Appendix)."""
from __future__ import annotations
import argparse, csv, math
from typing import List, Tuple
LOG2 = math.log(2.0)
EULER_GAMMA = 0.5772156649015329
XI_INF_REF = 3.5 - (math.pi ** 2) / 3.0
E_INF_REF  = 1.3525962171055381
KAPPA_REF  = -0.5267366742063
def sieve_totients(N: int) -> List[int]:
    if N < 0: raise ValueError("N must be nonnegative")
    phi = list(range(N + 1)); phi[0] = 0
    if N >= 1: phi[1] = 1
    for i in range(2, N + 1):
        if phi[i] == i:
            for j in range(i, N + 1, i):
                phi[j] -= phi[j] // i
    return phi
def sieve_prime_counts(N: int) -> List[int]:
    if N < 1: return [0] * (N + 1)
    is_prime = [True] * (N + 1); is_prime[0] = is_prime[1] = False
    limit = int(N ** 0.5)
    for p in range(2, limit + 1):
        if is_prime[p]:
            for m in range(p*p, N + 1, p): is_prime[m] = False
    pi = [0] * (N + 1); cnt = 0
    for i in range(2, N + 1):
        if is_prime[i]: cnt += 1
        pi[i] = cnt
    return pi
def Xi_terms(N: int) -> List[float]:
    t = [0.0] * (N + 1)
    for i in range(3, N + 1): t[i] = 2.0 / (i * i * (i - 1))
    return t
def E_terms(N: int) -> List[float]:
    t = [0.0] * (N + 1)
    for i in range(2, N + 1): t[i] = (math.log(i) / LOG2) / (i * i)
    return t
def Omega_terms(N: int, phi: List[int], pi: List[int]) -> List[float]:
    if len(phi) < N + 1 or len(pi) < N + 1: raise ValueError("phi and pi must be sized N+1")
    t = [0.0] * (N + 1)
    for i in range(2, N + 1):
        denom = phi[i] + pi[i]
        if denom <= 0: raise ArithmeticError(f"Unexpected denominator at i={i}")
        t[i] = 1.0 / (i * denom)
    return t
def kahan_prefix_sum(arr: List[float]) -> List[float]:
    s = 0.0; c = 0.0; out = [0.0] * len(arr)
    for i, x in enumerate(arr):
        y = float(x) - c; t = s + y
        c = (t - s) - y; s = t; out[i] = s
    return out
def triad_prefixes(N: int) -> Tuple[List[float], List[float], List[float]]:
    phi = sieve_totients(N); pi  = sieve_prime_counts(N)
    Tx  = Xi_terms(N); Te  = E_terms(N); To  = Omega_terms(N, phi, pi)
    return kahan_prefix_sum(Tx), kahan_prefix_sum(Te), kahan_prefix_sum(To)
def tail_bounds(n: int, D_omega: float = 10.0):
    if n < 2: return (0.0, 0.0, 0.0)
    R_xi = 2.0 / (n * n); R_e  = (math.log(n) + 1.0) / (n * LOG2)
    R_om = (math.exp(EULER_GAMMA) * math.log(max(math.log(max(n,3)), 1.0000001)) + D_omega) / max(n,1)
    return (R_xi, R_e, R_om)
def main():
    ap = argparse.ArgumentParser(description="Compute Xi, E, Omega, F up to N; optional CSV export.")
    ap.add_argument("--n", type=int, default=200_000)
    ap.add_argument("--csv", type=str, default="")
    ap.add_argument("--domega", type=float, default=10.0)
    args = ap.parse_args()
    N = int(args.n)
    Xi_pref, E_pref, Om_pref = triad_prefixes(N)
    XiN, EN, OmN = Xi_pref[N], E_pref[N], Om_pref[N]
    FN = XiN + OmN - EN
    print(f"N = {N}")
    print(f"Xi(N)     = {XiN:.15f}   (Xi_inf ≈ {XI_INF_REF:.15f}, tail ≤ {tail_bounds(N, args.domega)[0]:.3e})")
    print(f"E(N)      = {EN:.15f}   (E_inf ≈ {E_INF_REF:.15f}, tail ≤ {tail_bounds(N, args.domega)[1]:.3e})")
    print(f"Omega(N)  = {OmN:.15f}   (tail ≤ {tail_bounds(N, args.domega)[2]:.3e})")
    print(f"F(N)      = {FN:.15f}   (kappa_ref = {KAPPA_REF:.13f})")
    if args.csv:
        out = args.csv
        with open(out, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f); w.writerow(["n","Xi","Omega","E","F"])
            for n in range(1, N + 1):
                xi = Xi_pref[n]; om = Om_pref[n]; ee = E_pref[n]
                w.writerow([n, f"{xi:.15f}", f"{om:.15f}", f"{ee:.15f}", f"{(xi + om - ee):.15f}"])
        print(f"[done] wrote CSV to: {out}")
if __name__ == "__main__":
    main()
