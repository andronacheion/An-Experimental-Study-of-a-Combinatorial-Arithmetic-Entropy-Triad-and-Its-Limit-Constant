# src/tails.py
from __future__ import annotations
from math import log, e
from typing import Tuple

# Euler–Mascheroni constant (high-precision not needed for bounds)
EULER_GAMMA = 0.5772156649015329

def xi_tail_bound(n: int) -> float:
    """
    R_Xi(n) := sum_{i>n} 2/(i^2(i-1)) <= 2/n^2  (Lemma in paper).
    Valid for n >= 2.
    """
    n = max(int(n), 2)
    return 2.0 / (n * n)

def e_tail_bound(n: int) -> float:
    """
    R_E(n) := sum_{i>n} log2(i)/i^2 <= (ln n + 1)/(n ln 2).
    Valid for n >= 2.
    """
    n = max(int(n), 2)
    from math import log
    return (log(n) + 1.0) / (n * log(2.0))

def omega_tail_bound(n: int, D_omega: float = 10.0) -> float:
    """
    R_Omega(n) := sum_{i>n} 1/[i (phi(i)+pi(i))] <= (e^gamma ln ln n + D_omega)/n.
    Valid for n >= 3 (per Rosser–Schoenfeld style bound used in the paper).
    """
    n = max(int(n), 3)
    from math import log
    return ((pow(2.718281828459045, EULER_GAMMA) * log(log(n))) + D_omega) / n

def triad_tail_bound(n: int, D_omega: float = 10.0) -> Tuple[float, float, float, float]:
    """
    Returns (R_Xi, R_E, R_Omega, R_total).
    """
    rx = xi_tail_bound(n)
    re = e_tail_bound(n)
    ro = omega_tail_bound(n, D_omega=D_omega)
    return rx, re, ro, (rx + re + ro)
