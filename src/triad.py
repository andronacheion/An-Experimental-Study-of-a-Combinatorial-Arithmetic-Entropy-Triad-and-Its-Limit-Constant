# src/triad.py
# Definitions for Xi, Omega, E, F with Kahan/pairwise summation
from __future__ import annotations
from dataclasses import dataclass
from math import log, log2, isfinite
from typing import Iterable, Tuple, List

@dataclass
class Kahan:
    """Kahan summation (compensated) for running sums."""
    s: float = 0.0
    c: float = 0.0
    def add(self, x: float) -> None:
        y = x - self.c
        t = self.s + y
        self.c = (t - self.s) - y
        self.s = t
    def value(self) -> float:
        return self.s

def pairwise_sum(a: List[float]) -> float:
    """Pairwise (divide&conquer) sum: more stable than naive sum."""
    n = len(a)
    if n == 0:
        return 0.0
    if n == 1:
        return float(a[0])
    mid = n // 2
    return pairwise_sum(a[:mid]) + pairwise_sum(a[mid:])

# ---------- Arithmetic helpers (pure-Python sieves) ----------

def sieve_primes_bool(N: int) -> List[bool]:
    """Boolean prime table 0..N (True if prime)."""
    if N < 1:
        return [False] * (N + 1)
    is_prime = [True] * (N + 1)
    is_prime[0] = is_prime[1] = False
    p = 2
    while p * p <= N:
        if is_prime[p]:
            step = p
            start = p * p
            is_prime[start:N+1:step] = [False] * (((N - start) // step) + 1)
        p += 1
    return is_prime

def sieve_pi(N: int) -> List[int]:
    """pi(k) for k=0..N (cumulative prime counts)."""
    is_prime = sieve_primes_bool(N)
    pi = [0] * (N + 1)
    cnt = 0
    for k in range(N + 1):
        if is_prime[k]:
            cnt += 1
        pi[k] = cnt
    return pi

def sieve_phi(N: int) -> List[int]:
    """Euler totients phi(k) for k=0..N in O(N log log N)."""
    phi = list(range(N + 1))
    if N >= 1:
        phi[1] = 1
    for i in range(2, N + 1):
        if phi[i] == i:  # i is prime
            for j in range(i, N + 1, i):
                phi[j] -= phi[j] // i
    return phi

# ---------- Term generators (match the paper’s definitions) ----------

def xi_term(i: int) -> float:
    # 2 / (i^2 * (i-1))
    return 2.0 / (i * i * (i - 1.0))

def e_term(i: int) -> float:
    # log2(i) / i^2
    return log2(i) / (i * i)

def omega_term(i: int, phi_i: int, pi_i: int) -> float:
    # 1 / [ i * (phi(i) + pi(i)) ]
    return 1.0 / (i * (phi_i + pi_i))

# ---------- Per-n (scalar) evaluators with Kahan ----------

def Xi(n: int) -> float:
    """Xi(n) = sum_{i=3}^n 2/(i^2(i-1))."""
    if n < 3:
        return 0.0
    acc = Kahan()
    for i in range(3, n + 1):
        acc.add(xi_term(i))
    return acc.value()

def E(n: int) -> float:
    """E(n) = sum_{i=2}^n log2(i)/i^2."""
    if n < 2:
        return 0.0
    acc = Kahan()
    for i in range(2, n + 1):
        acc.add(e_term(i))
    return acc.value()

def Omega(n: int, phi: List[int], pi: List[int]) -> float:
    """Omega(n) = sum_{i=2}^n 1/[i (phi(i)+pi(i))]."""
    if n < 2:
        return 0.0
    acc = Kahan()
    for i in range(2, n + 1):
        acc.add(omega_term(i, phi[i], pi[i]))
    return acc.value()

def F(n: int, phi: List[int], pi: List[int]) -> float:
    """F(n) = Xi(n) + Omega(n) - E(n)."""
    return Xi(n) + Omega(n, phi, pi) - E(n)

# ---------- Streaming prefix (running) evaluators with Kahan ----------

def stream_prefix_values(N: int, phi: List[int], pi: List[int]):
    """
    Yields tuples (n, Xi(n), Omega(n), E(n), F(n)) for n=1..N using Kahan.
    Efficient single pass; matches the exact definitions and start indices.
    """
    xi_acc = Kahan()
    om_acc = Kahan()
    e_acc  = Kahan()

    xi_val = 0.0
    om_val = 0.0
    e_val  = 0.0

    for n in range(1, N + 1):
        if n >= 3:
            xi_acc.add(xi_term(n))
            xi_val = xi_acc.value()
        if n >= 2:
            e_acc.add(e_term(n))
            e_val = e_acc.value()
            om_acc.add(omega_term(n, phi[n], pi[n]))
            om_val = om_acc.value()
        F_val = xi_val + om_val - e_val
        yield (n, xi_val, om_val, e_val, F_val)

# --- Compatibility aliases (for older scripts) ---
# Some older scripts use these names:
sieve_totients = sieve_phi
sieve_prime_counts = sieve_pi

def Xi_terms(N: int):
    """Return list t where t[i] = 2 / (i^2 * (i-1)) for i>=3, else 0."""
    t = [0.0] * (N + 1)
    for i in range(3, N + 1):
        t[i] = 2.0 / (i * i * (i - 1.0))
    return t

def E_terms(N: int):
    """Return list t where t[i] = log2(i) / i^2 for i>=2, else 0."""
    t = [0.0] * (N + 1)
    for i in range(2, N + 1):
        t[i] = log2(i) / (i * i)
    return t

def Omega_terms(N: int, phi: List[int], pi: List[int]):
    """Return list t where t[i] = 1 / ( i * (phi[i] + pi[i]) ) for i>=2."""
    t = [0.0] * (N + 1)
    for i in range(2, N + 1):
        t[i] = 1.0 / (i * (phi[i] + pi[i]))
    return t

# Optional: control public API (nice to have)
__all__ = [
    "Kahan", "pairwise_sum",
    "sieve_primes_bool", "sieve_pi", "sieve_phi",
    "xi_term", "e_term", "omega_term",
    "Xi", "E", "Omega", "F", "stream_prefix_values",
    # compatibility exports
    "sieve_totients", "sieve_prime_counts",
    "Xi_terms", "E_terms", "Omega_terms",
]
