#!/usr/bin/env python3
"""Exact algebra checks for the Region C single-integral reduction.

Run with Python and SymPy installed:

    python3 verify_region_c_single_integrals.py

The script verifies the four elementary antiderivatives used in the
Fubini reduction and checks the coefficient of the first integral in r3
against the hyperbolic system.  An explicitly labelled obsolete coefficient
is also rejected as a regression test; it is not the coefficient in the paper.
"""

from __future__ import annotations

import sympy as sp


x, y, s, a = sp.symbols("x y s a", positive=True, real=True)
P, J1, J2 = sp.symbols("P J1 J2", real=True)


def clean(expr: sp.Expr) -> sp.Expr:
    return sp.factor(sp.cancel(sp.together(expr.rewrite(sp.exp))))


def verify_antiderivatives() -> None:
    ds = sp.exp(2 * s) - 1
    da = sp.exp(2 * a) - 1
    gamma1 = sp.exp(3 * s) / ds**3
    gamma2 = gamma1 * (ds - da) / ds

    i1_minus = 1 / (4 * da**2) - 1 / (4 * ds**2)
    i1_plus = 1 / (2 * da) + 1 / (4 * da**2) - 1 / (2 * ds) - 1 / (4 * ds**2)
    i2_minus = 1 / (12 * da**2) - 1 / (4 * ds**2) + da / (6 * ds**3)
    i2_plus = (
        1 / (4 * da)
        + 1 / (12 * da**2)
        - 1 / (2 * ds)
        - (1 - da) / (4 * ds**2)
        + da / (6 * ds**3)
    )

    checks = {
        "I1-minus": sp.diff(i1_minus, s) - gamma1 * sp.exp(-s),
        "I1-plus": sp.diff(i1_plus, s) - gamma1 * sp.exp(s),
        "I2-minus": sp.diff(i2_minus, s) - gamma2 * sp.exp(-s),
        "I2-plus": sp.diff(i2_plus, s) - gamma2 * sp.exp(s),
    }
    for name, residual in checks.items():
        residual = clean(residual)
        print(f"{name} residual: {residual}")
        assert residual == 0


def verify_r3_coefficient() -> None:
    dx = sp.exp(2 * x) - 1
    dy = sp.exp(2 * y) - 1
    da = sp.exp(2 * (x + y)) - 1

    c0 = sp.exp(2 * x) * dy**2 / da**2
    c1 = 8 * sp.exp(-x - y) * dx * dy
    c2 = 12 * sp.exp(-3 * x - y) * dx**2
    p = c0 * P + c1 * J1 + c2 * J2

    gamma1_at_a = sp.exp(3 * (x + y)) / da**3
    j1_a = -gamma1_at_a * P
    j2_a = -2 * sp.exp(2 * (x + y)) / da * (J1 - J2)
    d2 = 12 * sp.exp(-3 * x - 2 * y) * da * dx

    candidates = {
        "paper coefficient": 4 * sp.exp(-x - 2 * y) * da * dy,
        "obsolete coefficient (expected to fail)": (
            4 * sp.exp(-x - 2 * y) * (sp.exp(x + 2 * y) - 1) * dy
        ),
    }
    for name, d1 in candidates.items():
        r3 = d1 * J1 + d2 * J2
        first = clean(
            -sp.diff(p, x)
            + sp.diff(p, y)
            - 2 * sp.coth(y) * p
            + 2 * sp.csch(y) * r3
        )
        r3_y = (
            sp.diff(r3, y)
            + sp.diff(r3, J1) * j1_a
            + sp.diff(r3, J2) * j2_a
        )
        second = clean(r3_y + 2 * sp.csch(y) * p - 2 * sp.coth(y) * r3)
        print(f"{name} first-system residual: {first}")
        print(f"{name} second-system residual: {second}")
        if name == "paper coefficient":
            assert first == second == 0
        else:
            assert first != 0 and second != 0


if __name__ == "__main__":
    verify_antiderivatives()
    verify_r3_coefficient()
    print("PASS: Fubini kernels and the paper's r3 coefficient; obsolete variant rejected.")
