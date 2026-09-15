#!/usr/bin/env python3
"""Exact symbolic checks of Region C zero-mode HJB residuals."""

from __future__ import annotations

import sympy as sp


x, y, alpha, r, eta = sp.symbols("x y alpha r eta", positive=True, real=True)


def clean(expr: sp.Expr) -> sp.Expr:
    return sp.factor(sp.cancel(sp.together(expr.rewrite(sp.exp))))


def q(expr: sp.Expr) -> sp.Expr:
    """Unscaled derivative for gap direction (1,0,0,-1)."""
    return sp.diff(expr, x) - sp.diff(expr, alpha)


def q_second(expr: sp.Expr) -> sp.Expr:
    """Unscaled derivative for gap direction (1,-1,0,1)."""
    return sp.diff(expr, x) - sp.diff(expr, y) - sp.diff(expr, eta)


def objects() -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    a = x + y
    b = r - alpha
    d = lambda t: sp.exp(2 * t) - 1
    da, db = d(a), d(b)

    phi_b = sp.sinh(b) ** 3 / sp.sinh(r) ** 3
    gamma1 = sp.exp(3 * b) / db**3
    gamma2 = gamma1 * (db - da) / db

    i1m = 1 / (4 * da**2) - 1 / (4 * db**2)
    i1p = 1 / (2 * da) + 1 / (4 * da**2) - 1 / (2 * db) - 1 / (4 * db**2)
    i2m = 1 / (12 * da**2) - 1 / (4 * db**2) + da / (6 * db**3)
    i2p = (
        1 / (4 * da)
        + 1 / (12 * da**2)
        - 1 / (2 * db)
        - (1 - da) / (4 * db**2)
        + da / (6 * db**3)
    )

    e2alpha = sp.exp(2 * alpha)
    c1 = sp.Rational(3, 4) * sp.exp(-4 * alpha) * (e2alpha - 1) ** 2
    c2 = sp.Rational(3, 2) * sp.exp(-2 * alpha) * (e2alpha - 1)
    e2r = sp.exp(2 * r)
    f1 = sp.exp(-2 * r) * (e2alpha - e2r) ** 2 / sp.sinh(r) ** 5
    f2 = (
        sp.exp(-2 * r)
        * (e2alpha - e2r)
        * (e2alpha - (e2r + 1) / 2)
        / sp.sinh(r) ** 5
    )

    h1 = gamma1 * phi_b + c1 * f1 * i1m + c2 * f2 * (i1p - i1m) / 2
    h2 = gamma2 * phi_b + c1 * f1 * i2m + c2 * f2 * (i2p - i2m) / 2

    dx, dy = d(x), d(y)
    c0 = sp.exp(2 * x) * dy**2 / da**2
    c_one = 8 * sp.exp(-x - y) * dx * dy
    c_two = 12 * sp.exp(-3 * x - y) * dx**2
    b0 = c1 * sp.exp(-a) * f1 + c2 * sp.sinh(a) * f2
    kp = c0 * b0 + c_one * h1 + c_two * h2

    r_one = 4 * sp.exp(-x - 2 * y) * da * dy
    r_two = 12 * sp.exp(-3 * x - 2 * y) * da * dx
    kr = r_one * h1 + r_two * h2
    local = c0 * sp.sinh(a) ** 3 / sp.sinh(a + alpha) ** 3
    return local, kp, kr


def main() -> None:
    local, kp, kr = objects()
    for name, expr in zip(("local", "Kp", "Kr"), (local, kp, kr)):
        residual = clean(q(q(expr)) - expr)
        print(f"{name} residual: {residual}")
        assert residual == 0
    print("PASS: control 01110 is an exact Region C zero mode.")

    lambda_p = sp.sinh(y - eta) / sp.sinh(y)
    lambda_r = sp.sinh(eta) / sp.sinh(y)
    combined = (lambda_p * local, lambda_p * kp + lambda_r * kr)
    for name, expr in zip(("weighted-local", "weighted-kernel"), combined):
        residual = clean(q_second(q_second(expr)) - expr)
        print(f"{name} residual: {residual}")
        assert residual == 0
    print("PASS: control 01001 is an exact Region C zero mode.")


if __name__ == "__main__":
    main()
