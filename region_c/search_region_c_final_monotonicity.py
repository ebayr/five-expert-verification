#!/usr/bin/env python3
"""Exact formula library for the two final Region C traces.

Use q=exp(-2(y+alpha)), s=exp(-2y), with 0<q<=s<=1. After replacing
U(t)=cosh(t)*int_t^infinity h and reversing order, both targets are a local
h(rho) term plus an integral over 0<x<q. H(q,s)=s*G(q,s) is quadratic in s.

The historical filename is retained for imports, but this curated version
contains only exact formulas and moving-endpoint differentiation routines.
Numerical sampling and the exploratory driver are omitted. Run
verify_region_c_final_monotonicity.py for the assertion-based checks.
"""

from __future__ import annotations

import sympy as sp

A, Y, T = sp.symbols("A Y T", positive=True, real=True)
q, s, x = sp.symbols("q s x", positive=True, real=True)
S = x**6 - 3*x**5 - 3*x**4 + 12*x**3*sp.log(x) + 3*x**2 + 3*x - 1
h_x = sp.sqrt(2) * sp.sqrt(x) * S / (2 * (x - 1)**5 * (x + 1)**2)
local_h = 4*A*(Y - 1)**2*(A*Y + 1)/(9*(A*Y - 1)**3)

P_00010 = (
    4*A**4*T*Y**2 + 2*A**4*Y**2 - 6*A**3*T**2*Y**2
    + 3*A**3*T**2*Y - 4*A**3*T*Y**2 + 2*A**3*T*Y
    - 2*A**3*Y**2 + A**3*Y + A*T**4*Y - 2*A*T**4
    + 2*A*T**3*Y - 4*A*T**3 + 3*A*T**2*Y - 6*A*T**2
    + 2*T**4 + 4*T**3
)
W_00010 = -4*P_00010/(9*A**2*Y*(T - 1)**4)
P_01010 = (
    10*A**4*T*Y**2 + 5*A**4*Y**2 + 3*A**3*T**2*Y**2
    + 3*A**3*T**2*Y + 2*A**3*T*Y**2 + 2*A**3*T*Y
    + A**3*Y**2 + A**3*Y + 9*A**2*T**3*Y**2
    - 36*A**2*T**3*Y + 9*A**2*T**3 + 9*A**2*T**2*Y**2
    - 36*A**2*T**2*Y + 9*A**2*T**2 + 9*A**2*T*Y**2
    - 36*A**2*T*Y + 9*A**2*T + A*T**4*Y + A*T**4
    + 2*A*T**3*Y + 2*A*T**3 + 3*A*T**2*Y + 3*A*T**2
    + 5*T**4 + 10*T**3
)
W_01010 = -P_01010/(9*A**2*Y*(T - 1)**4)


def triangle_parts(weight: sp.Expr) -> tuple[sp.Expr, sp.Expr]:
    rules = {A: s/q, Y: 1/s, T: 1/x}
    local = sp.factor(sp.cancel(sp.together(local_h.subs(rules) * h_x.subs(x, q))))
    integrand = sp.factor(sp.cancel(sp.together(weight.subs(rules) * h_x / (2*x))))
    return local, integrand


def curvature(weight: sp.Expr) -> tuple[sp.Expr, sp.Expr]:
    local, integrand = triangle_parts(weight)
    return (
        sp.factor(sp.cancel(sp.together(sp.diff(s*local, s, 2)))),
        sp.factor(sp.cancel(sp.together(sp.diff(s*integrand, s, 2)))),
    )


def terminal_derivative(
    local: sp.Expr, integrand: sp.Expr, sign: int = -1,
    power: int = 2,
) -> tuple[int, sp.Expr]:
    """For M=sign*q^power*(local+integral), eliminate the integral in q."""
    scaled_local = sp.factor(sp.cancel(sp.together(sign*q**power*local)))
    scaled_integrand = sp.factor(sp.cancel(sp.together(sign*q**power*integrand)))
    numerator, denominator = sp.fraction(scaled_integrand)
    assert not denominator.has(q)
    degree = sp.degree(numerator, q)
    boundary = sp.Integer(0)
    for order in range(degree + 1):
        endpoint = sp.diff(scaled_integrand, q, order).subs(x, q)
        boundary = sp.factor(sp.cancel(sp.together(sp.diff(boundary, q) + endpoint)))
    terminal = sp.factor(sp.cancel(sp.together(
        sp.diff(scaled_local, q, degree + 1) + boundary
    )))
    return degree + 1, terminal


def differentiated_rep(
    local: sp.Expr, integrand: sp.Expr, order: int,
    sign: int = 1, power: int = 2,
) -> tuple[sp.Expr, sp.Expr]:
    """Return local and integral-kernel parts of the requested derivative."""
    scaled_local = sp.factor(sp.cancel(sp.together(sign*q**power*local)))
    scaled_integrand = sp.factor(sp.cancel(sp.together(sign*q**power*integrand)))
    boundary = sp.Integer(0)
    for previous_order in range(order):
        endpoint = sp.diff(scaled_integrand, q, previous_order).subs(x, q)
        boundary = sp.factor(sp.cancel(sp.together(sp.diff(boundary, q) + endpoint)))
    local_part = sp.factor(sp.cancel(sp.together(sp.diff(scaled_local, q, order) + boundary)))
    kernel = sp.factor(sp.cancel(sp.together(sp.diff(scaled_integrand, q, order))))
    return local_part, kernel
