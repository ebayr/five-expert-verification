#!/usr/bin/env python3
"""Exact checks for the final two Region C monotonicity arguments.

The large residual-to-single-integral calculation is checked in
``derive_region_c_remaining_traces.py``.  This file checks the subsequent
low-dimensional algebra: concavity, edge and interface terminal derivatives,
and the normalized derivative identities used to determine their signs.
"""

from __future__ import annotations

import sympy as sp

from search_region_c_final_monotonicity import (
    W_00010,
    W_01010,
    curvature,
    differentiated_rep,
    terminal_derivative,
    triangle_parts,
)


q, x = sp.symbols("q x", positive=True, real=True)


def clean(expr: sp.Expr) -> sp.Expr:
    return sp.factor(sp.cancel(sp.together(expr)))


def assert_zero(expr: sp.Expr) -> None:
    assert clean(expr) == 0


def main() -> None:
    B0 = (
        6*q**5 - 60*q**4*sp.log(q) + 125*q**4
        - 120*q**3*sp.log(q) - 80*q**3 - 60*q**2 + 10*q - 1
    )
    B1 = (
        22*q**6 - 300*q**5*sp.log(q) + 763*q**5
        - 1260*q**4*sp.log(q) + 615*q**4
        - 840*q**3*sp.log(q) - 1180*q**3 - 250*q**2 + 33*q - 3
    )
    E0 = (
        q**8 - 11*q**7 + 104*q**6 - 540*q**5*sp.log(q) + 1021*q**5
        - 1320*q**4*sp.log(q) - 540*q**3*sp.log(q) - 1021*q**3
        - 104*q**2 + 11*q - 1
    )
    E1 = (
        q**8 - 11*q**7 + 184*q**6 - 1500*q**5*sp.log(q) + 3261*q**5
        - 4200*q**4*sp.log(q) - 1500*q**3*sp.log(q) - 3261*q**3
        - 184*q**2 + 11*q - 1
    )
    N0 = (
        q**6 - 24*q**5 + 180*q**4*sp.log(q) - 375*q**4
        + 480*q**3*sp.log(q) + 180*q**2*sp.log(q) + 375*q**2
        + 24*q - 1
    )
    N1 = (
        q**6 - 34*q**5 + 300*q**4*sp.log(q) - 655*q**4
        + 840*q**3*sp.log(q) + 300*q**2*sp.log(q) + 655*q**2
        + 34*q - 1
    )

    # Curvature signs: B0<0 and B1<0 follow by backward monotonicity.
    assert_zero(sp.diff(B0, q, 5) - 720*(q - 1)**2/q**2)
    assert_zero(
        sp.diff(B1, q, 6)
        - 1440*(q - 1)*(11*q**2 - 14*q + 7)/q**3
    )
    for order in range(5):
        assert sp.limit(sp.diff(B0, q, order), q, 1) == 0
    for order in range(6):
        assert sp.limit(sp.diff(B1, q, order), q, 1) == 0

    expected_curvature_terminals = (
        -2*sp.sqrt(2)*q**sp.Rational(3, 2)*B0/(3*(q - 1)**10),
        sp.sqrt(2)*sp.sqrt(q)*B1/(2*(q - 1)**11),
    )

    edge_polynomials = (E0, E1)
    edge_normalizers = (
        60*q**3*(9*q**2 + 22*q + 9),
        300*q**3*(5*q**2 + 14*q + 5),
    )
    edge_derivatives = (
        (q - 1)**6*(27*q**4 + 52*q**3 + 162*q**2 + 52*q + 27)
        / (60*q**4*(9*q**2 + 22*q + 9)**2),
        (q - 1)**6*(5*q**4 + 12*q**3 + 158*q**2 + 12*q + 5)
        / (100*q**4*(5*q**2 + 14*q + 5)**2),
    )
    expected_edge_terminals = (
        -3*sp.sqrt(2)*E0/(4*q**sp.Rational(3, 2)*(q - 1)**11),
        -sp.sqrt(2)*E1/(4*q**sp.Rational(3, 2)*(q - 1)**11),
    )

    interface_polynomials = (N0, N1)
    interface_normalizers = (
        60*q**2*(3*q**2 + 8*q + 3),
        60*q**2*(5*q**2 + 14*q + 5),
    )
    interface_derivatives = (
        (q - 1)**8/(10*q**3*(3*q**2 + 8*q + 3)**2),
        (q - 1)**6*(5*q**2 - 34*q + 5)
        / (30*q**3*(5*q**2 + 14*q + 5)**2),
    )
    expected_interface_terminals = (
        sp.sqrt(2)*N0/(2*sp.sqrt(q)*(q - 1)**9),
        sp.sqrt(2)*N1/(4*sp.sqrt(q)*(q - 1)**9),
    )

    for index, weight in enumerate((W_00010, W_01010)):
        local, integrand = triangle_parts(weight)
        curv_local, curv_integrand = curvature(weight)
        order, terminal = terminal_derivative(curv_local, curv_integrand)
        assert order == index + 2
        assert_zero(terminal - expected_curvature_terminals[index])

        edge_local = clean(local.subs(sp.Symbol("s", positive=True, real=True), 1))
        # ``s`` from the source module is symbol-equal to a same-assumption
        # symbol constructed above; obtain it robustly from the expression.
        source_s = next(symbol for symbol in local.free_symbols if symbol.name == "s")
        assert_zero(sp.diff(source_s * local, source_s, 3))
        assert_zero(sp.diff(source_s * integrand, source_s, 3))
        edge_local = clean(local.subs(source_s, 1))
        edge_integrand = clean(integrand.subs(source_s, 1))
        order, terminal = terminal_derivative(
            edge_local, edge_integrand, sign=1, power=2
        )
        assert order == 5
        assert_zero(terminal - expected_edge_terminals[index])

        interface_local = clean(q*local.subs(source_s, q))
        interface_integrand = clean(q*integrand.subs(source_s, q))
        order, terminal = terminal_derivative(
            interface_local, interface_integrand, sign=1, power=0
        )
        assert order == 3
        assert_zero(terminal - expected_interface_terminals[index])

    for polynomial, normalizer, derivative in zip(
        edge_polynomials, edge_normalizers, edge_derivatives
    ):
        normalized = clean(polynomial/normalizer)
        assert_zero(sp.diff(normalized, q) - derivative)
        assert sp.limit(normalized, q, 1) == 0

    for polynomial, normalizer, derivative in zip(
        interface_polynomials, interface_normalizers, interface_derivatives
    ):
        normalized = clean(polynomial/normalizer)
        assert_zero(sp.diff(normalized, q) - derivative)
        assert sp.limit(normalized, q, 1) == 0

    # Endpoint asymptotics used in the backward/single-crossing arguments.
    for terminal in expected_edge_terminals:
        assert_zero(
            sp.limit((1 - q)**4*terminal, q, 1, dir="-")
            + 6*sp.sqrt(2)/7
        )
    assert_zero(
        sp.limit(sp.sqrt(q)*expected_interface_terminals[0], q, 0, dir="+")
        - sp.sqrt(2)/2
    )
    assert_zero(
        sp.limit(sp.sqrt(q)*expected_interface_terminals[1], q, 0, dir="+")
        - sp.sqrt(2)/4
    )
    assert_zero(
        sp.limit((1 - q)**2*expected_interface_terminals[1], q, 1, dir="-")
        + sp.sqrt(2)/14
    )

    print("PASS: both final Region C monotonicity chains are exact.")


if __name__ == "__main__":
    main()
