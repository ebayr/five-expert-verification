#!/usr/bin/env python3
"""Extract exact second-face residuals for the two final Region C cases.

This is an exact derivation and consistency check.  It reuses the representation
machinery from ``analyze_region_c_residual_pairs.py`` and prints each trace as

    sum_j c_j(y, alpha) U^(j)(y+alpha)
      + integral_(y+alpha)^infinity K(y, alpha, r) U(r) dr.

On the second face, ``alpha=z/3`` and ``y,z >= 0``.
"""

from __future__ import annotations

import sympy as sp

from analyze_region_c_residual_pairs import (
    base_rep,
    clean,
    hessian_reps,
    residual,
    substitute,
)
from verify_region_c_zero_modes import alpha, eta, r, x, y
from search_region_c_final_monotonicity import W_00010, W_01010
from search_region_c_final_monotonicity import local_h as expected_local_h


def extract(name: str):
    control = tuple(int(bit) for bit in name)
    base = base_rep()
    hessian = hessian_reps(base)
    full = residual(base, hessian, control)
    return substitute(full, {x: 0, eta: y})


def main() -> None:
    # Build the Hessian only once for the actual run.
    base = base_rep()
    hessian = hessian_reps(base)
    for name in ("00010", "01010"):
        control = tuple(int(bit) for bit in name)
        rep = substitute(residual(base, hessian, control), {x: 0, eta: y})
        print(f"\n===== G_{name} =====")
        for order, coefficient in enumerate(rep.coeff):
            normalized = clean(coefficient)
            if normalized != 0:
                print(f"c_{order} =")
                print(normalized)
        print("kernel =")
        normalized_kernel = clean(rep.kernel)
        print(normalized_kernel)

        # Replace the exponentials by algebraic variables and impose the
        # moving-endpoint constraint exp(2r)=exp(2alpha)exp(2y)T, T>=1.
        A, Y, R, T = sp.symbols("A Y R T", positive=True, real=True)
        algebraic = normalized_kernel.subs(
            {
                sp.exp(alpha): sp.sqrt(A),
                sp.exp(y): sp.sqrt(Y),
                sp.exp(r): sp.sqrt(R),
            },
            simultaneous=True,
        )
        # A second replacement catches even exponentials that SymPy kept as
        # exp(k*symbol) rather than powers of exp(symbol).
        algebraic = algebraic.subs(
            {
                sp.exp(2 * alpha): A,
                sp.exp(4 * alpha): A**2,
                sp.exp(6 * alpha): A**3,
                sp.exp(8 * alpha): A**4,
                sp.exp(2 * y): Y,
                sp.exp(4 * y): Y**2,
                sp.exp(2 * r): R,
                sp.exp(3 * r): R ** sp.Rational(3, 2),
                sp.exp(4 * r): R**2,
            },
            simultaneous=True,
        )
        endpoint_form = sp.factor(sp.cancel(sp.together(algebraic.subs(R, A * Y * T))))
        print("kernel with exp(2r)=A*Y*T =")
        print(endpoint_form)

        # Convert U=cosh(.)*int h to a single h-integral.  If T=exp(2t),
        # then K(r)cosh(r)dr is rational in T.  Its integral from AY to T,
        # added to c0*cosh(rho)+c1*sinh(rho), is the h-weight.
        c0 = clean(rep.coeff[0])
        c1 = clean(rep.coeff[1])
        assert all(clean(value) == 0 for value in rep.coeff[2:])
        c0_alg = c0.subs(
            {
                sp.exp(alpha): sp.sqrt(A),
                sp.exp(y): sp.sqrt(Y),
                sp.exp(2 * alpha): A,
                sp.exp(2 * y): Y,
                sp.exp(3 * alpha): A ** sp.Rational(3, 2),
                sp.exp(3 * y): Y ** sp.Rational(3, 2),
            },
            simultaneous=True,
        )
        c1_alg = c1.subs(
            {
                sp.exp(alpha): sp.sqrt(A),
                sp.exp(y): sp.sqrt(Y),
                sp.exp(2 * alpha): A,
                sp.exp(2 * y): Y,
                sp.exp(3 * alpha): A ** sp.Rational(3, 2),
                sp.exp(3 * y): Y ** sp.Rational(3, 2),
            },
            simultaneous=True,
        )
        c0_alg = sp.factor(sp.cancel(sp.together(c0_alg)))
        c1_alg = sp.factor(sp.cancel(sp.together(c1_alg)))
        cosh_rho = (A * Y + 1) / (2 * sp.sqrt(A * Y))
        sinh_rho = (A * Y - 1) / (2 * sp.sqrt(A * Y))
        base_weight = sp.factor(
            sp.cancel(sp.together(c0_alg * cosh_rho + c1_alg * sinh_rho))
        )
        local_h = sp.factor(sp.cancel(sp.together(-c1_alg * cosh_rho)))
        assert sp.factor(sp.cancel(sp.together(local_h - expected_local_h))) == 0

        # Recover the rational differential K(r)cosh(r)dr from the already
        # algebraized kernel.  Under R=exp(2r), dr=dR/(2R).
        k_alg = sp.factor(sp.cancel(sp.together(algebraic)))
        k_cosh_dr = sp.factor(
            sp.cancel(
                sp.together(
                    k_alg * (sp.sqrt(R) + 1 / sp.sqrt(R)) / 2 / (2 * R)
                )
            )
        )
        primitive = sp.integrate(k_cosh_dr, R, risch=True)
        assert not primitive.has(sp.Integral)
        h_weight = sp.factor(
            sp.cancel(
                sp.together(
                    base_weight + primitive.subs(R, T) - primitive.subs(R, A * Y)
                )
            )
        )
        expected_weight = {"00010": W_00010, "01010": W_01010}[name]
        assert sp.factor(sp.cancel(sp.together(h_weight - expected_weight))) == 0
        print("coefficient of h(rho) =")
        print(local_h)
        print("single-integral h-weight =")
        print(h_weight)

        q, s, xx = sp.symbols("q s xx", positive=True, real=True)
        hx = (
            sp.sqrt(2)
            * sp.sqrt(xx)
            * (
                xx**6
                - 3 * xx**5
                - 3 * xx**4
                + 12 * xx**3 * sp.log(xx)
                + 3 * xx**2
                + 3 * xx
                - 1
            )
            / (2 * (xx - 1) ** 5 * (xx + 1) ** 2)
        )
        triangle_rules = {A: s / q, Y: 1 / s, T: 1 / xx}
        local_triangle = sp.factor(
            sp.cancel(
                sp.together(
                    local_h.subs(triangle_rules)
                    * hx.subs(xx, q)
                )
            )
        )
        integrand_triangle = sp.factor(
            sp.cancel(
                sp.together(
                    h_weight.subs(triangle_rules) * hx / (2 * xx)
                )
            )
        )
        print("triangle local term (q=e^-2rho, s=e^-2y) =")
        print(local_triangle)
        int_num, int_den = sp.fraction(integrand_triangle)
        local_num, local_den = sp.fraction(local_triangle)
        print("triangle integral denominator =")
        print(sp.factor(int_den))
        print("triangle local denominator =")
        print(sp.factor(local_den))
        print("degrees in s (integrand numerator, local numerator) =")
        print(sp.degree(int_num, s), sp.degree(local_num, s))
    print("PASS: both residual-to-single-integral weights agree exactly.")


if __name__ == "__main__":
    main()
