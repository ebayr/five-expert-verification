#!/usr/bin/env python3
"""Exact verification of the Region C (p,r3) propagation system.

The calculation treats the boundary datum p0(a,w), and the two Volterra
integrals J1(a,w), J2(a,w), as arbitrary.  It uses only the fundamental
theorem identities

    d_a J1 = -Gamma1(a) p0(a,w),
    d_a J2 = -2 exp(2a)/(exp(2a)-1) (J1-J2).

Thus a zero residual verifies the propagation formulas independently of
the particular formula or normalization chosen for p0.
"""

from __future__ import annotations

import sympy as sp


x, y, z = sp.symbols("x y z", positive=True, real=True)
a, w = sp.symbols("a w", positive=True, real=True)
P, J1, J2 = sp.symbols("P J1 J2", real=True)


def coefficients() -> tuple[sp.Expr, ...]:
    dx = sp.exp(2 * x) - 1
    dy = sp.exp(2 * y) - 1
    da = sp.exp(2 * (x + y)) - 1

    c0 = sp.exp(2 * x) * dy**2 / da**2
    c1 = 8 * sp.exp(-x - y) * dx * dy
    c2 = 12 * sp.exp(-3 * x - y) * dx**2

    # The coefficient is exp(2x+2y)-1.  Replacing it by
    # exp(x+2y)-1 leaves nonzero hyperbolic-system residuals.
    d1 = 4 * sp.exp(-x - 2 * y) * da * dy
    d2 = 12 * sp.exp(-3 * x - 2 * y) * da * dx
    return c0, c1, c2, d1, d2


def propagation_formulas() -> tuple[sp.Expr, sp.Expr]:
    c0, c1, c2, d1, d2 = coefficients()
    p = c0 * P + c1 * J1 + c2 * J2
    r3 = d1 * J1 + d2 * J2
    return p, r3


def characteristic_derivative(expr: sp.Expr) -> sp.Expr:
    """Apply -d_x+d_y-2d_z while holding a=x+y,w=z-2x fixed."""
    # P,J1,J2 depend only on the invariants a,w, so this operator
    # differentiates only their displayed x,y coefficients.
    return -sp.diff(expr, x) + sp.diff(expr, y)


def y_derivative(expr: sp.Expr) -> sp.Expr:
    """Apply d_y, including a=x+y dependence of P,J1,J2."""
    gamma1_at_a = sp.exp(3 * (x + y)) / (sp.exp(2 * (x + y)) - 1) ** 3
    substitutions = {
        sp.Derivative(P, y): sp.Symbol("P_a"),
    }

    # Differentiate coefficients first, then add the a-derivatives of the
    # abstract quantities.  r3 has no direct P term, so P_a never appears.
    result = sp.diff(expr, y)
    result += sp.diff(expr, P) * sp.Symbol("P_a")
    result += sp.diff(expr, J1) * (-gamma1_at_a * P)
    result += sp.diff(expr, J2) * (
        -2
        * sp.exp(2 * (x + y))
        / (sp.exp(2 * (x + y)) - 1)
        * (J1 - J2)
    )
    # The direct sp.diff above sees P,J1,J2 as symbols, so only the added
    # chain-rule terms contribute their a-dependence.
    return result.xreplace(substitutions)


def verify() -> None:
    p, r3 = propagation_formulas()

    def exact_simplify(expr: sp.Expr) -> sp.Expr:
        # Rewriting coth/csch in exponentials is essential here; a plain
        # factor/cancel leaves a large expression that only looks nonzero.
        return sp.factor(sp.cancel(sp.together(expr.rewrite(sp.exp))))

    residual_1 = exact_simplify(
            characteristic_derivative(p)
            - 2 * sp.coth(y) * p
            + 2 * sp.csch(y) * r3
    )
    residual_2 = exact_simplify(
            y_derivative(r3)
            + 2 * sp.csch(y) * p
            - 2 * sp.coth(y) * r3
    )

    p_at_x_zero = sp.simplify(p.subs(x, 0))
    compatibility = sp.factor(sp.limit(p - r3, y, 0, dir="+"))

    print("first hyperbolic equation residual:", residual_1)
    print("second hyperbolic equation residual:", residual_2)
    print("initial trace p(0,y,z):", p_at_x_zero)
    print("compatibility limit p(x,0,z)-r3(x,0,z):", compatibility)
    if residual_1 != 0 or residual_2 != 0 or p_at_x_zero != P or compatibility != 0:
        raise SystemExit("hyperbolic-system verification failed")
    print("PASS: propagation equations, initial trace, and compatibility.")


if __name__ == "__main__":
    verify()
