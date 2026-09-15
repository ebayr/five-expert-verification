#!/usr/bin/env python3
"""Exact check of the paper's Region C boundary formula for p(0,y,z).

The two unevaluated integrals are represented by I1(y) and I2(y).  Their
derivatives are replaced by the corresponding lower-endpoint integrands, so
the final identity is checked for an arbitrary twice differentiable U.
"""

import sympy as sp


y, alpha = sp.symbols("y alpha", positive=True, real=True)
r = sp.symbols("r", positive=True, real=True)
U = sp.Function("U")
I1 = sp.Function("I1")
I2 = sp.Function("I2")

t = y + alpha
q = sp.exp(2 * alpha)
phi = sp.sinh(y) ** 3 / sp.sinh(t) ** 3

f1 = (
    sp.exp(-2 * r)
    * (q - sp.exp(2 * r)) ** 2
    * U(r)
    / sp.sinh(r) ** 5
)
f2 = (
    sp.exp(-2 * r)
    * (q - sp.exp(2 * r))
    * (q - (sp.exp(2 * r) + 1) / 2)
    * U(r)
    / sp.sinh(r) ** 5
)
f1_at_t = f1.subs(r, t)
f2_at_t = f2.subs(r, t)

p0 = (
    phi * U(t)
    + sp.Rational(3, 4)
    * sp.exp(-y - 4 * alpha)
    * (q - 1) ** 2
    * I1(y)
    + sp.Rational(3, 2)
    * sp.sinh(y)
    * sp.exp(-2 * alpha)
    * (q - 1)
    * I2(y)
)


def boundary_operator(expr: sp.Expr) -> sp.Expr:
    return sp.diff(sp.exp(-2 * y) * sp.diff(sp.exp(y) * expr, y), y)


residual = boundary_operator(p0) - phi * boundary_operator(U(t))
endpoint_rules = {
    sp.diff(I1(y), y): -f1_at_t,
    sp.diff(I1(y), y, 2): -sp.diff(f1_at_t, y),
    sp.diff(I2(y), y): -f2_at_t,
    sp.diff(I2(y), y, 2): -sp.diff(f2_at_t, y),
}
residual = residual.subs(endpoint_rules, simultaneous=True)
residual = residual.rewrite(sp.exp)
residual = sp.factor(sp.cancel(sp.together(residual)))

print("Boundary-ODE residual:")
print(residual)
if residual != 0:
    raise SystemExit("formula does not satisfy the boundary ODE identically")
print("PASS: the p(0,y,z) formula satisfies the boundary ODE for arbitrary U.")
