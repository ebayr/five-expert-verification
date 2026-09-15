#!/usr/bin/env python3
"""Exact identities for the diagonal, interface, and m--2 constructions.

Requires only SymPy. No files, notebooks, numerical samples, or external
services are read. See README.md for the analytic statements not checked.
"""

from time import perf_counter

import sympy as sp


y, p, r, d, w = sp.symbols("y p r d w", real=True)
m = sp.symbols("m", nonzero=True, real=True)
c = 1 / sp.sqrt(2)
VARIABLES = (y, p, r, d, w)
EXP_VARIABLES = sp.symbols("E_y E_p E_r E_d E_w", positive=True)
CHECK_COUNT = 0


def rational_form(expression):
    """Normalize hyperbolic identities by exact exponential substitution.

    Linear exponentials are written as monomials in independent positive
    symbols. Nonlinear exponentials, such as exp(-2*w*r), are left intact.
    Derivatives of unspecified functions remain symbolic algebraic atoms.
    No floating-point evaluation or sign sampling is performed.
    """
    expression = sp.sympify(expression)
    if expression.has(sp.Float):
        raise TypeError("Construction identities must not contain floating-point constants")
    expression = expression.rewrite(sp.exp)
    replacements = {}
    for exponential in expression.atoms(sp.exp):
        exponent = sp.expand(exponential.args[0])
        coefficients = [exponent.coeff(variable) for variable in VARIABLES]
        remainder = sp.expand(
            exponent - sum(a * v for a, v in zip(coefficients, VARIABLES))
        )
        if remainder.has(*VARIABLES) or any(
            coefficient.has(*VARIABLES) for coefficient in coefficients
        ):
            continue
        replacements[exponential] = sp.exp(remainder) * sp.prod(
            base**power for base, power in zip(EXP_VARIABLES, coefficients)
        )
    return sp.cancel(expression.xreplace(replacements))


def check(name, expression):
    """Require an exact zero; checks remain active under Python -O."""
    global CHECK_COUNT
    residual = rational_form(expression)
    if residual != 0:
        raise AssertionError(f"FAIL: {name}\nNonzero residual: {residual}")
    CHECK_COUNT += 1
    print(f"PASS: {name}", flush=True)


def diagonal_checks():
    """pc:eq:diagonalpair, pc:eq:Uode, pc:eq:Uintegral."""
    U = sp.Function("U")(y)
    S, C = sp.sinh(y), sp.cosh(y)
    ode = sp.diff(U, y, 2) + 5 * sp.coth(y) * sp.diff(U, y) - 6 * U
    ode += 3 * c * sp.coth(y)
    Q = C * U + S * sp.diff(U, y) / 6 - c * S / 2
    check(
        "diagonal companion elimination",
        sp.diff(Q, y) - 2 * sp.coth(y) * Q + 2 * U / S - S * ode / 6,
    )
    I = y / 16 - sp.sinh(2 * y) / 64 - sp.sinh(4 * y) / 64
    I += sp.sinh(6 * y) / 192
    check("primitive I'(y) = sinh(y)^4 cosh(y)^2", sp.diff(I, y) - S**4 * C**2)
    gamma = sp.Function("gamma")(y)
    check(
        "diagonal reduction-of-order integrating factor",
        sp.diff(sp.diff(gamma, y) * S**5 * C**2, y)
        + 3 * c * S**4 * C**2
        - S**5 * C * ode.subs(U, C * gamma).doit(),
    )
    H = 3 * c * I / (S**5 * C**2)
    D0 = c * (15 * I / S**6 - 3 * sp.coth(y) + sp.tanh(y) / 2)
    check(
        "initial curvature from the differentiated profile integrand",
        -2 * S * H - C * sp.diff(H, y) + c * sp.tanh(y) / 2 - D0,
    )
    check(
        "normal-compatibility identity modulo the profile ODE",
        (sp.diff(U, y, 3) - sp.diff(U, y)) / 3
        + 2 * sp.coth(y) * (sp.diff(U, y, 2) - U + c * sp.tanh(y) / 2)
        - (sp.diff(ode, y) + sp.coth(y) * ode) / 3,
    )

    # Algebraic substitution r=coth(y)>1, a=arccoth(r)=y. The same
    # a is retained as an indeterminate on both sides of the identity.
    a = sp.symbols("a", real=True)
    k = r**2 - 1
    Sr, Cr = 1 / sp.sqrt(k), r / sp.sqrt(k)
    Ir = a / 16 - 2 * Sr * Cr / 64
    Ir -= 4 * Sr * Cr * (Cr**2 + Sr**2) / 64
    Ir += 2 * (3 * Sr + 4 * Sr**3) * (4 * Cr**3 - 3 * Cr) / 192
    D0r = c * (15 * Ir / Sr**6 - 3 * r + sp.Rational(1, 2) / r)
    density = (8 - 25 * r**2 + 15 * r**4 - 15 * r * k**2 * a)
    density /= 16 * sp.sqrt(2) * r * sp.sqrt(k)
    check("density g(r) = -sinh(y)^3 D0(y), r=coth(y)", -Sr**3 * D0r - density)


def interface_checks():
    """id:eq:interfacepair through id:eq:companioninterface."""
    F = sp.Function("ell")(y, w)
    S, C = sp.sinh(y), sp.cosh(y)
    E = sp.diff(F, y, w) - sp.coth(y) * sp.diff(F, w)
    E += 2 * sp.coth(y) * sp.diff(F, y) - 2 * F + c * sp.coth(y)
    R = C * F + S * sp.diff(F, w) / 2 - c * S / 2
    check(
        "interface companion elimination",
        sp.diff(R, y) - 2 * sp.coth(y) * R + 2 * F / S - S * E / 2,
    )
    curvature = sp.diff(F, y, 2) - F
    check(
        "interface curvature factorization",
        sp.diff(E, y) + sp.coth(y) * E
        - sp.diff(curvature, w) - 2 * sp.coth(y) * curvature - c,
    )
    baseline = c * C * sp.atan(sp.exp(-y))
    baseline_companion = C * baseline - c * S / 2
    check(
        "arctangent baseline first-order equation",
        sp.diff(baseline, y) - sp.tanh(y) * baseline + c / 2,
    )
    check(
        "arctangent baseline curvature",
        sp.diff(baseline, y, 2) - baseline + c * sp.tanh(y) / 2,
    )
    check(
        "stationary baseline companion equation",
        sp.diff(baseline_companion, y) - 2 * sp.coth(y) * baseline_companion
        + 2 * baseline / S,
    )
    amplitude = sp.Function("D0")(y)
    transported = sp.exp(-2 * w * sp.coth(y)) * amplitude
    check("interface transported curvature", sp.diff(transported, w) + 2 * sp.coth(y) * transported)
    mode = sp.Function("C")(w)
    check(
        "interface remaining-mode ODE",
        (E - c * sp.coth(y)).subs(F, mode * sp.exp(-y)).doit() * S
        + sp.diff(mode, w) + 2 * mode,
    )

    # Pointwise kernels; differentiation of the integrals additionally
    # requires the analytic convergence arguments given in the paper.
    A = r * S - C
    first = A * sp.exp(-2 * w * r)
    second = -A**2 * sp.exp(-2 * w * r)
    check(
        "interface first equation for each density kernel",
        sp.diff(first, w) - 2 * second / S + 2 * sp.coth(y) * first,
    )
    check(
        "interface second equation for each density kernel",
        sp.diff(second, y) - 2 * sp.coth(y) * second + 2 * first / S,
    )
    check("interface moving-endpoint kernel vanishes", A.subs(r, sp.coth(y)))
    check(
        "interface second-derivative endpoint coefficient",
        sp.diff(A, y).subs(r, sp.coth(y)) * sp.diff(sp.coth(y), y) + S**-3,
    )
    # With reversed limits, dt contributes the positive factor 1/(r^2-1).
    k = r**2 - 1
    check(
        "interface Green-kernel change of variables",
        (C - r * S) / sp.sqrt(k) * (-k ** sp.Rational(3, 2)) / k - A,
    )


def propagation_checks():
    """pc:eq:Em through pc:eq:modeODE, plus the displayed K and P kernels."""
    F = sp.Function("F")(y, p)
    S = sp.sinh(y)
    E = sp.diff(F, y, 2) - m * sp.diff(F, y, p)
    E -= (m + 1) * sp.coth(y) * sp.diff(F, y)
    E += m * sp.coth(y) * sp.diff(F, p) + m * F
    curvature = sp.diff(F, y, 2) - F
    check(
        "general m--2 curvature factorization",
        sp.diff(E, y) + sp.coth(y) * E
        - sp.diff(curvature, y) + m * sp.diff(curvature, p)
        + m * sp.coth(y) * curvature,
    )
    R = sp.cosh(y) * F - S * (sp.diff(F, y) - m * sp.diff(F, p)) / m
    check(
        "general m--2 companion elimination",
        sp.diff(R, y) - 2 * sp.coth(y) * R + 2 * F / S + S * E / m,
    )
    Q = sp.Function("R")(y, p)
    check(
        "companion integrating-factor identity",
        sp.diff(Q / S**2, y) - (sp.diff(Q, y) - 2 * sp.coth(y) * Q) / S**2,
    )
    mode = sp.Function("C")(p)
    check(
        "general m--2 remaining-mode ODE",
        E.subs(F, mode * sp.exp(-y)).doit() * S
        - m * sp.diff(mode, p) - (m + 1) * mode,
    )
    t = y + d
    ratio = sp.sinh(r - d) / sp.sinh(r)
    check("integration-by-parts ratio derivative", sp.diff(ratio, r) - sp.sinh(d) / sp.sinh(r)**2)
    kernels = {}
    for order in (2, 3):
        endpoint = y + p / order
        datum = sp.Function("h")(endpoint)
        transported = (S / sp.sinh(endpoint)) ** order * datum
        check(
            f"m={order} characteristic curvature formula",
            sp.diff(transported, y) - order * sp.diff(transported, p)
            - order * sp.coth(y) * transported,
        )
        W = sp.sinh(r - t) * ratio**order
        V = (
            2 * order * sp.sinh(d) * sp.sinh(t) * sp.sinh(r - d) ** (order - 1)
            + order * (order - 1) * sp.sinh(d)**2 * sp.sinh(r - t)
            * sp.sinh(r - d) ** (order - 2)
        ) / sp.sinh(r) ** (order + 2)
        kernels[order] = V
        check(f"m={order} Green-kernel endpoint W(t)=0", W.subs(r, t))
        check(
            f"m={order} integration-by-parts endpoint W'(t)",
            sp.diff(W, r).subs(r, t) - (S / sp.sinh(t))**order,
        )
        check(f"m={order} twice-integrated kernel W''-W=V", sp.diff(W, r, 2) - W - V)
        check(f"m={order} zero-time integral kernel", V.subs(d, 0))
        check(
            f"m={order} generator local coefficient",
            (sp.diff((S / sp.sinh(t))**order, d) / order).subs(d, 0) + sp.coth(y),
        )
        check(
            f"m={order} generator integral kernel",
            (sp.diff(V, d) / order).subs(d, 0) - 2 * S / sp.sinh(r)**3,
        )
        modal_F = sp.exp(-y - sp.Rational(order + 1, order) * p)
        modal_R = sp.exp(-2 * y - sp.Rational(order + 1, order) * p)
        check(
            f"m={order} decaying mode first equation",
            sp.diff(modal_F, y) - order * sp.diff(modal_F, p)
            - order * sp.coth(y) * modal_F + order * modal_R / S,
        )
        check(
            f"m={order} decaying mode second equation",
            sp.diff(modal_R, y) - 2 * sp.coth(y) * modal_R + 2 * modal_F / S,
        )

    displayed_K = (4 * S * sp.sinh(d) + 6 * sp.sinh(d)**2 * sp.cosh(t)) / sp.sinh(r)**3
    displayed_K -= 6 * sp.sinh(d)**2 * sp.sinh(t) * sp.cosh(r) / sp.sinh(r)**4
    check("m=2 equals the displayed K integral kernel", kernels[2] - displayed_K)
    displayed_P = (
        sp.Rational(3, 4) * sp.exp(-y - 4 * d) * (sp.exp(2 * d) - 1)**2
        * sp.exp(-2 * r) * (sp.exp(2 * d) - sp.exp(2 * r))**2
        + sp.Rational(3, 2) * S * sp.exp(-2 * d) * (sp.exp(2 * d) - 1)
        * sp.exp(-2 * r) * (sp.exp(2 * d) - sp.exp(2 * r))
        * (sp.exp(2 * d) - (sp.exp(2 * r) + 1) / 2)
    ) / sp.sinh(r)**5
    check("m=3 equals the displayed exponential P integral kernel", kernels[3] - displayed_P)

    # Region B's elementary inhomogeneous pair, pc:eq:derivedBsystem.
    b0 = c * (1 + sp.exp(-2 * y) / 3) / 2
    b1 = 2 * c * sp.exp(-y) / 3
    check("Region B particular pair: first equation", sp.diff(b0, y) - 2 * sp.coth(y) * b0 + 2 * b1 / S + c)
    check("Region B particular pair: second equation", sp.diff(b1, y) - 2 * sp.coth(y) * b1 + 2 * b0 / S)


def main():
    start = perf_counter()
    print(f"Exact construction checks; SymPy {sp.__version__}", flush=True)
    diagonal_checks()
    interface_checks()
    propagation_checks()
    print(f"PASS: {CHECK_COUNT} exact identities in {perf_counter() - start:.3f}s", flush=True)


if __name__ == "__main__":
    main()
