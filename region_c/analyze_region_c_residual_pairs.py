#!/usr/bin/env python3
"""Classify the two face residuals for Region C controls.

The value is kept in the exact form

  sum_n c_n U^(n)(rho) + integral_rho^infinity K(r) U(r) dr.

Moving-endpoint differentiation is therefore exact and does not ask SymPy
to manipulate an unevaluated Integral.  The script checks whether the two
face residuals for a control satisfy the homogeneous 3-2 hyperbolic system.
"""

from __future__ import annotations

from dataclasses import dataclass
import itertools

import sympy as sp

from verify_region_c_zero_modes import alpha, eta, objects, r, x, y


rho = x + y + alpha
variables = (x, eta, y, alpha)


@dataclass
class Rep:
    coeff: list[sp.Expr]
    kernel: sp.Expr


def pad(values: list[sp.Expr], length: int) -> list[sp.Expr]:
    return values + [sp.Integer(0)] * (length - len(values))


def add(left: Rep, right: Rep, scale: sp.Expr = sp.Integer(1)) -> Rep:
    length = max(len(left.coeff), len(right.coeff))
    lc, rc = pad(left.coeff, length), pad(right.coeff, length)
    return Rep(
        [lc[index] + scale * rc[index] for index in range(length)],
        left.kernel + scale * right.kernel,
    )


def scale(rep: Rep, factor: sp.Expr) -> Rep:
    return Rep([factor * value for value in rep.coeff], factor * rep.kernel)


def derivative(rep: Rep, variable: sp.Symbol, rho_speed: sp.Expr | None = None) -> Rep:
    if rho_speed is None:
        rho_speed = sp.diff(rho, variable)
    out = [sp.Integer(0)] * (len(rep.coeff) + 1)
    for index, coefficient in enumerate(rep.coeff):
        out[index] += sp.diff(coefficient, variable)
        out[index + 1] += rho_speed * coefficient
    out[0] -= rho_speed * rep.kernel.subs(r, rho)
    return Rep(out, sp.diff(rep.kernel, variable))


def linear_combination(terms: list[tuple[sp.Expr, Rep]]) -> Rep:
    result = Rep([sp.Integer(0)], sp.Integer(0))
    for coefficient, rep in terms:
        result = add(result, rep, coefficient)
    return result


def substitute(rep: Rep, rules: dict[sp.Symbol, sp.Expr]) -> Rep:
    return Rep(
        [value.subs(rules, simultaneous=True) for value in rep.coeff],
        rep.kernel.subs(rules, simultaneous=True),
    )


def trim(rep: Rep) -> Rep:
    values = list(rep.coeff)
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    return Rep(values, rep.kernel)


def face_derivative(rep: Rep, variable: sp.Symbol) -> Rep:
    # After either face substitution rho=y+alpha.
    out = [sp.Integer(0)] * (len(rep.coeff) + 1)
    speed = sp.Integer(1) if variable in (y, alpha) else sp.Integer(0)
    face_rho = y + alpha
    for index, coefficient in enumerate(rep.coeff):
        out[index] += sp.diff(coefficient, variable)
        out[index + 1] += speed * coefficient
    out[0] -= speed * rep.kernel.subs(r, face_rho)
    return trim(Rep(out, sp.diff(rep.kernel, variable)))


def base_rep() -> Rep:
    local, kp, kr = objects()
    lambda_p = sp.sinh(y - eta) / sp.sinh(y)
    lambda_r = sp.sinh(eta) / sp.sinh(y)
    return Rep([lambda_p * local], lambda_p * kp + lambda_r * kr)


def hessian_reps(base: Rep) -> dict[tuple[int, int], Rep]:
    first = [derivative(base, variable) for variable in variables]
    result: dict[tuple[int, int], Rep] = {}
    for i in range(4):
        for j in range(i, 4):
            result[i, j] = derivative(first[i], variables[j])
    return result


def residual(base: Rep, hessian: dict[tuple[int, int], Rep], control: tuple[int, ...]) -> Rep:
    gaps = [control[index + 1] - control[index] for index in range(4)]
    # q-direction/sqrt(2) in (x,eta,y,alpha) coordinates.
    weights = [
        sp.Integer(gaps[0]),
        sp.Integer(gaps[1]),
        sp.Integer(gaps[1] + gaps[2]),
        sp.Rational(gaps[3] - gaps[1] - 2 * gaps[0], 3),
    ]
    square = Rep([sp.Integer(0)], sp.Integer(0))
    for i in range(4):
        for j in range(i, 4):
            multiplier = weights[i] * weights[j] * (1 if i == j else 2)
            if multiplier:
                square = add(square, hessian[i, j], multiplier)
    return add(base, square, -1)


def system_residual(first_face: Rep, second_face: Rep) -> tuple[Rep, Rep]:
    dy_f = face_derivative(first_face, y)
    da_f = face_derivative(first_face, alpha)
    transport = add(dy_f, da_f, -1)  # d_y-3d_z = d_y-d_alpha
    equation_one = linear_combination(
        [
            (1, transport),
            (-3 * sp.coth(y), first_face),
            (3 * sp.csch(y), second_face),
        ]
    )
    equation_two = linear_combination(
        [
            (1, face_derivative(second_face, y)),
            (2 * sp.csch(y), first_face),
            (-2 * sp.coth(y), second_face),
        ]
    )
    return trim(equation_one), trim(equation_two)


def clean(value: sp.Expr) -> sp.Expr:
    return sp.factor(sp.cancel(sp.together(value.rewrite(sp.exp))))


def zero_rep(rep: Rep) -> bool:
    return all(clean(value) == 0 for value in rep.coeff) and clean(rep.kernel) == 0


def equal_rep(left: Rep, right: Rep) -> bool:
    return zero_rep(add(left, right, -1))


def proportional_constant(left: Rep, right: Rep) -> sp.Expr | None:
    """Return c if left=c*right for a symbol-free constant c."""
    length = max(len(left.coeff), len(right.coeff))
    left_values = pad(left.coeff, length) + [left.kernel]
    right_values = pad(right.coeff, length) + [right.kernel]
    ratio: sp.Expr | None = None
    for left_value, right_value in zip(left_values, right_values):
        lc, rc = clean(left_value), clean(right_value)
        if lc == 0 and rc == 0:
            continue
        if lc == 0 or rc == 0:
            return None
        candidate = clean(lc / rc)
        if candidate.free_symbols:
            return None
        if ratio is None:
            ratio = candidate
        elif clean(candidate - ratio) != 0:
            return None
    if ratio is None:
        return None
    # Every coefficient and the kernel have already been normalized and
    # checked against this same constant ratio in the loop above.
    return ratio


def main() -> None:
    base = base_rep()
    print("building the ten coordinate Hessian representations")
    hessian = hessian_reps(base)
    controls = [(0,) + tail for tail in itertools.product((0, 1), repeat=4)]
    skip = {"00000", "00001", "00011", "00101", "01001", "01110"}
    # Include 00001 as a validation of the machinery.
    selected = [control for control in controls if "".join(map(str, control)) not in skip]
    selected.insert(0, (0, 0, 0, 0, 1))

    faces: dict[str, tuple[Rep, Rep]] = {}
    equations: dict[str, tuple[Rep, Rep]] = {}
    classification: dict[str, tuple[bool, bool]] = {}
    for control in selected:
        name = "".join(map(str, control))
        print(f"\ncontrol {name}")
        full = residual(base, hessian, control)
        first = substitute(full, {x: 0, eta: 0})
        second = substitute(full, {x: 0, eta: y})
        faces[name] = (first, second)
        equation_one, equation_two = system_residual(first, second)
        equations[name] = (equation_one, equation_two)
        print("  simplifying equation 1")
        first_zero = zero_rep(equation_one)
        print("  simplifying equation 2")
        second_zero = zero_rep(equation_two)
        classification[name] = (first_zero, second_zero)
        print(f"  homogeneous 3-2 equations: {first_zero}, {second_zero}")

    expected = {
        "00001": (True, True),
        "00010": (True, False),
        "00100": (False, False),
        "00110": (False, True),
        "00111": (False, True),
        "01000": (False, True),
        "01010": (False, False),
        "01011": (False, False),
        "01100": (False, False),
        "01101": (False, False),
        "01111": (False, True),
    }
    assert classification == expected

    print("\nexact collision classes")
    first_groups = [
        ("00100", "01000", "01111"),
        ("00110", "01010", "01101"),
        ("00111", "01011", "01100"),
    ]
    second_groups = [
        ("00010", "00100"),
        ("01000", "01111"),
        ("01010", "01011", "01100", "01101"),
    ]
    for group in first_groups:
        assert all(equal_rep(faces[group[0]][0], faces[name][0]) for name in group[1:])
        print("  first face:  " + " = ".join(group))
    for group in second_groups:
        assert all(equal_rep(faces[group[0]][1], faces[name][1]) for name in group[1:])
        print("  second face: " + " = ".join(group))

    print("\nopposite first-equation commutators inside equal-F classes")
    cancellations: list[tuple[str, str]] = []
    for group in first_groups:
        for left_index in range(len(group)):
            for right_index in range(left_index + 1, len(group)):
                left, right = group[left_index], group[right_index]
                if zero_rep(add(equations[left][0], equations[right][0])):
                    cancellations.append((left, right))
                    print(f"  F_{left}=F_{right} and E1_{left}=-E1_{right}")
    if not cancellations:
        print("  none")

    print("\nconstant weighted cancellations inside equal-F classes")
    weighted: list[tuple[str, str, sp.Expr]] = []
    for group in first_groups:
        left, right = group[0], group[1]
        ratio = proportional_constant(equations[left][0], equations[right][0])
        if ratio is not None and ratio < 0:
            weighted.append((left, right, ratio))
            print(f"  E1_{left}=({ratio}) E1_{right}")
        # In every collision class the last two controls also have equal G,
        # hence equal E1.  A ratio -2 between the first and either of them
        # proves the equal-weight three-control cancellation without asking
        # SymPy to expand the much larger three-term kernel sum.
        last_equal = equal_rep(equations[group[1]][0], equations[group[2]][0])
        assert last_equal and ratio == -2
        if last_equal and ratio == -2:
            print("  unweighted three-control cancellation: " + " + ".join(group))
    if not weighted:
        print("  none")
    print("PASS: hyperbolic-equation and collision classification is exact.")


if __name__ == "__main__":
    main()
