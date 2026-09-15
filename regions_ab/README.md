# Exact algebra checks for Regions A and B

This optional supplement checks selected symbolic identities and finite
certificates used in the manuscript. The proofs are written in the paper;
reading them does not require Mathematica or Python. These checks do not replace
the analytic arguments concerning integral convergence, endpoint limits,
monotonicity, or hyperbolic comparison principles. They are not a formal
verification of the whole paper or an independent all-controls PDE solver.

## Running

The only external dependency is SymPy. The scripts were tested with Python
3.12.14 and SymPy 1.14.0. After installing SymPy, run an entry point with ordinary
Python, for example:

```sh
python3 ab_a2a6_verify.py
```

An absolute path to a script works from any working directory. All auxiliary
modules are supplied alongside the entry points; no notebook, manuscript file,
external data, or network access is needed. Each entry point raises an exception
and exits nonzero if an exact check fails. Optimized Python (`-O` or `-OO`) is
explicitly rejected because it disables assertions. Successful checks print
`PASS` messages. No numerical sampling is used.

## Entry points and scope

Case numbers refer to the manuscript's numbering of controls.

| Script | Exact identities checked |
| --- | --- |
| `ab_a2a6_verify.py` | A2–A6: directional residual kernels and moving-endpoint terms; the A3 quintic and endpoint-domination identities; A4/A6 and A5/A6 boundary comparisons; hyperbolic trace identities for A2, A3, and A6. |
| `ab_a2a6_highorder.py` | The A2 endpoint-primitive derivative, the A6 differential reduction, its fourth parameter derivative and boundary derivatives, and the finite derivative/Taylor certificates for the four auxiliary one-variable functions. |
| `ab_a7a16_certificates.py` | Differential identities for the common A9/A10 one-variable reduction and transport certificates, the A16 reduction, and the integrated A7/A8 endpoint certificate, including the simplified identity `K'=192 sinh(t)^6`. Auxiliary labels such as `B3` in its output name polynomials, not Region B controls. |
| `ab_a7a16_polynomials.py` | Face-kernel factorizations and endpoint terms for A11, A14, and A16. |
| `ab_b34_exact_checks.py` | The B3/B7 and B4/B8 signed-kernel differences, four endpoint certificates, and the two diagonal differences. |
| `ab_bother_verify.py` | The paired-control and hyperbolic trace identities for B5–B10, their four-expert terms, and the four-expert reversal identities used in the reductions of B2, B11, B12, B13, and B16. |

The two supporting modules, `ab_bother_kernel.py` and
`ab_bother_fourrational.py`, contain exact definitions used by
`ab_bother_verify.py`; they are not additional tests.
`ab_a2a6_highorder.py` also imports the equality-checking helper from
`ab_a2a6_verify.py` without running that file's entry-point tests.

The current A2 endpoint argument uses the full-collision Hessian proved in the
smoothness section; the older primitive differentiation check remains a useful
additional algebra check. The shorter B8 proof directly invokes the comparison
lemma. Neither that invocation nor the analytic smoothness argument is proved
by these scripts. See [CHECKS.md](../CHECKS.md) for manuscript labels.

The scripts retain the algebraic substitutions used in the derivations. Their
symbol names are local to each script, and not a shared computational model of
the value function. In particular, zero polynomial identities are checked
exactly, while the signs and admissible ranges needed to turn those identities
into inequalities are justified in the manuscript.
