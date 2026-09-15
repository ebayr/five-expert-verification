# Optional exact checks for Region C

The paper contains the mathematical arguments. These scripts provide
optional, reproducible checks of finite symbolic identities; a reader does
not need Python or Mathematica to read the proof. They are not a formal
verification of the whole paper.

## Running the Python checks

Use Python 3.10 or later with SymPy installed. From this directory, run:

```sh
python3 verify_p_boundary.py
python3 verify_hyperbolic_system.py
python3 verify_region_c_single_integrals.py
python3 verify_region_c_zero_modes.py
python3 analyze_region_c_residual_pairs.py
python3 derive_region_c_remaining_traces.py
python3 verify_region_c_final_monotonicity.py
```

Each entry point prints a final `PASS` only after its checks finish. A failed
check raises an exception or exits nonzero. Do not use Python's `-O` option
or set `PYTHONOPTIMIZE`: these disable assertions. The two full residual
derivations are substantially slower than the small algebra checks.

No external data, private notebook, TeX parsing, network access, or repository
path is required. All formulas are embedded. SymPy is the only direct Python
package dependency; its ordinary dependency installation includes mpmath.

## Scope and relation to the manuscript

| Entry point | Exact identities checked | Manuscript topic |
| --- | --- | --- |
| `verify_p_boundary.py` | Boundary ODE for the displayed `p(0,y,z)` formula, with arbitrary `U` and endpoint differentiation | Construction of the value function |
| `verify_hyperbolic_system.py` | Two propagation equations, initial trace, and compatibility of `p` and `r3` | Region C propagation |
| `verify_region_c_single_integrals.py` | Four elementary Fubini antiderivatives and the `r3` coefficient | Single-integral representation |
| `verify_region_c_zero_modes.py` | Local and kernel identities for equality controls `01110` and `01001` | Equality modes and reduction to two faces |
| `analyze_region_c_residual_pairs.py` | Coordinate Hessian residuals, face collision identities, which equations vanish, and the exact three-control cancellations | Reduction of the remaining controls |
| `derive_region_c_remaining_traces.py` | Derivation of the two second-face residuals and their local and single-integral weights | Final two traces and quadratic reduction |
| `verify_region_c_final_monotonicity.py` | Quadratic dependence, moving-endpoint terminal derivatives, normalized sign identities, and selected endpoint limits | Concavity, the edge `s=1`, and monotonicity; additional `s=q` regression identities |

The last two scripts import `search_region_c_final_monotonicity.py`. Its
historical name is retained, but this copy is an exact formula library:
numerical sampling and the exploratory driver have been removed. The
residual derivation also imports `analyze_region_c_residual_pairs.py`, which
in turn imports `verify_region_c_zero_modes.py`. Keep these files together.

The single-integral check deliberately also tests an **obsolete coefficient**
and asserts that its residual is nonzero. That expected rejection is a
regression test, not a failure of the coefficient printed in the paper.

The checks do not replace the paper's analytic justifications of convergence,
differentiation under integrals, endpoint regularity, estimates at infinity,
minimum principles, or the reasoning that propagates derivative signs.
In particular, an exact identity is not by itself a proof of an inequality.

The current paper inherits the `s=q` edge signs from Region A. The Python and
Wolfram files also retain an older direct differentiation check along that edge;
it is supplementary and is not needed by the revised proof. See
[CHECKS.md](../CHECKS.md) for the equation-by-equation mapping.

## Optional Wolfram Language alternative

`five-expert-paper-checks.wl` independently expresses the final two-trace
moving-endpoint checks in Wolfram Language. In an installation with
`wolframscript`, run:

```sh
wolframscript -file five-expert-paper-checks.wl
```

It checks six terminal derivative formulas and six derivative identities
used for their signs, without `Integrate` or numerical integration. A failed
identity calls `Quit[1]`. This is a narrower alternative to the final Python
check, not a replacement for every Python entry point. The Wolfram script
has not been executed in the packaging environment; a Wolfram runtime was
not available there.
