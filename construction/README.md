# Exact construction identities

Run the standalone checker from any working directory:

```sh
python /path/to/verification-codes/construction/verify_construction.py
```

It requires SymPy (the package pins version 1.14.0), reads no other files,
and uses no private notebooks, local formula modules, numerical sampling,
network services, or machine-specific paths. It prints a named `PASS` for
each exact identity and exits with a nonzero status on failure. Checks remain
active under Python's `-O` option. Floating-point constants in the identities
are rejected.

## Coverage and manuscript crosswalk

The equations below are identified by their LaTeX labels, so the crosswalk
does not depend on equation numbering or pagination. The formulas are
transcribed into the checker rather than parsed from the manuscript.

| Part | Exact algebra checked | Manuscript labels |
| --- | --- | --- |
| Diagonal profile | Elimination of the companion trace; derivative of the primitive `I`; reduction-of-order integrating factor | `pc:eq:diagonalpair`, `pc:eq:Uode`, `pc:eq:Uintegral` |
| Initial curvature and density | The differentiated profile integrand gives `D0`; its substitution at `r=coth(s)` gives the displayed elementary `g(r)` with the correct sign | `id:eq:initialcurvature`, `id:eq:densityfromU`, `pc:eq:g` |
| Initial normal derivative | Differential identity reducing compatibility of `ell_w(s,0)=U'(s)/3` to the profile ODE | `id:eq:normalcompatibility` |
| Interface system | Companion elimination, curvature factorization, exponential transport equation, and ODE for the remaining homogeneous mode | `id:eq:interfacepair`, `id:eq:eliminated`, `id:eq:factorization`, `id:eq:curvaturetransport`, `id:eq:transportsolution` |
| Interface baseline and kernels | Arctangent baseline equations; both original system equations for the density kernels; zero moving-endpoint kernel and its second-derivative coefficient; algebraic Green-kernel change of variables | `id:eq:stationarypair`, `id:eq:greeninterface`, `id:eq:kernelinterface`, `id:eq:companioninterface`, `pc:eq:Lidentity` |
| Common propagation system | Companion integrating factor, elimination, curvature factorization, remaining-mode ODE; characteristic formula for `m=2,3` | `pc:eq:Tderivation`, `pc:eq:Em`, `pc:eq:curvaturefactorization`, `pc:eq:curvaturetransport`, `pc:eq:modeODE` |
| Integration by parts | `W(t)=0`, the explicit value of `W'(t)`, and `W''-W=V_m`, for each of `m=2,3` | `pc:eq:compactGreen`, `pc:eq:generalpropagator`, `pc:eq:Vm` |
| Displayed propagators | The `m=2` kernel equals the displayed `K` kernel; the `m=3` kernel equals the displayed exponential `P` kernel; zero-time and generator coefficients | `pc:eq:K`, `pc:eq:generator`, `pc:eq:Pderived`, `pc:eq:P` |
| Elementary pairs | The decaying mode solves both propagation equations for `m=2,3`; the Region B particular pair solves its two inhomogeneous equations | `pc:eq:m2system`, `pc:eq:modeODE`, `pc:eq:derivedBsystem` |

Hyperbolic expressions are rewritten exactly in exponential variables and
the resulting rational residuals are reduced to zero. Derivatives of
unspecified functions remain symbolic. In the density check, the common
quantity `arccoth(r)` is represented by an indeterminate; `r>1` is the
analytic domain, with the positive square root of `r^2-1`.

## What a successful run does not establish

This is an algebraic reproducibility check, not a proof assistant or an
independent proof of the construction theorem. In particular it does not:

- derive the reflected process, its boundary conditions, or scaling from
  the stochastic game;
- evaluate the improper integral giving `U(0)=45*pi^2/(512*sqrt(2))`;
- justify integration by parts, differentiation under improper integrals,
  the change of variables at singular endpoints, or passage from compactly
  supported data to bounded data;
- prove the required boundedness, decay, boundary limits, positivity of
  the kernels, or zero values of the boundary terms at infinity;
- prove that the bounded-mode selection yields uniqueness, the semigroup
  property, or equality of the integral solutions from differential
  identities alone;
- verify global regularity, all residual inequalities, or identification
  of the candidate with the game value.

The manuscript supplies these analytic arguments. For example, checking
the homogeneous-mode ODE is distinct from proving that boundedness and
the initial condition eliminate its free coefficient. Similarly, checking
the differential identity for the normal derivative does not check the
decay argument used to remove its two homogeneous modes.
