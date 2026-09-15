# Symbolic checks and manuscript crosswalk

This guide maps the executable checks to the LaTeX labels in
`ArXiv/five-expert-paper.tex`. Labels, rather than rendered equation numbers,
are used because the numbering changes as the paper is edited. A label can be
located by searching the manuscript source for `\label{...}`.

The scripts contain their own formulas: they do **not** read the manuscript or
automatically detect a change to its equations. A successful check means that
the asserted symbolic identities agree for the formulas embedded in that
script. It does not mean that every statement of the referenced proposition
has been checked computationally.

## Running and interpreting the checks

Install the dependencies in `requirements.txt`, then use `python run_checks.py`
from this directory. The runner and individual scripts also work when invoked
by an absolute path from another directory. The Region C classification is
optional and slower; include it with `--include-slow`. See [README.md](README.md)
for runner options. Do not use `-O`, `-OO`, or `PYTHONOPTIMIZE`, because some
checks rely on Python assertions.

The tables below describe source-level coverage, not a record of a fresh test
run. Execution results and the tested environment are recorded separately in
the package README. No numerical sampling is used as a proof certificate.

## Construction

[`construction/verify_construction.py`](construction/verify_construction.py)
checks the finite algebra behind the diagonal, interface and propagation
calculations. Its detailed inventory is in
[construction/README.md](construction/README.md). The following groups all
belong to this one entry point, not to separate test scripts.

| Calculation | Manuscript labels | Scope of the exact check |
| --- | --- | --- |
| Diagonal profile | `pc:eq:diagonalpair`, `pc:eq:Uode`, `pc:eq:Uintegral` | Companion elimination, the elementary primitive's derivative, and the reduction-of-order integrating factor. |
| Initial curvature and density | `id:eq:initialcurvature`, `id:eq:densityfromU`, `pc:eq:g`, `id:eq:normalcompatibility` | The curvature obtained from the profile integrand, its algebraic conversion to the displayed density, and the normal-compatibility differential identity modulo the profile ODE. |
| Interface equation and homogeneous mode | `id:eq:interfacepair`, `id:eq:eliminated`, `id:eq:factorization`, `id:eq:curvaturetransport`, `id:eq:transportsolution` | Elimination, curvature factorization, the exponential transport equation, and the ODE for the decaying homogeneous mode. |
| Interface baseline and integral kernels | `id:eq:stationarypair`, `id:eq:greeninterface`, `id:eq:kernelinterface`, `id:eq:companioninterface`, `pc:eq:Lidentity` | Baseline equations; the pointwise density-kernel equations; the vanishing moving-endpoint kernel and its second-derivative coefficient; and the algebraic Green-kernel change of variables. |
| Common propagation construction | `pc:eq:Tderivation`, `pc:eq:Em`, `pc:eq:curvaturefactorization`, `pc:eq:curvaturetransport`, `pc:eq:modeODE` | The companion integrating factor, elimination, factorization, remaining-mode ODE, and the characteristic curvature formula for `m=2,3`. |
| Green inversion and the displayed formulas | `pc:eq:compactGreen`, `pc:eq:generalpropagator`, `pc:eq:Vm`, `pc:eq:K`, `pc:eq:generator`, `pc:eq:Pderived`, `pc:eq:P` | The two lower-endpoint identities and `W''-W=V_m`, the zero-time integral kernels, generator coefficients, and equality with the displayed `K` and `P` integral kernels. |
| Elementary solution pairs | `pc:eq:m2system`, `pc:eq:modeODE`, `pc:eq:derivedBsystem` | Both equations for the decaying propagation modes and the Region B particular pair. |

These calculations do not evaluate `U(0)`, derive the stochastic reflection
model, justify Green inversion or the cutoff limit for bounded data, or prove
boundedness and uniqueness. For example, checking a mode's ODE does not prove
that the analytic boundary conditions exclude its free coefficient. Checking
the normal-compatibility differential identity likewise does not check the
decay argument removing its homogeneous modes. This standalone script uses
explicit failure exceptions even under `-O`; the rest of the package must
still be run without optimization as explained above.

## Regions A and B

Case numbers are those in the manuscript's control table
`ab:case-numbering`; complementary binary controls have the same residual.
The scripts use several local exponential coordinates and positive
normalizations. Their variable names are not one shared computational model.

| Entry point | Manuscript labels | What is asserted exactly |
| --- | --- | --- |
| [`regions_ab/ab_a2a6_verify.py`](regions_ab/ab_a2a6_verify.py) | `a26:coefficients`, `a26:residual`, `a26:homogeneous-system`; propositions `a26:case2`, `a26:case3`, `a26:case4`, `a26:case5`, `a26:case6`; `a26:A3-domination`, `a26:U6`, `a26:W6`, `a26:diagonal-comparison`, `a26:N4-diagonal` | Selected directional residual kernels and their moving-endpoint terms; the A3 polynomial and exact endpoint-domination identity; the A3/A4, A4/A6 and A5/A6 comparisons; the first hyperbolic equation for A2/A6 and the companion equation for A2/A3/A6, with integral, endpoint and four-expert contributions separated. |
| [`regions_ab/ab_a2a6_highorder.py`](regions_ab/ab_a2a6_highorder.py) | `a26:BJ`, `a26:B-expanded`, `a26:J-certificates` and the following derivative certificates in proposition `a26:case6` | The A6 differential reduction and the displayed parameter derivatives of its integral-free expression; derivative factorizations and initial derivatives of the four auxiliary functions `p0`–`p3`. It also retains a former A2 endpoint-primitive derivative as a supplemental check, not the current A2 endpoint proof. The inequalities deduced from these identities remain analytic arguments in the paper. |
| [`regions_ab/ab_a7a16_certificates.py`](regions_ab/ab_a7a16_certificates.py) | `a716:eq:F910`, `a716:eq:E-jets`, lemma `a716:lem:F-positive`; `a716:eq:H9`, lemma `a716:lem:G-positive`; `a716:eq:16-face`, proposition `a716:prop:A11A16`; `a716:eq:F78-zero`, `a716:eq:J78-second`, lemma `a716:lem:78-edges` | Integral elimination and differential identities for the common A9/A10 first-face certificate; selected A9 forcing derivatives and initial jets; the A16 derivative chain; and derivative identities used for the A7/A8 endpoint signs. Printed names such as `B3` designate auxiliary functions, not Region B controls. |
| [`regions_ab/ab_a7a16_polynomials.py`](regions_ab/ab_a7a16_polynomials.py) | `ad:eq:residualformula`, `a716:eq:11-polynomials`, propositions `a716:prop:A11A16` and `a716:prop:A14` | The two correction-kernel factorizations and moving-endpoint coefficients for A11, A14 and A16, obtained by directional differentiation of the common kernel. This does not independently check the four-expert residual inequalities used in the propositions. |
| [`regions_ab/ab_bother_verify.py`](regions_ab/ab_bother_verify.py) | `bo:eq:polynomial`, `bo:eq:residualkernel`, `bo:eq:traceidentities`, `bo:eq:firsterrors`, `bo:eq:seconderrors`; `bo:eq:Delta`, `bo:eq:Deltaidentities`; propositions and reductions `bo:prop:pairedcases`, `bo:eq:directcomparisons` | The B5/B6 and B9/B10 first-trace identities and error cancellations, the B8/B9 second-trace identity, the B7/B8 first-equation errors, selected companion equations and endpoint cancellations, and the four-expert reversal identities. The paired boundary principle and the sign arguments applying it are not executed by this script. |
| [`regions_ab/ab_b34_exact_checks.py`](regions_ab/ab_b34_exact_checks.py) | `b34:integraldifferences` and the polynomial endpoint and diagonal formulas in proposition `b34:comparison` | The B3/B7 and B4/B8 correction-kernel differences; two endpoint evaluations for `Q3`; two endpoint evaluations of the integrated `Q4` kernel; and the two diagonal kernel differences. The admissible parameter ranges, integral comparison and boundary propagation are supplied by the proof. |

Supporting modules `regions_ab/ab_bother_kernel.py` and
`regions_ab/ab_bother_fourrational.py` define the correction and four-expert
expressions imported by `ab_bother_verify.py`; they are not separate test
entry points. The high-order A2/A6 script imports its checking helper from
`ab_a2a6_verify.py` without executing that file's main test suite.

Two distinctions matter when reading the A/B output. The current A2 endpoint
identity `a26:A2-origin` follows from the collision jet
`sm:eq:fullcollisionjet`, not the old primitive calculation retained in the
high-order script. Also, the A7/A8 second-derivative check for `J78` starts
from an embedded explicit formula; it does not reconstruct that formula from
the defining endpoint integral. The paper gives those connections.

## Region C

| Entry point | Manuscript labels | What is asserted exactly |
| --- | --- | --- |
| [`region_c/verify_p_boundary.py`](region_c/verify_p_boundary.py) | `pc:eq:P`, `pc:eq:Pode` | The scalar identity for `P_yy-P`, including the differentiated lower endpoints of both tail integrals, for an arbitrary profile `U`. This is not a check of the equation defining `U`, nor by itself a verification of the full 3–2 system. |
| [`region_c/verify_hyperbolic_system.py`](region_c/verify_hyperbolic_system.py) | `pc:eq:ACsystem`, `pc:eq:nestedC`, `pc:eq:Cprefactors`, `pc:eq:Jnested` | Both propagation equations for the nested pair, its initial trace at `x=0`, and equality of the two traces as `y` tends to zero. The boundary datum and integrals are abstract quantities with their stipulated fundamental-theorem differentiation rules. |
| [`region_c/verify_region_c_single_integrals.py`](region_c/verify_region_c_single_integrals.py) | `pc:eq:Iprimitives`, `pc:eq:Cprefactors`; supporting the reduction `pc:eq:fubini` | Derivatives of the four elementary primitives and the consistency of the first companion prefactor with the hyperbolic system. An obsolete prefactor is deliberately checked to have a nonzero residual: that expected rejection is a regression test, not a failure of the paper's formula. Fubini's theorem and convergence are analytic obligations. |
| [`region_c/verify_region_c_zero_modes.py`](region_c/verify_region_c_zero_modes.py) | `pc:eq:kernelzeros`, proposition `rc:prop:zeromodes` | Local-term and kernel identities for the additional equality controls `01110` and `01001`. The former is complementary to `10001`, used in the A/B control table. The common equality control `00101` is not a third test in this script. |
| [`region_c/analyze_region_c_residual_pairs.py`](region_c/analyze_region_c_residual_pairs.py), optional slow entry point | `rc:eq:Fclasses`, `rc:eq:mixedclasses`, `rc:eq:Gclasses`, `rc:eq:Ecancel`; the exact equations stated before proposition `rc:prop:averaging`; proposition `rc:prop:00001` | Hessian-based residual representations for the ten averaging controls and the validation control `00001`; the classification of which two hyperbolic equations vanish; equality of collision traces; and the exact three-control cancellations. This is not an exhaustive run over all 48 regional control inequalities. |
| [`region_c/derive_region_c_remaining_traces.py`](region_c/derive_region_c_remaining_traces.py) | The two traces in `rc:eq:twotargets` and the conversion underlying `rc:eq:Hrepresentation`, using `rc:eq:Uprofile` | Derivation of the `00010` and `01010` second-face residuals from the coordinate Hessian; disappearance of derivatives of `U` of order two and higher; exact agreement of their local `h` coefficient and single-integral `h` weights with the companion formula library. It also prints the triangular-coordinate representation. |
| [`region_c/verify_region_c_final_monotonicity.py`](region_c/verify_region_c_final_monotonicity.py) | `rc:eq:Hrepresentation`, `rc:eq:J0`, `rc:eq:J1`, `rc:eq:K0`, `rc:eq:K1`; derivative identities in the proofs leading to `rc:eq:concavity` and `rc:eq:s1`; `paper:eq:pairderivative` | Quadratic dependence after multiplication by `s`; elimination of integrals in the terminal derivatives; the displayed derivative factorizations for their numerators; vanishing jets for the curvature numerators at `q=1`; and selected limits of the terminal expressions. It does not check all endpoint expansions or carry out the analytic integration of signs back to the original traces. |

The final monotonicity script also checks third-derivative identities along
`s=q`, using auxiliary numerators `N0` and `N1`. These are supplemental
identities retained from an alternative reduction, not the current proof of
`rc:eq:sq`: the manuscript obtains that edge from the Region A inequalities
and the established interface matching. In particular, the script's selected
limits of terminal derivatives must not be confused with a check of every
asymptotic expansion of `J_j`, `L_j` or the value function.

`region_c/search_region_c_final_monotonicity.py` is an exact formula library,
despite its historical filename. It supplies the two `h` weights and
moving-endpoint differentiation routines; it is not an additional numerical
search or test. The residual derivation imports the classification module,
which in turn imports the zero-mode formula definitions. Keep these modules
together when copying the supplement.

## Optional Wolfram Language checks

[`region_c/five-expert-paper-checks.wl`](region_c/five-expert-paper-checks.wl)
starts from the polynomial representation `rc:eq:Hrepresentation` and
implements the local-term/integrand update `paper:eq:pairderivative`. It checks
six terminal derivative identities and six identities for their numerator
derivatives. Four terminal identities correspond to `rc:eq:J0`, `rc:eq:J1`,
`rc:eq:K0` and `rc:eq:K1`; the remaining two concern the supplemental `s=q`
calculation described above. It checks identities, not endpoint estimates or
global inequalities.

The Wolfram file has not been executed in the preparation environment. Its
asserted formulas overlap the Python checks, but that is not a claim that this
particular Wolfram implementation has passed a runtime test. Mathematica is
not required for the Python supplement or for reading the paper.

## What these checks do not establish

The manuscript supplies the analytic arguments connecting the identities to
the results. In particular, this package is not a stand-alone verification of:

- convergence, dominated differentiation, or exchanges of improper integrals;
- the uniform tail and differentiated endpoint estimates, including
  `ad:eq:uniform-tail`, `rc:eq:facetails` and the smoothness estimates;
- the hypotheses and proofs of `ab:lem:two-two`, `bo:lem:paired`, or the
  equality-direction boundary reductions;
- every inequality in `ab:thm:complete` and `rc:thm:regionC`;
- closed-chamber regularity or the full global smoothness theorem
  `sm:thm:globalC2`;
- strictness throughout the exact COMB locus argument
  `paper:cor:comb-locus`, viscosity uniqueness, or identification with the
  discrete game's scaling limit in `paper:thm:verification`.

Some individual derivative factorizations and finite initial jets are checked
exactly. Their signs on the stated domains, the relevant limiting data, and
the comparison or Taylor arguments completing the proofs must still be read
in the manuscript. The supplement is not a formal proof-assistant certificate.
