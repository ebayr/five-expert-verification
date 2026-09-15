# Optional symbolic verification supplement

Companion to *Prediction with Five Experts and Geometric Stopping:
A Probabilistic Construction and Analytic Verification*, by Erhan Bayraktar, Ibrahim Ekren, and
Nikolaos Kolliopoulos.

Repository: <https://github.com/ebayr/five-expert-verification>.
See [CHECKS.md](CHECKS.md) for the manuscript-label mapping and
[SNAPSHOT.md](SNAPSHOT.md) for the manuscript source fingerprint.

## Purpose and limits

These files reproduce selected exact algebraic identities, residual kernels,
polynomial factorizations, and high-order differentiation certificates appearing
in the manuscript. The formulas are embedded in the scripts; no external data
or Mathematica notebooks are needed.

The paper supplies the mathematical arguments. In particular, successful code
checks do **not** by themselves establish the hypotheses of the hyperbolic
minimum principles, signs on whole domains, analytic endpoint and uniform
limit estimates, differentiation under improper integrals, global smoothness,
or identification with the game value. This is not a formal proof-assistant
certificate, and finite numerical sampling is not used as proof.

## Python checks

The supplied dependencies are SymPy 1.14.0 and mpmath 1.3.0. The package was
tested using Python 3.12.14. From this directory:

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python run_checks.py
```

The runner locates all inputs relative to its own file, so it can also be called
from a different working directory. It returns a nonzero status if a check
process fails or times out. Use `--verbose` to see each script's full output,
`--group construction`, `--group ab`, or `--group c` to select one part, and `--jobs 2` to run two check
processes at a time. The default timeout is 300 seconds per script; adjust it
with `--timeout` on a slower computer. Run without Python's `-O` option, because
individual checks use assertions. `--list` prints the selected entry points
without running them.

The default suite runs thirteen check scripts. The full Region C residual
classification is a fourteenth, slower entry point; include it with:

```sh
python run_checks.py --include-slow --jobs 2 --timeout 1800
```

The GitHub workflow runs this full command with Python 3.12 and the pinned
dependencies on each push or pull request. Its log reports a separate result
for every script; a green run records exact symbolic checks, not a certification
of the paper's analytic arguments. The workflow uses read-only repository
permissions and commit-pinned actions. Local test results are recorded in
[SNAPSHOT.md](SNAPSHOT.md).

| Directory | Manuscript topics |
| --- | --- |
| `construction/` | Propagator and interface differential identities, Green-kernel integration by parts, and interface density |
| `regions_ab/` | Case-specific residual identities and differentiation certificates for the A/B verification section |
| `region_c/` | Single-integral conversion, boundary datum, linear hyperbolic identities, equality modes, residual traces, and the final two monotonicity certificates |

Each directory contains a more detailed file guide. Supporting formula modules
are included because the check scripts import them. The residual classification
module also has the optional slow entry point described above.
Operational tests for the runner, distinct from mathematical checks, run with
`python -m unittest discover -s tests -v`.

## Mathematica / Wolfram Language alternative

`region_c/five-expert-paper-checks.wl` provides a compact alternative check of
twelve derivative identities for the final two Region C traces. In a fresh
Mathematica kernel, evaluate:

```wolfram
Get["/full/path/to/five-expert-verification/region_c/five-expert-paper-checks.wl"]
```

It is intended to print two lists of six `True` values and a final `PASS` line.
It avoids evaluating antiderivatives by repeatedly updating the local term and
integrand using Leibniz' rule, as described in the manuscript's appendix.

This Wolfram file was **not executed in a Wolfram kernel in the preparation
environment**. The corresponding differentiation formulas are also checked by
the Python supplement. Neither Mathematica nor a Wolfram license is required
to run the Python checks or compile the manuscript.

## Release and licensing

This directory is also maintained as `verification-codes/` alongside the paper
on Overleaf, but is independently runnable. The GitHub repository contains only
the curated code and documentation, not the Overleaf history, working notebooks,
correspondence, or manuscript drafts. No machine-specific paths or external
data are required.

No code license has been selected; licensing is deferred pending the authors'
decision. Public availability is not a grant of an MIT or other open-source
license. The manuscript citation is given above; its public paper link will be
added when available. The initial code release is designated `v1.0.0` and records
the matched manuscript source in `SNAPSHOT.md`.
