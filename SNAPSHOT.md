# Release v1.0.0

Prepared on September 15, 2026 for *Prediction with Five Experts and Geometric
Stopping: A Probabilistic Construction and Analytic Verification*, by Erhan
Bayraktar, Ibrahim Ekren, and Nikolaos Kolliopoulos.

## Matched manuscript

Source file: `ArXiv/five-expert-paper.tex` in the authors' paper project.
SHA-256 of that source file (not of the PDF):

```text
bc18d187c8f3fc83c1b1474ad16cf1920c8ba492431a93d9039db57979228083
```

This snapshot includes the shorter A2 collision-Hessian endpoint argument,
the A7/A8 identity `K'=192 sinh(t)^6`, the direct B8 comparison argument,
and the paper's code-availability paragraph. It does not claim an arXiv
submission date or identifier. The manuscript itself is not redistributed
in this code repository.

## Local validation

All 14 mathematical entry points passed with Python 3.12.14, SymPy 1.14.0,
and mpmath 1.3.0. They were run through `run_checks.py` in its three groups
from a working directory outside the package; the Region C run included
`--include-slow`. The construction entry point checks 45 exact identities.
The full Region C classification took approximately 154 seconds on the
preparation machine. Four operational runner tests also passed.

These are exact symbolic regression checks, not a formal verification of
the entire theorem. See [CHECKS.md](CHECKS.md) for coverage and exclusions.
The Wolfram alternative was not executed in a Wolfram kernel.

The GitHub workflow repeats the full mathematical suite and runner tests;
its actual status and logs are available in the repository's Actions tab.

## Reproducing this release

Select tag `v1.0.0`, install `requirements.txt`, then run:

```sh
python run_checks.py --include-slow --jobs 2 --timeout 1800
python -m unittest discover -s tests -v
```

Licensing is intentionally undecided pending the authors' agreement.
