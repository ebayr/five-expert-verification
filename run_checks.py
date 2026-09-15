#!/usr/bin/env python3
"""Run the curated symbolic checks; no notebooks or external data required."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import math
from pathlib import Path
import subprocess
import sys
import time


CHECKS = {
    "construction": (
        "construction/verify_construction.py",
    ),
    "ab": (
        "regions_ab/ab_a2a6_verify.py",
        "regions_ab/ab_a2a6_highorder.py",
        "regions_ab/ab_a7a16_certificates.py",
        "regions_ab/ab_a7a16_polynomials.py",
        "regions_ab/ab_bother_verify.py",
        "regions_ab/ab_b34_exact_checks.py",
    ),
    "c": (
        "region_c/verify_p_boundary.py",
        "region_c/verify_hyperbolic_system.py",
        "region_c/verify_region_c_single_integrals.py",
        "region_c/verify_region_c_zero_modes.py",
        "region_c/derive_region_c_remaining_traces.py",
        "region_c/verify_region_c_final_monotonicity.py",
    ),
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--group", choices=("all", *CHECKS), default="all")
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--timeout", type=float, default=300)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--list", action="store_true", help="list selected checks without running them")
    parser.add_argument("--include-slow", action="store_true",
                        help="also run the full Region C residual classification")
    args = parser.parse_args()
    if args.jobs < 1 or not math.isfinite(args.timeout) or args.timeout <= 0:
        parser.error("--jobs and --timeout must be positive")
    if not __debug__ or sys.flags.optimize:
        parser.error("run without -O: individual checks use assertions")
    base = Path(__file__).resolve().parent
    selected = [name for group, names in CHECKS.items()
                if args.group in ("all", group) for name in names]
    if args.include_slow and args.group in ("all", "c"):
        selected.append("region_c/analyze_region_c_residual_pairs.py")
    if args.list:
        print("\n".join(selected))
        return 0
    try:
        import sympy
        import mpmath
    except ImportError:
        parser.exit(2, "Install dependencies with: python -m pip install -r requirements.txt\n")

    print(f"Python {sys.version.split()[0]}; SymPy {sympy.__version__}; "
          f"mpmath {mpmath.__version__}", flush=True)
    def run(name):
        started = time.monotonic()
        try:
            result = subprocess.run(
                [sys.executable, "-B", str(base / name)],
                cwd=base, text=True, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, timeout=args.timeout, check=False,
            )
            success, output = result.returncode == 0, result.stdout
        except subprocess.TimeoutExpired:
            success, output = False, f"Timed out after {args.timeout:g} seconds.\n"
        except OSError as error:
            success, output = False, f"Could not run check: {error}\n"
        return name, success, output, time.monotonic() - started

    failures = []
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        for name, success, output, elapsed in pool.map(run, selected):
            print(f"{'PASS' if success else 'FAIL'} {name} ({elapsed:.2f}s)", flush=True)
            if args.verbose or not success:
                print(output, end="" if output.endswith("\n") else "\n", flush=True)
            if not success:
                failures.append(name)
    print(f"{len(selected) - len(failures)}/{len(selected)} check scripts passed.")
    print("Scope: symbolic checks only; see README.md for analytic obligations.")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
