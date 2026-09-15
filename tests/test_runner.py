"""Small operational tests; these are not mathematical verification checks."""

from pathlib import Path
import subprocess
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]


class RunnerTests(unittest.TestCase):
    def run_cli(self, *args, optimized=False):
        flags = ["-O"] if optimized else []
        return subprocess.run(
            [sys.executable, *flags, str(ROOT / "run_checks.py"), *args],
            cwd=ROOT.parent, text=True, capture_output=True, timeout=15,
        )

    def test_default_paths_exist(self):
        result = self.run_cli("--list")
        self.assertEqual(result.returncode, 0, result.stderr)
        names = result.stdout.splitlines()
        self.assertEqual(len(names), 13)
        self.assertEqual(len(set(names)), len(names))
        self.assertTrue(all((ROOT / name).is_file() for name in names))

    def test_slow_check_is_opt_in(self):
        result = self.run_cli("--list", "--group", "c", "--include-slow")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(result.stdout.splitlines()), 7)
        self.assertIn("analyze_region_c_residual_pairs.py", result.stdout)

    def test_invalid_resources_rejected(self):
        for args in (("--jobs", "0"), ("--timeout", "0"), ("--timeout", "nan")):
            with self.subTest(args=args):
                self.assertEqual(self.run_cli(*args).returncode, 2)

    def test_optimization_rejected(self):
        result = self.run_cli("--list", optimized=True)
        self.assertEqual(result.returncode, 2)
        self.assertIn("without -O", result.stderr)


if __name__ == "__main__":
    unittest.main()
