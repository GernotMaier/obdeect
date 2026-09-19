"""Regression tests for examples/run_examples.py CSV validation."""

import importlib.util
import unittest
from pathlib import Path

SCRIPT = Path(__file__).parents[2] / "examples" / "run_examples.py"
SPEC = importlib.util.spec_from_file_location("run_examples", SCRIPT)
RUN_EXAMPLES = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(RUN_EXAMPLES)


class TestRunExamplesValidation(unittest.TestCase):
    def test_allows_zero_detected_when_case_does_not_require_hits(self):
        rows = [{"status": "blocked_camera", "point_count": "2"} for _ in range(4)]
        counts = RUN_EXAMPLES.validate_rows("toy_flasher", rows, photons=4, require_detected=False)
        self.assertEqual(counts, {"blocked_camera": 4})

    def test_rejects_zero_detected_when_case_requires_hits(self):
        rows = [{"status": "blocked_camera", "point_count": "2"} for _ in range(4)]
        with self.assertRaisesRegex(ValueError, "no focal-plane hits"):
            RUN_EXAMPLES.validate_rows("toy_star", rows, photons=4, require_detected=True)


if __name__ == "__main__":
    unittest.main()
