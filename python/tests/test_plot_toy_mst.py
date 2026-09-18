"""Dependency-free checks for the optional path-visualization data reader."""

import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "obdeect" / "plotting.py"
SPEC = importlib.util.spec_from_file_location("plot_toy_mst", SCRIPT)
PLOT = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(PLOT)


class TestTracePathReader(unittest.TestCase):
    def test_reads_actual_ragged_trace_vertices(self):
        # T-VIS-009: no optical calculation is performed by the plot reader.
        csv_text = (
            "photon_id,wavelength_nm,status,point_count,path_length_m,x0_m,y0_m,z0_m,x1_m,y1_m,z1_m,x2_m,y2_m,z2_m\n"
            "4,400,detected,3,20,1,2,20,1,2,0,0.1,0.2,4.875\n"
            "5,400,blocked_mast,2,10,3,4,20,3,4,8,0,0,0\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "paths.csv"
            path.write_text(csv_text)
            paths = list(PLOT.read_paths(path))
        self.assertEqual(paths[0], ("detected", [(1.0, 2.0, 20.0), (1.0, 2.0, 0.0), (0.1, 0.2, 4.875)]))
        self.assertEqual(paths[1], ("blocked_mast", [(3.0, 4.0, 20.0), (3.0, 4.0, 8.0)]))


if __name__ == "__main__":
    unittest.main()
