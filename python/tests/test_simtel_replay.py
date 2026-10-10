"""Reference adapter boundary tests; the callback fixture is not simulator physics."""

import csv
import os
import runpy
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BUILDER = runpy.run_path(str(ROOT / "tools/build_simtel_replay.py"))


class ReferenceInstrumentationTest(unittest.TestCase):
    def test_instruments_copy_and_preserves_original_source_line_numbers(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder) / "simulator"
            (root / "common").mkdir(parents=True)
            main = '#include "sim_optical_backend.h"\n   if ( (options.save_pe_amp & 0x01) )\n'
            imaging = "".join(
                f"int trace_photon_{prefix}{name} (void)\n{{\n return -1;\n}}\n"
                for prefix, name in (
                    ("in_", "segmented"),
                    ("in_", "paraboloid"),
                    ("in_", "fresnel"),
                    ("with_", "secondary"),
                )
            )
            (root / "common/sim_telarray.c").write_text(main)
            (root / "common/sim_imaging.c").write_text(imaging)
            output = Path(folder) / "build"
            output.mkdir()
            copied_main, copied_imaging = BUILDER["instrument"](root, output)
            self.assertEqual((root / "common/sim_telarray.c").read_text(), main)
            self.assertEqual((root / "common/sim_imaging.c").read_text(), imaging)
            self.assertIn(
                "obdeect_reference_replay(array.optics,array.max_tel)", copied_main.read_text()
            )
            self.assertEqual(copied_imaging.read_text().count("__LINE__"), 4)
            self.assertIn('#line 1 "common/sim_imaging.c"', copied_imaging.read_text())
            self.assertEqual(copied_imaging.read_text().count("\n"), imaging.count("\n") + 2)

    def test_rejects_changed_setup_boundary_and_unknown_recipe(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder) / "source"
            (root / "common").mkdir(parents=True)
            (root / "common/sim_telarray.c").write_text("changed setup")
            with self.assertRaisesRegex(ValueError, "unsupported"):
                BUILDER["instrument"](root, Path(folder))
            with self.assertRaisesRegex(ValueError, "exactly"):
                BUILDER["build"](root, Path(folder) / "build", {})


@unittest.skipUnless(shutil.which("c++"), "C++ compiler is unavailable")
class ReferenceReplayBoundaryTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.folder = tempfile.TemporaryDirectory()
        cls.root = Path(cls.folder.name)
        stub = cls.root / "callback.cpp"
        stub.write_text(r"""
#include <cmath>
extern "C" {
struct ReferenceMeasurement {
 int status,mirror,loss_line;
 double position[3],direction[3],time_ns,relative_efficiency;
 double primary[3],secondary[3],primary_cosine,secondary_cosine;
 double propagation_group_index,primary_envelope,secondary_envelope;
 double upstream_optical_efficiency;
 int pixel_status,pixel_id;
 double pixel_time_ns,camera_absolute_efficiency,pixel_x_cm,pixel_y_cm;
};
int obdeect_reference_replay(void *,unsigned);
void obdeect_reference_measure(void *,unsigned index,const double *p,const double *d,
 double wavelength,double distance,double weight,ReferenceMeasurement *r) {
 r->status=p[0]<0?1:0; r->mirror=3; r->loss_line=r->status?123:0;
 for(int i=0;i<3;++i){r->position[i]=p[i];r->direction[i]=d[i];
 r->primary[i]=NAN;r->secondary[i]=NAN;}
 r->time_ns=7; r->relative_efficiency=weight;
 r->propagation_group_index=1.00023;r->primary_envelope=0.9;r->secondary_envelope=0.8;
 r->upstream_optical_efficiency=0.72;
 r->pixel_status=3;r->pixel_id=1;r->pixel_time_ns=7;r->camera_absolute_efficiency=0.36;
 r->pixel_x_cm=0.;r->pixel_y_cm=0.;
 if(index!=0||wavelength!=400||distance!=1e30)r->status=-1;
}
}
int main(){return obdeect_reference_replay(nullptr,1);}
""")
        cls.executable = cls.root / "replay"
        subprocess.run(
            [
                "c++",
                "-std=c++20",
                "-I" + str(ROOT / "cpp/include"),
                str(ROOT / "tools/simtel_replay.cpp"),
                str(stub),
                "-o",
                str(cls.executable),
            ],
            check=True,
            capture_output=True,
        )

    @classmethod
    def tearDownClass(cls):
        cls.folder.cleanup()

    def replay(self, rows, *, alias=False):
        source = self.root / "source.csv"
        output = source if alias else self.root / "arrivals.csv"
        with source.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow([
                "run_id",
                "event_id",
                "array_id",
                "telescope_id",
                "photon_id",
                "x_m",
                "y_m",
                "z_m",
                "dx",
                "dy",
                "dz",
                "wavelength_nm",
                "time_ns",
                "weight",
                "bunch_id",
            ])
            writer.writerows(rows)
        original = source.read_bytes()
        result = subprocess.run(
            [str(self.executable)],
            capture_output=True,
            text=True,
            env={
                **os.environ,
                "OBDEECT_REPLAY_INPUT": str(source),
                "OBDEECT_REPLAY_OUTPUT": str(output),
                "OBDEECT_REPLAY_TELESCOPE_INDEX": "0",
            },
        )
        self.assertEqual(source.read_bytes(), original)
        if result.returncode:
            return result, []
        with output.open(newline="") as handle:
            return result, list(csv.DictReader(handle))

    def test_preserves_identity_units_time_and_unavailable_loss_fields(self):
        rows = [
            [9, 2, 4, 7, 100, 1.234567890123, 0, 20, 0, 0, -1, 400, 8, 0.5, 12],
            [9, 2, 4, 7, 101, -1, 0, 20, 0, 0, -1, 400, 8, 0.5, 12],
            [9, 2, 5, 7, 100, 2, 0, 20, 0, 0, -1, 400, 8, 0.5, 12],
        ]
        result, output = self.replay(rows)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(output), 3)
        self.assertEqual(output[0]["photon_id"], "100")
        self.assertEqual(output[0]["telescope_id"], "7")
        self.assertEqual(output[0]["bunch_id"], "12")
        self.assertEqual(output[2]["array_id"], "5")
        self.assertAlmostEqual(float(output[0]["camera_x_m"]), rows[0][5], places=14)
        self.assertEqual(float(output[0]["arrival_time_ns"]), 15)
        self.assertEqual(output[1]["status"], "lost_unclassified")
        self.assertEqual(output[1]["loss_source_line"], "123")
        self.assertEqual(float(output[0]["propagation_group_index"]), 1.00023)
        self.assertEqual(float(output[0]["upstream_optical_efficiency"]), 0.72)
        for name in (
            "camera_x_m",
            "camera_dx",
            "arrival_time_ns",
            "relative_efficiency",
            "mirror_index",
            "primary_x_m",
        ):
            self.assertEqual(output[1][name], "")

    def test_rejects_duplicates_unresolved_spectrum_and_source_alias(self):
        row = [9, 2, 4, 7, 100, 1, 0, 20, 0, 0, -1, 400, 8, 1, 12]
        result, _ = self.replay([row, row])
        self.assertIn("duplicate", result.stderr)
        row[11] = 0
        result, _ = self.replay([row])
        self.assertIn("resolved positive wavelength", result.stderr)
        row[11] = 400
        result, _ = self.replay([row], alias=True)
        self.assertIn("aliases", result.stderr)


if __name__ == "__main__":
    unittest.main()
