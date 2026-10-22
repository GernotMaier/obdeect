"""Synthetic telescope parameter records for tests (not CTAO models).

The dual-mirror prescription reuses the published SST design polynomials as a
realistic numerical test case; nothing here is used outside the tests.
"""

from __future__ import annotations

import json
from pathlib import Path


def record(name: str, value, unit=None, *, file: bool = False) -> dict:
    return {
        "parameter": name,
        "parameter_version": "1.0.0",
        "value": value,
        "unit": unit,
        "file": file,
    }


# Seven circular 1 m mirrors on a Davies-Cotton dish; focal lengths automatic.
MIRROR_LIST = """# x y diameter focal shape
   0.0     0.0  100.0  0.0  0
 125.0     0.0  100.0  0.0  0
-125.0     0.0  100.0  0.0  0
  62.5   108.3  100.0  0.0  0
 -62.5   108.3  100.0  0.0  0
  62.5  -108.3  100.0  0.0  0
 -62.5  -108.3  100.0  0.0  0
"""

CYLINDERS = """# %ECSV 1.0
# ---
# datatype:
# - {name: id, datatype: string}
# - {name: group, datatype: string}
# - {name: x1, unit: m, datatype: float64}
# - {name: y1, unit: m, datatype: float64}
# - {name: z1, unit: m, datatype: float64}
# - {name: x2, unit: m, datatype: float64}
# - {name: y2, unit: m, datatype: float64}
# - {name: z2, unit: m, datatype: float64}
# - {name: diameter, unit: m, datatype: float64}
# meta:
#   coordinate_frame: primary_mirror_frame
# schema: astropy-2.0
id group x1 y1 z1 x2 y2 z2 diameter
mast-0 camera_support 0.0 1.6 0.1 0.0 0.5 7.9 0.1
"""

POLYHEDRA = """# %ECSV 1.0
# ---
# datatype:
# - {name: id, datatype: string}
# - {name: group, datatype: string}
# - {name: x1, unit: m, datatype: float64}
# - {name: y1, unit: m, datatype: float64}
# - {name: z1, unit: m, datatype: float64}
# - {name: x2, unit: m, datatype: float64}
# - {name: y2, unit: m, datatype: float64}
# - {name: z2, unit: m, datatype: float64}
# - {name: x3, unit: m, datatype: float64}
# - {name: y3, unit: m, datatype: float64}
# - {name: z3, unit: m, datatype: float64}
# - {name: x4, unit: m, datatype: float64}
# - {name: y4, unit: m, datatype: float64}
# - {name: z4, unit: m, datatype: float64}
# meta:
#   coordinate_frame: primary_mirror_frame
# schema: astropy-2.0
id group x1 y1 z1 x2 y2 z2 x3 y3 z3 x4 y4 z4
plate-0 dish_structure 1.0 -1.0 3.0 1.2 -1.0 3.0 1.0 -0.8 3.0 1.0 -1.0 3.2
"""


def segmented(**changes) -> tuple[dict, dict]:
    """A small single-mirror telescope with f = 8 m: records and file contents."""
    records = {
        "mirror_class": record("mirror_class", 0),
        "focal_length": record("focal_length", 800.0, "cm"),
        "mirror_offset": record("mirror_offset", -50.0, "cm"),
        "dish_shape_length": record("dish_shape_length", 800.0, "cm"),
        "parabolic_dish": record("parabolic_dish", False),
        "mirror_focal_length": record("mirror_focal_length", 0.0, "cm"),
        "mirror_list": record("mirror_list", "mirrors.dat", file=True),
        "focus_offset": record(
            "focus_offset", [0.0, 20.0, 0.0, 0.0], ["cm", "deg", "null", "null"]
        ),
        "camera_body_diameter": record("camera_body_diameter", 60.0, "cm"),
        "camera_body_shape": record("camera_body_shape", 2),
        "camera_depth": record("camera_depth", 0.0, "cm"),
        "telescope_transmission": record("telescope_transmission", [0.9, 1, 0.1, 5.0, 2.0, 0.0]),
    }
    files = {"mirror_list": MIRROR_LIST}
    for name, value in changes.items():
        if isinstance(value, str) and name in (
            "telescope_obscuration_cylinders",
            "telescope_obscuration_quadrilaterals",
        ):
            records[name] = record(name, f"{name}.ecsv", file=True)
            files[name] = value
        elif value is None:
            records.pop(name, None)
        else:
            records[name] = value
    return records, files


def dual(**changes) -> tuple[dict, dict]:
    """An SST-like Schwarzschild-Couder telescope (published SST design polynomials)."""
    cm = ["cm"] * 11
    records = {
        "mirror_class": record("mirror_class", 2),
        "focal_length": record("focal_length", 215.0, "cm"),
        "mirror_offset": record("mirror_offset", 0.0, "cm"),
        "primary_mirror_parameters": record(
            "primary_mirror_parameters",
            [
                0.0,
                0.000608051,
                -7.36248e-10,
                5.82125e-15,
                -6.76447e-20,
                -3.89542e-24,
                -5.28036e-29,
                2.99107e-34,
                4.39153e-38,
                6.17433e-43,
                -2.73586e-47,
            ],
            cm,
        ),
        "primary_mirror_diameter": record("primary_mirror_diameter", 424.1, "cm"),
        "primary_mirror_hole_diameter": record("primary_mirror_hole_diameter", 0.0, "cm"),
        "secondary_mirror_parameters": record(
            "secondary_mirror_parameters",
            [
                310.84,
                0.00229358,
                2.8273e-08,
                -2.7689e-12,
                8.80066e-17,
                3.37315e-21,
                -1.02973e-25,
                -6.72882e-30,
                -3.06437e-34,
                3.1718e-38,
                -3.71217e-43,
            ],
            cm,
        ),
        "secondary_mirror_diameter": record("secondary_mirror_diameter", 180.0, "cm"),
        "secondary_mirror_hole_diameter": record("secondary_mirror_hole_diameter", 0.0, "cm"),
        "focal_surface_parameters": record(
            "focal_surface_parameters",
            [
                258.88,
                -0.00471698,
                -1.04952e-07,
                -4.67036e-12,
                -2.59788e-16,
                -1.61847e-20,
                -1.08033e-24,
                -7.55454e-29,
                -5.46285e-33,
                -4.0516e-37,
                -3.06502e-41,
            ],
            cm,
        ),
        "camera_body_diameter": record("camera_body_diameter", 57.0, "cm"),
        "camera_body_shape": record("camera_body_shape", 2),
        "camera_depth": record("camera_depth", 50.0, "cm"),
        "secondary_mirror_shadow_diameter": record("secondary_mirror_shadow_diameter", 214.0, "cm"),
        "secondary_mirror_shadow_offset": record("secondary_mirror_shadow_offset", 0.0, "cm"),
        "secondary_mirror_baffle": record(
            "secondary_mirror_baffle", [294.8, 273.6, 97.0, 10.0, 0.0], "cm"
        ),
        "telescope_transmission": record("telescope_transmission", [0.96, 1, 0.0368, 6.795, 0.45]),
    }
    records.update(changes)
    return records, {}


def write_simulation_models(
    root: Path, model: str, records: dict, files: dict, *, version: str = "9.9.9"
) -> None:
    """A minimal simulation-models checkout holding one telescope model."""
    (root / "productions" / version).mkdir(parents=True)
    table = {name: rec["parameter_version"] for name, rec in records.items()}
    manifest = {
        "model_version": version,
        "production_table_name": model,
        "parameters": {model: table},
    }
    (root / "productions" / version / f"{model}.json").write_text(json.dumps(manifest))
    for name, rec in records.items():
        directory = root / "model_parameters" / model / name
        directory.mkdir(parents=True)
        (directory / f"{name}-{rec['parameter_version']}.json").write_text(
            json.dumps({**rec, "instrument": model})
        )
        if rec.get("file"):
            (directory / rec["value"]).write_text(files[name])
