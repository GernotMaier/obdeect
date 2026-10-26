"""Compile simulation-models parameter records into an obdeect optical model.

The optical model reproduces sim_telarray's telescope set-up (``tel_setup_*``
in ``sim_imaging.c`` and ``sim_config.c``) up to the camera entrance:
mirrors, their random errors, shadowing and optical responses. Camera filter,
light guides, quantum efficiency and pixels belong to the camera simulation.

The input is a mapping from parameter name to its simulation-models record
plus a function returning the contents of a parameter's file, so the same
compiler serves a simulation-models checkout and simtools.
"""

from __future__ import annotations

import functools
import hashlib
import json
import math
import random
from collections.abc import Callable, Mapping
from statistics import NormalDist
from typing import Any

from obdeect.formats import (
    SHAPES,
    ModelError,
    components,
    length_m,
    length_unit_m,
    lengths_m,
    number,
    read_mirror_list,
    read_response,
    read_segmentation,
    read_spatial_map,
    read_structure_table,
    secondary_segment_frame,
    segmentation_from_record,
    signed_length_m,
)
from obdeect.simtel_tables import evaluate

MODEL_FORMAT = "obdeect.optical-model.v2"
# Structure table (masts, quadrupods, dish frame); see docs/optical-elements.md#telescope-structure.
STRUCTURE_PARAMETER = "telescope_structure"
ALL_LEGS = ["incoming", "after_primary", "after_secondary"]


class _Parameters:
    """Parameter records; remembers which ones the compiler used."""

    def __init__(self, records: Mapping[str, dict[str, Any]], read_file: Callable[[str], str]):
        self._records = records
        self._read_file = read_file
        self.used: set[str] = set()
        self.files: dict[str, str] = {}

    def __contains__(self, name: str) -> bool:
        return name in self._records

    def record(self, name: str) -> dict[str, Any]:
        if name not in self._records:
            raise ModelError(f"missing parameter {name}")
        self.used.add(name)
        return self._records[name]

    def value(self, name: str, default: Any = None) -> Any:
        return self.record(name).get("value") if name in self else default

    def file(self, name: str) -> tuple[str, str]:
        """Contents of the parameter's file and its RPOL options (``file#rpol:options``)."""
        value = self.value(name)
        if not isinstance(value, str) or not value or value.lower() == "none":
            raise ModelError(f"parameter {name} names no file")
        contents = self._read_file(name)
        self.files[name] = hashlib.sha256(contents.encode()).hexdigest()
        return contents, value.partition("#rpol:")[2]

    def has_file(self, name: str) -> bool:
        value = self.value(name) if name in self else None
        return isinstance(value, str) and bool(value) and value.lower() != "none"


def compile_model(
    records: Mapping[str, dict[str, Any]],
    read_file: Callable[[str], str],
    *,
    configuration_seed: int,
    zenith_deg: float | None = None,
    propagation_group_index: float = 1.0,
    provenance: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Return the optical model as a JSON-serialisable dictionary.

    Parameters
    ----------
    records
        simulation-models parameter records by name (``value``, ``unit``, ...).
    read_file
        Returns the text of the file named by a parameter (its ``value``).
    configuration_seed
        Seed of the one-time random mirror errors (focal lengths, placement,
        alignment). The realised geometry is stored in the model.
    zenith_deg
        Telescope zenith angle for zenith-dependent focus and alignment errors;
        None uses the reference angle of each parameter.
    propagation_group_index
        Group refractive index of air at the telescope, for arrival times.
    provenance
        Identity of the model (name, version, site), stored unchanged.
    """
    if isinstance(configuration_seed, bool) or not isinstance(configuration_seed, int):
        raise ModelError("configuration seed must be an integer")
    if not 0 <= configuration_seed < 2**64:
        raise ModelError("configuration seed must fit in 64 bits")
    p = _Parameters(records, read_file)
    notes: list[str] = []
    mirror_class = p.value("mirror_class", 0)
    if mirror_class == 0:
        optics = _segmented_telescope(p, configuration_seed, zenith_deg)
    elif mirror_class == 2:
        optics = _dual_mirror_telescope(p, zenith_deg)
    else:
        raise ModelError(f"mirror_class {mirror_class} is not supported")

    next_id = 1 + max(
        _surface_ids(optics)
    )  # Obscurers follow the optical surfaces in the id space.
    obscurers = optics.pop("obscurers")
    structure = _structure(p, next_id + len(obscurers))
    obscurers += structure
    for index, obscurer in enumerate(obscurers[: len(obscurers) - len(structure)]):
        obscurer["id"] = next_id + index
    optics["obscurers"] = obscurers

    transmission = None
    if "telescope_transmission" in p:
        if structure:
            p.used.add("telescope_transmission")
            notes.append(
                "telescope_transmission is not applied: "
                "the model's structure geometry shadows instead"
            )
        else:
            transmission = _telescope_transmission(p, dual=mirror_class == 2)
    optics["telescope_transmission"] = transmission
    if not propagation_group_index >= 1.0:
        raise ModelError("propagation group index must be >= 1")
    optics["propagation_group_index"] = propagation_group_index
    model = {
        "format": MODEL_FORMAT,
        "sampling": _sampling(p, optics),
        "provenance": {
            **(provenance or {}),
            "configuration_seed": configuration_seed,
            "zenith_deg": zenith_deg,
            "structure": bool(structure),
            "parameters": {name: records[name].get("parameter_version") for name in sorted(p.used)},
            "files_sha256": dict(sorted(p.files.items())),
        },
        "notes": notes,
        "optics": optics,
    }
    model["sha256"] = model_hash(model)
    return model


def _sampling(p: _Parameters, optics: dict[str, Any]) -> dict[str, float]:
    """sim_telarray star-light sampling (sim_config.c setup_starlight): a disk of
    1.2 times the dish radius, perpendicular to the optical axis at the dish
    reference height; rays start above every surface and obscurer."""
    overall = (
        signed_length_m(p.record("mirror_offset"), "mirror_offset") if "mirror_offset" in p else 0.0
    )
    primary = optics["primary"]
    if "facets" in primary:
        dish = max(
            math.hypot(*f["centre_m"][:2]) + 0.5 * f["diameter_m"] for f in primary["facets"]
        )
    else:
        dish = primary["asphere"]["outer_radius_m"]
    heights = [optics["focal_surface"].get("plane_z_m", 0.0)]
    for role in ("primary", "secondary"):
        mirror = optics.get(role) or {}
        heights += [f["centre_m"][2] for f in mirror.get("facets", [])]
        if "asphere" in mirror:
            heights.append(mirror["asphere"]["vertex_z_m"])
    if "asphere" in optics["focal_surface"]:
        heights.append(optics["focal_surface"]["asphere"]["vertex_z_m"])
    for obscurer in optics["obscurers"]:
        for key in ("start_m", "end_m", "centre_m"):
            if key in obscurer:
                heights.append(obscurer[key][2])
        heights += [v[2] for v in obscurer.get("vertices_m", [])]
    return {
        "disk_radius_m": 1.2 * dish,
        "disk_z_m": -overall,
        "start_z_m": max(heights) + 1.0,
    }


def model_hash(model: dict[str, Any]) -> str:
    """SHA-256 of the canonical JSON of the model without its ``sha256`` entry."""
    content = {key: value for key, value in model.items() if key != "sha256"}
    text = json.dumps(content, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(text.encode()).hexdigest()


def _surface_ids(optics: dict[str, Any]) -> list[int]:
    ids = [optics["focal_surface"]["id"]]
    for role in ("primary", "secondary"):
        mirror = optics.get(role)
        if not mirror:
            continue
        ids += [f["id"] for f in mirror.get("facets", [])]
        ids += [s["id"] for s in mirror.get("segments", [])]
        if "asphere_id" in mirror:
            ids.append(mirror["asphere_id"])
    return ids


# --- Segmented single-mirror telescopes (LST, MST): sim_imaging.c tel_setup_primary ---


def _segmented_telescope(p: _Parameters, seed: int, zenith_deg: float | None) -> dict[str, Any]:
    focal_length = length_m(p.record("focal_length"), "focal_length")
    overall_offset = signed_length_m(p.record("mirror_offset"), "mirror_offset")
    facets = read_mirror_list(p.file("mirror_list")[0])
    fallback = (
        length_m(p.record("mirror_focal_length"), "mirror_focal_length", allow_zero=True)
        if ("mirror_focal_length" in p)
        else 0.0
    )
    for facet in facets:
        facet["focal_length_m"] = abs(facet["catalogue_focal_length_m"]) or fallback
    _panel_prescription(facets, p, seed)
    _nominal_placement(facets, p, focal_length, overall_offset)
    _panel_alignment(facets, p, seed, zenith_deg)
    primary: dict[str, Any] = {
        "facets": [
            {
                "id": f["id"],
                "centre_m": f["centre_m"],
                "normal": f["normal"],
                "tangent": f["tangent"],
                "map_tangent": f["map_tangent"],
                "shape": f["shape"],
                "diameter_m": f["diameter_m"],
                "radius_of_curvature_m": 2.0 * f["focal_length_m"],
            }
            for f in facets
        ]
    }
    _mirror_response(p, primary, "mirror_reflectivity", "mirror_degraded_reflection")
    _scatter(p, primary, surface_slopes=False)
    for name in ("primary_mirror_degraded_map", "primary_degraded_map"):
        if p.has_file(name):
            contents, filename, options = _file_args(p, name)
            primary["degradation"] = read_spatial_map(contents, filename, filename_options=options)
    focus = _focus_offset(p, zenith_deg)
    obscurers = _camera_body(p, focal_length - overall_offset, after_primary=False)
    return {
        "primary": primary,
        "secondary": None,
        "secondary_blocks_incoming": False,
        "focal_surface": {"id": len(facets), "plane_z_m": focal_length + focus - overall_offset},
        "obscurers": obscurers,
    }


def _file_args(p: _Parameters, name: str) -> tuple[str, str, str]:
    contents, options = p.file(name)
    return contents, name, options


def _panel_prescription(facets: list[dict[str, Any]], p: _Parameters, seed: int) -> None:
    """Focal lengths and distances to the focus, with sim_telarray's random errors."""
    if p.value("flip_mirrors", False) not in (False, True, 0, 1):
        raise ModelError("flip_mirrors must be boolean")
    if p.value("flip_mirrors", False):
        for facet in facets:
            facet["centre_m"][:2] = facet["centre_m"][:2][::-1]
            facet["shape"] = {
                "hexagon_flat_x": "hexagon_flat_y",
                "hexagon_flat_y": "hexagon_flat_x",
            }.get(facet["shape"], facet["shape"])
    focal = length_m(p.record("focal_length"), "focal_length")
    dish = length_m(p.record("dish_shape_length"), "dish_shape_length", allow_zero=True) or focal
    parabolic = p.value("parabolic_dish")
    if not isinstance(parabolic, bool):
        raise ModelError("parabolic_dish must be boolean")
    radii = [math.hypot(*f["centre_m"][:2]) for f in facets]
    smallest, largest = min(radii), max(radii)
    grading = (
        signed_length_m(p.record("grading_of_focal_length"), "grading_of_focal_length")
        if "grading_of_focal_length" in p
        else 0.0
    )
    distance_sigma = (
        length_m(
            p.record("mirror_align_random_distance"),
            "mirror_align_random_distance",
            allow_zero=True,
        )
        if "mirror_align_random_distance" in p
        else 0.0
    )
    sigma, limit = (
        lengths_m(p.record("random_focal_length"), "random_focal_length", 2)
        if "random_focal_length" in p
        else (0.0, 0.0)
    )
    if sigma < 0 or limit < 0:
        raise ModelError("random_focal_length widths must be non-negative")
    scales = p.value("mirror_f_scale", [1.0] * len(facets))
    if not isinstance(scales, list):
        scales = [scales] * len(facets)
    if len(scales) != len(facets) or any(number(v, "mirror_f_scale") <= 0 for v in scales):
        raise ModelError("mirror_f_scale needs one positive scale per mirror")
    opt = [number(v, "mirror_opt") for v in p.value("mirror_opt", [0.0, 0.0, 0.0])]
    if len(opt) != 3:
        raise ModelError("mirror_opt needs three coefficients")
    for facet, radius, scale in zip(facets, radii, scales, strict=True):
        z_file = facet["centre_m"][2]
        if z_file:
            height = abs(z_file)
        elif parabolic:
            height = radius**2 / (4 * dish)
        else:
            if radius > dish:
                raise ModelError("a mirror lies outside the dish sphere")
            height = dish - math.sqrt(dish * dish - radius * radius)
        nominal_distance = math.hypot(focal - height, radius)
        distance = nominal_distance * (1 + opt[2] + opt[1] * radius**2 / focal**2)
        rng = random.Random(f"obdeect-panel-distance-v1:{seed}:{facet['id']}")
        if z_file <= 0 and distance_sigma:
            distance += _gauss(rng) * distance_sigma
        if not distance > radius:
            raise ModelError("a mirror's distance to the focus must exceed its radial offset")
        facet["distance_m"] = distance
        catalogue = facet["catalogue_focal_length_m"]
        if catalogue > 0:
            facet["focal_length_m"] = catalogue * scale
            continue
        nominal = (
            abs(catalogue) * scale if catalogue < 0 else facet["focal_length_m"] or nominal_distance
        )
        if catalogue == 0 and largest > smallest:
            nominal += grading * (radius - 0.5 * (smallest + largest)) / (largest - smallest)
        rng = random.Random(f"obdeect-panel-focal-v1:{seed}:{facet['id']}")
        if limit > 100 * sigma:
            error = limit * (rng.random() - 0.5)
        elif sigma:
            # Clipped Gaussian by inverse CDF; sim_telarray rejects and redraws.
            lo = NormalDist().cdf(-limit / sigma) if limit else 0.0
            hi = NormalDist().cdf(limit / sigma) if limit else 1.0
            u = min(max(lo + (hi - lo) * rng.random(), 1e-300), 1 - 1e-16)
            error = sigma * NormalDist().inv_cdf(u)
        else:
            error = 0.0
        facet["focal_length_m"] = nominal + error
        if not facet["focal_length_m"] > 0:
            raise ModelError("a mirror's focal length must be positive")


def _gauss(rng: random.Random) -> float:
    """Standard normal deviate from uniforms only (stable across Python versions)."""
    return math.sqrt(-2.0 * math.log(1.0 - rng.random())) * math.cos(2.0 * math.pi * rng.random())


def _nominal_placement(
    facets: list[dict[str, Any]], p: _Parameters, focal: float, overall: float
) -> None:
    """Mirror centres on the dish (sim_imaging.c tel_setup_primary)."""
    opt = [number(v, "mirror_opt") for v in p.value("mirror_opt", [0.0, 0.0, 0.0])]
    for facet in facets:
        x, y, z_file = facet["centre_m"]
        radius = math.hypot(x, y)
        distance = facet["distance_m"]
        z = z_file if z_file > 0 else focal - math.sqrt(distance * distance - radius * radius)
        facet["centre_m"] = [x, y, z - overall]
        facet["inclination"] = (
            0.5 * math.asin(radius / distance) * (1 + opt[0] * radius**2 / focal**2)
        )


def _panel_alignment(
    facets: list[dict[str, Any]], p: _Parameters, seed: int, zenith_deg: float | None
) -> None:
    """Random mirror tilts and the mirror frame (Euler angles of sim_telarray's make_trafo)."""
    widths = []
    for name in ("mirror_align_random_horizontal", "mirror_align_random_vertical"):
        if name not in p:
            widths.append(0.0)
            continue
        values, units = components(p.record(name), name, 4)
        if units[:2] != ["deg", "deg"]:
            raise ModelError(f"{name}: the first two components must be in degrees")
        rms, reference, cosine, sine = (math.radians(number(v, name)) for v in values)
        theta = math.radians(zenith_deg) if zenith_deg is not None else reference
        widths.append(
            math.sqrt(
                rms**2
                + cosine**2 * (math.cos(theta) - math.cos(reference)) ** 2
                + sine**2 * (math.sin(theta) - math.sin(reference)) ** 2
            )
        )
    for facet in facets:
        x, y, _ = facet["centre_m"]
        phi = math.atan2(y, x)
        rng = random.Random(f"obdeect-panel-alignment-v1:{seed}:{facet['id']}")
        horizontal, vertical = (_gauss(rng) * width for width in widths)
        inclination = (
            facet.pop("inclination") - horizontal * math.sin(phi) + vertical * math.cos(phi)
        )
        gamma = horizontal * math.cos(phi) + vertical * math.sin(phi)
        c, t = math.cos(phi), math.sin(phi)
        ci, si, cg, sg = (
            math.cos(inclination),
            math.sin(inclination),
            math.cos(gamma),
            math.sin(gamma),
        )
        facet["normal"] = [-c * si * cg - t * sg, -t * si * cg + c * sg, ci * cg]
        u = [c * ci, t * ci, si]
        v = [-t * cg + c * si * sg, c * cg + t * si * sg, -ci * sg]
        # Aperture frame rotated by the mirror azimuth (sim_telarray's phi_m).
        facet["tangent"] = [c * a - t * b for a, b in zip(u, v, strict=True)]
        facet["map_tangent"] = u


def _focus_offset(p: _Parameters, zenith_deg: float | None) -> float:
    """Camera focus offset; the zenith terms apply only when both are non-zero,
    as in sim_config.c."""
    if "focus_offset" not in p:
        return 0.0
    values, units = components(p.record("focus_offset"), "focus_offset", 4)
    unit = length_unit_m(units[0], "focus_offset")
    base = number(values[0], "focus_offset") * unit
    if not (values[2] and values[3]):
        return base
    if zenith_deg is None:
        raise ModelError("a zenith-dependent focus offset needs the telescope zenith angle")
    theta, reference = math.radians(zenith_deg), math.radians(number(values[1], "focus_offset"))
    return base + unit * (
        number(values[2], "focus_offset") * (math.cos(theta) - math.cos(reference))
        + number(values[3], "focus_offset") * (math.sin(theta) - math.sin(reference))
    )


def _camera_body(p: _Parameters, front_z: float, *, after_primary: bool) -> list[dict[str, Any]]:
    """Camera body shadow as sim_telarray tests it.

    Single-mirror telescopes test the projection of the incoming ray onto the
    camera front plane and, with a camera depth, onto its back plane. Dual-mirror
    telescopes test the rays between the mirrors against the camera back plane
    and an open circular side wall.
    """
    if "camera_body_diameter" not in p:
        return []
    diameter = length_m(p.record("camera_body_diameter"), "camera_body_diameter", allow_zero=True)
    if diameter == 0:
        return []
    shape = SHAPES.get(p.value("camera_body_shape", 0))
    if shape is None:
        raise ModelError("unsupported camera_body_shape")
    depth = (
        length_m(p.record("camera_depth"), "camera_depth", allow_zero=True)
        if "camera_depth" in p
        else 0.0
    )

    def plate(name: str, z: float, legs: list[str]) -> dict[str, Any]:
        return {
            "name": name,
            "group": "camera",
            "legs": legs,
            "shape": "plate",
            "centre_m": [0.0, 0.0, z],
            "normal": [0.0, 0.0, 1.0],
            "tangent": [1.0, 0.0, 0.0],
            "aperture": shape,
            "diameter_m": diameter,
        }

    if not after_primary:
        planes = [plate("camera-front", front_z, ["incoming"])]
        if depth > 0:
            planes.append(plate("camera-back", front_z + depth, ["incoming"]))
        return planes
    back = front_z - depth
    bodies = [plate("camera-back", back, ["after_primary"])]
    if depth > 0:
        bodies.append({
            "name": "camera-side-wall",
            "group": "camera",
            "legs": ["after_primary"],
            "shape": "frustum",
            "start_m": [0.0, 0.0, back],
            "end_m": [0.0, 0.0, front_z],
            "start_radius_m": 0.5 * diameter,
            "end_radius_m": 0.5 * diameter,
            "thickness_m": 0.0,
            "solid": False,
        })
    return bodies


def _mirror_response(
    p: _Parameters, mirror: dict[str, Any], table: str, degraded: str, same_as: str | None = None
) -> None:
    """Reflectivity resampled on a 1 nm grid, times the degradation factor.

    A secondary mirror without its own table (value "same" or no parameter)
    reflects like the primary, as in sim_telarray.
    """
    if same_as and (table not in p or str(p.value(table)).lower() == "same"):
        table = same_as
    if not p.has_file(table):
        return
    contents, options = p.file(table)
    factor = number(p.value(degraded, 1.0), degraded) if degraded in p else 1.0
    if not 0 <= factor <= 1:
        raise ModelError(f"{degraded} must be a fraction")
    mirror["reflectivity"] = json.loads(_response_json(contents, table, options, factor))


@functools.lru_cache(maxsize=16)
def _response_json(contents: str, name: str, options: str, factor: float) -> str:
    """Resampled response, cached: models differing only in mirror errors share it."""
    return json.dumps(
        response_grid(read_response(contents, name, filename_options=options), factor)
    )


def response_grid(table: dict[str, Any], factor: float = 1.0) -> dict[str, Any]:
    """Resample a response table onto a 1 nm wavelength grid.

    1-D tables use their declared RPOL interpolation; 2-D tables are bilinear,
    as in sim_telarray.
    """
    x, y, values = table["x"], table["y"], table["values"]
    start, stop = math.floor(x[0]), math.ceil(x[-1])
    grid = [float(w) for w in range(start, stop + 1)]
    if not y:
        rows = [[factor * evaluate(x, values, table["interpolation"], w)] for w in grid]
    else:
        ny = len(y)
        columns = [[values[i * ny + j] for i in range(len(x))] for j in range(ny)]
        linear = {"boundary": "clamp", "scheme": 1}
        rows = [[factor * evaluate(x, column, linear, w) for column in columns] for w in grid]
    return {
        "wavelength_min_nm": grid[0],
        "wavelength_step_nm": 1.0,
        "incidence_deg": list(y) or [0.0],
        "values": [[min(1.0, max(0.0, v)) for v in row] for row in rows],
    }


def _scatter(p: _Parameters, mirror: dict[str, Any], *, surface_slopes: bool) -> None:
    if "mirror_reflection_random_angle" not in p:
        return
    name = "mirror_reflection_random_angle"
    values, units = components(p.record(name), name, 3)
    if units[0] != "deg" or units[2] not in ("deg", None):
        raise ModelError(f"{name}: the scattering angles must be in degrees")
    first, fraction, second = (number(v, "mirror_reflection_random_angle") for v in values)
    if not (0 <= first < 90 and 0 <= second < 90 and 0 <= fraction < 1):
        raise ModelError("invalid mirror_reflection_random_angle")
    mirror["scatter"] = {
        "sigma1_rad": math.radians(first),
        "fraction2": fraction,
        "sigma2_rad": math.radians(second),
        "method": "surface_slopes" if surface_slopes else "outgoing_angles",
    }


def _telescope_transmission(p: _Parameters, *, dual: bool) -> dict[str, float]:
    """sim_config.c conventions: method 1 enables the angular term; a zero scale
    means the camera half-angle."""
    values = p.value("telescope_transmission")
    values = values if isinstance(values, list) else [values]
    if not 1 <= len(values) <= 6:
        raise ModelError("telescope_transmission requires 1 to 6 components")
    t0, method, amplitude, scale, power, outer = [
        number(v, "telescope_transmission") for v in values
    ] + [0.0] * (6 - len(values))
    if not 0 <= t0 <= 1 or method not in (0, 1) or amplitude < 0:
        raise ModelError("invalid telescope_transmission")
    amplitude = amplitude if method == 1 else 0.0
    if not scale:
        body = (
            length_m(p.record("camera_body_diameter"), "camera_body_diameter", allow_zero=True)
            if "camera_body_diameter" in p
            else 0.0
        )
        scale_rad = 0.5 * body / length_m(p.record("focal_length"), "focal_length")
    else:
        scale_rad = math.radians(scale)
    if amplitude and not scale_rad > 0:
        raise ModelError("an angular telescope_transmission needs a positive scale")
    return {
        "on_axis": t0,
        "amplitude": amplitude,
        "angular_scale_rad": scale_rad or 1.0,
        "power": power or 2.0,
        # The dual-mirror tracer uses the original form without outer power.
        "outer_power": 1.0 if dual else (outer or 1.0),
    }


# --- Dual-mirror telescopes (SST, SCT): sim_imaging.c tel_setup_secondary ---


def _dual_mirror_telescope(p: _Parameters, zenith_deg: float | None) -> dict[str, Any]:
    overall = (
        signed_length_m(p.record("mirror_offset"), "mirror_offset") if "mirror_offset" in p else 0.0
    )
    primary_surface = _asphere(p, "primary_mirror")
    secondary_surface = _asphere(p, "secondary_mirror")
    focal_surface = _asphere(p, "focal_surface", aperture=False)
    # sim_telarray rotates the secondary frame by pi about y: in the telescope
    # frame its sag changes sign.
    coefficients = secondary_surface["coefficients_m"]
    if coefficients[0] < primary_surface["coefficients_m"][0]:
        coefficients[0] = -coefficients[0]
    coefficients[1:] = [-c for c in coefficients[1:]]
    for surface in (primary_surface, secondary_surface, focal_surface):
        surface["coefficients_m"][0] -= overall
    camera_front = focal_surface["coefficients_m"][0]
    focal_surface["coefficients_m"][0] -= _focus_offset(p, zenith_deg)
    focal_surface["outer_radius_m"] = secondary_surface["outer_radius_m"]

    primary = {"asphere_id": 0, "asphere": _vertex_form(primary_surface), "segments": []}
    secondary = {"asphere_id": 1, "asphere": _vertex_form(secondary_surface), "segments": []}
    next_id = 3
    for mirror, name, transform in (
        (primary, "primary_mirror_segmentation", False),
        (secondary, "secondary_mirror_segmentation", True),
    ):
        if name not in p or p.value(name) is None:
            continue
        if p.has_file(name) and p.record(name).get("file"):
            segments = read_segmentation(p.file(name)[0])
        else:
            segments = segmentation_from_record(p.record(name))
        if transform:
            segments = secondary_segment_frame(segments)
        for segment in segments:
            segment["id"] = next_id
            next_id += 1
        mirror["segments"] = segments
    _mirror_response(p, primary, "mirror_reflectivity", "mirror_degraded_reflection")
    _mirror_response(
        p,
        secondary,
        "secondary_mirror_reflectivity",
        "secondary_mirror_degraded_reflection",
        same_as="mirror_reflectivity",
    )
    _scatter(p, primary, surface_slopes=True)
    _scatter(p, secondary, surface_slopes=True)
    for mirror, names, flip in (
        (primary, ("primary_mirror_degraded_map", "primary_degraded_map"), False),
        (secondary, ("secondary_mirror_degraded_map", "secondary_degraded_map"), True),
    ):
        for name in names:
            if p.has_file(name):
                contents, filename, options = _file_args(p, name)
                degradation = read_spatial_map(contents, filename, filename_options=options)
                mirror["degradation"] = _flip_x(degradation) if flip else degradation

    obscurers = _camera_body(p, camera_front, after_primary=True)
    blocks_incoming = True
    shadow = ("secondary_mirror_shadow_diameter", "secondary_mirror_shadow_offset")
    if all(name in p for name in shadow):
        diameter = signed_length_m(p.record(shadow[0]), shadow[0])
        if diameter < 0:
            diameter = 2 * secondary_surface["outer_radius_m"]
        offset = signed_length_m(p.record(shadow[1]), shadow[1])
        if offset <= 0:
            # Concave secondary: shadow plane at its rim; otherwise at its vertex.
            offset = secondary_surface["coefficients_m"][0]
            if secondary_surface["coefficients_m"][1] < 0:
                rho = 0.5 * diameter / secondary_surface["radial_scale_m"]
                offset = sum(
                    c * rho ** (2 * i) for i, c in enumerate(secondary_surface["coefficients_m"])
                )
        else:
            offset -= overall
        if diameter > 0:
            blocks_incoming = False
            obscurers.insert(
                0,
                {
                    "name": "secondary-shadow",
                    "group": "secondary",
                    "legs": ["incoming"],
                    "shape": "plate",
                    "centre_m": [0.0, 0.0, offset],
                    "normal": [0.0, 0.0, 1.0],
                    "tangent": [1.0, 0.0, 0.0],
                    "aperture": "circle",
                    "diameter_m": diameter,
                },
            )
    if "secondary_mirror_baffle" in p:
        baffle = _baffle(p.record("secondary_mirror_baffle"), overall)
        if baffle:
            obscurers.append(baffle)
    return {
        "primary": primary,
        "secondary": secondary,
        "secondary_blocks_incoming": blocks_incoming,
        "focal_surface": {"id": 2, "asphere": _vertex_form(focal_surface)},
        "obscurers": obscurers,
    }


def _flip_x(table: dict[str, Any]) -> dict[str, Any]:
    """Map in the secondary's own frame -> telescope frame (x -> -x)."""
    return {
        **table,
        "x_m": [-v for v in reversed(table["x_m"])],
        "values": list(reversed(table["values"])),
    }


def _vertex_form(surface: dict[str, Any]) -> dict[str, Any]:
    coefficients = list(surface["coefficients_m"])
    vertex = coefficients[0]
    coefficients[0] = 0.0
    while len(coefficients) > 1 and coefficients[-1] == 0.0:
        coefficients.pop()
    return {
        "vertex_z_m": vertex,
        "inner_radius_m": surface["inner_radius_m"],
        "outer_radius_m": surface["outer_radius_m"],
        "radial_scale_m": surface["radial_scale_m"],
        "coefficients_m": coefficients,
    }


def _asphere(p: _Parameters, prefix: str, *, aperture: bool = True) -> dict[str, Any]:
    """Even polynomial z(r) in metres, as in sim_config.c.

    Without a reference radius the coefficients are physical-centimetre powers;
    with one, z = R sum a_i (r/R)^(2i).
    """
    values, units = components(p.record(f"{prefix}_parameters"), f"{prefix}_parameters")
    if any(u != "cm" for u in units):
        raise ModelError(f"{prefix}_parameters must be in centimetres")
    if f"{prefix}_ref_radius" in p:
        scale = length_m(p.record(f"{prefix}_ref_radius"), f"{prefix}_ref_radius")
        coefficients = [number(v, prefix) * scale for v in values]
    else:
        scale = 1.0
        coefficients = [number(v, prefix) * 0.01 ** (1 - 2 * i) for i, v in enumerate(values)]
    surface = {
        "coefficients_m": coefficients,
        "radial_scale_m": scale,
        "inner_radius_m": 0.0,
        "outer_radius_m": 1.0,
    }
    if aperture:
        surface["outer_radius_m"] = 0.5 * length_m(
            p.record(f"{prefix}_diameter"), f"{prefix}_diameter"
        )
        surface["inner_radius_m"] = 0.5 * length_m(
            p.record(f"{prefix}_hole_diameter"), f"{prefix}_hole_diameter", allow_zero=True
        )
        if surface["inner_radius_m"] >= surface["outer_radius_m"]:
            raise ModelError(f"{prefix} hole must be smaller than the mirror")
    return surface


def _baffle(record: dict[str, Any], overall: float) -> dict[str, Any] | None:
    """Secondary baffle (sim_telarray manual 12.3.3): an open cone wall.

    The values are z1, z2, R1, the wall thickness and R2 (zero: R2 = R1). As in
    sim_config.c, a thickness not below R1 is the outer radius at z1.
    """
    values = record.get("value")
    if not isinstance(values, list) or not 3 <= len(values) <= 5:
        raise ModelError("secondary_mirror_baffle requires 3 to 5 lengths")
    z1, z2, r1, thickness, r2 = (
        lengths_m(record, "secondary_mirror_baffle", len(values)) + [0.0, 0.0]
    )[:5]
    if r1 == 0:
        return None
    r2 = r2 or r1
    if thickness >= r1 > 0:
        thickness -= r1
    if r1 < 0 or r2 < 0 or thickness < 0 or z1 == z2:
        raise ModelError("invalid secondary_mirror_baffle")
    return {
        "name": "secondary-baffle",
        "group": "secondary",
        "legs": ["incoming", "after_primary"],
        "shape": "frustum",
        "start_m": [0.0, 0.0, z1 - overall],
        "end_m": [0.0, 0.0, z2 - overall],
        "start_radius_m": r1,
        "end_radius_m": r2,
        "thickness_m": thickness,
        "solid": False,
    }


# --- Telescope structure (beyond sim_telarray): masts, quadrupods, dish frame ---


def _structure(p: _Parameters, first_id: int) -> list[dict[str, Any]]:
    """Opaque structure in the dish frame (z = 0 at the dish vertex).

    One row per element, by its ``shape``: ``cylinder`` (axis end points 1 and
    2, ``diameter``) or ``convex_polyhedron`` (the convex hull of its corner
    points 1, 2, ...). Columns a row does not use are nan.
    """
    if not p.has_file(STRUCTURE_PARAMETER):
        return []
    name = STRUCTURE_PARAMETER
    table = read_structure_table(p.file(name)[0], name)
    if table["frame"] != "primary_mirror_frame":
        raise ModelError(f"{name}: unsupported coordinate_frame {table['frame']}")
    if not {"id", "shape", "diameter", "x1", "x2"} <= set(table["columns"]):
        raise ModelError(f"{name} needs the columns id, shape, diameter and x1, y1, z1, ...")
    corners = sum(1 for column in table["columns"] if column[0] == "x" and column[1:].isdigit())
    overall = (
        signed_length_m(p.record("mirror_offset"), "mirror_offset") if "mirror_offset" in p else 0.0
    )
    obscurers: list[dict[str, Any]] = []
    for row in table["rows"]:
        where = f"{name} {row['id']}"
        points = [[row[f"{a}{i}"] for a in "xyz"] for i in range(1, corners + 1)]
        if any(math.isnan(v) != math.isnan(point[0]) for point in points for v in point):
            raise ModelError(f"{where}: a point has nan and numeric coordinates")
        used = [[x, y, z - overall] for x, y, z in points if not math.isnan(x)]
        if any(math.isnan(point[0]) for point in points[: len(used)]):
            raise ModelError(f"{where}: unused points must come last")
        base = {
            "id": first_id + len(obscurers),
            "name": row["id"],
            "group": row.get("group", "structure"),
            "legs": ALL_LEGS,
        }
        if row["shape"] == "cylinder":
            if len(used) != 2 or not row["diameter"] > 0.0:
                raise ModelError(f"{where}: a cylinder needs two end points and a diameter")
            obscurers.append({
                **base,
                "shape": "cylinder",
                "start_m": used[0],
                "end_m": used[1],
                "diameter_m": row["diameter"],
            })
        elif row["shape"] == "convex_polyhedron":
            if not math.isnan(row["diameter"]):
                raise ModelError(f"{where}: a convex_polyhedron has no diameter")
            obscurers.append({
                **base,
                "shape": "polyhedron",
                "vertices_m": used,
                "planes": convex_hull(used),
            })
        else:
            raise ModelError(f"{where}: unsupported shape {row['shape']}")
    return obscurers


def convex_hull(vertices: list[list[float]]) -> list[list[float]]:
    """Face planes [nx, ny, nz, d] (outward unit normal, n.x <= d inside) of
    the convex hull of a few points."""

    def sub(a, b):
        return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]

    def dot(a, b):
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]

    scale = max(math.dist(a, b) for a in vertices for b in vertices)
    tolerance = 1e-9 * scale
    planes: list[list[float]] = []
    n = len(vertices)
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                a, b, c = vertices[i], vertices[j], vertices[k]
                e, f = sub(b, a), sub(c, a)
                normal = [
                    e[1] * f[2] - e[2] * f[1],
                    e[2] * f[0] - e[0] * f[2],
                    e[0] * f[1] - e[1] * f[0],
                ]
                length = math.sqrt(dot(normal, normal))
                if length <= 1e-12 * scale * scale:
                    continue
                normal = [v / length for v in normal]
                offset = dot(normal, a)
                sides = [dot(normal, v) - offset for v in vertices]
                if all(s <= tolerance for s in sides):
                    plane = [*normal, offset]
                elif all(s >= -tolerance for s in sides):
                    plane = [-normal[0], -normal[1], -normal[2], -offset]
                else:
                    continue
                if not any(
                    dot(plane[:3], q[:3]) > 1 - 1e-12 and abs(plane[3] - q[3]) <= tolerance
                    for q in planes
                ):
                    planes.append(plane)
    if len(planes) < 4:
        raise ModelError("a structure element has no volume (its vertices are coplanar)")
    return planes
