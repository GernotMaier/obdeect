#!/usr/bin/env python3
"""Optional headless plots of compiled telescope geometry and recorded optical diagnostics."""

import argparse
import csv
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Literal
from xml.sax.saxutils import escape

from obdeect.camera_config import PIXEL_APERTURE_SHAPES


class _FileRenderer:
    """Create file-only figures without importing pyplot or a GUI backend."""

    def __init__(self):
        from matplotlib import cm, colormaps
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        from matplotlib.figure import Figure

        self.cm = cm
        self.get_cmap = colormaps.get_cmap
        self._figure = Figure
        self._canvas = FigureCanvasAgg

    def figure(self, **options):
        figure = self._figure(**options)
        self._canvas(figure)
        return figure

    def subplots(self, *args, figsize=None, **options):
        figure = self.figure(figsize=figsize)
        return figure, figure.subplots(*args, **options)

    @staticmethod
    def close(figure):
        figure.clear()


def _save_figure(figure, output: Path, **options):
    """Make vector exports reproducible by removing date and random ID metadata."""
    import matplotlib

    metadata = {}
    if output.suffix.lower() == ".svg":
        metadata["Date"] = None
    elif output.suffix.lower() == ".pdf":
        metadata = {"CreationDate": None, "ModDate": None}
    with matplotlib.rc_context({"svg.hashsalt": "obdeect-plot-v1"}):
        figure.savefig(output, metadata=metadata, **options)


def _parse_path_row(path: Path, row: dict[str, str]) -> list[tuple[float, float, float]]:
    """Validate and return the ragged path vertices from one CSV record."""
    points = []
    try:
        count = int(row["point_count"])
        if not 1 <= count <= 65:
            raise ValueError("point_count must be between 1 and 65")
        for index in range(count):
            point = tuple(float(row[f"{axis}{index}_m"]) for axis in "xyz")
            if not all(math.isfinite(value) for value in point):
                raise ValueError("non-finite path vertex")
            points.append(point)
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"{path}: invalid path vertices: {error}") from error
    return points


def _native_arrivals(path: Path):
    """Return the validated native iterator, or None for legacy diagnostics."""
    with path.open(newline="", encoding="utf-8") as handle:
        fields = csv.DictReader(handle).fieldnames or ()
    if "contract_version" in fields:
        from obdeect.result_contract import iter_arrivals

        return iter_arrivals(path)
    return None


def read_paths(path: Path):
    """Yield terminal status and validated vertices from a trace CSV."""
    arrivals = _native_arrivals(path)
    if arrivals is not None:
        for arrival in arrivals:
            yield arrival.status, list(arrival.interaction_points_m)
        return
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            yield row["status"], _parse_path_row(path, row)


def read_trace_rows(path: Path):
    """Yield validated path records with source and spectral metadata."""
    arrivals = _native_arrivals(path)
    if arrivals is not None:
        # Plot filters need original CSV columns beyond the arrival contract.
        # Exhaust validation without retaining a second copy of every photon.
        for _ in arrivals:
            pass
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            yield row, row["status"], _parse_path_row(path, row)


def focal_plane_hits(path: Path):
    """Yield detected focal-plane ``(x, y, weight)`` samples from a trace CSV.

    The final ragged path vertex is the detector-surface intersection.  Modern
    CSV output supplies source weights and throughput; older files remain
    usable with unit weight.
    """
    arrivals = _native_arrivals(path)
    if arrivals is not None:
        for arrival in arrivals:
            if arrival.detected:
                yield arrival.focal_x_m, arrival.focal_y_m, arrival.optical_weight
        return
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row.get("status") != "detected":
                continue
            try:
                point_index = int(row["point_count"]) - 1
            except (KeyError, TypeError, ValueError) as error:
                raise ValueError(f"{path}: detected row has an invalid point_count") from error
            if point_index < 0:
                raise ValueError(f"{path}: detected row has no path vertices")
            try:
                source_weight = float(row.get("source_weight", "1") or "1")
                throughput = float(row.get("throughput", "1") or "1")
                x = float(row[f"x{point_index}_m"])
                y = float(row[f"y{point_index}_m"])
            except (KeyError, TypeError, ValueError) as error:
                raise ValueError(f"{path}: detected row has invalid focal-plane data") from error
            if not all(math.isfinite(value) for value in (source_weight, throughput, x, y)):
                raise ValueError(f"{path}: detected row has non-finite focal-plane data")
            if source_weight < 0.0 or throughput < 0.0:
                raise ValueError(f"{path}: detected row has negative optical weight")
            yield (
                x,
                y,
                source_weight * throughput,
            )


def compiled_focal_plane_extent(path: Path) -> float:
    """Return the active focal-plane radius carried by a compiled model."""
    try:
        optical_model = json.loads(path.read_text(encoding="utf-8"))
        trace_model = optical_model["trace_model"]
        if trace_model["kind"] == "segmented" or (
            trace_model["kind"] == "axisymmetric" and trace_model.get("detector_surfaces")
        ):
            surfaces = trace_model["detector_surfaces"]
            if not isinstance(surfaces, list) or not surfaces:
                raise ValueError("segmented model requires detector surfaces")
            extents = []
            for index, surface in enumerate(surfaces):
                polygon = _polygon_from_surface(
                    surface,
                    identifier=surface.get("id", index),
                    role="detector",
                    description=f"detector surface {index}",
                )
                if surface["shape"] == "circle":
                    centre = _vector3(surface["centre_m"], f"detector surface {index} centre_m")
                    extents.append(
                        math.hypot(centre[0], centre[1]) + float(surface["diameter_m"]) * 0.5
                    )
                else:
                    extents.extend(
                        math.hypot(vertex[0], vertex[1]) for vertex in polygon.vertices_m
                    )
            extent = max(extents)
        elif trace_model["kind"] == "axisymmetric":
            extent = float(trace_model["detector"]["outer_radius_m"])
        else:
            raise ValueError("unsupported trace-model kind")
    except (KeyError, TypeError, ValueError, json.JSONDecodeError) as error:
        raise ValueError(f"{path}: invalid compiled focal-plane geometry") from error
    if not math.isfinite(extent) or extent <= 0.0:
        raise ValueError(f"{path}: invalid compiled focal-plane extent")
    return extent


TELESCOPE_NAMES = ("reference-mst", "LST", "MST", "SST", "SCT")


@dataclass(frozen=True)
class PlotPolygon:
    """One finite compiled surface rendered as an exact boundary polygon."""

    identifier: int
    role: Literal["primary", "secondary", "detector", "obscurer"]
    vertices_m: tuple[tuple[float, float, float], ...]


@dataclass(frozen=True)
class AxisymmetricSurface:
    """One finite rotationally symmetric compiled surface."""

    role: Literal["primary", "secondary", "detector"]
    vertex_z_m: float
    inner_radius_m: float
    outer_radius_m: float
    radial_scale_m: float
    coefficient_m: tuple[float, ...]


@dataclass(frozen=True)
class PlotObscurer:
    """A finite opaque cylinder supplied by the trace model."""

    identifier: int
    first_endpoint_m: tuple[float, float, float]
    second_endpoint_m: tuple[float, float, float]
    diameter_m: float


@dataclass(frozen=True)
class PlotPixel:
    """Finite camera-layout footprint, used only in the focal-plane panel."""

    identifier: int
    vertices_xy_m: tuple[tuple[float, float], ...]


@dataclass(frozen=True)
class PlotOpticalModel:
    """Validated renderer-facing view of a compiled optical model."""

    kind: Literal["segmented", "axisymmetric"]
    model_label: str
    model_sha256: str | None
    polygons: tuple[PlotPolygon, ...]
    axisymmetric_surfaces: tuple[AxisymmetricSurface, ...]
    obscurers: tuple[PlotObscurer, ...]
    unavailable_roles: tuple[str, ...]
    legacy_geometry: bool = True
    frame_origin: str = "unspecified"
    camera_pixels: tuple[PlotPixel, ...] = ()


def _vector3(value, description: str) -> tuple[float, float, float]:
    """Return a finite three-vector from compiled-model JSON."""
    if not isinstance(value, list) or len(value) != 3:
        raise ValueError(f"{description} must be a three-vector")
    result = tuple(float(component) for component in value)
    if not all(math.isfinite(component) for component in result):
        raise ValueError(f"{description} must be finite")
    return result


def _cross(left: tuple[float, float, float], right: tuple[float, float, float]):
    return (
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    )


def _normalised(vector: tuple[float, float, float], description: str):
    length = math.sqrt(sum(component * component for component in vector))
    if not math.isfinite(length) or length <= 0.0:
        raise ValueError(f"{description} must have non-zero length")
    return tuple(component / length for component in vector)


def _aperture_local_vertices(shape: str, diameter_m: float):
    """Return aperture-boundary points in a facet's tangent plane.

    The model's diameter is the circle diameter, square side length, or regular
    hexagon flat-to-flat distance, respectively.
    """
    if not math.isfinite(diameter_m) or diameter_m <= 0.0:
        raise ValueError("facet diameter_m must be finite and positive")
    if shape == "square":
        half = diameter_m * 0.5
        return ((-half, -half), (half, -half), (half, half), (-half, half))
    if shape in {"hexagon_flat_x", "hexagon_flat_y"}:
        radius = diameter_m / math.sqrt(3.0)
        offset_deg = 30.0 if shape == "hexagon_flat_x" else 0.0
        return tuple(
            (
                radius * math.cos(math.radians(offset_deg + 60.0 * index)),
                radius * math.sin(math.radians(offset_deg + 60.0 * index)),
            )
            for index in range(6)
        )
    if shape == "circle":
        return tuple(
            (
                diameter_m * 0.5 * math.cos(2.0 * math.pi * index / 32.0),
                diameter_m * 0.5 * math.sin(2.0 * math.pi * index / 32.0),
            )
            for index in range(32)
        )
    raise ValueError(f"unsupported facet shape {shape!r}")


def _polygon_from_surface(
    surface: dict,
    *,
    identifier: int,
    role: Literal["primary", "secondary", "detector", "obscurer"],
    description: str,
) -> PlotPolygon:
    """Build one finite aperture polygon from a serialised tangent frame."""
    if (
        not isinstance(identifier, int)
        or isinstance(identifier, bool)
        or not 0 <= identifier < 4294967295
    ):
        raise ValueError(f"{description} has invalid component id")
    centre = _vector3(surface["centre_m"], f"{description} centre_m")
    normal = _normalised(_vector3(surface["normal"], f"{description} normal"), description)
    tangent = _normalised(_vector3(surface["tangent"], f"{description} tangent"), description)
    if abs(sum(left * right for left, right in zip(normal, tangent, strict=True))) > 1.0e-9:
        raise ValueError(f"{description} normal and tangent must be perpendicular")
    binormal = _normalised(_cross(normal, tangent), description)
    diameter = float(surface["diameter_m"])
    vertices = tuple(
        tuple(
            centre[axis] + tangent[axis] * local_x + binormal[axis] * local_y for axis in range(3)
        )
        for local_x, local_y in _aperture_local_vertices(str(surface["shape"]), diameter)
    )
    if not all(math.isfinite(value) for vertex in vertices for value in vertex):
        raise ValueError(f"{description} has non-finite aperture vertices")
    return PlotPolygon(identifier, role, vertices)


def _axisymmetric_surface(
    surface: dict, role: Literal["primary", "secondary", "detector"]
) -> AxisymmetricSurface:
    """Validate one C++ ``EvenPolynomialSurface`` serialisation."""
    try:
        coefficients = tuple(float(value) for value in surface["coefficient_m"])
        vertex = float(surface["vertex_z_m"])
        inner = float(surface["inner_radius_m"])
        outer = float(surface["outer_radius_m"])
        scale = float(surface["radial_scale_m"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"invalid axisymmetric {role} surface") from error
    if (
        len(coefficients) != 13
        or not all(math.isfinite(value) for value in (vertex, inner, outer, scale, *coefficients))
        or inner < 0.0
        or outer <= inner
        or scale <= 0.0
    ):
        raise ValueError(f"invalid axisymmetric {role} surface")
    return AxisymmetricSurface(role, vertex, inner, outer, scale, coefficients)


def _model_label(optical_model: dict) -> str:
    """Return a compact provenance label without requiring CTAO-specific fields."""
    provenance = optical_model.get("provenance", {})
    if not isinstance(provenance, dict):
        return "compiled optical model"
    model = provenance.get("model")
    version = provenance.get("model_version")
    if isinstance(model, str) and isinstance(version, str):
        return f"{model}/{version}"
    if isinstance(model, str):
        return model
    return "compiled optical model"


def _surface_height_and_slope(surface: AxisymmetricSurface, radius: float):
    """Compiled even-polynomial sag and derivative at a radius."""
    q = (radius / surface.radial_scale_m) ** 2
    height, derivative_q = 0.0, 0.0
    for coefficient in reversed(surface.coefficient_m):
        derivative_q = derivative_q * q + height
        height = height * q + coefficient
    return surface.vertex_z_m + height, derivative_q * 2 * radius / surface.radial_scale_m**2


def _aspheric_mask_polygon(mask: dict, surface: AxisymmetricSurface) -> PlotPolygon:
    """Boundary from a finite mask; hex outlines use the centre tangent plane.

    Annular-sector edges sample the actual radial gap boundary at 128 points.
    Hex masks use the tracer's slope-flattened centre frame; curvature remains
    explicitly unrendered between these planar footprint vertices.
    """
    identifier = int(mask["id"])
    if mask["shape"] == "hexagon":
        x, y = (float(value) for value in mask["centre_xy_m"])
        diameter, rotation = float(mask["diameter_m"]), math.radians(float(mask["rotation_deg"]))
        if not all(math.isfinite(value) for value in (x, y, diameter, rotation)):
            raise ValueError("non-finite hex mask")
        radius = math.hypot(x, y)
        z, slope = _surface_height_and_slope(surface, radius)
        cx, cy = (x / radius, y / radius) if radius else (1.0, 0.0)
        cosine = 1 / math.hypot(1, slope)
        sine = slope * cosine
        vertices = []
        for u, v in _aperture_local_vertices("hexagon_flat_x", diameter):
            local_x = u * math.cos(rotation) - v * math.sin(rotation)
            local_y = u * math.sin(rotation) + v * math.cos(rotation)
            radial, azimuth = local_x * cx + local_y * cy, -local_x * cy + local_y * cx
            vertices.append((
                x + radial * cosine * cx - azimuth * cy,
                y + radial * cosine * cy + azimuth * cx,
                z + radial * sine,
            ))
        return PlotPolygon(identifier, surface.role, tuple(vertices))
    if mask["shape"] != "annular_sector":
        raise ValueError("unsupported finite aspheric mask")
    inner, outer = float(mask["inner_radius_m"]), float(mask["outer_radius_m"])
    start, span = math.radians(float(mask["start_deg"])), math.radians(float(mask["span_deg"]))
    gap = float(mask["gap_m"])
    if (
        not all(math.isfinite(value) for value in (inner, outer, start, span, gap))
        or inner < 0
        or outer <= inner
        or not 0 < span <= 2 * math.pi
        or gap < 0
    ):
        raise ValueError("invalid annular-sector mask")
    inner = max(inner, gap, gap / span, surface.inner_radius_m)
    outer = min(outer, surface.outer_radius_m)
    if inner >= outer:
        raise ValueError("annular-sector gap removes the complete footprint")

    def end(radius):
        return start + span - (gap / radius if radius else 0)

    samples = 128
    outline = [(inner + (outer - inner) * i / (samples - 1), start) for i in range(samples)]
    outline.extend(
        (outer, start + (end(outer) - start) * i / (samples - 1)) for i in range(samples)
    )
    outline.extend(
        (r, end(r)) for r in (outer - (outer - inner) * i / (samples - 1) for i in range(samples))
    )
    outline.extend(
        (inner, end(inner) - (end(inner) - start) * i / (samples - 1)) for i in range(samples)
    )
    vertices = tuple(
        (
            radius * math.cos(angle),
            radius * math.sin(angle),
            _surface_height_and_slope(surface, radius)[0],
        )
        for radius, angle in outline
    )
    return PlotPolygon(identifier, surface.role, vertices)


def _camera_pixel_polygons(optical_model: dict) -> tuple[PlotPixel, ...]:
    """Use finite parsed camera entrance shapes; never infer footprints from centres.

    Shape codes follow sim_telarray/common/mc_aux.h Pix_Type and are boundary
    import conventions, not telescope-dependent drawing rules.
    """
    camera = optical_model.get("camera", {})
    if camera.get("entrance_surfaces"):
        return tuple(
            PlotPixel(
                plane["source_pixel_id"],
                tuple(
                    (point[0], point[1])
                    for point in _polygon_from_surface(
                        plane,
                        identifier=plane["source_pixel_id"],
                        role="detector",
                        description="compiled pixel entrance",
                    ).vertices_m
                ),
            )
            for plane in camera["entrance_surfaces"]
        )
    types = {item["id"]: item for item in camera.get("pixel_types", [])}
    if not types:
        return ()
    angle = math.radians(float(camera.get("rotation_deg", 0)))
    if not math.isfinite(angle):
        raise ValueError("camera rotation must be finite")
    cosine, sine = math.cos(angle), math.sin(angle)
    result = []
    for pixel in camera.get("pixels", []):
        pixel_type = types[pixel["type_id"]]
        shape = PIXEL_APERTURE_SHAPES[pixel_type["funnel_shape_code"]]
        diameter = float(pixel_type["funnel_diameter_m"])
        x, y = (float(value) for value in pixel["centre_xy_m"])
        if not math.isfinite(x) or not math.isfinite(y):
            raise ValueError("camera pixel centre must be finite")
        vertices = tuple(
            ((x + dx) * cosine - (y + dy) * sine, (x + dx) * sine + (y + dy) * cosine)
            for dx, dy in _aperture_local_vertices(shape, diameter)
        )
        result.append(PlotPixel(int(pixel["id"]), vertices))
    if len({item.identifier for item in result}) != len(result):
        raise ValueError("camera pixel IDs must be unique")
    return tuple(result)


def _plot_obscurer(obscurer: dict, index: int) -> PlotObscurer:
    """Validate one finite opaque cylinder for plotting."""
    try:
        identifier = int(obscurer["id"])
        first = _vector3(obscurer["first_endpoint_m"], f"obscurer {index} first endpoint")
        second = _vector3(obscurer["second_endpoint_m"], f"obscurer {index} second endpoint")
        diameter = float(obscurer["diameter_m"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"invalid obscurer {index}") from error
    if (
        not 0 <= identifier < 4294967295
        or not math.isfinite(diameter)
        or diameter <= 0.0
        or first == second
    ):
        raise ValueError(f"invalid obscurer {index}")
    return PlotObscurer(identifier, first, second, diameter)


def load_plot_optical_model(path: Path) -> PlotOpticalModel:
    """Load only finite geometry explicitly represented by a compiled model.

    The renderer intentionally uses the same finite surfaces and cylinders as
    the transport model.  It neither consults a telescope name nor constructs
    customary camera/support hardware that is absent from the JSON.
    """
    try:
        optical_model = json.loads(path.read_text(encoding="utf-8"))
        trace_model = optical_model["trace_model"]
        kind = trace_model["kind"]
    except (KeyError, TypeError, OSError, json.JSONDecodeError) as error:
        raise ValueError(f"{path}: invalid compiled optical model") from error
    if (
        optical_model.get("format", "obdeect.compiled-optical-model.v1")
        != "obdeect.compiled-optical-model.v1"
    ):
        raise ValueError(f"{path}: unsupported compiled optical-model format")
    if kind not in {"segmented", "axisymmetric"}:
        raise ValueError(f"{path}: unsupported trace-model kind {kind!r}")
    model_sha = optical_model.get("optical_model_sha256")
    if model_sha is not None and (
        not isinstance(model_sha, str)
        or len(model_sha) != 64
        or any(character not in "0123456789abcdef" for character in model_sha)
    ):
        raise ValueError(f"{path}: invalid optical_model_sha256")
    declared = optical_model.get("plot_geometry", {})
    legacy = "plot_geometry" not in optical_model
    origin = "unspecified"
    if not legacy:
        if not isinstance(declared, dict) or declared.get("schema_version") != 1:
            raise ValueError(f"{path}: unsupported plot_geometry schema_version")
        frame = declared.get("frame", {})
        if (
            not isinstance(frame, dict)
            or frame.get("unit") != "m"
            or frame.get("axes") != "+x,+y,+z"
            or not isinstance(frame.get("origin"), str)
        ):
            raise ValueError(f"{path}: invalid plot_geometry frame")
        origin = frame["origin"]
        components = declared.get("components")
        if not isinstance(components, list):
            raise ValueError(f"{path}: invalid plot_geometry components")
        for component in components:
            if not isinstance(component, dict):
                raise ValueError(f"{path}: invalid plot_geometry component")
            role = component.get("role")
            expected = {
                "primary": "primary_facets" if kind == "segmented" else "primary",
                "secondary": "secondary",
                "detector": "detector_surfaces"
                if kind == "segmented" or trace_model.get("detector_surfaces")
                else "detector",
                "opaque_cylinder": "cylinder_obscurers",
            }.get(role)
            allowed_sources = (
                {"trace_model.primary_to_secondary_planes", "trace_model.incoming_obscurer_planes"}
                if role == "obscurer"
                else {"trace_model.primary_to_secondary_cylinders"}
                if role == "opaque_cylinder" and kind == "axisymmetric"
                else {f"trace_model.{expected}"}
                if expected is not None
                else set()
            )
            if component.get("source") not in allowed_sources:
                raise ValueError(f"{path}: invalid plot_geometry component source")
            source = trace_model.get(component["source"].removeprefix("trace_model."))
            identifier = component.get("id")
            if (
                not isinstance(identifier, int)
                or isinstance(identifier, bool)
                or not 0 <= identifier < 4294967295
            ):
                raise ValueError(f"{path}: invalid plot_geometry component id")
            if isinstance(source, list):
                if not any(
                    isinstance(item, dict) and item.get("id", index) == identifier
                    for index, item in enumerate(source)
                ):
                    raise ValueError(f"{path}: plot_geometry references absent component")
            elif not isinstance(source, dict):
                raise ValueError(f"{path}: plot_geometry references absent geometry")
    unavailable: tuple[str, ...] = ()
    if isinstance(declared, dict) and "unavailable_roles" in declared:
        roles = declared["unavailable_roles"]
        if not isinstance(roles, list) or not all(isinstance(role, str) for role in roles):
            raise ValueError(f"{path}: invalid plot_geometry unavailable_roles")
        unavailable = tuple(sorted(set(roles)))

    try:
        camera_pixels = _camera_pixel_polygons(optical_model)
        if kind == "segmented":
            facets = trace_model["primary_facets"]
            detectors = trace_model.get("detector_surfaces", [])
            obscurers = trace_model.get("cylinder_obscurers", [])
            if (
                not isinstance(facets, list)
                or not facets
                or not isinstance(detectors, list)
                or not isinstance(obscurers, list)
            ):
                raise ValueError("segmented model requires primary facets")
            polygons = tuple(
                _polygon_from_surface(
                    facet,
                    identifier=facet.get("id", index),
                    role="primary",
                    description=f"primary facet {index}",
                )
                for index, facet in enumerate(facets)
            ) + tuple(
                _polygon_from_surface(
                    detector,
                    identifier=detector.get("id", len(facets) + index),
                    role="detector",
                    description=f"detector surface {index}",
                )
                for index, detector in enumerate(detectors)
            )
            polygons += tuple(
                _polygon_from_surface(
                    plane,
                    identifier=plane["id"],
                    role="obscurer",
                    description="incoming camera housing",
                )
                for plane in trace_model.get("incoming_obscurer_planes", [])
            )
            obscurer_rows = tuple(
                _plot_obscurer(obscurer, index) for index, obscurer in enumerate(obscurers)
            )
            identifiers = [polygon.identifier for polygon in polygons]
            identifiers.extend(obscurer.identifier for obscurer in obscurer_rows)
            if len(identifiers) != len(set(identifiers)):
                raise ValueError("compiled component IDs must be unique")
            return PlotOpticalModel(
                "segmented",
                _model_label(optical_model),
                model_sha,
                polygons,
                (),
                tuple(obscurer_rows),
                unavailable,
                legacy,
                origin,
                camera_pixels,
            )
        obscurer_planes = tuple(
            _polygon_from_surface(
                plane,
                identifier=plane["id"],
                role="obscurer",
                description="compiled obscurer plane",
            )
            for field in ("primary_to_secondary_planes", "incoming_obscurer_planes")
            for plane in trace_model.get(field, [])
        )
        surfaces = tuple(
            _axisymmetric_surface(trace_model[source], role)
            for source, role in (
                ("primary", "primary"),
                ("secondary", "secondary"),
                ("detector", "detector"),
            )
            if role != "detector" or not trace_model.get("detector_surfaces")
        )
        mask_polygons = tuple(
            _aspheric_mask_polygon(mask, surface)
            for surface in surfaces
            if surface.role != "detector"
            for mask in trace_model.get(f"{surface.role}_segments", [])
        ) + tuple(
            _polygon_from_surface(
                plane,
                identifier=plane["id"],
                role="detector",
                description="compiled entrance",
            )
            for plane in trace_model.get("detector_surfaces", [])
        )
        mask_polygons += obscurer_planes
        obscurer_cylinders = tuple(
            _plot_obscurer(cylinder, index)
            for index, cylinder in enumerate(trace_model.get("primary_to_secondary_cylinders", []))
        )
        mask_ids = [polygon.identifier for polygon in mask_polygons]
        mask_ids.extend(obscurer.identifier for obscurer in obscurer_cylinders)
        if len(mask_ids) != len(set(mask_ids)):
            raise ValueError("finite mask component IDs must be unique")
        return PlotOpticalModel(
            "axisymmetric",
            _model_label(optical_model),
            model_sha,
            mask_polygons,
            surfaces,
            obscurer_cylinders,
            unavailable,
            legacy,
            origin,
            camera_pixels,
        )
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"{path}: invalid compiled {kind} geometry: {error}") from error


def _axisymmetric_profile(surface: AxisymmetricSurface, samples: int = 256):
    """Return the exact even-polynomial section profile in telescope metres."""
    if samples < 2:
        raise ValueError("axisymmetric profile needs at least two samples")
    radii = tuple(
        surface.inner_radius_m
        + (surface.outer_radius_m - surface.inner_radius_m) * index / (samples - 1)
        for index in range(samples)
    )

    def sag(radius: float) -> float:
        radius_squared = (radius / surface.radial_scale_m) ** 2
        result = 0.0
        for coefficient in reversed(surface.coefficient_m):
            result = result * radius_squared + coefficient
        return surface.vertex_z_m + result

    result = tuple((radius, sag(radius)) for radius in radii)
    if not all(math.isfinite(height) for _, height in result):
        raise ValueError("axisymmetric profile has non-finite sampled height")
    return result


def compiled_facet_polygons(path: Path):
    """Return exact primary-aperture polygons from a segmented trace model.

    Each result is ``(vertices, centre_z)`` where vertices are global 3-D
    positions.  The function deliberately uses the serialised trace geometry,
    so the face-on and 3-D views show what the tracer actually consumes.
    """
    try:
        trace_model = json.loads(path.read_text(encoding="utf-8"))["trace_model"]
        if trace_model["kind"] != "segmented":
            raise ValueError("mirror-face views require a segmented trace model")
        optical_model = load_plot_optical_model(path)
        return [
            (
                polygon.vertices_m,
                sum(vertex[2] for vertex in polygon.vertices_m) / len(polygon.vertices_m),
            )
            for polygon in optical_model.polygons
            if polygon.role == "primary"
        ]
    except (KeyError, TypeError, ValueError, json.JSONDecodeError) as error:
        raise ValueError(f"{path}: invalid compiled primary-facet geometry: {error}") from error


_ROLE_STYLE = {
    "primary": {"color": "#607d8b", "label": "M1"},
    "secondary": {"color": "#455a64", "label": "M2"},
    "detector": {"color": "#e69f00", "label": "detector"},
    "obscurer": {"color": "#a14b3b", "label": "obscurer"},
}


def _optical_model_bounds(optical_model: PlotOpticalModel):
    """Return physical 3-D bounds covering every finite optical-model primitive."""
    values = [[], [], []]
    for polygon in optical_model.polygons:
        for vertex in polygon.vertices_m:
            for coordinate, value in enumerate(vertex):
                values[coordinate].append(value)
    for surface in optical_model.axisymmetric_surfaces:
        if any(p.role == surface.role for p in optical_model.polygons):
            continue
        profile = _axisymmetric_profile(surface)
        values[0].extend((-surface.outer_radius_m, surface.outer_radius_m))
        values[1].extend((-surface.outer_radius_m, surface.outer_radius_m))
        values[2].extend(height for _, height in profile)
    for obscurer in optical_model.obscurers:
        radius = obscurer.diameter_m * 0.5
        for endpoint in (obscurer.first_endpoint_m, obscurer.second_endpoint_m):
            for coordinate, value in enumerate(endpoint):
                values[coordinate].extend((value - radius, value + radius))
    if not all(values):
        raise ValueError("compiled model has no finite renderable geometry")
    result = []
    for coordinate_values in values:
        lower, upper = min(coordinate_values), max(coordinate_values)
        result.append((lower, upper))
    return tuple(result)


def _set_equal_3d_bounds(axis, optical_model: PlotOpticalModel) -> None:
    """Apply equal physical scales with a small readable margin."""
    bounds = _optical_model_bounds(optical_model)
    centres = tuple((lower + upper) * 0.5 for lower, upper in bounds)
    span = max(upper - lower for lower, upper in bounds)
    span = max(span, 1.0) * 1.05
    for coordinate, centre in enumerate(centres):
        getattr(axis, f"set_{'xyz'[coordinate]}lim")(centre - span * 0.5, centre + span * 0.5)
    axis.set_box_aspect((1.0, 1.0, 1.0))


def _cylinder_mesh(obscurer: PlotObscurer, sides: int = 12):
    """Return finite cylinder side faces from exactly the trace-model endpoints."""
    first, second = obscurer.first_endpoint_m, obscurer.second_endpoint_m
    direction = _normalised(
        tuple(right - left for left, right in zip(first, second, strict=True)), "cylinder"
    )
    reference = (0.0, 0.0, 1.0) if abs(direction[2]) < 0.9 else (0.0, 1.0, 0.0)
    first_basis = _normalised(_cross(direction, reference), "cylinder basis")
    second_basis = _cross(direction, first_basis)
    radius = obscurer.diameter_m * 0.5

    def ring(centre):
        return tuple(
            tuple(
                centre[axis]
                + radius
                * (
                    first_basis[axis] * math.cos(2.0 * math.pi * index / sides)
                    + second_basis[axis] * math.sin(2.0 * math.pi * index / sides)
                )
                for axis in range(3)
            )
            for index in range(sides)
        )

    first_ring, second_ring = ring(first), ring(second)
    faces = [first_ring, second_ring]
    faces.extend(
        (
            first_ring[index],
            first_ring[(index + 1) % sides],
            second_ring[(index + 1) % sides],
            second_ring[index],
        )
        for index in range(sides)
    )
    return faces


def _axisymmetric_mesh(surface: AxisymmetricSurface, sides: int = 32):
    """Tessellate only the already-defined polynomial surface for display."""
    profile = _axisymmetric_profile(surface)
    rings = [
        tuple(
            (
                radius * math.cos(2.0 * math.pi * side / sides),
                radius * math.sin(2.0 * math.pi * side / sides),
                height,
            )
            for side in range(sides)
        )
        for radius, height in profile
    ]
    return [
        (
            rings[row][side],
            rings[row][(side + 1) % sides],
            rings[row + 1][(side + 1) % sides],
            rings[row + 1][side],
        )
        for row in range(len(rings) - 1)
        for side in range(sides)
    ]


def _draw_optical_model_assembly(axis, optical_model: PlotOpticalModel) -> None:
    """Draw a restrained orthographic 3-D optical model without synthetic hardware."""
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    for role in ("primary", "secondary", "detector", "obscurer"):
        polygons = [
            polygon.vertices_m for polygon in optical_model.polygons if polygon.role == role
        ]
        if polygons:
            style = _ROLE_STYLE[role]
            axis.add_collection3d(
                Poly3DCollection(
                    polygons,
                    facecolors=style["color"],
                    edgecolors="#263238",
                    linewidths=0.22,
                    alpha=0.9 if role == "primary" else 0.78,
                    label=f"{style['label']} ({len(polygons)})",
                )
            )
    for surface in optical_model.axisymmetric_surfaces:
        if any(p.role == surface.role for p in optical_model.polygons):
            continue
        style = _ROLE_STYLE[surface.role]
        axis.add_collection3d(
            Poly3DCollection(
                _axisymmetric_mesh(surface),
                facecolors=style["color"],
                edgecolors="#263238",
                linewidths=0.08,
                alpha=0.82,
                label=style["label"],
            )
        )
    for index, obscurer in enumerate(optical_model.obscurers):
        axis.add_collection3d(
            Poly3DCollection(
                _cylinder_mesh(obscurer),
                facecolors=_ROLE_STYLE["obscurer"]["color"],
                edgecolors="#57261e",
                linewidths=0.15,
                alpha=0.62,
                label=f"opaque cylinders ({len(optical_model.obscurers)})" if index == 0 else None,
            )
        )
    axis.legend(fontsize=8, loc="upper left")
    _set_equal_3d_bounds(axis, optical_model)
    axis.set_proj_type("ortho")
    axis.view_init(elev=24, azim=-55)
    axis.set(xlabel="x [m]", ylabel="y [m]", zlabel="z [m]")


def _convex_hull(points):
    """Deterministic finite silhouette in a projection."""
    points = sorted(set(points))
    if len(points) <= 2:
        return points

    def turn(a, b, c):
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])

    lower, upper = [], []
    for target, sequence in ((lower, points), (upper, reversed(points))):
        for point in sequence:
            while len(target) >= 2 and turn(target[-2], target[-1], point) <= 0:
                target.pop()
            target.append(point)
    return lower[:-1] + upper[:-1]


def _cylinder_footprint(obscurer: PlotObscurer, horizontal: int, vertical: int):
    """Project a capped finite cylinder using 64 documented angular samples."""
    return _convex_hull([
        (vertex[horizontal], vertex[vertical])
        for face in _cylinder_mesh(obscurer, sides=64)
        for vertex in face
    ])


def _draw_optical_model_section(axis, optical_model: PlotOpticalModel, coordinate: int = 0) -> None:
    """Draw finite surface boundaries in one telescope-frame axial projection."""
    horizontal = "xy"[coordinate]
    for role in ("primary", "secondary", "detector", "obscurer"):
        style = _ROLE_STYLE[role]
        polygons = [item for item in optical_model.polygons if item.role == role]
        for index, polygon in enumerate(polygons):
            points = [(vertex[coordinate], vertex[2]) for vertex in polygon.vertices_m]
            points.append(points[0])
            axis.plot(
                *zip(*points, strict=True),
                color=style["color"],
                linewidth=1.0,
                label=f"{style['label']} ({len(polygons)})" if index == 0 else None,
            )
    for surface in optical_model.axisymmetric_surfaces:
        if any(p.role == surface.role for p in optical_model.polygons):
            continue
        profile = _axisymmetric_profile(surface)
        positive = [(radius, height) for radius, height in profile]
        negative = [(-radius, height) for radius, height in reversed(profile)]
        for branch in (negative, positive):
            axis.plot(
                *zip(*branch, strict=True),
                color=_ROLE_STYLE[surface.role]["color"],
                linewidth=1.4,
                label=_ROLE_STYLE[surface.role]["label"] if branch is positive else None,
            )
    from matplotlib.patches import Polygon

    for obscurer in optical_model.obscurers:
        # Only cylinders reaching the section plane are included. The polygon
        # is their finite orthographic footprint, not an inferred ray shadow.
        other = 1 - coordinate
        lower = min(obscurer.first_endpoint_m[other], obscurer.second_endpoint_m[other])
        upper = max(obscurer.first_endpoint_m[other], obscurer.second_endpoint_m[other])
        direction = _normalised(
            tuple(
                b - a
                for a, b in zip(obscurer.first_endpoint_m, obscurer.second_endpoint_m, strict=True)
            ),
            "cylinder",
        )
        radial_extent = obscurer.diameter_m / 2 * math.sqrt(max(0, 1 - direction[other] ** 2))
        if lower > radial_extent or upper < -radial_extent:
            continue
        footprint = _cylinder_footprint(obscurer, coordinate, 2)
        axis.add_patch(
            Polygon(
                footprint,
                color=_ROLE_STYLE["obscurer"]["color"],
                alpha=0.6,
                label=f"opaque cylinder #{obscurer.identifier}",
            )
        )
    if axis.get_legend_handles_labels()[0]:
        axis.legend(fontsize=8, loc="upper right")
    axis.axvline(0.0, color="#546e7a", linewidth=0.7, linestyle="--", zorder=-1)
    axis.set(
        xlabel=f"telescope {horizontal} [m]",
        ylabel="telescope z [m]",
        title="B. axial projection (in-plane cylinders)",
    )
    axis.set_aspect("equal", adjustable="box")
    from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

    span = axis.get_xlim()[1] - axis.get_xlim()[0]
    decade = 10 ** math.floor(math.log10(span / 5))
    length = max(value * decade for value in (1, 2, 5) if value * decade <= span / 5)
    axis.add_artist(
        AnchoredSizeBar(
            axis.transData,
            length,
            f"{length:g} m",
            "lower right",
            pad=0.4,
            color="#263238",
            frameon=False,
        )
    )


def _draw_optical_model_pupil(
    axis, optical_model: PlotOpticalModel, *, title: str = "C. entrance pupil", colour_by="status"
) -> None:
    """Draw exact finite aperture boundaries looking along the optical axis."""
    from matplotlib.collections import PolyCollection
    from matplotlib.patches import Circle

    for role in ("primary", "secondary", "detector", "obscurer"):
        polygons = [
            [(vertex[0], vertex[1]) for vertex in polygon.vertices_m]
            for polygon in optical_model.polygons
            if polygon.role == role
        ]
        if polygons:
            style = _ROLE_STYLE[role]
            collection = PolyCollection(
                polygons,
                facecolors=style["color"],
                edgecolors="#263238",
                linewidths=0.3,
                alpha=0.82 if role == "primary" else 0.72,
                label=f"{style['label']} ({len(polygons)})",
            )
            if colour_by == "facet-z" and role == "primary":
                import numpy as np

                heights = [
                    sum(v[2] for v in p.vertices_m) / len(p.vertices_m)
                    for p in optical_model.polygons
                    if p.role == "primary"
                ]
                collection.set_array(np.asarray(heights))
                collection.set_cmap("viridis")
                axis.figure.colorbar(collection, ax=axis, label="primary aperture-centre z [m]")
            axis.add_collection(collection)
    for surface in optical_model.axisymmetric_surfaces:
        if any(p.role == surface.role for p in optical_model.polygons):
            continue
        style = _ROLE_STYLE[surface.role]
        axis.add_patch(
            Circle(
                (0.0, 0.0),
                surface.outer_radius_m,
                fill=False,
                color=style["color"],
                linewidth=1.2,
                label=style["label"],
            )
        )
        if surface.inner_radius_m > 0.0:
            axis.add_patch(
                Circle(
                    (0.0, 0.0),
                    surface.inner_radius_m,
                    fill=False,
                    color=style["color"],
                    linewidth=1.0,
                )
            )
    from matplotlib.patches import Polygon

    for obscurer in optical_model.obscurers:
        axis.add_patch(
            Polygon(
                _cylinder_footprint(obscurer, 0, 1),
                color=_ROLE_STYLE["obscurer"]["color"],
                alpha=0.6,
            )
        )
    axis.autoscale_view()
    axis.legend(fontsize=8, loc="best")
    axis.set(xlabel="telescope x [m]", ylabel="telescope y [m]", title=title)
    axis.set_aspect("equal", adjustable="box")


def _draw_focal_geometry(axis, optical_model: PlotOpticalModel) -> None:
    """Draw only compiled detector geometry in the telescope x/y plane."""
    from matplotlib.collections import PolyCollection
    from matplotlib.patches import Circle

    polygons = [
        [(vertex[0], vertex[1]) for vertex in polygon.vertices_m]
        for polygon in optical_model.polygons
        if polygon.role == "detector"
    ]
    if polygons:
        axis.add_collection(
            PolyCollection(
                polygons, facecolors="#e69f00", edgecolors="#5f4600", linewidths=0.4, alpha=0.7
            )
        )
        axis.autoscale_view()
    for surface in optical_model.axisymmetric_surfaces:
        if surface.role == "detector":
            axis.add_patch(Circle((0.0, 0.0), surface.outer_radius_m, color="#e69f00", alpha=0.7))
            if surface.inner_radius_m > 0.0:
                axis.add_patch(Circle((0.0, 0.0), surface.inner_radius_m, color="white"))
    if optical_model.camera_pixels:
        axis.add_collection(
            PolyCollection(
                [pixel.vertices_xy_m for pixel in optical_model.camera_pixels],
                facecolors="none",
                edgecolors="#263238",
                linewidths=0.35,
            )
        )
        axis.autoscale_view()
    if not polygons and not any(s.role == "detector" for s in optical_model.axisymmetric_surfaces):
        axis.text(
            0.5,
            0.5,
            "Finite detector geometry unavailable",
            transform=axis.transAxes,
            ha="center",
            va="center",
            fontsize=11,
        )
    axis.autoscale_view()
    axis.set(xlabel="focal-plane x [m]", ylabel="focal-plane y [m]", title="D. focal geometry")
    axis.set_aspect("equal", adjustable="box")


def _optical_model_footer(optical_model: PlotOpticalModel) -> str:
    """Describe rendered and known-unavailable geometry in one compact line."""
    primary_count = sum(polygon.role == "primary" for polygon in optical_model.polygons)
    detector_count = sum(polygon.role == "detector" for polygon in optical_model.polygons)
    rendered = []
    if primary_count:
        rendered.append(f"{primary_count} M1 finite panel footprints")
    secondary_count = sum(p.role == "secondary" for p in optical_model.polygons)
    if secondary_count:
        rendered.append(f"{secondary_count} M2 finite panel footprints")
    if optical_model.axisymmetric_surfaces:
        rendered.extend(
            _ROLE_STYLE[surface.role]["label"] for surface in optical_model.axisymmetric_surfaces
        )
    if detector_count:
        rendered.append(f"{detector_count} detector surfaces")
    if optical_model.camera_pixels:
        rendered.append(
            f"{len(optical_model.camera_pixels)} 2D camera footprints (focal panel only)"
        )
    if optical_model.obscurers:
        rendered.append(f"{len(optical_model.obscurers)} opaque cylinders")
    coverage = (
        "not modelled: " + ", ".join(optical_model.unavailable_roles)
        if optical_model.unavailable_roles
        else "mechanical geometry coverage not declared"
    )
    sha = optical_model.model_sha256[:12] if optical_model.model_sha256 else "unavailable"
    return (
        f"Telescope frame: +z optical axis; origin {optical_model.frame_origin}; metres | "
        f"SHA256: {sha} | "
        f"rendered: {', '.join(rendered)}\n{coverage}"
        + (" | plot_geometry unavailable (legacy model)" if optical_model.legacy_geometry else "")
        + (
            " | planar aperture boundaries; curvature not rendered"
            if optical_model.polygons
            else " | 256 radial samples"
        )
    )


def draw_telescope(
    plt,
    optical_model_path: Path,
    output: Path,
    section_plane: str = "xz",
    *,
    input_path: Path | None = None,
    show_paths: bool = True,
    panel: str | None = None,
    bins: int = 64,
    max_paths: int = 300,
    selection="first",
    seed=0,
    statuses=(),
    component_id=None,
    context_radius_m: float | None = None,
    colour_by="status",
) -> None:
    """Render the four-panel plate, or an explicitly selected geometry panel."""
    optical_model = load_plot_optical_model(optical_model_path)
    coordinate = {"xz": 0, "yz": 1}[section_plane]
    full_plate = panel in (None, "focal-plane-hits")
    if not full_plate and panel not in {"assembly", "section", "pupil", "focal-geometry"}:
        raise ValueError("unsupported telescope panel")
    figure = plt.figure(figsize=(14, 11) if full_plate else (10, 8))
    names = ("assembly", "section", "pupil", "focal-geometry") if full_plate else (panel,)
    axes = {}
    for index, name in enumerate(names):
        axis = (
            figure.add_subplot(2, 2, index + 1, projection="3d" if name == "assembly" else None)
            if full_plate
            else figure.add_subplot(projection="3d" if name == "assembly" else None)
        )
        axes[name] = axis
        if name == "assembly":
            _draw_optical_model_assembly(axis, optical_model)
        elif name == "section":
            _draw_optical_model_section(axis, optical_model, coordinate)
        elif name == "pupil":
            _draw_optical_model_pupil(axis, optical_model)
        else:
            _draw_focal_geometry(axis, optical_model)
    if panel == "focal-plane-hits":
        if input_path is None:
            raise ValueError("focal-plane-hits requires trace input")
        samples = list(focal_plane_hits(input_path))
        if not samples:
            raise ValueError("no detected focal-plane hits")
        x, y, weights = zip(*samples, strict=True)
        extent = compiled_focal_plane_extent(optical_model_path)
        focal = axes["focal-geometry"]
        focal.clear()
        density = focal.hist2d(
            x,
            y,
            bins=bins,
            weights=weights,
            cmap="viridis",
            range=[[-extent, extent], [-extent, extent]],
        )
        figure.colorbar(density[3], ax=focal, label="weighted detected photons / bin")
        focal.set(
            xlabel="focal-plane x [m]", ylabel="focal-plane y [m]", title="D. focal-plane hits"
        )
        focal.set_aspect("equal")
    if input_path is not None and show_paths:
        records = select_trace_rows(
            input_path,
            max_paths,
            selection=selection,
            seed=seed,
            statuses=statuses,
            component_id=component_id,
        )
        context_radius_m = (
            context_radius_m
            or max(abs(value) for pair in _optical_model_bounds(optical_model) for value in pair)
            * 1.2
        )
        targets = {
            name: axis
            for name, axis in axes.items()
            if not full_plate or name in {"assembly", "section"}
        }
        _overlay_paths(plt, targets, records, coordinate, context_radius_m, colour_by)
    figure.suptitle(
        f"{optical_model.model_label} compiled telescope geometry [{optical_model.kind}]"
    )
    figure.subplots_adjust(bottom=0.11, top=0.93, wspace=0.3, hspace=0.26)
    figure.text(0.5, 0.025, _optical_model_footer(optical_model), ha="center", fontsize=8)
    _save_figure(figure, output, dpi=180, bbox_inches="tight")
    plt.close(figure)


def _overlay_paths(plt, axes, records, coordinate, radius, colour_by):
    """Recorded flights in 3D/axial projections; entrance or focal markers in 2D."""
    colour_map, normaliser, colour_label = _trace_colour_scale(plt, records, colour_by)
    if colour_map is not None:
        axis = next(iter(axes.values()))
        axis.figure.colorbar(
            plt.cm.ScalarMappable(norm=normaliser, cmap=colour_map),
            ax=axis,
            label=colour_label,
            shrink=0.65,
        )
    field = {
        "wavelength": "wavelength_nm",
        "arrival-time": "arrival_time_ns",
        "incidence-angle": "incidence_focal_deg",
    }.get(colour_by)
    seen = {name: set() for name in axes}
    clipped = False
    for row, status, points in records:
        colour = (
            colour_map(normaliser(float(row[field])))
            if field
            else _STATUS_COLOURS.get(status, "black")
        )
        segments = clip_recorded_path(points, radius)
        clipped |= any(abs(value) > radius for point in points for value in point)
        for name, axis in axes.items():
            axial = name.startswith("section")
            section_coordinate = 1 if name == "section-yz" else coordinate
            label = status if colour_by == "status" and status not in seen[name] else None
            if name == "assembly":
                for first, second in segments:
                    axis.plot(
                        *zip(first, second, strict=True),
                        color=colour,
                        alpha=0.5,
                        linewidth=0.7,
                        label=label,
                    )
                    label = None
                in_bounds = [p for p in points if all(abs(v) <= radius for v in p)]
                if in_bounds:
                    axis.scatter(*zip(*in_bounds, strict=True), color=colour, s=7)
                if all(abs(v) <= radius for v in points[-1]):
                    axis.scatter(*points[-1], color=colour, marker="x", s=15)
            elif axial:
                for first, second in segments:
                    axis.plot(
                        (first[section_coordinate], second[section_coordinate]),
                        (first[2], second[2]),
                        color=colour,
                        alpha=0.5,
                        linewidth=0.7,
                        label=label,
                    )
                    label = None
                for point in points:
                    if all(abs(v) <= radius for v in point):
                        axis.plot(
                            point[section_coordinate],
                            point[2],
                            marker=".",
                            color=colour,
                            markersize=3,
                        )
                terminal = points[-1]
                if all(abs(v) <= radius for v in terminal):
                    axis.plot(
                        terminal[section_coordinate],
                        terminal[2],
                        marker="x",
                        color=colour,
                        markersize=4,
                    )
            elif name == "pupil" and all(abs(v) <= radius for v in points[0][:2]):
                axis.scatter(points[0][0], points[0][1], s=9, color=colour, label=label)
            elif status == "detected":
                axis.scatter(points[-1][0], points[-1][1], s=9, color=colour, label=label)
            seen[name].add(status)
    for name, axis in axes.items():
        if clipped and (name.startswith("section") or name == "assembly"):
            if name == "assembly":
                axis.plot([], [], [], color="0.4", linestyle="--", label="incoming path clipped")
            else:
                axis.plot([], [], color="0.4", linestyle="--", label="incoming path clipped")
        if records and (colour_by == "status" or clipped):
            handles, labels = axis.get_legend_handles_labels()
            if handles:
                axis.legend(handles, labels, fontsize=8)
        if name == "assembly":
            axis.set(xlim=(-radius, radius), ylim=(-radius, radius), zlim=(-radius, radius))
            axis.set_box_aspect((1, 1, 1))
        if name.startswith("section"):
            axis.set_xlim(-radius, radius)
            axis.set_ylim(-radius, radius)


def draw_pupil(
    plt, optical_model_path: Path, output: Path, *, label_panels=False, colour_by="status"
) -> None:
    """Render a full-size exact entrance-pupil geometry view."""
    optical_model = load_plot_optical_model(optical_model_path)
    figure, axis = plt.subplots(figsize=(9, 8))
    if colour_by == "facet-z" and optical_model.kind != "segmented":
        raise ValueError("facet-z requires a segmented optical model")
    _draw_optical_model_pupil(
        axis, optical_model, title="Compiled entrance pupil", colour_by=colour_by
    )
    if label_panels:
        if optical_model.kind != "segmented":
            raise ValueError("panel labels require a segmented optical model")
        primary = [p for p in optical_model.polygons if p.role == "primary"]
        if len(primary) > 300:
            raise ValueError("panel labels are limited to 300 facets")
        for polygon in primary:
            centre = tuple(
                sum(v[i] for v in polygon.vertices_m) / len(polygon.vertices_m) for i in (0, 1)
            )
            axis.text(*centre, str(polygon.identifier), ha="center", va="center", fontsize=6)
    figure.suptitle(f"{optical_model.model_label} [{optical_model.kind}]")
    figure.subplots_adjust(bottom=0.1, top=0.93)
    figure.text(0.5, 0.025, _optical_model_footer(optical_model), ha="center", fontsize=8)
    _save_figure(figure, output, dpi=180, bbox_inches="tight")
    plt.close(figure)


def draw_section(plt, optical_model_path: Path, output: Path, section_plane: str = "xz") -> None:
    """Render one full-size x-z or y-z compiled optical section."""
    optical_model = load_plot_optical_model(optical_model_path)
    figure, axis = plt.subplots(figsize=(10, 7))
    _draw_optical_model_section(axis, optical_model, {"xz": 0, "yz": 1}[section_plane])
    figure.suptitle(f"{optical_model.model_label} compiled optical section [{optical_model.kind}]")
    figure.subplots_adjust(bottom=0.12, top=0.92)
    figure.text(0.5, 0.025, _optical_model_footer(optical_model), ha="center", fontsize=8)
    _save_figure(figure, output, dpi=180, bbox_inches="tight")
    plt.close(figure)


def draw_assembly_3d(plt, optical_model_path: Path, output: Path) -> None:
    """Render one full-size orthographic compiled assembly view."""
    optical_model = load_plot_optical_model(optical_model_path)
    figure = plt.figure(figsize=(11, 8))
    axis = figure.add_subplot(projection="3d")
    _draw_optical_model_assembly(axis, optical_model)
    figure.suptitle(f"{optical_model.model_label} compiled assembly [{optical_model.kind}]")
    figure.subplots_adjust(bottom=0.1, top=0.92)
    figure.text(0.5, 0.025, _optical_model_footer(optical_model), ha="center", fontsize=8)
    _save_figure(figure, output, dpi=180, bbox_inches="tight")
    plt.close(figure)


def draw_reference_telescope(axis, telescope: str):
    """Draw a side projection; reference models show optical surfaces only."""
    from matplotlib.patches import Rectangle

    if telescope == "reference-mst":
        mirror_x = [(-6.0 + 12.0 * index / 200.0) for index in range(201)]
        mirror_z = [9.75 - math.sqrt(9.75**2 - x**2) for x in mirror_x]
        axis.plot(mirror_x, mirror_z, color="0.15", linewidth=2.0, label="spherical primary")
        axis.add_patch(
            Rectangle(
                (-0.55, 4.625), 1.1, 0.5, color="black", alpha=0.25, label="camera body projection"
            )
        )
        for base_x, top_x in ((-4.2, -0.7), (0.0, 0.0), (4.2, 0.7)):
            axis.plot([base_x, top_x], [0.30, 4.625], color="0.3", linewidth=1)
        axis.plot([-2.0, 2.0], [4.875, 4.875], color="black", linewidth=0.7, label="focal screen")
        return
    # Only validated analytic reference outlines are drawn here.
    # SC geometry must come from a compiled-optical-model export.
    if telescope in ("SST", "SCT"):
        axis.text(
            0.02,
            0.98,
            "Unvalidated SC paths; geometry outline unavailable",
            transform=axis.transAxes,
            va="top",
        )
        return
    radius, camera_z, focal_radius = {
        "LST": (11.5, 28.0, 1.5),
        "MST": (6.0, 16.0, 1.0),
    }[telescope]
    primary_x = [(-radius + 2.0 * radius * index / 200.0) for index in range(201)]
    # An outline deliberately avoids reimplementing C++ prescription math.
    primary_z = (
        [x * x / (4 * camera_z) for x in primary_x]
        if telescope == "LST"
        else [2 * camera_z - math.sqrt((2 * camera_z) ** 2 - x * x) for x in primary_x]
    )
    axis.plot(primary_x, primary_z, color="0.15", linewidth=2.0, label="M1 reference aperture")
    axis.plot(
        [-focal_radius, focal_radius],
        [camera_z, camera_z],
        color="black",
        linewidth=0.9,
        label="focal plane",
    )


def draw_structure(plt, telescope: str, output: Path):
    """Render both side projections of the available telescope geometry."""
    figure, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    for axis, horizontal in zip(axes, "xy", strict=True):
        draw_reference_telescope(axis, telescope)
        axis.set(xlabel=f"telescope {horizontal} [m]", ylabel="telescope z [m]")
        axis.set_aspect("equal", adjustable="box")
        axis.set_ylim(-0.5, {"reference-mst": 6.0, "LST": 30.0, "MST": 18.0}.get(telescope, 8.0))
    axes[0].legend(loc="best")
    description = (
        "reference mirror, camera and supports"
        if telescope == "reference-mst"
        else "reference optical outline"
    )
    figure.suptitle(f"{telescope} {description}")
    figure.tight_layout()
    _save_figure(figure, output, dpi=160, bbox_inches="tight")
    plt.close(figure)


def draw_compiled_structure(plt, optical_model_path: Path, output: Path):
    """Compatibility alias for the compiled section renderer."""
    draw_section(plt, optical_model_path, output)


def draw_compiled_mirror(plt, optical_model_path: Path, output: Path):
    """Compatibility alias for the pupil renderer."""
    draw_pupil(plt, optical_model_path, output)


def _detector_polygons(optical_model_path: Path):
    """Return detector apertures through the common validated geometry loader."""
    return [
        p.vertices_m
        for p in load_plot_optical_model(optical_model_path).polygons
        if p.role == "detector"
    ]


def draw_compiled_3d(plt, optical_model_path: Path, output: Path):
    """Compatibility alias for the assembly renderer."""
    draw_assembly_3d(plt, optical_model_path, output)


def clip_recorded_path(points, radius: float):
    """Clip recorded straight flight segments to a telescope-frame cube.

    This only crops stored vertices for display; it performs no optical solve.
    """
    if not math.isfinite(radius) or radius <= 0:
        raise ValueError("context radius must be finite and positive")
    segments = []
    for first, second in zip(points, points[1:]):
        direction = tuple(b - a for a, b in zip(first, second, strict=True))
        lower, upper = 0.0, 1.0
        for origin, delta in zip(first, direction, strict=True):
            if delta == 0:
                if abs(origin) > radius:
                    upper = -1
                    break
            else:
                a, b = sorted(((-radius - origin) / delta, (radius - origin) / delta))
                lower, upper = max(lower, a), min(upper, b)
        if lower <= upper:
            segments.append(
                tuple(
                    tuple(a + t * d for a, d in zip(first, direction, strict=True))
                    for t in (lower, upper)
                )
            )
    return segments


_STATUS_COLOURS = {
    "detected": "#0072b2",
    "blocked_obscurer": "#d55e00",
    "absorbed_material": "#d55e00",
    "escaped_material": "#999999",
    "missed_secondary": "#cc79a7",
    "intersection_failure": "#000000",
    "blocked_camera": "#d55e00",
    "blocked_mast": "#d55e00",
    "missed_primary": "#777777",
    "missed_screen": "#cc79a7",
    "no_detector": "#cc79a7",
    "escaped_optical_model": "#999999",
    "interaction_limit": "#e69f00",
    "invalid_input": "#000000",
}


def select_trace_rows(
    path: Path,
    max_paths: int = 300,
    *,
    selection: str = "first",
    seed: int = 0,
    statuses=(),
    component_id: int | None = None,
):
    """Select recorded paths without changing vertices or using global randomness."""
    records = list(read_trace_rows(path))
    if component_id is not None and not all("terminal_surface_id" in row for row, _, _ in records):
        raise ValueError("component filter requires recorded terminal_surface_id")
    selected = [
        record
        for record in records
        if (not statuses or record[1] in statuses)
        and (component_id is None or record[0].get("terminal_surface_id") == str(component_id))
    ]
    if selection == "first":
        return selected[:max_paths]
    if selection != "stratified":
        raise ValueError("unknown path selection")
    groups = {}
    for index, record in enumerate(selected):
        row, status, _ = record
        identity = row.get("photon_id", str(index))
        rank = hashlib.sha256(f"{seed}:{identity}".encode()).hexdigest()
        groups.setdefault((status, row.get("source_kind", "unknown")), []).append((rank, record))
    queues = [
        iter(record for _, record in sorted(groups[key], key=lambda item: item[0]))
        for key in sorted(groups)
    ]
    result = []
    while queues and len(result) < max_paths:
        pending = []
        for queue in queues:
            record = next(queue, None)
            if record is not None:
                result.append(record)
                pending.append(queue)
                if len(result) == max_paths:
                    break
        queues = pending
    return result


def entrance_samples(path: Path):
    """Entrance point 0 and input weights; identities come only from diagnostics."""
    for row, status, points in read_trace_rows(path):
        x = float(row.get("entrance_x_m", points[0][0]))
        y = float(row.get("entrance_y_m", points[0][1]))
        weight = float(row.get("source_weight", "1") or "1")
        if not all(math.isfinite(value) for value in (x, y, weight)) or weight < 0:
            raise ValueError("invalid entrance point or source weight")
        yield x, y, weight, status, row.get("terminal_surface_id")


def weighted_loss_bins(samples, bins: int, extent: float, *, statuses=(), component_id=None):
    """Lost input weight / incident weight; empty bins remain None."""
    if bins < 1 or not math.isfinite(extent) or extent <= 0:
        raise ValueError("positive bins and finite positive extent are required")
    incident = [[0.0] * bins for _ in range(bins)]
    lost = [[0.0] * bins for _ in range(bins)]
    for x, y, weight, status, identifier in samples:
        if abs(x) > extent or abs(y) > extent:
            continue
        ix = min(bins - 1, int((x + extent) * bins / (2 * extent)))
        iy = min(bins - 1, int((y + extent) * bins / (2 * extent)))
        incident[iy][ix] += weight
        if (
            status != "detected"
            and (not statuses or status in statuses)
            and (component_id is None or identifier == str(component_id))
        ):
            lost[iy][ix] += weight
    return [
        [lost[y][x] / incident[y][x] if incident[y][x] else None for x in range(bins)]
        for y in range(bins)
    ]


def draw_loss_map(
    plt,
    path: Path,
    optical_model_path: Path,
    output: Path,
    bins: int = 64,
    statuses=(),
    component_id: int | None = None,
):
    """Draw factual terminal outcomes and weighted loss fractions at entrance."""
    import numpy as np

    optical_model = load_plot_optical_model(optical_model_path)
    samples = list(entrance_samples(path))
    if not samples:
        raise ValueError("loss map requires entrance records")
    if component_id is not None and any(sample[4] is None for sample in samples):
        raise ValueError("component filter requires recorded terminal_surface_id")
    bounds = _optical_model_bounds(optical_model)
    extent = max(abs(value) for bound in bounds[:2] for value in bound)
    extent = max(extent, max(max(abs(sample[0]), abs(sample[1])) for sample in samples), 1e-6)
    figure, axes = plt.subplots(1, 2, figsize=(13, 6))
    _draw_optical_model_pupil(axes[0], optical_model, title="Recorded entrance outcomes")
    for status in sorted({sample[3] for sample in samples}):
        group = [
            sample
            for sample in samples
            if sample[3] == status
            and (not statuses or status in statuses)
            and (component_id is None or sample[4] == str(component_id))
        ]
        if group:
            axes[0].scatter(
                [sample[0] for sample in group],
                [sample[1] for sample in group],
                s=9,
                color=_STATUS_COLOURS.get(status, "black"),
                label=status,
            )
    axes[0].legend(fontsize=8)
    fractions = np.asarray(
        weighted_loss_bins(samples, bins, extent, statuses=statuses, component_id=component_id),
        dtype=float,
    )
    image = axes[1].imshow(
        np.ma.masked_invalid(fractions),
        origin="lower",
        extent=(-extent, extent, -extent, extent),
        vmin=0,
        vmax=1,
        cmap="magma",
        interpolation="nearest",
    )
    figure.colorbar(image, ax=axes[1], label="lost input weight / incident-bin weight")
    axes[1].set(
        xlabel="entrance x [m]", ylabel="entrance y [m]", title="Weighted terminal loss fraction"
    )
    axes[1].set_aspect("equal")
    identities = (
        "recorded terminal component IDs"
        if any(sample[4] not in (None, "", "4294967295") for sample in samples)
        else "component unknown"
    )
    figure.suptitle(f"{optical_model.model_label}: {identities}; entrance = recorded point 0")
    figure.subplots_adjust(bottom=0.15, wspace=0.3)
    figure.text(0.5, 0.02, _optical_model_footer(optical_model), ha="center", fontsize=7)
    _save_figure(figure, output, dpi=160, bbox_inches="tight")
    plt.close(figure)


def _trace_colour_scale(plt, records, colour_by):
    """Validate an explicitly requested recorded numeric colour quantity."""
    if colour_by == "status":
        return None, None, None
    from matplotlib.colors import Normalize

    field, label = {
        "wavelength": ("wavelength_nm", "wavelength [nm]"),
        "arrival-time": ("arrival_time_ns", "recorded arrival time [ns]"),
        "incidence-angle": ("incidence_focal_deg", "focal incidence angle [deg]"),
    }[colour_by]
    try:
        values = [float(row[field]) for row, _, _ in records]
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"numeric path colour requires recorded {field}") from error
    if not values or not all(math.isfinite(value) for value in values):
        raise ValueError(f"numeric path colour requires finite {field} values")
    return plt.get_cmap("viridis"), Normalize(min(values), max(values)), label


def draw_rays(
    plt,
    path: Path,
    telescope: str,
    output: Path,
    max_paths: int,
    optical_model_path: Path | None = None,
    *,
    selection="first",
    seed=0,
    statuses=(),
    component_id=None,
    colour_by="status",
    context_radius_m=None,
):
    """Render recorded paths in the two telescope side projections."""
    figure, axes = plt.subplots(1, 2, figsize=(12, 7), sharey=True)
    optical_model = (
        load_plot_optical_model(optical_model_path) if optical_model_path is not None else None
    )
    if context_radius_m is None:
        context_radius_m = (
            max(abs(value) for pair in _optical_model_bounds(optical_model) for value in pair) * 1.2
            if optical_model
            else {"reference-mst": 8.0, "LST": 32.0, "MST": 20.0}.get(telescope, 8.0)
        )
    records = select_trace_rows(
        path,
        max_paths,
        selection=selection,
        seed=seed,
        statuses=statuses,
        component_id=component_id,
    )
    for axis, coordinate, horizontal in zip(axes, (0, 1), "xy", strict=True):
        if optical_model is not None:
            _draw_optical_model_section(axis, optical_model, coordinate)
        else:
            draw_reference_telescope(axis, telescope)
            axis.set(xlabel=f"telescope {horizontal} [m]", ylabel="telescope z [m]")
            axis.set_aspect("equal", adjustable="box")
            axis.set_ylim(
                -0.5, {"reference-mst": 8.0, "LST": 32.0, "MST": 20.0}.get(telescope, 8.0)
            )
    _overlay_paths(
        plt,
        dict(zip(("section-xz", "section-yz"), axes, strict=True)),
        records,
        0,
        context_radius_m,
        colour_by,
    )
    labels = (
        ", ".join(sorted({row.get("source_kind", "unknown") for row, _, _ in records}))
        or "unknown source"
    )
    incidence = [
        float(row["incidence_focal_deg"]) for row, _, _ in records if row.get("incidence_focal_deg")
    ]
    angle = f"; focal incidence mean {sum(incidence) / len(incidence):.3g}°" if incidence else ""
    title = optical_model.model_label if optical_model is not None else telescope
    figure.suptitle(f"{title} recorded ray paths — {labels}{angle}")
    if optical_model is not None:
        figure.text(0.5, 0.02, _optical_model_footer(optical_model), ha="center", fontsize=7)
    else:
        figure.text(0.5, 0.02, "Telescope frame: +z optical axis; metres", ha="center", fontsize=8)
    figure.tight_layout(rect=(0, 0.13, 1, 0.94))
    _save_figure(figure, output, dpi=160, bbox_inches="tight")
    plt.close(figure)


def draw_focal_plane(
    plt, path: Path, output: Path, bins: int, telescope: str, optical_model_path: Path | None = None
):
    """Render surviving-photon intensity and its Cartesian projections."""
    import numpy as np

    optical_model = load_plot_optical_model(optical_model_path) if optical_model_path else None
    model_label = optical_model.model_label if optical_model else telescope
    samples = list(focal_plane_hits(path))
    if not samples:
        raise SystemExit("No surviving photons reached the focal plane in this CSV")
    values = np.asarray(samples, dtype=float)
    x, y, weights = values.T
    extent = (
        compiled_focal_plane_extent(optical_model_path)
        if optical_model_path is not None
        else max(float(np.max(np.abs(x))), float(np.max(np.abs(y))), 1.0e-6) * 1.05
    )

    figure = plt.figure(figsize=(9, 8))
    grid = figure.add_gridspec(
        2, 2, width_ratios=(4, 1.25), height_ratios=(1.25, 4), hspace=0.06, wspace=0.06
    )
    top = figure.add_subplot(grid[0, 0])
    image = figure.add_subplot(grid[1, 0])
    right = figure.add_subplot(grid[1, 1], sharey=image)
    histogram = image.hist2d(
        x,
        y,
        bins=bins,
        range=[[-extent, extent], [-extent, extent]],
        weights=weights,
        cmap="viridis",
    )
    figure.colorbar(histogram[3], ax=image, label="weighted surviving photons / bin")
    bin_edges = np.linspace(-extent, extent, bins + 1)
    top.hist(x, bins=bin_edges, weights=weights, color="tab:blue")
    right.hist(y, bins=bin_edges, weights=weights, orientation="horizontal", color="tab:blue")
    image.set(
        xlabel="focal-plane x [m]",
        ylabel="focal-plane y [m]",
        title=f"{model_label} focal plane",
    )
    image.set_aspect("equal", adjustable="box")
    top.set(ylabel="weighted surviving photons")
    right.set(xlabel="weighted surviving photons")
    top.tick_params(labelbottom=False)
    right.tick_params(labelleft=False)
    if optical_model is not None:
        figure.subplots_adjust(bottom=0.14)
        figure.text(0.5, 0.02, _optical_model_footer(optical_model), ha="center", fontsize=7)
    _save_figure(figure, output, dpi=160, bbox_inches="tight")
    plt.close(figure)


def draw_focal_plane_svg(
    path: Path, output: Path, bins: int, telescope: str, optical_model_path: Path | None = None
):
    """Write a weighted PSF map and measured summary without plotting dependencies."""
    from obdeect.analysis import analyse_trace_csv

    optical_model = load_plot_optical_model(optical_model_path) if optical_model_path else None
    if optical_model is not None:
        telescope = optical_model.model_label
    samples = list(focal_plane_hits(path))
    result = analyse_trace_csv(path)
    extent = (
        compiled_focal_plane_extent(optical_model_path)
        if optical_model_path is not None
        else max(max(abs(x), abs(y)) for x, y, _ in samples) * 1.05
    )
    extent = max(extent, 1e-6)
    histogram = [[0.0] * bins for _ in range(bins)]
    for x, y, weight in samples:
        if abs(x) > extent or abs(y) > extent:
            continue
        ix = min(bins - 1, max(0, int((x + extent) * bins / (2 * extent))))
        iy = min(bins - 1, max(0, int((y + extent) * bins / (2 * extent))))
        histogram[iy][ix] += weight
    maximum = max(max(row) for row in histogram)
    left = 95
    top = 90
    size = 520
    cell = size / bins
    elements = [
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 720 760">',
        '<rect width="720" height="760" fill="white"/>',
        f"<title>{escape(telescope)} weighted focal-plane PSF</title>",
        f'<text x="95" y="38" font-size="24">{escape(telescope)} focal-plane PSF</text>',
        '<rect x="95" y="90" width="520" height="520" fill="#101827"/>',
    ]
    for iy, row in enumerate(histogram):
        for ix, weight in enumerate(row):
            if weight <= 0:
                continue
            fraction = math.sqrt(weight / maximum)
            red = round(25 + 230 * fraction)
            green = round(35 + 190 * fraction)
            blue = round(70 + 85 * (1 - fraction))
            elements.append(
                f'<rect x="{left + ix * cell:.4f}" y="{top + (bins - 1 - iy) * cell:.4f}" '
                f'width="{cell:.4f}" height="{cell:.4f}" fill="#{red:02x}{green:02x}{blue:02x}"/>'
            )
    centroid_x = left + (result.centroid_x_m + extent) * size / (2 * extent)
    centroid_y = top + (extent - result.centroid_y_m) * size / (2 * extent)
    elements.extend([
        f'<circle cx="{centroid_x:.4f}" cy="{centroid_y:.4f}" r="5" '
        'fill="none" stroke="white" stroke-width="2"/>',
        '<text x="95" y="645" font-size="16">x and y [m]; '
        "colour: weighted detected photons per bin</text>",
        f'<text x="95" y="676" font-size="16">Centroid: '
        f"({result.centroid_x_m:.6g}, {result.centroid_y_m:.6g}) m</text>",
        f'<text x="95" y="702" font-size="16">D80: {result.d80_m:.6g} m; '
        f"throughput: {result.optical_throughput:.4g}</text>",
        f'<text x="95" y="728" font-size="14">{bins} x {bins} bins; '
        f"axis extent: ±{extent:.6g} m; input: {escape(path.name)}</text>",
        "</svg>",
    ])
    if optical_model is not None:
        elements[0] = '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 720 815">'
        elements[1] = '<rect width="720" height="815" fill="white"/>'
        elements.insert(
            -1,
            f'<text x="95" y="755" font-size="11">'
            f"Compiled model SHA256: {escape(optical_model.model_sha256 or 'unavailable')}</text>",
        )
        elements.insert(
            -1,
            '<text x="95" y="780" font-size="11">Telescope frame: '
            f"+z optical axis; origin {escape(optical_model.frame_origin)}; metres</text>",
        )
    output.write_text("\n".join(elements) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(
        description="Render compiled telescope geometry and recorded optical diagnostics."
    )
    parser.add_argument("--input", type=Path, help="trace CSV (required for rays and focal plane)")
    parser.add_argument("--output", type=Path, default=Path("optical-diagnostics.png"))
    parser.add_argument("--max-paths", type=int, default=300)
    parser.add_argument(
        "--view",
        choices=(
            "structure",
            "telescope",
            "pupil",
            "section",
            "assembly-3d",
            "compiled-structure",
            "compiled-mirror",
            "compiled-3d",
            "rays",
            "focal-plane",
            "loss-map",
        ),
        default="rays",
        help=(
            "compiled-model geometry, recorded ray paths, or focal-plane image; "
            "structure uses a reference outline only when no optical model is supplied"
        ),
    )
    parser.add_argument(
        "--focal-plane",
        action="store_true",
        help=(
            "plot surviving focal-plane photons and weighted x/y projections "
            "instead of photon paths"
        ),
    )
    parser.add_argument("--bins", type=int, default=64, help="focal-plane histogram bins per axis")
    parser.add_argument(
        "--optical-model-json",
        type=Path,
        help="compiled optical model JSON for compiled views or focal-plane extent",
    )
    parser.add_argument(
        "--section-plane",
        choices=("xz", "yz"),
        default="xz",
        help="telescope-frame plane for --view section and telescope",
    )
    parser.add_argument(
        "--telescope",
        choices=TELESCOPE_NAMES,
        default="reference-mst",
        help="outline only; does not modify trace coordinates",
    )
    parser.add_argument("--path-selection", choices=("first", "stratified"))
    parser.add_argument("--selection-seed", type=int)
    parser.add_argument("--status", nargs="+", default=[])
    parser.add_argument("--component-id", type=int)
    parser.add_argument("--label-panels", action="store_true")
    parser.add_argument(
        "--panel", choices=("assembly", "section", "pupil", "focal-geometry", "focal-plane-hits")
    )
    parser.add_argument("--show-paths", action=argparse.BooleanOptionalAction)
    parser.add_argument("--context-radius-m", type=float)
    parser.add_argument(
        "--colour-by",
        choices=("status", "wavelength", "arrival-time", "incidence-angle", "facet-z"),
        default="status",
    )
    args = parser.parse_args()
    if args.context_radius_m is not None and (
        not math.isfinite(args.context_radius_m) or args.context_radius_m <= 0
    ):
        parser.error("--context-radius-m must be finite and positive")
    if args.component_id is not None and not 0 <= args.component_id < 4294967295:
        parser.error("--component-id must be a finite compiled component identifier")
    if args.panel is not None and args.view != "telescope":
        parser.error("--panel requires --view telescope")
    if args.panel == "focal-plane-hits" and (args.input is None or args.view != "telescope"):
        parser.error("--panel focal-plane-hits requires telescope view and --input")
    from obdeect.result_contract import TERMINAL_STATUSES

    if any(status not in TERMINAL_STATUSES for status in args.status):
        parser.error("unsupported terminal status filter")
    if args.max_paths < 1 or args.bins < 1:
        parser.error("--max-paths and --bins must be positive")
    if (args.status or args.component_id is not None) and args.view not in {
        "rays",
        "telescope",
        "loss-map",
    }:
        parser.error("trace filters require rays, telescope, or loss-map view")
    if args.colour_by == "facet-z" and args.view not in {"pupil", "compiled-mirror"}:
        parser.error("facet-z requires --view pupil")
    if args.colour_by not in {"status", "facet-z"} and (
        args.input is None or args.view not in {"rays", "telescope"}
    ):
        parser.error("numeric path colours require --input with rays or telescope view")
    if args.output.suffix.lower() not in {".png", ".pdf", ".svg"}:
        parser.error("supported outputs are PNG, PDF, and SVG")
    if args.input is None and (
        args.status
        or args.component_id is not None
        or args.path_selection is not None
        or args.selection_seed is not None
        or args.show_paths is not None
        or args.context_radius_m is not None
    ):
        parser.error("trace selection requires --input")
    args.path_selection = args.path_selection or "first"
    args.selection_seed = args.selection_seed or 0
    args.show_paths = args.show_paths is not False
    if args.label_panels and args.view not in {"pupil", "compiled-mirror"}:
        parser.error("--label-panels requires --view pupil")
    if args.view == "loss-map" and args.optical_model_json is None:
        parser.error("--optical-model-json is required for --view loss-map")
    compiled_views = {
        "telescope",
        "pupil",
        "section",
        "assembly-3d",
        "compiled-structure",
        "compiled-mirror",
        "compiled-3d",
    }
    if args.input is None and args.view not in {"structure", *compiled_views}:
        parser.error("--input is required for rays and focal plane")
    if args.focal_plane and args.view != "rays":
        parser.error("--focal-plane cannot be combined with --view")
    if args.view in compiled_views and args.optical_model_json is None:
        parser.error(f"--optical-model-json is required for --view {args.view}")

    if (args.focal_plane or args.view == "focal-plane") and args.output.suffix.lower() == ".svg":
        draw_focal_plane_svg(
            args.input, args.output, args.bins, args.telescope, args.optical_model_json
        )
        return

    try:
        plt = _FileRenderer()
    except ImportError as error:
        raise SystemExit(
            "Install optional visualization dependency: python -m pip install matplotlib"
        ) from error

    if args.view == "loss-map":
        draw_loss_map(
            plt,
            args.input,
            args.optical_model_json,
            args.output,
            args.bins,
            args.status,
            args.component_id,
        )
        return
    if args.view == "structure":
        if args.optical_model_json is not None:
            draw_section(plt, args.optical_model_json, args.output, args.section_plane)
        else:
            draw_structure(plt, args.telescope, args.output)
        return
    if args.view == "telescope":
        draw_telescope(
            plt,
            args.optical_model_json,
            args.output,
            args.section_plane,
            input_path=args.input,
            show_paths=args.show_paths,
            panel=args.panel,
            bins=args.bins,
            max_paths=args.max_paths,
            selection=args.path_selection,
            seed=args.selection_seed,
            statuses=args.status,
            component_id=args.component_id,
            context_radius_m=args.context_radius_m,
            colour_by=args.colour_by,
        )
        return
    if args.view in {"pupil", "compiled-mirror"}:
        draw_pupil(
            plt,
            args.optical_model_json,
            args.output,
            label_panels=args.label_panels,
            colour_by=args.colour_by,
        )
        return
    if args.view in {"section", "compiled-structure"}:
        draw_section(plt, args.optical_model_json, args.output, args.section_plane)
        return
    if args.view in {"assembly-3d", "compiled-3d"}:
        draw_assembly_3d(plt, args.optical_model_json, args.output)
        return
    if args.focal_plane or args.view == "focal-plane":
        draw_focal_plane(
            plt, args.input, args.output, args.bins, args.telescope, args.optical_model_json
        )
        return
    draw_rays(
        plt,
        args.input,
        args.telescope,
        args.output,
        args.max_paths,
        args.optical_model_json,
        selection=args.path_selection,
        seed=args.selection_seed,
        statuses=args.status,
        component_id=args.component_id,
        colour_by=args.colour_by,
        context_radius_m=args.context_radius_m,
    )


if __name__ == "__main__":
    main()
