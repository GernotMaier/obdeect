#!/usr/bin/env python3
"""Optional visualization for CSV paths from the reference and CTAO reference tracers."""

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Literal
from xml.sax.saxutils import escape


def _parse_path_row(path: Path, row: dict[str, str]) -> list[tuple[float, float, float]]:
    """Validate and return the ragged path vertices from one CSV record."""
    points = []
    try:
        count = int(row["point_count"])
        if not 1 <= count <= 4:
            raise ValueError("point_count must be between 1 and 4")
        for index in range(count):
            point = tuple(float(row[f"{axis}{index}_m"]) for axis in "xyz")
            if not all(math.isfinite(value) for value in point):
                raise ValueError("non-finite path vertex")
            points.append(point)
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"{path}: invalid path vertices: {error}") from error
    return points


def read_paths(path: Path):
    """Yield terminal status and validated vertices from a trace CSV."""
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            yield row["status"], _parse_path_row(path, row)


def read_trace_rows(path: Path):
    """Yield validated path records with source and spectral metadata."""
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            yield row, row["status"], _parse_path_row(path, row)


def focal_plane_hits(path: Path):
    """Yield detected focal-plane ``(x, y, weight)`` samples from a trace CSV.

    The final ragged path vertex is the detector-surface intersection.  Modern
    CSV output supplies source weights and throughput; older files remain
    usable with unit weight.
    """
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
        if trace_model["kind"] == "segmented":
            surfaces = trace_model["detector_surfaces"]
            extent = max(
                math.hypot(float(surface["centre_m"][0]), float(surface["centre_m"][1]))
                + float(surface["diameter_m"]) * 0.5
                for surface in surfaces
            )
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
    role: Literal["primary", "detector"]
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
class PlotScene:
    """Validated renderer-facing view of a compiled optical model."""

    kind: Literal["segmented", "axisymmetric"]
    model_label: str
    model_sha256: str | None
    polygons: tuple[PlotPolygon, ...]
    axisymmetric_surfaces: tuple[AxisymmetricSurface, ...]
    obscurers: tuple[PlotObscurer, ...]
    unavailable_roles: tuple[str, ...]


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
    surface: dict, *, identifier: int, role: Literal["primary", "detector"], description: str
) -> PlotPolygon:
    """Build one finite aperture polygon from a serialised tangent frame."""
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


def load_plot_scene(path: Path) -> PlotScene:
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
    unavailable: tuple[str, ...] = ()
    if isinstance(declared, dict) and "unavailable_roles" in declared:
        roles = declared["unavailable_roles"]
        if not isinstance(roles, list) or not all(isinstance(role, str) for role in roles):
            raise ValueError(f"{path}: invalid plot_geometry unavailable_roles")
        unavailable = tuple(sorted(set(roles)))

    try:
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
                    identifier=int(facet.get("id", index)),
                    role="primary",
                    description=f"primary facet {index}",
                )
                for index, facet in enumerate(facets)
            ) + tuple(
                _polygon_from_surface(
                    detector,
                    identifier=int(detector.get("id", len(facets) + index)),
                    role="detector",
                    description=f"detector surface {index}",
                )
                for index, detector in enumerate(detectors)
            )
            obscurer_rows = []
            for index, obscurer in enumerate(obscurers):
                first = _vector3(obscurer["first_endpoint_m"], f"obscurer {index} first endpoint")
                second = _vector3(
                    obscurer["second_endpoint_m"], f"obscurer {index} second endpoint"
                )
                diameter = float(obscurer["diameter_m"])
                if not math.isfinite(diameter) or diameter <= 0.0 or first == second:
                    raise ValueError(f"invalid obscurer {index}")
                obscurer_rows.append(PlotObscurer(int(obscurer["id"]), first, second, diameter))
            identifiers = [polygon.identifier for polygon in polygons]
            identifiers.extend(obscurer.identifier for obscurer in obscurer_rows)
            if len(identifiers) != len(set(identifiers)):
                raise ValueError("compiled component IDs must be unique")
            return PlotScene(
                "segmented",
                _model_label(optical_model),
                model_sha,
                polygons,
                (),
                tuple(obscurer_rows),
                unavailable,
            )
        surfaces = tuple(
            _axisymmetric_surface(trace_model[source], role)
            for source, role in (
                ("primary", "primary"),
                ("secondary", "secondary"),
                ("detector", "detector"),
            )
        )
        return PlotScene(
            "axisymmetric", _model_label(optical_model), model_sha, (), surfaces, (), unavailable
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

    return tuple((radius, sag(radius)) for radius in radii)


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
        scene = load_plot_scene(path)
        return [
            (
                polygon.vertices_m,
                sum(vertex[2] for vertex in polygon.vertices_m) / len(polygon.vertices_m),
            )
            for polygon in scene.polygons
            if polygon.role == "primary"
        ]
    except (KeyError, TypeError, ValueError, json.JSONDecodeError) as error:
        raise ValueError(f"{path}: invalid compiled primary-facet geometry: {error}") from error


_ROLE_STYLE = {
    "primary": {"color": "#607d8b", "label": "M1"},
    "secondary": {"color": "#455a64", "label": "M2"},
    "detector": {"color": "#e69f00", "label": "detector"},
    "obscurer": {"color": "#a14b3b", "label": "opaque cylinder"},
}


def _scene_bounds(scene: PlotScene):
    """Return physical 3-D bounds covering every finite scene primitive."""
    values = [[], [], []]
    for polygon in scene.polygons:
        for vertex in polygon.vertices_m:
            for coordinate, value in enumerate(vertex):
                values[coordinate].append(value)
    for surface in scene.axisymmetric_surfaces:
        profile = _axisymmetric_profile(surface)
        values[0].extend((-surface.outer_radius_m, surface.outer_radius_m))
        values[1].extend((-surface.outer_radius_m, surface.outer_radius_m))
        values[2].extend(height for _, height in profile)
    for obscurer in scene.obscurers:
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


def _set_equal_3d_bounds(axis, scene: PlotScene) -> None:
    """Apply equal physical scales with a small readable margin."""
    bounds = _scene_bounds(scene)
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
    profile = _axisymmetric_profile(surface, samples=64)
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


def _draw_scene_assembly(axis, scene: PlotScene) -> None:
    """Draw a restrained orthographic 3-D scene without synthetic hardware."""
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    for role in ("primary", "detector"):
        polygons = [polygon.vertices_m for polygon in scene.polygons if polygon.role == role]
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
    for surface in scene.axisymmetric_surfaces:
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
    for obscurer in scene.obscurers:
        axis.add_collection3d(
            Poly3DCollection(
                _cylinder_mesh(obscurer),
                facecolors=_ROLE_STYLE["obscurer"]["color"],
                edgecolors="#57261e",
                linewidths=0.15,
                alpha=0.62,
                label=f"opaque cylinders ({len(scene.obscurers)})",
            )
        )
        break
    _set_equal_3d_bounds(axis, scene)
    axis.set_proj_type("ortho")
    axis.view_init(elev=24, azim=-55)
    axis.set(xlabel="x [m]", ylabel="y [m]", zlabel="z [m]", title="A. orthographic assembly")


def _draw_scene_section(axis, scene: PlotScene, coordinate: int = 0) -> None:
    """Draw the exact finite geometry in one telescope-frame axial projection."""
    horizontal = "xy"[coordinate]
    for role in ("primary", "detector"):
        style = _ROLE_STYLE[role]
        for polygon in (item for item in scene.polygons if item.role == role):
            points = [(vertex[coordinate], vertex[2]) for vertex in polygon.vertices_m]
            points.append(points[0])
            axis.plot(*zip(*points, strict=True), color=style["color"], linewidth=1.0)
    for surface in scene.axisymmetric_surfaces:
        profile = _axisymmetric_profile(surface)
        positive = [(radius, height) for radius, height in profile]
        negative = [(-radius, height) for radius, height in reversed(profile)]
        axis.plot(
            *zip(*(negative + positive), strict=True),
            color=_ROLE_STYLE[surface.role]["color"],
            linewidth=1.4,
        )
    for obscurer in scene.obscurers:
        first = obscurer.first_endpoint_m
        second = obscurer.second_endpoint_m
        axis.plot(
            (first[coordinate], second[coordinate]),
            (first[2], second[2]),
            color=_ROLE_STYLE["obscurer"]["color"],
            linewidth=max(1.0, obscurer.diameter_m * 4.0),
            alpha=0.75,
        )
    axis.axvline(0.0, color="#546e7a", linewidth=0.7, linestyle="--", zorder=-1)
    axis.set(
        xlabel=f"telescope {horizontal} [m]", ylabel="telescope z [m]", title="B. axial section"
    )
    axis.set_aspect("equal", adjustable="box")


def _draw_scene_pupil(axis, scene: PlotScene, *, title: str = "C. entrance pupil") -> None:
    """Draw exact finite aperture boundaries looking along the optical axis."""
    from matplotlib.collections import PolyCollection
    from matplotlib.patches import Circle

    for role in ("primary", "detector"):
        polygons = [
            [(vertex[0], vertex[1]) for vertex in polygon.vertices_m]
            for polygon in scene.polygons
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
            axis.add_collection(collection)
    for surface in scene.axisymmetric_surfaces:
        style = _ROLE_STYLE[surface.role]
        axis.add_patch(
            Circle(
                (0.0, 0.0), surface.outer_radius_m, fill=False, color=style["color"], linewidth=1.2
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
    axis.autoscale_view()
    axis.set(xlabel="telescope x [m]", ylabel="telescope y [m]", title=title)
    axis.set_aspect("equal", adjustable="box")


def _draw_focal_geometry(axis, scene: PlotScene) -> None:
    """Draw only compiled detector geometry in the telescope x/y plane."""
    from matplotlib.collections import PolyCollection
    from matplotlib.patches import Circle

    polygons = [
        [(vertex[0], vertex[1]) for vertex in polygon.vertices_m]
        for polygon in scene.polygons
        if polygon.role == "detector"
    ]
    if polygons:
        axis.add_collection(
            PolyCollection(
                polygons, facecolors="#e69f00", edgecolors="#5f4600", linewidths=0.4, alpha=0.7
            )
        )
        axis.autoscale_view()
    for surface in scene.axisymmetric_surfaces:
        if surface.role == "detector":
            axis.add_patch(Circle((0.0, 0.0), surface.outer_radius_m, color="#e69f00", alpha=0.7))
            if surface.inner_radius_m > 0.0:
                axis.add_patch(Circle((0.0, 0.0), surface.inner_radius_m, color="white"))
    axis.set(xlabel="focal-plane x [m]", ylabel="focal-plane y [m]", title="D. focal geometry")
    axis.set_aspect("equal", adjustable="box")


def _scene_footer(scene: PlotScene) -> str:
    """Describe rendered and known-unavailable geometry in one compact line."""
    primary_count = sum(polygon.role == "primary" for polygon in scene.polygons)
    detector_count = sum(polygon.role == "detector" for polygon in scene.polygons)
    rendered = []
    if primary_count:
        rendered.append(f"{primary_count} M1 facets")
    if scene.axisymmetric_surfaces:
        rendered.extend(
            _ROLE_STYLE[surface.role]["label"] for surface in scene.axisymmetric_surfaces
        )
    if detector_count:
        rendered.append(f"{detector_count} detector surfaces")
    if scene.obscurers:
        rendered.append(f"{len(scene.obscurers)} opaque cylinders")
    coverage = (
        "not modelled: " + ", ".join(scene.unavailable_roles)
        if scene.unavailable_roles
        else "mechanical geometry coverage not declared"
    )
    sha = scene.model_sha256[:12] if scene.model_sha256 else "unavailable"
    return (
        f"Telescope frame (x, y, z), metres | model SHA256: {sha} | "
        f"rendered: {', '.join(rendered)} | {coverage}"
    )


def draw_telescope(plt, optical_model_path: Path, output: Path, section_plane: str = "xz") -> None:
    """Render the four-panel model-faithful telescope diagnostic plate."""
    scene = load_plot_scene(optical_model_path)
    coordinate = {"xz": 0, "yz": 1}[section_plane]
    figure = plt.figure(figsize=(14, 11))
    grid = figure.add_gridspec(2, 2)
    assembly = figure.add_subplot(grid[0, 0], projection="3d")
    section = figure.add_subplot(grid[0, 1])
    pupil = figure.add_subplot(grid[1, 0])
    focal = figure.add_subplot(grid[1, 1])
    _draw_scene_assembly(assembly, scene)
    _draw_scene_section(section, scene, coordinate)
    _draw_scene_pupil(pupil, scene)
    _draw_focal_geometry(focal, scene)
    figure.suptitle(f"{scene.model_label} compiled telescope geometry [{scene.kind}]")
    figure.subplots_adjust(bottom=0.09, top=0.93, wspace=0.3, hspace=0.26)
    figure.text(0.5, 0.025, _scene_footer(scene), ha="center", fontsize=8)
    figure.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(figure)


def draw_pupil(plt, optical_model_path: Path, output: Path) -> None:
    """Render a full-size exact entrance-pupil geometry view."""
    scene = load_plot_scene(optical_model_path)
    figure, axis = plt.subplots(figsize=(9, 8))
    _draw_scene_pupil(axis, scene, title="Compiled entrance pupil")
    figure.suptitle(f"{scene.model_label} [{scene.kind}]")
    figure.subplots_adjust(bottom=0.1, top=0.93)
    figure.text(0.5, 0.025, _scene_footer(scene), ha="center", fontsize=8)
    figure.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(figure)


def draw_section(plt, optical_model_path: Path, output: Path, section_plane: str = "xz") -> None:
    """Render one full-size x-z or y-z compiled optical section."""
    scene = load_plot_scene(optical_model_path)
    figure, axis = plt.subplots(figsize=(10, 7))
    _draw_scene_section(axis, scene, {"xz": 0, "yz": 1}[section_plane])
    figure.suptitle(f"{scene.model_label} compiled optical section [{scene.kind}]")
    figure.subplots_adjust(bottom=0.12, top=0.92)
    figure.text(0.5, 0.025, _scene_footer(scene), ha="center", fontsize=8)
    figure.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(figure)


def draw_assembly_3d(plt, optical_model_path: Path, output: Path) -> None:
    """Render one full-size orthographic compiled assembly view."""
    scene = load_plot_scene(optical_model_path)
    figure = plt.figure(figsize=(11, 8))
    axis = figure.add_subplot(projection="3d")
    _draw_scene_assembly(axis, scene)
    figure.suptitle(f"{scene.model_label} compiled assembly [{scene.kind}]")
    figure.subplots_adjust(bottom=0.1, top=0.92)
    figure.text(0.5, 0.025, _scene_footer(scene), ha="center", fontsize=8)
    figure.savefig(output, dpi=180, bbox_inches="tight")
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
    figure.savefig(output, dpi=160, bbox_inches="tight")
    plt.close(figure)


def draw_compiled_structure(plt, optical_model_path: Path, output: Path):
    """Plot the explicit compiled primary facets and focal boundary.

    The plot consumes only geometry present in the compiled optical_model JSON. Any
    unavailable secondary, support, or obscuration geometry is listed in the
    figure annotation instead of being replaced by an illustrative outline.
    """
    optical_model = json.loads(optical_model_path.read_text(encoding="utf-8"))
    facets = optical_model.get("primary", {}).get("facets", [])
    if not facets:
        raise SystemExit("compiled optical model contains no primary facets")
    figure, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    for axis, coordinate in zip(axes, (0, 1), strict=True):
        for facet in facets:
            centre = facet.get("nominal_centre_m")
            if not isinstance(centre, list) or len(centre) != 3:
                continue
            axis.scatter(centre[coordinate], centre[2], s=8, c="tab:blue", alpha=0.65)
        camera = optical_model.get("camera", {})
        elements = camera.get("pixels", [])
        if elements:
            extent = max(
                math.hypot(float(item["centre_xy_m"][0]), float(item["centre_xy_m"][1]))
                for item in elements
            )
            focal = float(optical_model.get("focal_length_m", facets[0].get("focal_length_m", 0.0)))
            axis.plot(
                [-extent, extent],
                [focal, focal],
                color="tab:orange",
                linewidth=2,
                label="focal boundary",
            )
        axis.set(xlabel=f"compiled telescope {'xy'[coordinate]} [m]", ylabel="compiled z [m]")
        axis.set_aspect("equal", adjustable="box")
    axes[0].scatter([], [], s=8, c="tab:blue", label="mirror facet centres")
    axes[0].legend(loc="best")
    blockers = optical_model.get("report", {}).get("trace_blockers", [])
    subtitle = "compiled model geometry"
    if blockers:
        subtitle += "; unresolved: " + ", ".join(str(item) for item in blockers[:2])
    figure.suptitle(subtitle)
    figure.tight_layout()
    figure.savefig(output, dpi=160, bbox_inches="tight")
    plt.close(figure)


def draw_compiled_mirror(plt, optical_model_path: Path, output: Path):
    """Render the primary as seen from above, coloured by facet-centre height."""
    from matplotlib.collections import PolyCollection

    facets = compiled_facet_polygons(optical_model_path)
    polygons = [[(vertex[0], vertex[1]) for vertex in vertices] for vertices, _ in facets]
    heights = [height for _, height in facets]
    figure, axis = plt.subplots(figsize=(9, 8))
    collection = PolyCollection(
        polygons,
        array=heights,
        cmap="viridis",
        edgecolors="0.15",
        linewidths=0.32,
    )
    axis.add_collection(collection)
    axis.autoscale_view()
    axis.set_aspect("equal", adjustable="box")
    axis.set(xlabel="telescope x [m]", ylabel="telescope y [m]", title="Primary mirror face")
    colourbar = figure.colorbar(collection, ax=axis, pad=0.02)
    colourbar.set_label("facet-centre z [m]")
    figure.suptitle("Compiled mirror panels, viewed along the optical axis")
    figure.tight_layout()
    figure.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(figure)


def _detector_polygons(optical_model_path: Path):
    """Return detector-surface boundary polygons available in a segmented model."""
    try:
        trace_model = json.loads(optical_model_path.read_text(encoding="utf-8"))["trace_model"]
        detectors = trace_model["detector_surfaces"]
        result = []
        for index, detector in enumerate(detectors):
            centre = _vector3(detector["centre_m"], f"detector {index} centre_m")
            normal = _normalised(
                _vector3(detector["normal"], f"detector {index} normal"), "detector normal"
            )
            tangent = _normalised(
                _vector3(detector["tangent"], f"detector {index} tangent"), "detector tangent"
            )
            binormal = _normalised(_cross(normal, tangent), "detector tangent and normal")
            vertices = tuple(
                tuple(
                    centre[axis] + tangent[axis] * local_x + binormal[axis] * local_y
                    for axis in range(3)
                )
                for local_x, local_y in _aperture_local_vertices(
                    str(detector["shape"]), float(detector["diameter_m"])
                )
            )
            result.append(vertices)
        return result
    except (KeyError, TypeError, ValueError, json.JSONDecodeError) as error:
        raise ValueError(
            f"{optical_model_path}: invalid compiled detector geometry: {error}"
        ) from error


def draw_compiled_3d(plt, optical_model_path: Path, output: Path):
    """Render recorded optical surfaces in a CAD-like orthographic 3-D view."""
    from matplotlib import colors
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    facets = compiled_facet_polygons(optical_model_path)
    detectors = _detector_polygons(optical_model_path)
    heights = [height for _, height in facets]
    colour_map = plt.get_cmap("viridis")
    normaliser = colors.Normalize(vmin=min(heights), vmax=max(heights) or 1.0)
    figure = plt.figure(figsize=(11, 8), layout="constrained")
    axis = figure.add_subplot(projection="3d")
    mirror = Poly3DCollection(
        [vertices for vertices, _ in facets],
        facecolors=[colour_map(normaliser(height)) for height in heights],
        edgecolors="0.12",
        linewidths=0.18,
        alpha=0.98,
    )
    axis.add_collection3d(mirror)
    if detectors:
        camera = Poly3DCollection(
            detectors,
            facecolors="tab:orange",
            edgecolors="0.15",
            linewidths=0.4,
            alpha=0.72,
        )
        axis.add_collection3d(camera)

    all_vertices = [vertex for vertices, _ in facets for vertex in vertices]
    all_vertices.extend(vertex for vertices in detectors for vertex in vertices)
    limits = []
    for coordinate in range(3):
        values = [vertex[coordinate] for vertex in all_vertices]
        centre = (min(values) + max(values)) * 0.5
        limits.append((centre, max(max(values) - min(values), 1.0)))
    largest_span = max(span for _, span in limits)
    for coordinate, (centre, _) in enumerate(limits):
        getattr(axis, f"set_{'xyz'[coordinate]}lim")(
            centre - largest_span * 0.55, centre + largest_span * 0.55
        )
    axis.set_box_aspect((1, 1, max(limits[2][1] / largest_span, 0.25)))
    axis.set_proj_type("ortho")
    axis.view_init(elev=23, azim=-57)
    axis.set(xlabel="x [m]", ylabel="y [m]", zlabel="z [m]")
    axis.set_title("Compiled optical geometry")
    colourbar = figure.colorbar(
        plt.cm.ScalarMappable(norm=normaliser, cmap=colour_map), ax=axis, shrink=0.62, pad=0.08
    )
    colourbar.set_label("primary facet-centre z [m]")
    figure.text(
        0.5,
        0.02,
        "Primary panels and focal detector surfaces present in the optical-model JSON",
        ha="center",
        fontsize=9,
    )
    figure.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(figure)


def draw_rays(
    plt,
    path: Path,
    telescope: str,
    output: Path,
    max_paths: int,
    optical_model_path: Path | None = None,
):
    """Render recorded paths in the two telescope side projections."""
    figure, axes = plt.subplots(1, 2, figsize=(12, 7), sharey=True)
    colours = {
        "detected": "tab:blue",
        "blocked_camera": "tab:red",
        "blocked_mast": "tab:orange",
        "missed_primary": "0.5",
        "missed_screen": "tab:purple",
    }
    source_colours = {"star": "tab:blue", "illuminator": "tab:green", "laser": "tab:red"}
    source_seen = set()
    incidence = []
    for count, (row, status, points) in enumerate(read_trace_rows(path)):
        if count >= max_paths:
            break
        source = row.get("source_kind", "unknown")
        source_seen.add(source)
        colour = source_colours.get(source, colours.get(status, "black"))
        if row.get("incidence_focal_deg"):
            incidence.append(float(row["incidence_focal_deg"]))
        for axis, coordinate in zip(axes, (0, 1), strict=True):
            axis.plot(
                [point[coordinate] for point in points],
                [point[2] for point in points],
                color=colour,
                alpha=0.25,
                linewidth=0.7,
            )
    scene = load_plot_scene(optical_model_path) if optical_model_path is not None else None
    for axis, coordinate, horizontal in zip(axes, (0, 1), "xy", strict=True):
        if scene is not None:
            _draw_scene_section(axis, scene, coordinate)
        else:
            draw_reference_telescope(axis, telescope)
            axis.set(xlabel=f"telescope {horizontal} [m]", ylabel="telescope z [m]")
            axis.set_aspect("equal", adjustable="box")
            axis.set_ylim(
                -0.5, {"reference-mst": 8.0, "LST": 32.0, "MST": 20.0}.get(telescope, 8.0)
            )
    axes[0].legend(loc="best")
    labels = ", ".join(sorted(source_seen)) or "unknown source"
    angle = f"; focal incidence mean {sum(incidence) / len(incidence):.3g}°" if incidence else ""
    title = scene.model_label if scene is not None else telescope
    figure.suptitle(f"{title} recorded ray paths — {labels}{angle}")
    figure.tight_layout()
    figure.savefig(output, dpi=160, bbox_inches="tight")
    plt.close(figure)


def draw_focal_plane(
    plt, path: Path, output: Path, bins: int, telescope: str, optical_model_path: Path | None = None
):
    """Render surviving-photon intensity and its Cartesian projections."""
    import numpy as np

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
        title=f"{telescope} focal plane",
    )
    image.set_aspect("equal", adjustable="box")
    top.set(ylabel="weighted surviving photons")
    right.set(xlabel="weighted surviving photons")
    top.tick_params(labelbottom=False)
    right.tick_params(labelleft=False)
    figure.savefig(output, dpi=160, bbox_inches="tight")
    plt.close(figure)


def draw_focal_plane_svg(
    path: Path, output: Path, bins: int, telescope: str, optical_model_path: Path | None = None
):
    """Write a weighted PSF map and measured summary without plotting dependencies."""
    from obdeect.analysis import analyse_trace_csv

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
    output.write_text("\n".join(elements) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(
        description="Plot an obdeect telescope reference outline and traced paths."
    )
    parser.add_argument("--input", type=Path, help="trace CSV (required for rays and focal plane)")
    parser.add_argument("--output", type=Path, default=Path("artificial_mst_paths.png"))
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
        ),
        default="rays",
        help="telescope outline, compiled-model geometry, recorded ray paths, or focal-plane image",
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
    args = parser.parse_args()
    if args.max_paths < 1 or args.bins < 1:
        parser.error("--max-paths and --bins must be positive")
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

    if args.focal_plane and args.output.suffix.lower() == ".svg":
        draw_focal_plane_svg(
            args.input, args.output, args.bins, args.telescope, args.optical_model_json
        )
        return

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise SystemExit(
            "Install optional visualization dependency: python -m pip install matplotlib"
        ) from error

    if args.view == "structure":
        draw_structure(plt, args.telescope, args.output)
        return
    if args.view == "telescope":
        draw_telescope(plt, args.optical_model_json, args.output, args.section_plane)
        return
    if args.view in {"pupil", "compiled-mirror"}:
        draw_pupil(plt, args.optical_model_json, args.output)
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
    )


if __name__ == "__main__":
    main()
