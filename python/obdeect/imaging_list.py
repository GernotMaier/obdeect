"""simtools imaging-list interchange at the optical backend boundary."""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path

from obdeect.result_contract import ArrivalContractError, OpticalArrival, read_arrivals


@dataclass(frozen=True)
class ImagingListMetadata:
    """Simulation inputs required by the common imaging-list analysis."""

    sampling_radius_m: float
    entrance_z_m: float
    focal_length_m: float
    camera_rotation_deg: float
    prime_focus: bool
    optical_model_sha256: str
    source_distance_m: float
    off_axis_x_deg: float
    off_axis_y_deg: float
    zenith_angle_deg: float = 0.0
    detector_configuration_seed: int = 0


def load_imaging_metadata(
    path: Path,
    source_distance_m: float,
    off_axis_x_deg: float,
    off_axis_y_deg: float,
    zenith_angle_deg: float = 0.0,
) -> ImagingListMetadata:
    """Resolve the sim_telarray star-launch convention from an explicit optical model.

    sim_config.c::star_light samples a disk with radius 1.2 times the outer
    mirror catalogue extent. This adapter convention does not change the core.
    """
    model = json.loads(path.read_text(encoding="utf-8"))
    blockers = model.get("report", {}).get("trace_blockers", [])
    if any("panel alignment and distance" in item for item in blockers):
        raise ValueError(
            "optical model omits configured mirror effects; recompile the optical model"
        )
    alignment = model.get("primary", {}).get("alignment", {})
    alignment_zenith = alignment.get("zenith_angle_deg")
    if alignment_zenith is not None and not math.isclose(
        alignment_zenith, zenith_angle_deg, abs_tol=1e-12
    ):
        raise ValueError(
            "optical model panel alignment was compiled for a different telescope zenith"
        )
    trace = model["trace_model"]
    kind = trace["kind"]
    if kind == "segmented":
        facets = trace["primary_facets"]

        def facet_circumradius(facet: dict) -> float:
            factor = {
                "circle": 1.0,
                "square": math.sqrt(2),
                "hexagon": 2 / math.sqrt(3),
                "hexagon_flat_x": 2 / math.sqrt(3),
                "hexagon_flat_y": 2 / math.sqrt(3),
            }.get(facet["shape"])
            if factor is None:
                raise ValueError(f"unsupported primary facet shape: {facet['shape']}")
            return facet["diameter_m"] * factor / 2

        radius = max(
            math.hypot(*facet["centre_m"][:2]) + facet_circumradius(facet) for facet in facets
        )
        top = max(
            facet["centre_m"][2] + facet_circumradius(facet) * math.hypot(*facet["normal"][:2])
            for facet in facets
        )
    elif kind == "axisymmetric":
        radius = trace["primary"]["outer_radius_m"]
        top = max(trace[role]["vertex_z_m"] for role in ("primary", "secondary", "detector"))
    else:
        raise ValueError("imaging lists require segmented or axisymmetric optical models")
    for surface in trace.get("detector_surfaces", []):
        top = max(top, surface["centre_m"][2] + surface["diameter_m"] / 2)
    for field in ("cylinder_obscurers", "primary_to_secondary_cylinders"):
        for cylinder in trace.get(field, []):
            top = max(
                top,
                cylinder["first_endpoint_m"][2] + cylinder["diameter_m"] / 2,
                cylinder["second_endpoint_m"][2] + cylinder["diameter_m"] / 2,
            )
    for field in ("incoming_obscurer_planes", "primary_to_secondary_planes"):
        for plane in trace.get(field, []):
            top = max(top, plane["centre_m"][2] + plane["diameter_m"] / 2)
    for surface in trace.get("opaque_obscurers", []):
        if surface["shape"] == "quadrilateral":
            top = max(top, *(point[2] for point in surface["vertices_m"]))
        elif surface["shape"] in ("hollow_frustum", "solid_frustum"):
            top = max(
                top,
                surface["first_endpoint_m"][2] + surface["first_radius_m"] + surface["thickness_m"],
                surface["second_endpoint_m"][2]
                + surface["second_radius_m"]
                + surface["thickness_m"],
            )
    values = (
        radius,
        top,
        model["focal_length_m"],
        model["camera"]["rotation_deg"],
        source_distance_m,
        off_axis_x_deg,
        off_axis_y_deg,
        zenith_angle_deg,
    )
    if not all(math.isfinite(value) for value in values):
        raise ValueError("imaging-list metadata must be finite")
    if radius <= 0 or model["focal_length_m"] <= 0 or source_distance_m <= max(0.0, top) + radius:
        raise ValueError("invalid imaging-list sampling radius, focal length, or source distance")
    return ImagingListMetadata(
        1.2 * radius,
        max(0.0, top) + radius,
        model["focal_length_m"],
        model["camera"]["rotation_deg"],
        kind == "segmented",
        model["optical_model_sha256"],
        source_distance_m,
        off_axis_x_deg,
        off_axis_y_deg,
        zenith_angle_deg,
        model.get("random_seeds", {}).get("detector_configuration_seed", 0),
    )


_COLUMNS = (
    "Telescope ID",
    "Pixel number (-1: pixel acceptance not evaluated)",
    "X coordinate in camera imaging frame [cm]",
    "Y coordinate in camera imaging frame [cm]",
    "X initial position in telescope frame [cm]",
    "Y initial position in telescope frame [cm]",
    "Z initial position in telescope frame [cm]",
    "Arrival time at imaging surface [ns]",
    "Primary surface ID",
    "X coordinate at imaging surface [cm]",
    "Y coordinate at imaging surface [cm]",
    "x direction slope at imaging surface",
    "y direction slope at imaging surface",
    "Angle to telescope axis at imaging surface [deg]",
    "Star number",
    "X coordinate at 2f distance (unavailable)",
    "Y coordinate at 2f distance (unavailable)",
    "Relative optical efficiency factor (before pixel)",
    "Efficiency including pixel (unavailable)",
    "x direction slope at pixel (unavailable)",
    "y direction slope at pixel (unavailable)",
    "Angle to pixel normal (unavailable)",
    "Arrival time at pixel (unavailable)",
    "Angle of incidence at focal surface, w.r.t. optical axis [deg]",
    "Angle of incidence onto focal surface, w.r.t. normal [deg]",
    "X at reflection point on primary mirror [cm]",
    "Y at reflection point on primary mirror [cm]",
    "Z at reflection point on primary mirror [cm]",
    "Angle of incidence onto primary mirror [deg]",
    "X at reflection point on secondary mirror [cm]",
    "Y at reflection point on secondary mirror [cm]",
    "Z at reflection point on secondary mirror [cm]",
    "Angle of incidence onto secondary mirror [deg]",
)


def _imaging_row(arrival: OpticalArrival, metadata: ImagingListMetadata) -> list[float | int]:
    """Convert one geometric arrival without weighting or applying pixel acceptance."""
    expected_points = 3 if metadata.prime_focus else 4
    if len(arrival.interaction_points_m) != expected_points:
        raise ArrivalContractError("imaging arrival has an unsupported optical path")
    direction = arrival.final_direction
    if direction is None or direction[2] == 0:
        raise ArrivalContractError("imaging arrival requires a final direction")
    x, y = arrival.focal_x_m, arrival.focal_y_m
    if x is None or y is None:
        raise ArrivalContractError("imaging arrival requires focal coordinates")
    theta = math.radians(metadata.camera_rotation_deg)
    # The common reader rotates camera coordinates back into the analysis frame.
    cx = x * math.cos(theta) + y * math.sin(theta)
    cy = -x * math.sin(theta) + y * math.cos(theta)
    sx = direction[0] / direction[2]
    sy = direction[1] / direction[2]
    axis_angle = math.degrees(math.atan(math.hypot(sx, sy)))
    primary = arrival.interaction_points_m[1]
    secondary = arrival.interaction_points_m[2] if not metadata.prime_focus else (math.nan,) * 3
    ids = arrival.interaction_surface_ids
    return [
        arrival.telescope_id,
        -1,
        cx * 100,
        cy * 100,
        *(value * 100 for value in arrival.interaction_points_m[0]),
        arrival.arrival_time_ns if arrival.arrival_time_ns is not None else math.nan,
        ids[0] if ids else -1,
        cx * 100,
        cy * 100,
        sx,
        sy,
        axis_angle,
        1,
        math.nan,
        math.nan,
        arrival.throughput,
        *(math.nan for _ in range(5)),
        axis_angle,
        arrival.incidence_focal_deg if arrival.incidence_focal_deg is not None else math.nan,
        *(value * 100 for value in primary),
        arrival.incidence_primary_deg if arrival.incidence_primary_deg is not None else math.nan,
        *(value * 100 for value in secondary),
        arrival.incidence_secondary_deg if not metadata.prime_focus else math.nan,
    ]


def write_imaging_list(
    input_path: Path,
    output_path: Path,
    metadata: ImagingListMetadata,
    emitted_photons: int,
    *,
    ray_tracing_seed: int = 0,
) -> int:
    """Write the shared imaging-list contract; loss rows count only in the header."""
    if input_path.resolve() == output_path.resolve() or (
        output_path.exists() and input_path.samefile(output_path)
    ):
        raise ValueError("imaging-list output must not overwrite native arrivals")
    arrivals = read_arrivals(input_path)
    if emitted_photons <= 0 or len(arrivals) != emitted_photons:
        raise ArrivalContractError("native arrival count disagrees with emitted_photons")
    if (
        isinstance(ray_tracing_seed, bool)
        or not isinstance(ray_tracing_seed, int)
        or not 0 <= ray_tracing_seed < 2**64
    ):
        raise ArrivalContractError("ray-tracing seed must be uint64")
    if any(arrival.source_kind != "star" or arrival.source_weight != 1 for arrival in arrivals):
        raise ArrivalContractError("imaging lists require unit-weight star photons")
    area = math.pi * metadata.sampling_radius_m**2
    if not math.isfinite(area) or metadata.sampling_radius_m <= 0:
        raise ArrivalContractError("imaging sampling area must be finite and positive")
    for arrival in arrivals:
        if arrival.sampling_area_m2 is not None and not math.isclose(
            arrival.sampling_area_m2, area, rel_tol=1e-12
        ):
            raise ArrivalContractError("native sampling area disagrees with imaging-list metadata")
    count = 0
    with output_path.open("w", encoding="utf-8") as output:
        output.write(
            "# Optical imaging list\n"
            "# imaging_list_contract = simtools-imaging-list-v1\n"
            "# backend = obdeect\n"
            "# source = finite_distance_star\n"
            "# source_sampling = keyed_uniform_disk\n"
            f"# ray_tracing_seed = {ray_tracing_seed}\n"
            f"# detector_configuration_seed = {metadata.detector_configuration_seed}\n"
            f"# optical_model_sha256 = {metadata.optical_model_sha256}\n"
            f"# Focal_length = {metadata.focal_length_m * 100:.17g} cm\n"
            f"# Camera rotation angle = {metadata.camera_rotation_deg:.17g} deg\n"
            f"# source_distance [km] = {metadata.source_distance_m / 1000:.17g}\n"
            f"# off_axis_x [deg] = {metadata.off_axis_x_deg:.17g}\n"
            f"# off_axis_y [deg] = {metadata.off_axis_y_deg:.17g}\n"
            f"# zenith_angle [deg] = {metadata.zenith_angle_deg:.17g}\n"
            f"# sampling_radius [m] = {metadata.sampling_radius_m:.17g}\n"
            f"# launch_plane_z [m] = {metadata.entrance_z_m:.17g}\n"
            "# boundary = continuous_focal_surface; pixel acceptance not applied\n"
            "# unavailable physical observables are represented by nan\n"
            f"# Telescope 0 with {emitted_photons} photons from 1 star(s) "
            f"falling on an area of {area:.17g} m^2\n"
        )
        wavelengths = sorted({arrival.wavelength_nm for arrival in arrivals})
        output.write(
            "# wavelengths [nm] = " + ",".join(f"{value:.17g}" for value in wavelengths) + "\n"
        )
        for index, description in enumerate(_COLUMNS, start=1):
            output.write(f"# Column {index}: {description}\n")
        for arrival in arrivals:
            if arrival.detected:
                output.write(" ".join(f"{value:.17g}" for value in _imaging_row(arrival, metadata)))
                output.write("\n")
                count += 1
    return count
