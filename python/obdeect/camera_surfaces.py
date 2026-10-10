"""Compile explicit pixel entrance and cathode planes in the telescope frame."""

import math
from collections import defaultdict
from typing import Any

from obdeect.camera_config import PIXEL_APERTURE_SHAPES, CameraConfigError


def compile_camera_surfaces(
    camera: dict[str, Any], focal: dict[str, Any], orientation_mode: int, *, reflected: bool
) -> dict[str, list[dict[str, Any]]]:
    """Follow simtel camera_setup_inclined_pixels, retaining finite physical planes.

    Modes 0/1 place pixels individually; modes 2/3 place their module fronts in
    one plane. Modes 1/3 keep all normals parallel to the optical axis. A prime
    focus camera reverses z and preserves both transverse coordinates.
    The focal polynomial already includes its telescope-frame vertex placement.
    """
    if isinstance(orientation_mode, bool) or orientation_mode not in (0, 1, 2, 3):
        raise CameraConfigError("unsupported pixel orientation mode")
    types = {entry["id"]: entry for entry in camera["pixel_types"]}
    pixels = camera["pixels"]
    required = {"z_offset_m", "rotation_deg", "normal_slopes", "module", "enabled"}
    if any(not required.issubset(pixel) for pixel in pixels):
        raise CameraConfigError("physical pixel placement requires offsets, normals and modules")
    coefficients = focal["coefficient_m"]
    scale = focal["radial_scale_m"]

    def sag(x: float, y: float) -> float:
        q = math.hypot(x, y) / scale
        return sum(coefficient * q ** (2 * index) for index, coefficient in enumerate(coefficients))

    def slopes(x: float, y: float) -> tuple[float, float]:
        radius = math.hypot(x, y)
        if radius == 0:
            return 0.0, 0.0
        q = radius / scale
        derivative = sum(
            2 * index * coefficient * q ** (2 * index - 1) / scale
            for index, coefficient in enumerate(coefficients)
            if index
        )
        return -derivative * x / radius, -derivative * y / radius

    modules = defaultdict(list)
    for pixel in pixels:
        modules[pixel["module"]].append(pixel)
    module_planes = {}
    if orientation_mode in (2, 3):
        for identifier, members in modules.items():
            x = sum(pixel["centre_xy_m"][0] for pixel in members) / len(members)
            y = sum(pixel["centre_xy_m"][1] for pixel in members) / len(members)
            z = sum(sag(*pixel["centre_xy_m"]) for pixel in members) / len(members)
            nx, ny = slopes(x, y) if orientation_mode == 2 else (0.0, 0.0)
            module_planes[identifier] = nx, ny, x * nx + y * ny + z
    entrances, cathodes = [], []
    for pixel in pixels:
        if not pixel["enabled"]:
            continue
        x, y = pixel["centre_xy_m"]
        nx, ny = pixel["normal_slopes"]
        rotation = math.radians(pixel["rotation_deg"])
        if orientation_mode in (1, 3) and (nx != 0 or ny != 0 or rotation != 0):
            raise CameraConfigError("manual pixel alignment conflicts with parallel orientation")
        if orientation_mode in (2, 3):
            nx, ny, plane = module_planes[pixel["module"]]
            z = plane - x * nx - y * ny
        else:
            z = sag(x, y)
            if orientation_mode == 0 and nx == 0 and ny == 0:
                nx, ny = slopes(x, y)
        z += (-1 if reflected else 1) * pixel["z_offset_m"]
        norm_xy = math.hypot(nx, ny)
        norm = math.sqrt(1 + norm_xy * norm_xy)
        normal = [nx / norm, ny / norm, 1 / norm]
        if norm_xy:
            squared = norm_xy * norm_xy
            u = [
                nx * nx / (squared * norm) + ny * ny / squared,
                nx * ny / (squared * norm) - nx * ny / squared,
                -nx / norm,
            ]
            v = [u[1], ny * ny / (squared * norm) + nx * nx / squared, -ny / norm]
        else:
            u, v = [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]
        tangent = [
            math.cos(rotation) * a + math.sin(rotation) * b for a, b in zip(u, v, strict=True)
        ]
        if reflected:
            # The prime-focus tracer reverses camera z and explicitly restores
            # camera x before pixel assignment. Preserve both transverse axes.
            normal = [normal[0], normal[1], -normal[2]]
            tangent = [tangent[0], tangent[1], -tangent[2]]
        camera_angle = math.radians(camera.get("rotation_deg", 0.0))
        cosine, sine = math.cos(camera_angle), math.sin(camera_angle)
        x, y = cosine * x - sine * y, sine * x + cosine * y
        normal = [
            cosine * normal[0] - sine * normal[1],
            sine * normal[0] + cosine * normal[1],
            normal[2],
        ]
        tangent = [
            cosine * tangent[0] - sine * tangent[1],
            sine * tangent[0] + cosine * tangent[1],
            tangent[2],
        ]
        pixel_type = types[pixel["type_id"]]
        centre = [x, y, z]
        metadata = {
            "source_pixel_id": pixel["id"],
            "source_type_id": pixel["type_id"],
            "enabled": pixel["enabled"],
            "response_y_sign": -1 if reflected else 1,
            "assignment_radius_m": pixel_type["funnel_diameter_m"]
            * 0.5
            * (2 / math.sqrt(3) if pixel_type["funnel_shape_code"] in (1, 3) else 1),
            "normal": normal,
            "tangent": tangent,
        }
        entrances.append({
            **metadata,
            "centre_m": centre,
            "shape": PIXEL_APERTURE_SHAPES[pixel_type["funnel_shape_code"]],
            "diameter_m": pixel_type["funnel_diameter_m"],
        })
        cathodes.append({
            **metadata,
            "centre_m": [
                value - pixel_type["funnel_depth_m"] * component
                for value, component in zip(centre, normal, strict=True)
            ],
            "shape": PIXEL_APERTURE_SHAPES[pixel_type["cathode_shape_code"]],
            "diameter_m": pixel_type["cathode_diameter_m"],
        })
    if not entrances:
        raise CameraConfigError("physical camera has no enabled pixel entrances")
    radius_by_type = {
        identifier: entry["funnel_diameter_m"]
        * 0.5
        * (2 / math.sqrt(3) if entry["funnel_shape_code"] in (1, 3) else 1)
        for identifier, entry in types.items()
    }
    x_low = min(
        0, *(pixel["centre_xy_m"][0] - radius_by_type[pixel["type_id"]] for pixel in pixels)
    )
    x_high = max(
        0, *(pixel["centre_xy_m"][0] + radius_by_type[pixel["type_id"]] for pixel in pixels)
    )
    y_low = min(
        0, *(pixel["centre_xy_m"][1] - radius_by_type[pixel["type_id"]] for pixel in pixels)
    )
    y_high = max(
        0, *(pixel["centre_xy_m"][1] + radius_by_type[pixel["type_id"]] for pixel in pixels)
    )
    angle = math.radians(camera.get("rotation_deg", 0))
    assignment = dict(
        nx=int(math.sqrt(4 * len(pixels)) + 2),
        ny=int(math.sqrt(4 * len(pixels)) + 2),
        x_low_m=x_low,
        x_high_m=x_high,
        y_low_m=y_low,
        y_high_m=y_high,
        x_basis=[math.cos(angle), math.sin(angle), 0.0],
        y_basis=[-math.sin(angle), math.cos(angle), 0.0],
    )
    if reflected:
        assignment["reference_plane_z_m"] = coefficients[0]
    return {
        "entrance_surfaces": entrances,
        "cathode_surfaces": cathodes,
        "assignment_grid": assignment,
    }
