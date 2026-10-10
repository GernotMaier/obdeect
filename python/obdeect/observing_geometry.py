"""Resolve an immutable optical model for one telescope pointing at the boundary."""

from __future__ import annotations

import copy
import hashlib
import json
import math
from pathlib import Path

from obdeect.camera_surfaces import compile_camera_surfaces
from obdeect.optical_model_compiler import (
    OpticalModelCompileError,
    apply_panel_alignment,
    build_plot_geometry,
    build_trace_model,
    derive_nominal_single_reflector,
    resolve_focus_offset,
)
from obdeect.pixel_response_compiler import bind_pixel_responses


def resolve_observing_geometry(model: dict, zenith_deg: float) -> dict:
    """Re-evaluate deformation without accumulating rotations or re-drawing errors."""
    if not math.isfinite(zenith_deg) or not -180 <= zenith_deg <= 180:
        raise OpticalModelCompileError("telescope zenith must be in [-180, 180] degrees")
    context = model.get("observing_geometry")
    if context is None:
        alignment = model.get("primary", {}).get("alignment", {}).get("zenith_angle_deg")
        if alignment is not None and not math.isclose(alignment, zenith_deg, abs_tol=1e-12):
            raise OpticalModelCompileError("recompile optical model to retain observing geometry")
        return model
    if context["zenith_angle_deg"] == zenith_deg:
        return model
    result = copy.deepcopy(model)
    parameters = context["parameters"]
    primary, camera = result["primary"], result.get("camera")
    prime = result["trace_model"]["kind"] == "segmented"
    if prime:
        derive_nominal_single_reflector(primary["facets"], parameters)
        apply_panel_alignment(
            primary["facets"], parameters, context["detector_configuration_seed"], zenith_deg
        )
        primary["alignment"]["zenith_angle_deg"] = zenith_deg
    focus = parameters.get("focus_offset")
    if focus is not None and camera is not None and camera.get("entrance_surfaces"):
        offset = resolve_focus_offset(focus, zenith_deg)
        if prime:
            units = {"m": 1, "cm": 0.01, "mm": 0.001}
            overall = parameters["mirror_offset"]
            focal = {
                "coefficient_m": [
                    result["focal_length_m"] + offset - overall["value"] * units[overall["unit"]]
                ],
                "radial_scale_m": 1,
            }
            result["detector_vertex_z_m"] = focal["coefficient_m"][0]
            mode = 1
        else:
            focal = copy.deepcopy(context["nominal_focal_surface"])
            focal["coefficient_m"][0] -= offset
            result["focal_surface"] = focal
            mode = parameters["pixels_parallel"]["value"]
        camera.update(compile_camera_surfaces(camera, focal, mode, reflected=prime))
    trace = build_trace_model(result)
    for name in (
        "telescope_transmission",
        "propagation_group_index",
        "opaque_obscurers",
        "primary_degradation",
        "secondary_degradation",
        "camera_degradation",
        "camera_degradation_in_detector_frame",
        "primary_degradation_in_facet_frame",
    ):
        if name in result["trace_model"]:
            trace[name] = result["trace_model"][name]
    if camera is not None and camera.get("pixel_optical_response_tables"):
        trace["pixel_responses"] = bind_pixel_responses(camera, trace)
    result["trace_model"] = trace
    result["observing_geometry"]["zenith_angle_deg"] = zenith_deg
    result["report"]["production_trace_ready"] = False
    result["plot_geometry"] = build_plot_geometry(trace, result["report"])
    result.pop("optical_model_sha256", None)
    result["optical_model_sha256"] = hashlib.sha256(
        json.dumps(result, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    return result


def prepare_observing_model(source: Path, output: Path, zenith_deg: float) -> Path:
    """Write a provenance-bound model only when the observing geometry changes."""
    model = json.loads(source.read_text())
    resolved = resolve_observing_geometry(model, zenith_deg)
    if resolved is model:
        return source
    if source.resolve() == output.resolve():
        raise OpticalModelCompileError(
            "observing output must not overwrite the input optical model"
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(resolved, indent=2) + "\n")
    return output
