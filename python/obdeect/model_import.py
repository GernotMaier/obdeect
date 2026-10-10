"""Resolve and record a selected simulation-models production for optical-model compilation."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from obdeect.parameter_roles import NON_TRANSPORT_PARAMETER_ROLES


class ImportError(ValueError):
    """An input violates the simulation-models manifest contract."""


# These are independently derived reference products, not prescriptions for
# individual ray propagation. Retain their records/assets for comparison.
REFERENCE_DIAGNOSTIC_PARAMETERS = frozenset({
    "effective_focal_length",
    "optics_properties",
    "camera_filter_photon_incident_angle",
    "primary_mirror_incidence_angle",
    "secondary_mirror_incidence_angle",
})


# Optical transport ends at the physical detector surface.  This explicit
# allow-list prevents camera electronics, trigger, gain, and calibration
# settings from leaking into an optical model.
RAY_TRACING_PARAMETERS = frozenset({
    "array_element_position_ground",
    "axes_offsets",
    "camera_body_diameter",
    "camera_body_shape",
    "camera_depth",
    "camera_filter",
    "camera_filter_photon_incident_angle",
    "camera_pixel_layout",
    "camera_pixel_types",
    "camera_pixels",
    "camera_rotate",
    "camera_transmission",
    "dish_shape_length",
    "effective_focal_length",
    "focal_length",
    "focus_offset",
    "lightguide_efficiency_vs_incidence_angle",
    "lightguide_efficiency_vs_wavelength",
    "mirror_align_random_distance",
    "mirror_align_random_horizontal",
    "mirror_align_random_vertical",
    "mirror_class",
    "mirror_degraded_reflection",
    "mirror_focal_length",
    "mirror_list",
    "mirror_offset",
    "mirror_reflection_random_angle",
    "mirror_reflectivity",
    "optics_properties",
    "parabolic_dish",
    "random_focal_length",
    "grading_of_focal_length",
    "mirror_f_scale",
    "mirror_opt",
    "flip_mirrors",
    "camera_scale_factor",
    "camera_degraded_efficiency",
    "camera_degraded_map",
    "primary_degraded_map",
    "secondary_degraded_map",
    "pixels_parallel",
    "primary_mirror_segmentation",
    "primary_mirror_degraded_map",
    "primary_mirror_diameter",
    "primary_mirror_hole_diameter",
    "primary_mirror_incidence_angle",
    "primary_mirror_parameters",
    "primary_mirror_ref_radius",
    "secondary_mirror_segmentation",
    "secondary_mirror_baffle",
    "secondary_mirror_degraded_map",
    "secondary_mirror_degraded_reflection",
    "secondary_mirror_diameter",
    "secondary_mirror_hole_diameter",
    "secondary_mirror_incidence_angle",
    "secondary_mirror_parameters",
    "secondary_mirror_ref_radius",
    "secondary_mirror_reflectivity",
    "secondary_mirror_shadow_diameter",
    "secondary_mirror_shadow_offset",
    "focal_surface_parameters",
    "focal_surface_ref_radius",
    "camera_config_file",
    "telescope_axis_height",
    "telescope_obscuration_cylinders",
    "telescope_obscuration_quadrilaterals",
    "telescope_random_angle",
    "telescope_random_error",
    "telescope_sphere_radius",
    "telescope_transmission",
})


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def component(value: str) -> bool:
    return (
        isinstance(value, str)
        and value not in ("", ".", "..")
        and "/" not in value
        and "\\" not in value
    )


def within_root(path: Path, root: Path) -> Path:
    if not path.resolve().is_relative_to(root):
        raise ImportError(f"source path escapes checkout: {path}")
    return path


def load_json(path: Path) -> dict[str, Any]:
    try:
        content = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ImportError(f"cannot read JSON {path}: {error}") from error
    if not isinstance(content, dict):
        raise ImportError(f"JSON root must be an object: {path}")
    return content


def record(path: Path, root: Path) -> dict[str, str]:
    return {"path": path.relative_to(root).as_posix(), "sha256": sha256(path)}


def resolve_model(root: Path, model: str, version: str) -> dict[str, Any]:
    """Return complete JSON-serializable production data or raise ImportError."""
    root = root.resolve()
    # Released simulation-models checkouts commonly contain the data package
    # in a nested ``simulation-models/`` directory, while exported data trees
    # use that directory itself as their root.  Accept both layouts without
    # searching arbitrary descendants.
    nested_root = root / "simulation-models"
    if not (root / "productions").is_dir() and (nested_root / "productions").is_dir():
        root = nested_root
    if not component(model) or not component(version):
        raise ImportError("model and version must be simple path components")
    records: dict[str, dict[str, str]] = {}
    resolved: dict[str, tuple[str, str]] = {}

    def inherit(instrument: str, stack: tuple[str, ...] = ()) -> None:
        if not component(instrument) or instrument in stack:
            raise ImportError("invalid or cyclic design-model inheritance")
        manifest_path = root / "productions" / version / f"{instrument}.json"
        manifest = load_json(within_root(manifest_path, root))
        if manifest.get("model_version") != version:
            raise ImportError(f"manifest version mismatch in {manifest_path}")
        if manifest.get("production_table_name") != instrument:
            raise ImportError(f"manifest table name mismatch in {manifest_path}")
        tables = manifest.get("parameters")
        if (
            not isinstance(tables, dict)
            or set(tables) != {instrument}
            or not isinstance(tables[instrument], dict)
        ):
            raise ImportError(
                "production manifest must contain exactly the requested parameter table"
            )
        key = "production_manifest" if instrument == model else f"design_manifest:{instrument}"
        records[key] = record(manifest_path, root)
        designs = manifest.get("design_model", {})
        if not isinstance(designs, dict) or set(designs) - {instrument}:
            raise ImportError("invalid design-model mapping")
        parent = designs.get(instrument)
        if parent is not None:
            inherit(parent, (*stack, instrument))
        for name, parameter_version in sorted(tables[instrument].items()):
            if not component(name) or not component(parameter_version):
                raise ImportError("parameter names and versions must be strings")
            resolved[name] = (instrument, parameter_version)
            parameter_path = (
                root / "model_parameters" / instrument / name / f"{name}-{parameter_version}.json"
            )
            records[f"source_parameter:{instrument}:{name}"] = record(
                within_root(parameter_path, root), root
            )

    inherit(model)
    if not any(key.startswith("design_manifest:") for key in records):
        records = {
            key: value for key, value in records.items() if not key.startswith("source_parameter:")
        }
    parameters: dict[str, Any] = {}
    source_parameter_coverage: dict[str, Any] = {}
    assets: dict[str, dict[str, str]] = {}
    sites = set()
    for name, (instrument, parameter_version) in sorted(resolved.items()):
        if not component(name) or not component(parameter_version):
            raise ImportError("parameter names and versions must be strings")
        parameter_path = (
            root / "model_parameters" / instrument / name / f"{name}-{parameter_version}.json"
        )
        parameter = load_json(within_root(parameter_path, root))
        if isinstance(parameter.get("site"), str):
            sites.add(parameter["site"])
        if parameter.get("instrument") != instrument or parameter.get("parameter") != name:
            raise ImportError(f"parameter identity mismatch in {parameter_path}")
        if parameter.get("parameter_version") != parameter_version:
            raise ImportError(f"parameter version mismatch in {parameter_path}")
        if not all(key in parameter for key in ("type", "value", "unit", "file")):
            raise ImportError(f"parameter record is incomplete: {parameter_path}")
        if not isinstance(parameter["file"], bool):
            raise ImportError(f"file flag must be boolean: {parameter_path}")
        outside_scope = name in NON_TRANSPORT_PARAMETER_ROLES
        source_parameter_coverage[name] = {
            "disposition": "outside_optical_scope" if outside_scope else "unsupported",
            "record": record(parameter_path, root),
            "reason": NON_TRANSPORT_PARAMETER_ROLES[name]
            if outside_scope
            else "semantics not reviewed",
        }
        if name not in RAY_TRACING_PARAMETERS:
            if parameter.get("required_for_trace") is True:
                raise ImportError(f"required optical parameter is unsupported: {name}")
            continue
        source_parameter_coverage[name]["disposition"] = "selected"
        source_parameter_coverage[name]["reason"] = "passed to optical-model compiler"
        parameters[name] = parameter
        records[f"parameter:{name}"] = record(parameter_path, root)
        if parameter["file"] and parameter["value"] is not None:
            value = parameter["value"]
            asset_name = value.partition("#rpol:")[0] if isinstance(value, str) else value
            if asset_name != value and name not in {
                "mirror_reflectivity",
                "secondary_mirror_reflectivity",
                "camera_filter",
                "lightguide_efficiency_vs_incidence_angle",
                "lightguide_efficiency_vs_wavelength",
                "camera_degraded_map",
                "primary_degraded_map",
                "secondary_degraded_map",
                "primary_mirror_degraded_map",
                "secondary_mirror_degraded_map",
            }:
                raise ImportError(f"RPOL options are unsupported for asset {name}")
            if not component(asset_name):
                raise ImportError(f"unsafe or invalid asset name in {parameter_path}")
            # Historical exports use ``model_parameters/Files`` while current
            # simulation-models records commonly keep the asset beside the
            # parameter JSON. Support only these two documented locations;
            # never search arbitrary descendants or silently choose a file.
            candidates = (
                root / "model_parameters" / "Files" / asset_name,
                parameter_path.parent / asset_name,
            )
            existing = [
                within_root(candidate, root) for candidate in candidates if candidate.is_file()
            ]
            if len(existing) != 1:
                locations = ", ".join(str(candidate) for candidate in candidates)
                raise ImportError(f"declared model asset is missing or ambiguous: {locations}")
            asset_path = existing[0]
            assets[name] = record(asset_path, root)
    environment = None
    if len(sites) > 1:
        raise ImportError("telescope parameters disagree on the site")
    if sites:
        site = next(iter(sites))
        if not component(site):
            raise ImportError("site must be a simple path component")
        instrument = f"OBS-{site}"
        manifest_path = root / "productions" / version / f"{instrument}.json"
        if manifest_path.is_file():
            from obdeect.propagation import ambient_group_index

            manifest = load_json(within_root(manifest_path, root))
            if manifest.get("model_version") != version:
                raise ImportError("site environment model version mismatch")
            table = manifest.get("parameters", {}).get(instrument, {})
            records["environment_manifest"] = record(manifest_path, root)
            selected = {}
            for name in ("atmospheric_profile", "corsika_observation_level"):
                parameter_version = table.get(name)
                if not component(parameter_version):
                    raise ImportError(f"missing site environment parameter: {name}")
                path = within_root(
                    root
                    / "model_parameters"
                    / instrument
                    / name
                    / f"{name}-{parameter_version}.json",
                    root,
                )
                selected[name] = load_json(path)
                if (
                    selected[name].get("instrument") != instrument
                    or selected[name].get("parameter") != name
                    or selected[name].get("parameter_version") != parameter_version
                ):
                    raise ImportError("site environment parameter identity mismatch")
                records[f"environment_parameter:{name}"] = record(path, root)
            profile = selected["atmospheric_profile"]
            profile_path = within_root(
                root / "model_parameters" / instrument / "atmospheric_profile" / profile["value"],
                root,
            )
            records["environment_atmospheric_profile"] = record(profile_path, root)
            level = selected["corsika_observation_level"]
            if level.get("unit") != "m" or profile.get("file") is not True:
                raise ImportError("unsupported site environment units or profile")
            try:
                index = ambient_group_index(profile_path.read_text(), float(level["value"]))
            except (ValueError, KeyError, IndexError) as error:
                raise ImportError(f"invalid site environment: {error}") from error
            environment = dict(propagation_group_index=index, observation_level_m=level["value"])
    return {
        "format": "obdeect.simulation-models-optical-model-ir.v1",
        "model": model,
        "model_version": version,
        "source_root": root.name,
        "input_records": records,
        "assets": assets,
        "parameters": parameters,
        "source_parameter_coverage": source_parameter_coverage,
        "environment": environment,
    }
