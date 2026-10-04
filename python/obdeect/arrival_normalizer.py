"""Convert versioned obdeect arrival records into production-comparison rows."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from obdeect.result_contract import ArrivalContractError, OpticalArrival, read_arrivals

_FIELDS = (
    "photon_id",
    "run_id",
    "event_id",
    "array_id",
    "telescope_id",
    "bunch_id",
    "status",
    "focal_x_m",
    "focal_y_m",
    "path_length_m",
    "arrival_time_ns",
    "incidence_primary_deg",
    "incidence_secondary_deg",
    "incidence_focal_deg",
    "wavelength_nm",
    "source_weight",
    "throughput",
    "terminal_surface_id",
    "final_dx",
    "final_dy",
    "final_dz",
    "interaction_surface_ids",
    "response_loss_fraction",
    "terminal_loss_fraction",
)


def comparison_row(arrival: OpticalArrival) -> dict[str, float | int | str]:
    """Return one finite comparison row from an optical-boundary arrival.

    Loss rows retain their terminal status and last interaction's transverse
    position. Their unavailable incidence values are represented as zero only
    because the production comparator checks category equality first and does
    not compare optical residuals for matched losses.
    """
    last_point = arrival.interaction_points_m[-1]
    direction = arrival.final_direction
    return {
        "photon_id": arrival.photon_id,
        "run_id": arrival.run_id,
        "event_id": arrival.event_id,
        "array_id": arrival.array_id,
        "telescope_id": arrival.telescope_id,
        "bunch_id": arrival.bunch_id,
        "status": arrival.status,
        "focal_x_m": arrival.focal_x_m if arrival.detected else last_point[0],
        "focal_y_m": arrival.focal_y_m if arrival.detected else last_point[1],
        "path_length_m": arrival.path_length_m,
        "arrival_time_ns": arrival.arrival_time_ns
        if arrival.arrival_time_ns is not None
        else arrival.emission_time_ns + arrival.path_length_m / 0.299792458,
        "incidence_primary_deg": arrival.incidence_primary_deg or 0.0,
        "incidence_secondary_deg": arrival.incidence_secondary_deg or 0.0,
        "incidence_focal_deg": arrival.incidence_focal_deg or 0.0,
        "wavelength_nm": arrival.wavelength_nm,
        "source_weight": arrival.source_weight,
        "throughput": arrival.throughput,
        **(
            {
                "response_loss_fraction": arrival.response_loss_fraction,
                "terminal_loss_fraction": arrival.terminal_loss_fraction,
            }
            if arrival.response_loss_fraction is not None
            else {}
        ),
        "interaction_surface_ids": ";".join(map(str, arrival.interaction_surface_ids))
        if arrival.interaction_surface_ids is not None
        else "",
        "terminal_surface_id": str(arrival.terminal_surface_id)
        if arrival.terminal_surface_id is not None
        else "unrecorded",
        **(
            {
                name: direction[index]
                for index, name in enumerate(("final_dx", "final_dy", "final_dz"))
            }
            if direction is not None
            else {}
        ),
    }


def normalize_arrivals(
    input_path: Path,
    output_path: Path,
    *,
    optical_model_sha256: str | None = None,
    source_sha256: str | None = None,
) -> int:
    """Write one normalized row per obdeect arrival and return its count."""
    for name, value in (
        ("optical_model_sha256", optical_model_sha256),
        ("source_sha256", source_sha256),
    ):
        if value is not None and (
            not isinstance(value, str)
            or len(value) != 64
            or any(character not in "0123456789abcdef" for character in value)
        ):
            raise ArrivalContractError(f"invalid {name}")
    arrivals = read_arrivals(input_path)
    if (optical_model_sha256 is not None or source_sha256 is not None) and any(
        arrival.arrival_time_ns is None for arrival in arrivals
    ):
        raise ArrivalContractError("provenance-bound normalization requires recorded arrival time")
    try:
        with output_path.open("w", newline="", encoding="utf-8") as handle:
            fields = list(_FIELDS)
            if not all(arrival.interaction_surface_ids is not None for arrival in arrivals):
                fields.remove("interaction_surface_ids")
            if not all(arrival.final_direction is not None for arrival in arrivals):
                fields = [field for field in fields if not field.startswith("final_d")]
            if not all(arrival.response_loss_fraction is not None for arrival in arrivals):
                fields = [
                    field
                    for field in fields
                    if field not in {"response_loss_fraction", "terminal_loss_fraction"}
                ]
            provenance = {
                name: value
                for name, value in (
                    ("optical_model_sha256", optical_model_sha256),
                    ("source_sha256", source_sha256),
                )
                if value is not None
            }
            writer = csv.DictWriter(handle, fieldnames=[*fields, *provenance], lineterminator="\n")
            writer.writeheader()
            for arrival in arrivals:
                writer.writerow({
                    key: value
                    for key, value in (comparison_row(arrival) | provenance).items()
                    if key in writer.fieldnames
                })
    except OSError as error:
        raise ArrivalContractError(f"cannot write {output_path}: {error}") from error
    return len(arrivals)


def main() -> None:
    """Normalize a versioned obdeect arrival CSV for production validation."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="obdeect-arrival-v1 CSV")
    parser.add_argument("--output", type=Path, required=True, help="normalized comparison CSV")
    parser.add_argument("--optical-model-sha256", help="declared compiled model content hash")
    parser.add_argument("--source-sha256", help="declared frozen input photon file hash")
    args = parser.parse_args()
    try:
        count = normalize_arrivals(
            args.input,
            args.output,
            optical_model_sha256=args.optical_model_sha256,
            source_sha256=args.source_sha256,
        )
    except ArrivalContractError as error:
        raise SystemExit(f"arrival normalization failed: {error}") from error
    print(f"Normalized {count} obdeect arrivals into {args.output}")


if __name__ == "__main__":
    main()
