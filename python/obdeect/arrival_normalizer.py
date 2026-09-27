"""Convert versioned obdeect arrival records into production-comparison rows."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from obdeect.result_contract import ArrivalContractError, OpticalArrival, read_arrivals

_FIELDS = (
    "photon_id",
    "status",
    "focal_x_m",
    "focal_y_m",
    "path_length_m",
    "arrival_time_ns",
    "incidence_primary_deg",
    "incidence_secondary_deg",
    "incidence_focal_deg",
)


def comparison_row(arrival: OpticalArrival) -> dict[str, float | int | str]:
    """Return one finite comparison row from an optical-boundary arrival.

    Loss rows retain their terminal status and last interaction's transverse
    position. Their unavailable incidence values are represented as zero only
    because the production comparator checks category equality first and does
    not compare optical residuals for matched losses.
    """
    last_point = arrival.interaction_points_m[-1]
    return {
        "photon_id": arrival.photon_id,
        "status": arrival.status,
        "focal_x_m": arrival.focal_x_m if arrival.detected else last_point[0],
        "focal_y_m": arrival.focal_y_m if arrival.detected else last_point[1],
        "path_length_m": arrival.path_length_m,
        "arrival_time_ns": arrival.emission_time_ns + arrival.path_length_m / 0.299792458,
        "incidence_primary_deg": arrival.incidence_primary_deg or 0.0,
        "incidence_secondary_deg": arrival.incidence_secondary_deg or 0.0,
        "incidence_focal_deg": arrival.incidence_focal_deg or 0.0,
    }


def normalize_arrivals(input_path: Path, output_path: Path) -> int:
    """Write one normalized row per obdeect arrival and return its count."""
    arrivals = read_arrivals(input_path)
    try:
        with output_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=_FIELDS, lineterminator="\n")
            writer.writeheader()
            for arrival in arrivals:
                writer.writerow(comparison_row(arrival))
    except OSError as error:
        raise ArrivalContractError(f"cannot write {output_path}: {error}") from error
    return len(arrivals)


def main() -> None:
    """Normalize a versioned obdeect arrival CSV for production validation."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="obdeect-arrival-v1 CSV")
    parser.add_argument("--output", type=Path, required=True, help="normalized comparison CSV")
    args = parser.parse_args()
    try:
        count = normalize_arrivals(args.input, args.output)
    except ArrivalContractError as error:
        raise SystemExit(f"arrival normalization failed: {error}") from error
    print(f"Normalized {count} obdeect arrivals into {args.output}")


if __name__ == "__main__":
    main()
