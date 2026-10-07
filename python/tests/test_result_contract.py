import csv
from pathlib import Path

import pytest
from obdeect.result_contract import ArrivalContractError, read_arrivals

HEADER = (
    "contract_version,photon_id,source_kind,wavelength_nm,emission_time_ns,source_weight,throughput,status,"
    "point_count,path_length_m,x0_m,y0_m,z0_m,x1_m,y1_m,z1_m\n"
)


def write_trace(tmp_path: Path, row: str) -> Path:
    path = tmp_path / "trace.csv"
    path.write_text(HEADER + row, encoding="utf-8")
    return path


def test_read_arrival_exposes_optical_boundary(tmp_path: Path) -> None:
    arrivals = read_arrivals(
        write_trace(
            tmp_path, "obdeect-arrival-v1,0,star,400,2,3,0.5,detected,2,4,0,0,1,0.1,0.2,2\n"
        )
    )

    assert len(arrivals) == 1
    assert arrivals[0].focal_x_m == 0.1
    assert arrivals[0].focal_y_m == 0.2
    assert arrivals[0].optical_weight == 1.5


def test_read_arrival_accepts_cpp_escaped_optical_model_status(tmp_path: Path) -> None:
    arrivals = read_arrivals(
        write_trace(
            tmp_path,
            "obdeect-arrival-v1,0,star,400,0,1,0,escaped_optical_model,1,1,0,0,1\n",
        )
    )

    assert arrivals[0].status == "escaped_optical_model"


def test_read_arrival_rejects_unknown_contract_version(tmp_path: Path) -> None:
    row = "obdeect-arrival-v2,0,star,400,0,1,1,detected,1,1,0,0,1\n"
    with pytest.raises(ArrivalContractError, match="contract_version"):
        read_arrivals(write_trace(tmp_path, row))


def test_read_arrival_exposes_optional_incidence_angles(tmp_path: Path) -> None:
    path = tmp_path / "angles.csv"
    path.write_text(
        HEADER.replace("path_length_m,", "path_length_m,incidence_primary_deg,incidence_focal_deg,")
        + "obdeect-arrival-v1,0,star,400,2,1,1,detected,2,4,12.5,0.3,0,0,1,0.1,0.2,2\n",
        encoding="utf-8",
    )
    arrival = read_arrivals(path)[0]
    assert arrival.incidence_primary_deg == 12.5
    assert arrival.incidence_secondary_deg is None
    assert arrival.incidence_focal_deg == 0.3


@pytest.mark.parametrize(
    "row",
    [
        "obdeect-arrival-v1,0,star,400,0,-1,1,detected,1,1,0,0,1\n",
        "obdeect-arrival-v1,0,star,400,0,1,1,detected,1,-1,0,0,1\n",
        "obdeect-arrival-v1,0,star,400,0,1,1.1,detected,1,1,0,0,1\n",
    ],
)
def test_read_arrival_rejects_invalid_optical_scalars(tmp_path: Path, row: str) -> None:
    with pytest.raises(ArrivalContractError):
        read_arrivals(write_trace(tmp_path, row))


@pytest.mark.parametrize(
    "row",
    [
        "obdeect-arrival-v1,0,moon,400,0,1,1,detected,1,1,0,0,1\n",
        "obdeect-arrival-v1,0,star,400,0,1,1,unknown,1,1,0,0,1\n",
    ],
)
def test_read_arrival_rejects_unknown_enums(tmp_path: Path, row: str) -> None:
    with pytest.raises(ArrivalContractError):
        read_arrivals(write_trace(tmp_path, row))


@pytest.mark.parametrize(
    "row,match",
    [
        ("obdeect-arrival-v1,0,star,400,0,1,1,detected,1,1,0,0,1\n", "detector vertex"),
        ("obdeect-arrival-v1,0,star,400,0,1,0.8,missed_primary,1,1,0,0,1\n", "nonzero throughput"),
        ("obdeect-arrival-v1,0,star,400,0,1,0,escaped_optical model,1,1,0,0,1\n", "invalid status"),
    ],
)
def test_read_arrival_rejects_status_invariants(tmp_path, row, match):
    with pytest.raises(ArrivalContractError, match=match):
        read_arrivals(write_trace(tmp_path, row))


def test_read_arrival_rejects_duplicate_batch_identity(tmp_path):
    row = "obdeect-arrival-v1,0,star,400,0,1,0,missed_primary,1,1,0,0,1\n"
    with pytest.raises(ArrivalContractError, match="duplicate photon_id"):
        read_arrivals(write_trace(tmp_path, row + row))


def test_zero_response_detector_arrival_is_legitimate(tmp_path):
    row = "obdeect-arrival-v1,0,star,400,0,1,0,detected,2,1,0,0,1,0,0,2\n"
    arrival = read_arrivals(write_trace(tmp_path, row))[0]
    assert arrival.detected
    assert arrival.optical_weight == 0


def test_context_identity_and_recorded_detector_time(tmp_path):
    path = tmp_path / "context.csv"
    fields = (
        ",run_id,event_id,array_id,telescope_id,bunch_id,arrival_time_ns,"
        "terminal_surface_id,final_dx,final_dy,final_dz"
    )
    header = HEADER.rstrip("\n") + fields + "\n"
    base = "obdeect-arrival-v1,7,replay,400,0,1,0,detected,2,1,0,0,1,0,0,2"
    path.write_text(header + base + ",1,2,3,4,5,99,8,0,0,1\n" + base + ",1,9,3,4,5,99,8,0,0,1\n")
    arrivals = read_arrivals(path)
    assert len(arrivals) == 2
    assert arrivals[0].arrival_time_ns == 99
    assert arrivals[0].terminal_surface_id == 8
    assert arrivals[0].final_direction == (0, 0, 1)
    with path.open("a") as handle:
        handle.write(base + ",1,2,3,4,5,99,8,0,0,1\n")
    with pytest.raises(ArrivalContractError, match="duplicate"):
        read_arrivals(path)


@pytest.mark.parametrize(
    "response,terminal,valid",
    [(0.2, 0, True), (0.1, 0, False), (-0.2, 0.4, False), (0.1, 0.1, False)],
)
def test_optical_fraction_ledger_closes_at_detector(tmp_path, response, terminal, valid):
    path = tmp_path / "ledger.csv"
    header = HEADER.rstrip("\n") + ",response_loss_fraction,terminal_loss_fraction\n"
    row = f"obdeect-arrival-v1,0,star,400,0,1,0.8,detected,2,1,0,0,1,0,0,2,{response},{terminal}\n"
    path.write_text(header + row)
    if valid:
        arrival = read_arrivals(path)[0]
        assert arrival.response_loss_fraction == 0.2
        assert arrival.terminal_loss_fraction == 0
    else:
        with pytest.raises(ArrivalContractError):
            read_arrivals(path)


def test_optical_loss_ledger_retains_response_loss_before_geometric_loss(tmp_path):
    path = tmp_path / "lost.csv"
    header = (
        HEADER.rstrip("\n")
        + ",response_loss_fraction,terminal_loss_fraction,interaction_surface_ids\n"
    )
    path.write_text(
        header + "obdeect-arrival-v1,0,star,400,0,2,0,missed_screen,2,1,0,0,1,0,0,2,0.2,0.8,3\n"
    )
    arrival = read_arrivals(path)[0]
    assert arrival.response_loss_fraction * arrival.source_weight == 0.4
    assert arrival.terminal_loss_fraction * arrival.source_weight == 1.6
    assert arrival.interaction_surface_ids == (3,)


def test_bounded_material_paths_retain_phase_and_group_observables(tmp_path):
    path = tmp_path / "material.csv"
    row = dict(
        contract_version="obdeect-arrival-v1",
        photon_id=1,
        source_kind="replay",
        wavelength_nm=400,
        emission_time_ns=0,
        source_weight=1,
        throughput=1,
        status="detected",
        point_count=65,
        path_length_m=64,
        optical_path_m=70,
        arrival_time_ns=300,
    )
    for index in range(65):
        row.update({f"x{index}_m": 0, f"y{index}_m": 0, f"z{index}_m": index})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=row)
        writer.writeheader()
        writer.writerow(row)
    arrival = read_arrivals(path)[0]
    assert len(arrival.interaction_points_m) == 65
    assert arrival.optical_path_m == 70
    assert arrival.arrival_time_ns == 300
    assert arrival.path_length_m == 64
