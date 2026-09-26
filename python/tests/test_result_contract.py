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
