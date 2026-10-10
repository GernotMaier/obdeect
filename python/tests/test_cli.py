from pathlib import Path

import obdeect
import pytest
from obdeect import cli


def test_executable_path_resolves_packaged_native_file(tmp_path: Path, monkeypatch) -> None:
    native = tmp_path / "_native"
    native.mkdir()
    executable = native / "obdeect-simtools-raytrace"
    executable.write_text("#!/bin/sh\n", encoding="utf-8")
    executable.chmod(0o755)
    monkeypatch.setattr(obdeect, "__path__", [str(tmp_path)])

    assert cli.executable_path() == executable


def test_executable_path_reports_missing_file(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(obdeect, "__path__", [str(tmp_path)])

    with pytest.raises(FileNotFoundError, match="Reinstall obdeect-dev"):
        cli.executable_path()


def test_executable_path_resolves_split_editable_install(tmp_path: Path, monkeypatch) -> None:
    source = tmp_path / "source" / "obdeect"
    installed = tmp_path / "site-packages" / "obdeect"
    source.mkdir(parents=True)
    native = installed / "_native"
    native.mkdir(parents=True)
    executable = native / "obdeect-simtools-raytrace"
    executable.write_text("#!/bin/sh\n", encoding="utf-8")
    executable.chmod(0o755)
    monkeypatch.setattr(cli, "__file__", str(source / "cli.py"))
    monkeypatch.setattr(obdeect, "__path__", [str(source), str(installed)])

    assert cli.executable_path() == executable


def test_executable_path_rejects_nonexecutable_file(tmp_path: Path, monkeypatch) -> None:
    if cli.os.name == "nt":
        pytest.skip("Windows does not use Unix executable permissions")
    native = tmp_path / "_native"
    native.mkdir()
    executable = native / "obdeect-simtools-raytrace"
    executable.write_text("#!/bin/sh\n", encoding="utf-8")
    executable.chmod(0o644)
    monkeypatch.setattr(obdeect, "__path__", [str(tmp_path)])

    with pytest.raises(PermissionError, match="not executable"):
        cli.executable_path()


def test_ground_replay_preserves_identity_and_transforms_at_boundary(tmp_path, monkeypatch):
    import csv
    import json

    model = tmp_path / "model.json"
    model.write_text(
        json.dumps({
            "telescope_frame": dict(
                axis_origin_m=[1, 2, 3],
                axes_offsets_m=[0, 0],
                known_error_deg=0,
                unknown_error_deg=0,
            )
        })
    )
    photons = tmp_path / "photons.csv"
    photons.write_text("photon_id,x_m,y_m,z_m,dx,dy,dz,time_ns,weight\n17,1,2,4,0,0,-1,12,0.7\n")
    original = photons.read_text()
    monkeypatch.setattr(
        cli.sys,
        "argv",
        [
            "obdeect",
            "--photon-input-frame",
            "ground",
            "--pointing-zenith-deg",
            "0",
            "--pointing-azimuth-deg",
            "0",
            "--optical-model",
            str(model),
            "--photon-input",
            str(photons),
            "--output",
            str(tmp_path / "out.csv"),
        ],
    )
    monkeypatch.setattr(cli, "executable_path", lambda: Path("/native"))

    def run(command):
        assert "--photon-input-frame" not in command
        prepared = Path(command[command.index("--photon-input") + 1])
        with prepared.open() as stream:
            row = next(csv.DictReader(stream))
        assert [float(row[key]) for key in ("x_m", "y_m", "z_m")] == [0, 0, 1]
        assert [row[key] for key in ("photon_id", "time_ns", "weight")] == ["17", "12", "0.7"]
        return 0

    monkeypatch.setattr(cli.subprocess, "call", run)
    assert cli.simtools_raytrace_main() == 0
    assert photons.read_text() == original


@pytest.mark.parametrize(
    "arguments,message",
    [
        (["--photon-input-frame", "invalid"], "ground or telescope"),
        (["--photon-input-frame", "ground"], "explicit pointing"),
        (["--photon-input-frame"], "missing value"),
    ],
)
def test_ground_replay_rejects_missing_boundary_context(monkeypatch, arguments, message):
    monkeypatch.setattr(cli.sys, "argv", ["obdeect", *arguments])
    with pytest.raises(ValueError, match=message):
        cli.simtools_raytrace_main()
