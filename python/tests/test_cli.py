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
