from pathlib import Path

import pytest
from obdeect import cli


def test_executable_path_resolves_packaged_native_file(tmp_path: Path, monkeypatch) -> None:
    native = tmp_path / "_native"
    native.mkdir()
    executable = native / "obdeect-simtools-raytrace"
    executable.write_text("#!/bin/sh\n", encoding="utf-8")
    executable.chmod(0o755)
    monkeypatch.setattr(cli, "__file__", str(tmp_path / "cli.py"))

    assert cli.executable_path() == executable


def test_executable_path_reports_missing_file(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(cli, "__file__", str(tmp_path / "cli.py"))

    with pytest.raises(FileNotFoundError, match="Reinstall obdeect-dev"):
        cli.executable_path()
