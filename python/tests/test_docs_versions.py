"""Version selection and publishing layout for the documentation site."""

import importlib.util
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "docs" / "build_versions.py"
spec = importlib.util.spec_from_file_location("build_versions", SCRIPT)
versions = importlib.util.module_from_spec(spec)
spec.loader.exec_module(versions)


def test_release_tags_are_stable_and_sorted_numerically():
    assert versions.release_tags([
        "v0.9.0",
        "v0.10.0",
        "v1.0.0rc1",
        "v1.0.0-prep",
        "other",
        "v0.2.1",
    ]) == ["v0.10.0", "v0.9.0", "v0.2.1"]


def test_manifest_marks_only_the_latest_release_as_stable():
    entries = versions.switcher_entries("https://example.org/project/", ["v0.10.0", "v0.9.0"])
    assert entries == [
        {"name": "dev", "version": "dev", "url": "https://example.org/project/dev/"},
        {
            "name": "v0.10.0 (stable)",
            "version": "v0.10.0",
            "url": "https://example.org/project/v0.10.0/",
            "preferred": True,
        },
        {"name": "v0.9.0", "version": "v0.9.0", "url": "https://example.org/project/v0.9.0/"},
    ]
    assert versions.switcher_entries("https://example.org/project", []) == entries[:1]


def test_legacy_index_keeps_historical_content_and_adds_navigation(tmp_path):
    (tmp_path / "INDEX.md").write_text("# Historical docs\n\nRelease-specific text.\n")
    for name in ["USER_GUIDE", "PHYSICS", "VALIDATION"]:
        (tmp_path / f"{name}.md").write_text(f"# {name}\n")
    versions.prepare_legacy_docs(tmp_path)
    index = (tmp_path / "index.md").read_text()
    assert "Release-specific text." in index
    assert "USER_GUIDE\nPHYSICS\nVALIDATION" in index
    assert "INDEX.md" not in {path.name for path in tmp_path.iterdir()}


def test_publish_preserves_root_pages_and_provides_version_aliases(tmp_path):
    dev = tmp_path / "build"
    dev.mkdir()
    (dev / "index.html").write_text("development")
    (dev / "usage.html").write_text("existing link")
    site = tmp_path / "site"
    versions.publish_dev(dev, site)
    release = site / "v0.3.0"
    release.mkdir()
    (release / "index.html").write_text("release")
    versions.publish_manifest(site, "https://example.org/project", ["v0.3.0"])
    assert (site / "usage.html").read_text() == "existing link"
    assert (site / "dev" / "index.html").read_text() == "development"
    assert (site / "stable" / "index.html").read_text() == "release"
    assert (site / "switcher.json").exists()


def test_refuse_to_publish_into_an_existing_directory(tmp_path):
    with pytest.raises(FileExistsError):
        versions.publish_dev(tmp_path / "build", tmp_path)
