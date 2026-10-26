"""Build dev and tagged documentation as one GitHub Pages artifact."""

import argparse
import io
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def release_tags(tags):
    """Select final vMAJOR.MINOR.PATCH tags, newest first."""
    return sorted(
        (tag for tag in tags if re.fullmatch(r"v\d+\.\d+\.\d+", tag)),
        key=lambda tag: tuple(map(int, tag[1:].split("."))),
        reverse=True,
    )


def switcher_entries(base_url, tags):
    base_url = base_url.rstrip("/")
    entries = [{"name": "dev", "version": "dev", "url": f"{base_url}/dev/"}]
    for index, tag in enumerate(tags):
        entry = {"name": tag, "version": tag, "url": f"{base_url}/{tag}/"}
        if index == 0:
            entry.update(name=f"{tag} (stable)", preferred=True)
        entries.append(entry)
    return entries


def prepare_legacy_docs(source):
    """Give pre-Sphinx releases a root and navigation without replacing their text."""
    legacy = source / "INDEX.md"
    if "index.md" not in {path.name for path in source.iterdir()} and legacy.exists():
        text = legacy.read_text()
        legacy.unlink()
        pages = [
            name
            for name in ("USER_GUIDE", "PHYSICS", "VALIDATION")
            if (source / f"{name}.md").exists()
        ]
        (source / "index.md").write_text(
            text.rstrip() + "\n\n```{toctree}\n:hidden:\n\n" + "\n".join(pages) + "\n```\n"
        )


def publish_dev(build, site):
    """Preserve root URLs and provide an explicit dev URL; refuse stale output."""
    if site.exists():
        raise FileExistsError(f"Use a fresh output directory: {site}")
    shutil.copytree(build, site)
    shutil.copytree(build, site / "dev")


def publish_manifest(site, base_url, tags):
    if tags:
        shutil.copytree(site / tags[0], site / "stable")
    (site / "switcher.json").write_text(
        json.dumps(switcher_entries(base_url, tags), indent=2) + "\n"
    )
    (site / ".nojekyll").touch()


def build_docs(source, destination, base_url, version, package_path=None):
    env = dict(
        os.environ,
        OBDEECT_DOCS_VERSION=version,
        OBDEECT_DOCS_BASE_URL=base_url,
        OBDEECT_DOCS_SOURCE=str(source),
    )
    if package_path is not None:
        env["PYTHONPATH"] = str(package_path)
    subprocess.run(
        [
            sys.executable,
            "-m",
            "sphinx",
            "-W",
            "--keep-going",
            "-b",
            "html",
            "-c",
            str(ROOT / "docs"),
            str(source),
            str(destination),
        ],
        env=env,
        check=True,
    )


def build_site(site, base_url, dev_ref=None):
    """Rebuild releases from git, never from copies of the development pages."""
    if site.exists():
        raise FileExistsError(f"Use a fresh output directory: {site}")
    tags = release_tags(
        subprocess.check_output(["git", "tag", "--list"], cwd=ROOT, text=True).splitlines()
    )
    published = []
    with tempfile.TemporaryDirectory(prefix="obdeect-docs-") as directory:
        work = Path(directory)
        dev_source = ROOT / "docs"
        dev_package = None
        if dev_ref is not None:
            checkout = export_ref(dev_ref, work / "dev-source")
            dev_source = checkout / "docs"
            dev_package = install_package(checkout, work / "dev-package")
        build_docs(dev_source, work / "dev-html", base_url, "dev", dev_package)
        publish_dev(work / "dev-html", site)
        for tag in tags:
            checkout = export_ref(tag, work / tag)
            source = checkout / "docs"
            if not any((source / name).exists() for name in ("index.md", "INDEX.md")):
                print(f"Omitting {tag}: no public documentation index", flush=True)
                continue
            prepare_legacy_docs(source)
            package = install_package(checkout, work / f"{tag}-package", tag.removeprefix("v"))
            build_docs(source, site / tag, base_url, tag, package)
            published.append(tag)
        publish_manifest(site, base_url, published)


def export_ref(ref, destination):
    archive = subprocess.check_output(["git", "archive", ref], cwd=ROOT)
    with tarfile.open(fileobj=io.BytesIO(archive)) as tar:
        tar.extractall(destination, filter="data")
    return destination


def install_package(checkout, destination, version=None):
    # --target keeps historical APIs separate from the development installation.
    env = dict(os.environ)
    if version is not None:
        env["SETUPTOOLS_SCM_PRETEND_VERSION"] = version
    subprocess.run(
        [
            sys.executable,
            "-m",
            "pip",
            "install",
            "--no-deps",
            "--target",
            str(destination),
            str(checkout),
        ],
        check=True,
        env=env,
    )
    return destination


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--base-url", required=True)
    parser.add_argument("--dev-ref", help="Git ref for dev pages (default: working tree)")
    args = parser.parse_args()
    build_site(args.output.resolve(), args.base_url, args.dev_ref)
