"""Sphinx configuration of the obdeect documentation (layout as in simtools)."""

import tomllib
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]

project = "obdeect"
copyright = "2026, obdeect developers"  # noqa: A001
with open(ROOT / "CITATION.cff", encoding="utf-8") as file:
    author = ", ".join(
        f"{a['given-names']} {a['family-names']}" for a in yaml.safe_load(file)["authors"]
    )
with open(ROOT / "pyproject.toml", "rb") as file:
    distribution = tomllib.load(file)["project"]["name"]
try:
    release = version(distribution)
except PackageNotFoundError:
    release = "unknown"
version = release

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.githubpages",
]
source_suffix = {".md": "markdown"}
root_doc = "index"
exclude_patterns = ["_build", "changes", "figures"]
myst_heading_anchors = 3
myst_enable_extensions = ["dollarmath"]
nitpicky = False
autodoc_typehints = "description"
autodoc_member_order = "bysource"

html_theme = "pydata_sphinx_theme"
html_title = project
html_logo = "_static/images/obdeect_small.jpg"
html_static_path = ["_static"]
html_theme_options = {
    "navbar_start": ["navbar-logo"],
    "navbar_end": ["navbar-icon-links"],
    "show_toc_level": 1,
    "navigation_with_keys": False,
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/GernotMaier/obdeect",
            "icon": "fa-brands fa-github",
        },
    ],
    "logo": {"alt_text": "obdeect"},
}
