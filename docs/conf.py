"""Sphinx configuration of the obdeect documentation (layout as in simtools)."""

import os
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
    "sphinx_design",
]
source_suffix = {".md": "markdown"}
root_doc = "index"
exclude_patterns = ["_build", "changes", "figures"]
myst_heading_anchors = 3
myst_enable_extensions = ["colon_fence", "dollarmath"]
nitpicky = False
autodoc_typehints = "description"
autodoc_member_order = "bysource"

html_theme = "pydata_sphinx_theme"
html_title = project
html_logo = "_static/images/obdeect_21ratio_small.jpg"
html_css_files = ["obdeect.css"]
DOCS_SOURCE = Path(os.environ.get("OBDEECT_DOCS_SOURCE", Path(__file__).parent))
DOCS_VERSION = os.environ.get("OBDEECT_DOCS_VERSION", "dev")
DOCS_BASE_URL = os.environ.get("OBDEECT_DOCS_BASE_URL", "https://gammasim.github.io/obdeect")
html_static_path = [str(DOCS_SOURCE / "_static")]
# The current theme/CSS also applies to releases predating the documentation site.
html_static_path.append(str(Path(__file__).parent / "_static"))
if DOCS_VERSION != "dev":
    release = version = DOCS_VERSION.removeprefix("v")
html_theme_options = {
    "navbar_start": ["navbar-logo", "version-switcher"],
    "switcher": {
        "json_url": f"{DOCS_BASE_URL.rstrip('/')}/switcher.json",
        "version_match": DOCS_VERSION,
    },
    # The manifest is generated and tested with the site, before it is published.
    "check_switcher": False,
    "navbar_end": ["navbar-icon-links"],
    "show_toc_level": 1,
    "navigation_with_keys": False,
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/gammasim/obdeect",
            "icon": "fa-brands fa-github",
        },
    ],
    "logo": {"alt_text": "obdeect"},
}
