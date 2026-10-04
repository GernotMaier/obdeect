"""Require a new Towncrier fragment for changes submitted in a pull request."""

from __future__ import annotations

import argparse
import re
import subprocess
from pathlib import Path

_FRAGMENT = re.compile(r"docs/changes/[a-zA-Z0-9_-]+\.(feature|bugfix)\.md$")


def check(base: str) -> None:
    """Reject pull requests with no nonempty newly added release note."""
    files = subprocess.check_output(
        ["git", "diff", "--name-only", "--diff-filter=A", f"{base}...HEAD"], text=True
    ).splitlines()
    if not any(_FRAGMENT.fullmatch(name) and Path(name).read_text().strip() for name in files):
        raise ValueError(
            "Add a nonempty docs/changes/<issue>.feature.md or <issue>.bugfix.md "
            "Towncrier fragment describing the change for gamma-ray astronomers."
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", required=True, help="Pull request base commit")
    args = parser.parse_args()
    try:
        check(args.base)
    except (ValueError, subprocess.CalledProcessError) as error:
        raise SystemExit(str(error)) from error


if __name__ == "__main__":
    main()
