"""Figures of the documentation (written to docs/_static/figures).

They need inputs that are not part of this repository, given as arguments:

    python docs/figures/make_figures.py \
        --simulation-models <simulation-models checkout>/simulation-models \
        --validate-optics <results of run_psf.sh: TEL/BACKEND/results/*.ecsv> \
        [--baffle <results of run_baffle.sh: TEL/BACKEND-with|without/results/*.ecsv>]

The validation scripts are in the obdeect-validation repository.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from astropy.table import QTable  # noqa: E402
from obdeect import OpticalModel  # noqa: E402

OUTPUT = Path(__file__).resolve().parents[1] / "_static" / "figures"
# Reference palette (dataviz skill): series 1 and 2, text and grid inks.
BLUE, ORANGE = "#2a78d6", "#eb6834"
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
STRUCTURE_MODEL = ("MSTx-FlashCam", "7.0.0")

plt.rcParams.update({
    "axes.edgecolor": MUTED,
    "axes.labelcolor": INK,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.6,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "font.size": 9,
    "savefig.dpi": 150,
})


def _results(directory: Path) -> QTable | None:
    files = sorted((directory / "results").glob("*.ecsv"))
    return QTable.read(files[0]) if files else None


def validation(root: Path) -> None:
    """Deviation of obdeect from sim_telarray: D80 and effective area, all offsets."""
    telescopes, d80, area = [], [], []
    for tel in sorted(p.name for p in root.iterdir() if p.is_dir()):
        s, o = _results(root / tel / "sim_telarray"), _results(root / tel / "obdeect")
        if s is None or o is None:
            continue
        telescopes.append(tel)
        d80.append(100 * (o["psf_cm"].value / s["psf_cm"].value - 1))
        area.append(100 * (o["eff_area"].value / s["eff_area"].value - 1))
    figure, axes = plt.subplots(1, 2, figsize=(9, 4.2), sharey=True)
    for ax, values, title in zip(axes, (d80, area), ("D80", "Effective area"), strict=True):
        for row, v in enumerate(values):
            ax.plot(v, np.full(len(v), row), "o", ms=4, color=BLUE, alpha=0.8)
        ax.axvline(0, color=MUTED, lw=1)
        ax.set(title=f"{title}: obdeect / sim_telarray - 1", xlabel="deviation [%]")
        ax.set_yticks(range(len(telescopes)), telescopes)
        ax.invert_yaxis()
        limit = max(1.2, max(np.abs(np.concatenate(values))) * 1.1)
        ax.set_xlim(-limit, limit)
    figure.suptitle(
        "simtools-validate-optics, both backends: star at 10 km, zenith 20 deg, "
        "0-4 deg off axis (one dot per offset)",
        color=MUTED,
    )
    figure.tight_layout()
    figure.savefig(OUTPUT / "validation_psf_area.png")
    plt.close(figure)


def _hull(points: np.ndarray) -> np.ndarray:
    """Convex hull of 2-D points (monotone chain), closed."""
    p = sorted(map(tuple, np.round(points, 9)))
    if len(p) < 3:
        return np.array(p)

    def half(sequence):
        out = []
        for q in sequence:
            while (
                len(out) >= 2
                and (out[-1][0] - out[-2][0]) * (q[1] - out[-2][1])
                - (out[-1][1] - out[-2][1]) * (q[0] - out[-2][0])
                <= 0
            ):
                out.pop()
            out.append(q)
        return out

    lower, upper = half(p), half(reversed(p))
    ring = lower[:-1] + upper[:-1]
    return np.array([*ring, ring[0]])


def structure(simulation_models: Path) -> None:
    """The MST structure seen from the sky and from the side, by the light each element stops."""
    name, version = STRUCTURE_MODEL
    model = OpticalModel.from_simulation_models(
        simulation_models, name, version, configuration_seed=1, zenith_deg=20.0
    )
    photons = model.star_photons(1_000_000, seed=11)
    hits = model.obscurer_hits(photons, seed=12)
    body = hits["group"] != "camera"
    sole = dict(zip(hits["name"][body], 100 * hits["sole"][body] / hits["unobscured"], strict=True))
    crossed = dict(zip(hits["name"][body], hits["crossed"][body], strict=True))
    optics = model.document["optics"]
    colormap = plt.get_cmap("Blues")
    top = max(sole.values())
    figure, axes = plt.subplots(1, 2, figsize=(10, 5.2))
    for ax, (i, j), title in zip(
        axes, ((0, 1), (0, 2)), ("seen from the sky (x, y)", "side view (x, z)"), strict=True
    ):
        for facet in optics["primary"].get("facets", []):
            c = facet["centre_m"]
            ax.add_patch(
                plt.Circle(
                    (c[i], c[j]), 0.6 * facet["diameter_m"], fc="none", ec=GRID, lw=0.5, zorder=0
                )
            )
        for o in optics["obscurers"]:
            if o["group"] == "camera":
                continue
            share = sole.get(o["name"], 0.0)
            colour = colormap(0.25 + 0.75 * min(share / top, 1.0)) if share > 0 else "#c3c2b7"
            style = {"color": colour, "zorder": 1 + share, "alpha": 1.0 if share > 0 else 0.6}
            if crossed.get(o["name"], 0) == 0:
                style |= {"ls": ":"}
            if o["shape"] == "cylinder":
                a, b = o["start_m"], o["end_m"]
                ax.plot([a[i], b[i]], [a[j], b[j]], lw=1.6, solid_capstyle="round", **style)
            elif o["shape"] == "polyhedron":
                ring = _hull(np.array(o["vertices_m"])[:, [i, j]])
                ax.fill(ring[:, 0], ring[:, 1], ec=style["color"], fc=style["color"], lw=0.5)
        ax.set(aspect="equal", xlabel="x [m]", ylabel="y [m]" if j == 1 else "z [m]", title=title)
    mappable = plt.cm.ScalarMappable(cmap=colormap, norm=plt.Normalize(0, top))
    bar = figure.colorbar(mappable, ax=axes, shrink=0.7, pad=0.02)
    bar.set_label("light stopped by this element alone [% of the unobscured light]")
    never = int(sum(v == 0 for v in crossed.values()))
    figure.suptitle(
        f"{name} {version} telescope_structure, star on axis: grey = stops nothing alone, "
        f"dotted = never met ({never} of {len(crossed)} elements)",
        color=MUTED,
    )
    figure.savefig(OUTPUT / "structure_census.png", bbox_inches="tight")
    plt.close(figure)


def baffle(root: Path) -> None:
    """Effective area with and without the secondary baffle, both backends."""
    telescopes = sorted(p.name for p in root.iterdir() if p.is_dir())
    figure, axes = plt.subplots(1, len(telescopes), figsize=(4.5 * len(telescopes), 3.8))
    for ax, tel in zip(np.atleast_1d(axes), telescopes, strict=True):
        for backend, colour, marker in (("sim_telarray", BLUE, "o"), ("obdeect", ORANGE, "s")):
            for state, line in (("with", "-"), ("without", "--")):
                table = _results(root / tel / f"{backend}-{state}")
                if table is None:
                    continue
                offset = np.hypot(table["off_x"].value, table["off_y"].value)
                ax.plot(
                    offset,
                    table["eff_area"].value,
                    line,
                    marker=marker,
                    ms=4,
                    lw=1.5,
                    color=colour,
                    label=f"{backend}, {state} baffle",
                )
        ax.set(title=tel, xlabel="off-axis angle [deg]", ylabel="effective area [m$^2$]")
    np.atleast_1d(axes)[0].legend(frameon=False)
    figure.tight_layout()
    figure.savefig(OUTPUT / "baffle_area.png")
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--simulation-models", type=Path, required=True)
    parser.add_argument("--validate-optics", type=Path, required=True)
    parser.add_argument("--baffle", type=Path)
    args = parser.parse_args()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    validation(args.validate_optics)
    structure(args.simulation_models)
    if args.baffle:
        baffle(args.baffle)
    print(f"figures written to {OUTPUT}")


if __name__ == "__main__":
    main()
