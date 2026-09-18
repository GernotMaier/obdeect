# obdeect

`obdeect` is an IACT-only ray-tracing prototype for blue Cherenkov photons.
The first runnable slice is deliberately small: it traces a 400-nm parallel
artificial source through one spherical, MST-inspired primary mirror, a camera
shadow and four mast supports to a focal screen. It is not yet a CTAO model or
a CORSIKA7 reader; those remain later compatibility milestones.

## Repository layout

```text
cpp/include/obdeect/  C++20 trace kernels and model types
cpp/src/              executable entry points
cpp/tests/            dependency-free C++ tests
python/obdeect/       Python diagnostic package
python/tests/         dependency-free Python tests
.github/workflows/    formatting, lint, build and test CI
```

## Build and run

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
./build/obdeect_toy --photons 100000 --output toy_mst_paths.csv
python -m pip install .
obdeect-plot-toy toy_mst_paths.csv --output toy_mst_paths.png
python -m unittest python/tests/test_plot_toy_mst.py
```

The executable needs only a C++20 compiler and the standard library. The
Python package declares Matplotlib as its only runtime dependency and installs
the `obdeect-plot-toy` command. `--no-structure` removes the camera and masts
for the analytic one-mirror focus baseline.

For contributor workflow, acceptance rules and local lint commands, see
[CONTRIBUTING.md](CONTRIBUTING.md). The project is BSD-3-Clause licensed;
citation metadata is in [CITATION.cff](CITATION.cff).

## Validation status

The analytic one-sphere path has an executable IACTrace cross-check and the
published MST central-facet `f=R/2` relation has a C++ test. Full 86-facet MST
area, PSF and timing validation is intentionally not claimed yet; see
[VALIDATION.md](VALIDATION.md) for measured residuals, literature references
and the remaining eligibility criteria.

## Coordinate and model convention

The mirror vertex is at `z=0 m`; incoming artificial Cherenkov photons begin
at `z=20 m` and travel in `-z`. The spherical centre is at `z=9.75 m`, making
the paraxial screen location `z=4.875 m`. The simple structure is four
finite-cylinder support legs plus a circular camera face evaluated before the
primary reflection. CSV records source, termination/mirror, and focal-screen
vertices so the plot displays the actual traced path.
