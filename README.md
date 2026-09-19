# obdeect

`obdeect` is a portable C++20 ray-tracing prototype for imaging atmospheric
Cherenkov telescopes (IACTs). It provides deterministic artificial sources,
analytic optical-reference traces, CSV diagnostics, and focal-plane PSF
analysis. The C++ trace core depends only on the standard library; the Python
package supplies import, plotting, and analysis commands.

## Status

This is **not** a CTAO production simulator, a CORSIKA/EventIO reader, or a
replacement for sim_telarray or ROBAST. `obdeect_ctao` traces analytic
reference prescriptions for LST, MST, SST, and SCT; it does not yet trace the
corresponding imported telescope models. In particular, production facet
alignment, camera geometry, structures, curved focal surfaces, and scene-wide
material transport remain incomplete. See the [production review](docs/PRODUCTION_REVIEW.md),
the [replacement gate](docs/SIMTELARRAY_REPLACEMENT.md), and [TODO.md](TODO.md)
for the evidence required to make a production claim.

## Build and test

Requirements: a C++20 compiler, CMake 3.20+, Ninja, and Python 3.14+.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'

cmake --preset debug
cmake --build --preset debug
ctest --test-dir build/debug --output-on-failure
python -m pytest python/tests

# Optional contributor checks (also run in CI where available)
ruff format --check python tools examples
ruff check python tools examples
pre-commit run --all-files
```

The C++ tests cover optical kernels, source generation, input parsing, tracing,
materials, and scene contracts. The Python tests cover model import, scene
compilation, plotting, and weighted PSF analysis. CI runs these tests, the
examples, formatting/lint checks, and a native address/undefined-behaviour
sanitizer build.

## Quick start

Trace the MST-inspired reference scene and render its paths:

```bash
./build/debug/obdeect_toy --photons 10000 --output toy_paths.csv
obdeect-plot-toy toy_paths.csv --output toy_paths.png
```

Trace an analytic LST prescription:

```bash
./build/debug/obdeect_ctao --telescope LST --photons 10000 --output lst_paths.csv
obdeect-plot-toy lst_paths.csv --telescope LST --focal-plane --output lst_focal_plane.png
```

`obdeect_toy` supports deterministic `star`, `illuminator`, and `laser`
sources, field offsets, finite source distance/divergence, and a configurable
monochromatic wavelength. It is a simple spherical MST-inspired scene with a
camera shadow and four mast supports, not an MST model.

```bash
./build/debug/obdeect_toy --source star --field-x-deg 0.5 --wavelength-nm 400
./build/debug/obdeect_toy --source illuminator --distance-m 50
./build/debug/obdeect_toy --source laser --distance-m 50 --divergence-deg 0.1
```

For a tested set of runnable source and analytic-reference examples, run:

```bash
python examples/run_examples.py --build build/debug --output out/examples --photons 1000
```

## PSF and field-angle scans

`obdeect-psf` derives weighted centroid, R80/D80, optical throughput, and
terminal-loss closure from trace CSV output. D80 is twice the radius of the
centroid-centred circle containing 80% of detected optical weight.

```bash
obdeect-psf derive lst_paths.csv --field-x-deg 0 --output lst_psf.json

obdeect-psf scan --executable ./build/debug/obdeect_toy --photons 10000 \
  --field-x-deg 0 0.5 1.0 --output-dir out/toy_psf --plot out/toy_psf.png
```

A scan writes one trace CSV per offset and `psf_scan.csv`/`psf_scan.json`.
These metrics describe the selected deterministic reference configuration;
they are not validation of a production telescope model.

## Model provenance and scene compilation

The Python adapters can select a pinned `simulation-models` record, hash its
declared assets, and compile a deliberately incomplete, provenance-checked
scene hand-off. Compilation retains mirror-list centres, shapes, diameters,
and focal lengths, then reports normals, alignment, camera, structures, and
materials as trace blockers rather than guessing them. It does not connect a
production scene to a native trace executable yet.

```bash
obdeect-import-simulation-models /path/to/simulation-models LSTN-design \
  --version 6.3.0 --output lstn.ir.json
obdeect-compile-scene lstn.ir.json --source-root /path/to/simulation-models \
  --output lstn.scene.json
```

CSV trace output includes source weight, wavelength, emission time, terminal
status, path length, and path vertices. The current executable-level
throughput is one for a detector-surface hit and zero for a loss; material and
coating kernels exist but are not yet integrated into these analytic scenes.

## Project map

| Path | Purpose |
| --- | --- |
| `cpp/include/obdeect/` | Standard-library-only C++ geometry, optics, source, and trace kernels. |
| `cpp/tests/` | Native unit and contract tests. |
| `python/obdeect/` | Import, scene-compilation, plotting, and PSF commands. |
| `python/tests/` | Python unit tests. |
| `examples/` | End-to-end executable and rendering smoke tests. |
| `docs/` | Production scope and replacement evidence. |

The project is BSD-3-Clause licensed; citation metadata is in
[CITATION.cff](CITATION.cff).
