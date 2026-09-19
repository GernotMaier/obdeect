# obdeect

`obdeect` is a portable C++20 ray-tracing prototype for imaging atmospheric
Cherenkov telescopes (IACTs). It provides deterministic artificial sources,
analytic optical-reference traces, CSV diagnostics, and focal-plane PSF
analysis. The C++ trace core depends only on the standard library; the Python
package supplies import, plotting, and analysis commands.

## Build and test

Requirements: a C++20 compiler, CMake 3.20+, Ninja, and Python 3.14+.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'

Compile and run the C++ executables:

```bash
# Configure the build system
cmake -S . -B build
# Build the project
cmake --build build
# Run the tests and executables
ctest --test-dir build --output-on-failure
# Run the Python unit tests
python -m unittest discover -s python/tests
```

## Plotting and testing

MST (Medium-Sized Telescope) example run:

```bash
./build/obdeect_toy --photons 100000 --output toy_mst_paths.csv
obdeect-plot-toy toy_mst_paths.csv --output toy_mst_paths.png
```

LST (Large-Sized Telescope) example run:

```bash
./build/obdeect_ctao --telescope LST --photons 100000 --output lst_paths.csv
obdeect-plot-toy lst_paths.csv --telescope LST --output lst_paths.png
```

Create a focal-plane intensity image with weighted 1-D projections:

```bash
obdeect-plot-toy lst_paths.csv --telescope LST --focal-plane --bins 96 \
  --output lst_focal_plane.png
```

## Artificial calibration sources

The baseline executable uses 400-nm photons and supports three deterministic
source models over the same entrance pupil:

```bash
# Plane wave from an on/off-axis star; angles are telescope-frame degrees.
./build/obdeect_toy --source star --field-x-deg 0.5 --field-y-deg 0.0

# Finite-distance point flasher, with per-photon inverse-square weights.
./build/obdeect_toy --source illuminator --distance-m 50

# Collimated or finite-divergence calibration laser.
./build/obdeect_toy --source laser --distance-m 50 --divergence-deg 0.1
```

CSV output stores every traced path vertex plus wavelength, emission time,
source weight and the trace's wavelength-independent throughput. The latter is
currently one for a detector-surface hit and zero for every loss; a compiled
coating/material scene will replace it with wavelength-dependent transport.

## CTAO reference models

`obdeect_ctao --telescope LST|MST|SST|SCT` writes ragged photon paths for the
same Python plotter. The catalogue is pinned to public `simulation-models`
6.3.0 identifiers and its import API requires explicit provenance. It contains
optical prescriptions only: it does not load model JSON, facet positions,
camera pixels, alignment, structures, throughput, or the SCT coordinate
transform. LST ideal-paraboloid and MST central-sphere baselines are executable
and tested; SST/SCT are two-mirror prescription scaffolding awaiting their
model-specific geometry validation.

The core now also has a finite-facet and tabulated-coating kernel for use by a
future model importer. A facet list is not currently embedded in any reference
model, so this must not be interpreted as segmented CTAO telescope support.
The explicit evidence required before claiming sim_telarray/ROBAST-level
coverage is in [docs/SIMTELARRAY_REPLACEMENT.md](docs/SIMTELARRAY_REPLACEMENT.md).

### Import a pinned simulation-models production

The standard-library importer records every selected parameter record and all
declared model-file assets with SHA-256 hashes. It neither downloads nor copies
external model data.

```bash
obdeect-import-simulation-models /path/to/simulation-models LSTN-design \
  --version 6.3.0 --output lstn-design.ir.json
```

This command selects and hashes the source parameter records and their declared
assets. It is provenance input for scene compilation, not a ray-tracing command
and does not download, copy, or silently interpret model files.

The emitted `obdeect.simulation-models-ir.v1` JSON is the auditable hand-off
from model selection to the future C++ scene compiler. An unresolved parameter
file, missing declared asset, identity mismatch, or unsafe path fails the
import; no field is silently discarded.

## Coordinate and model convention

The mirror vertex is at `z=0 m`; incoming artificial Cherenkov photons begin
at `z=20 m` and travel in `-z`. The spherical centre is at `z=9.75 m`, making
the paraxial screen location `z=4.875 m`. The simple structure is four
finite-cylinder support legs plus a circular camera face evaluated before the
primary reflection. CSV records source, termination/mirror, and focal-screen
vertices so the plot displays the actual traced path.

## Core architecture status

The C++20 core is split by responsibility so a photon kernel never needs Python, YAML, EventIO or a plotting dependency.

| Component | Current implementation |
| --- | --- |
| `math.hpp` | `Vec3`, dot/cross/norm and checked direction normalisation. |
| `photon_buffer.hpp`, `abi.hpp` | SoA input contract, result buffers, status values and validation. |
| `source.hpp` | Deterministic 400-nm star, point illuminator and laser sources. |
| `tables.hpp` | Immutable no-extrapolation 1-D response interpolation. |
| `intersections.hpp`, `geometry.hpp` | Plane, sphere/cap, disk, finite cylinder and axisymmetric-surface dispatch. |
| `interactions.hpp` | Checked specular reflection. |
| `axisymmetric_optics.hpp` | Paraboloid/even-polynomial SC surfaces and forward Newton intersection. |
| `facets.hpp` | Finite circular facet intersection, nearest-hit selection and tabulated coating response. |
| `scene.hpp`, `trace.hpp` | Immutable directed toy-scene compilation and scalar SoA block tracing. |
| `diagnostics.hpp` | Status/weight closure summary. |
| `model_import.hpp` | Canonical CTAO model provenance/import target. |

## Linting

Install CMake, Ninja, a C++20 compiler, Python 3.14+, and the development
tools:

```bash
python -m pip install --upgrade pip ruff
python -m pip install .
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

## Citation and License

The project is BSD-3-Clause licensed; citation metadata is in
[CITATION.cff](CITATION.cff).
