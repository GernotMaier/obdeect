# obdeect

`obdeect` is fast ray-tracing code for imaging atmospheric Cherenkov telescopes (IACTs).

Dependencies area C++20 compiler and the standard library for the core ray-tracing code.

## Repository layout

```text
cpp/include/obdeect/  C++20 trace kernels and model types
cpp/src/              executable entry points
python/obdeect/       Python diagnostic package
```

## Installation and testing

Set up the environment:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install .
```

Compile and run the C++ executables:

```bash
# Configure the build system
cmake -S . -B build
# Build the project
cmake --build build
# Run the tests and executables
ctest --test-dir build --output-on-failure
# Run the Python unit tests
pytest python/tests/test_plot_toy_mst.py
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

## Artificial calibration sources

The baseline executable uses 400-nm photons and supports three deterministic
source models over the same entrance pupil:

```bash
# Plane wave from an on/off-axis star; angles are telescope-frame degrees.
./build/obdeect_toy --source star --field-x-deg 0.5 --field-y-deg 0.0

# Finite-distance point flasher, with per-ray inverse-square weights.
./build/obdeect_toy --source illuminator --distance-m 50

# Collimated or finite-divergence calibration laser.
./build/obdeect_toy --source laser --distance-m 50 --divergence-deg 0.1
```

CSV output stores every traced path vertex plus wavelength, emission time,
source weight and the trace's wavelength-independent throughput. The latter is
currently one for a detector-surface hit and zero for every loss; a compiled
coating/material scene will replace it with wavelength-dependent transport.

## CTAO reference models

`obdeect_ctao --telescope LST|MST|SST|SCT` writes ragged ray paths for the
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
python tools/import_simulation_models.py /path/to/simulation-models LSTN-design \
  --version 6.3.0 --output lstn-design.ir.json
```

After installation, use the equivalent stable command without relying on the
checkout layout:

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
python -m unittest python/tests/test_plot_toy_mst.py
ruff format --check python
ruff check python
clang-format --dry-run --Werror cpp/include/obdeect/toy_mst.hpp cpp/src/toy_main.cpp cpp/tests/test_core.cpp
```
