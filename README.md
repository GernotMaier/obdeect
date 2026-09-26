# obdeect

`obdeect` is a C++20 optical ray tracing prototype for imaging atmospheric
Cherenkov telescopes. The C++ core uses only the standard library. Python
provides plotting, PSF analysis, and model import tools.

## Start here

Requirements: a C++20 compiler, CMake 3.20+, Ninja or Make, and Python 3.10+.

```sh
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'
cmake --preset debug
cmake --build --preset debug
ctest --test-dir build/debug --output-on-failure
python -m unittest discover -s python/tests
```

Trace the simple spherical telescope and make three plots:

```sh
./build/debug/obdeect_toy --photons 10000 --output toy.csv
obdeect-plot-toy --telescope toy-mst --view structure --output structure.png
obdeect-plot-toy toy.csv --telescope toy-mst --view rays --output rays.png
obdeect-plot-toy toy.csv --telescope toy-mst --view focal-plane --output focal.png
```

The structure plot shows the toy mirror, camera, and four supports in two
side projections. The ray plot uses actual CSV vertices. The focal plot shows
detected weight per bin and x/y projections. `--max-paths` limits displayed
rays; `--bins` controls the focal histogram. Coordinates are telescope frame
metres. The older `--focal-plane` flag remains supported.

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
same Python plotter. The analytic catalogue uses public `simulation-models`
6.3.0 identifiers as a tested baseline; the production importer accepts an
explicitly selected version from the supplied checkout and records provenance.
It contains optical prescriptions only: it does not load model JSON, facet positions,
camera pixels, alignment, structures, throughput, or the SCT coordinate
transform. LST ideal-paraboloid and MST central-sphere baselines are executable
and tested; SST/SCT are two-mirror prescription scaffolding awaiting their
model-specific geometry validation.

The core now also has a finite-facet and tabulated-coating kernel for use by a
future model importer. A facet list is not currently embedded in any reference
model, so this must not be interpreted as segmented CTAO telescope support.
The explicit status, remaining work, and evidence required before claiming
sim_telarray/ROBAST-level coverage are in [docs/STATUS.md](docs/STATUS.md).
Reference manifest and comparison schemas are described in
[docs/REFERENCE_MANIFEST.md](docs/REFERENCE_MANIFEST.md) and
[docs/COMPARISON_CONTRACT.md](docs/COMPARISON_CONTRACT.md).

### Import a selected simulation-models production

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

### CORSIKA 7 EventIO input

The optional C++ EventIO adapter reads CORSIKA IACT `TELFIL` files directly.
In this workspace it builds against the adjacent EventIO C source checkout;
an installed library can instead be selected with `-DEventio_ROOT=/path/to/install`.

```bash
cmake -S . -B build/eventio -DOBDEECT_BUILD_EVENTIO_INPUT=ON
cmake --build build/eventio --target obdeect_eventio_summary
./build/eventio/obdeect_eventio_summary /path/to/CORSIKA_TELFIL
```

`EventioPhotonReader` yields one weighted `OpticalPhoton` per CORSIKA bunch,
with a `PhotonBatchContext` for run, event, reused array and telescope. It
supports full, compact and 3D records in both array layouts. A zero wavelength
stays unresolved; CEFFIC photoelectron-like bunches are rejected. Use
`resolve_eventio_spectrum()` to create deterministic spectral children before
any wavelength-dependent response. `AtmosphereTransmissionTable` reads the
simulation transmission table and `attenuate_eventio_direct_beam()` applies
direct-beam extinction between emission and telescope arrival. Check that the
table's observation altitude and CORSIKA production configuration match before
using it. These stages are separate from the file decoder and optical kernel.
The reader's CORSIKA-local arrival rays must be transformed into the selected
telescope scene frame before tracing. See the
[input and atmosphere plan](docs/CORSIKA7_EVENTIO_INPUT_PLAN.md) for conventions
and remaining production-file validation.

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

Install CMake, Ninja, a C++20 compiler, Python 3.10+, and the development
tools:

```bash
python -m pip install --upgrade pip ruff
python -m pip install .
```

## Examples and tutorials

[Tutorials](docs/TUTORIALS.md) cover sources, optical reference models,
field-angle scans, and model provenance. Run the examples with:

```sh
python examples/run_examples.py --build build/debug --output out/examples --photons 1000
```

Each case writes a CSV, structure/ray/focal plots where applicable, and status
counts. See [examples/README.md](examples/README.md) for the case list.

## Scope

The toy scene is MST inspired; it is not a CTAO production telescope. The
`obdeect_ctao` command traces analytic LST/MST optical baselines and contains
experimental SST/SCT prescriptions. CTAO plots show optical outlines only:
structures, segmented mirrors, camera pixels, and materials are not built into
those traces. Importing a `simulation-models` record records provenance and
partial geometry but does not make it trace ready. Current limitations and
validation evidence are in [status](docs/STATUS.md). The current
[code review](docs/CODE_REVIEW.md) records fixes and remaining duplication.

## PSF scan

```bash
obdeect-psf scan --executable ./build/debug/obdeect_toy --photons 10000 \
  --field-x-deg 0 0.5 1.0 --output-dir out/toy_psf --plot out/toy_psf.png
```

A scan writes one trace CSV per offset and `psf_scan.csv`/`psf_scan.json`.
These metrics describe the selected deterministic reference configuration;
they are not validation of a production telescope model.

The [production 7.0.0 sim_telarray reference](docs/reference/7.0.0/README.md)
includes focal-plane PSFs, integration radii, cumulative distributions, and
archived photon lists for LST, both MST cameras, and SST across their selected
North/South sites. These are inputs for future obdeect production-scene
comparisons.

## Model provenance and scene compilation

The Python adapters can select a `simulation-models` record, hash its
declared assets, and compile a deliberately incomplete, provenance-checked
scene hand-off. Compilation retains mirror-list geometry, derives nominal
single-reflector panel normals, and reports run-specific alignment, detector
surfaces, structures, and materials as trace blockers. It does not connect a
production scene to a native trace executable yet.

```bash
obdeect-import-simulation-models /path/to/simulation-models LSTN-design \
  --version 7.0.0 --output lstn.ir.json
obdeect-compile-scene lstn.ir.json --source-root /path/to/simulation-models \
  --simtel-root /path/to/sim_telarray --output lstn.scene.json
```

`--simtel-root` is needed when a camera response table is found in
sim_telarray's `cfg/CTA` search path rather than simulation-models `Files`.
The selected file path and SHA-256 hash are recorded in the scene provenance.

CSV trace output includes source weight, wavelength, emission time, terminal
status, path length, and path vertices. The current executable-level
throughput is one for a detector-surface hit and zero for a loss; material and
coating kernels exist but are not yet integrated into these analytic scenes.

## Project map

| Path | Purpose |
| --- | --- |
| `cpp/include/obdeect/` | Geometry, optics, sources, and tracing |
| `cpp/tests/` | Native tests |
| `python/obdeect/` | Import, plotting, and PSF commands |
| `python/tests/` | Python tests |
| `examples/` | Runnable end to end examples |
| `docs/` | Tutorials, status, and comparison contracts |

BSD-3-Clause license. Citation metadata is in [CITATION.cff](CITATION.cff).
## Generative AI disclosure

Generative AI tools were used to write much of this project; outputs were
reviewed and validated by the authors.

Generative AI tools (mostly ChatGPT 5.6) were used to write the entire code of this project. All AI-assisted outputs were reviewed, validated, and, where necessary, modified by the authors to ensure accuracy and reliability.
