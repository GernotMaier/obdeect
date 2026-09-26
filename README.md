# obdeect

`obdeect` is a C++20 optical ray tracing prototype for imaging atmospheric
Cherenkov telescopes. The C++ core uses only the standard library. Python
provides plotting, PSF analysis, and model import tools.

## Start here

You need a C++20 compiler, CMake 3.20+, Ninja, and Python 3.10+.

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

The mirror vertex is `z=0`; incoming artificial photons travel toward `-z`.
CSV rows hold source and interaction vertices, status, source weight, and
throughput. A detector hit has throughput one in executable reference scenes;
a loss has zero. These are geometry examples, not photoelectron or effective
area predictions.

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
Generative AI tools were used to write much of this project; outputs were
reviewed and validated by the authors.
