# obdeect

`obdeect` traces optical photons from a defined source to a physical detector
surface. It is intended for CTAO optical studies: field-of-view stars, nearby
illuminators, and calibration/alignment lasers. It does not simulate camera
electronics, triggering, or sensor conversion.

## Install

```sh
python3 -m pip install obdeect-dev
```

| Command | Purpose |
| --- | --- |
| `obdeect-import-simulation-models` | Select optical ray-tracing inputs from a `simulation-models` production. |
| `obdeect-compile-scene` | Compile those inputs into an auditable optical scene and native surface table. |
| `obdeect-simtools-raytrace` | Trace stars, illuminators, or lasers through a compiled scene. |
| `obdeect-plot-reference` | Plot paths, focal-plane distributions, or compiled geometry. |
| `obdeect-psf` | Derive a PSF from arrival records or run a field-angle scan. |

Run any command with `--help` for its input and output contract. Native
programs also accept `-h`.

## CTAO telescope simulations

Choose an explicit `simulation-models` version and model, then compile before
tracing. The importer keeps only parameters that affect optical transport
(geometry and optical response); camera readout, trigger, and calibration
parameters are deliberately excluded.

```sh
MODEL_ROOT=/path/to/simulation-models

obdeect-import-simulation-models \
  --source-root "$MODEL_ROOT" --model LSTN-design --version 7.0.0 \
  --output lst.ir.json
obdeect-compile-scene \
  --input lst.ir.json --source-root "$MODEL_ROOT" \
  --output lst.scene.json --native-output lst.scene.csv
obdeect-simtools-raytrace \
  --scene-file lst.scene.csv --source star --photons 10000 \
  --field-x-deg 0.0 --field-y-deg 0.0 --wavelength-nm 300,400,500 \
  --output lst-arrivals.csv
```

Use the same workflow with the selected production model for each CTAO
telescope configuration:

```sh
# LST
obdeect-import-simulation-models --source-root "$MODEL_ROOT" --model LSTN-design --version 7.0.0 --output lst.ir.json
# MST FlashCam or NectarCam
obdeect-import-simulation-models --source-root "$MODEL_ROOT" --model MSTx-FlashCam --version 7.0.0 --output mst-flashcam.ir.json
obdeect-import-simulation-models --source-root "$MODEL_ROOT" --model MSTx-NectarCam --version 7.0.0 --output mst-nectarcam.ir.json
# SST
obdeect-import-simulation-models --source-root "$MODEL_ROOT" --model SSTS-design --version 7.0.0 --output sst.ir.json
```

Compile each IR and pass the resulting `--scene-file` to
`obdeect-simtools-raytrace`. A compiled surface table is a CSV representation
of finite optical surfaces (panel centres, normals, boundaries, and detector
surface), not an invented “scene file”. The native tracer accepts the compiled
general panel geometry; it does not assume axisymmetry.

## Sources

Stars use a plane wave with two field offsets. For a nearby illuminator, give
its physical location in telescope coordinates; rays are emitted from that
point and are therefore non-parallel across the dish. For a laser, give its
origin, beam axis, and divergence. `--wavelength-nm` accepts a comma-separated
discrete spectrum; samples are distributed deterministically across it.

```sh
# Off-axis star
obdeect-simtools-raytrace --scene-file lst.scene.csv --source star \
  --field-x-deg 0.5 --field-y-deg -0.2 --wavelength-nm 350,400,450 \
  --photons 30000 --output star.csv

# Nearby dish illuminator (50–1000 m is typical)
obdeect-simtools-raytrace --scene-file lst.scene.csv --source illuminator \
  --source-x-m 0 --source-y-m 0 --source-z-m 100 --wavelength-nm 400 \
  --photons 30000 --output illuminator.csv

# Calibration/alignment laser, for example near the dish centre
obdeect-simtools-raytrace --scene-file lst.scene.csv --source laser \
  --source-x-m 0 --source-y-m 0 --source-z-m 0 \
  --direction-x 0 --direction-y 0 --direction-z -1 --divergence-deg 0.05 \
  --wavelength-nm 355 --photons 30000 --output laser.csv
```

The source options describe the ray launch. Pulse shape, atmospheric
attenuation, and scattered light must be supplied by a qualified source and
atmosphere model before optical transport; they are not silently inferred.

## Results and plots

The trace command reports the source, number of photons, detected photons, and
output path. Its `obdeect-arrival-v1` CSV records the source kind, wavelength,
emission time, weight, terminal status, path length, interaction coordinates,
and primary/secondary/focal incidence angles.

```sh
obdeect-plot-reference --input star.csv --view focal-plane --output star-psf.png
obdeect-photon-distributions --input star.csv --output star-distributions.svg
obdeect-psf derive --input star.csv --output star-psf.json
```

## Development

```sh
cmake --preset debug
cmake --build --preset debug
ctest --test-dir build/debug --output-on-failure
python3 -m unittest discover -s python/tests
```

`obdeect-demo-mst` and `obdeect-analytic-optics` are developer diagnostics.
They are not CTAO telescope models and are not part of the production workflow.
The test executables are written to `build/<preset>/tests/`, separate from
user-facing native programs.

## Validation scope

The current evidence and open validation gates are in [docs/STATUS.md](docs/STATUS.md).
Do not interpret a successful trace as equivalence to sim_telarray: production
equivalence requires the documented LST, MST, and SST geometry, source, PSF,
effective-area, focal-length, and incidence-angle comparison matrix.

Use the fail-closed [production validation procedure](docs/PRODUCTION_VALIDATION.md)
to compile a selected model, compare the identical resolved photon block with
sim_telarray, and retain the machine-readable residual summary.
