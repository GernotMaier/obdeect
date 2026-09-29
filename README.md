# obdeect

`obdeect` traces optical photons through model-derived telescope geometry and
writes photon paths and detector arrivals. It currently supports **nominal
geometry studies**, not full CTAO production simulations.

For panel tests and more plots, see [optical-study recipes](docs/USE_CASES.md).

## Install from this checkout

```sh
python3 -m pip install -e .
```

The package installs Python commands and the native tracer. For C++ development,
run `cmake --preset debug && cmake --build --preset debug`. Executables are in
`build/debug/` and tests in `build/debug/tests/`.

## Trace a star with a CTAO model

Set `MODEL_ROOT` to a `simulation-models` checkout. The example uses production
7.0.0; select the version of your own model checkout. Run from the directory
where you want output files.

```sh
MODEL_ROOT=/path/to/simulation-models
MODEL=LSTN-design
obdeect-import-simulation-models \
  --source-root "$MODEL_ROOT" --model "$MODEL" --version 7.0.0 \
  --output "$MODEL.ir.json"
obdeect-compile-scene \
  --input "$MODEL.ir.json" --source-root "$MODEL_ROOT" \
  --output "$MODEL.geometry.json" --native-output "$MODEL.surfaces.csv"
obdeect-simtools-raytrace \
  --scene-file "$MODEL.surfaces.csv" --source star \
  --field-x-deg 0.5 --field-y-deg 0 --wavelength-nm 300,400,500 \
  --photons 100000 --output "$MODEL.arrivals.csv"
obdeect-psf derive --input "$MODEL.arrivals.csv" --output "$MODEL.psf.json"
```

Repeat those steps with the following production table names. SCT import and
export work, but its native trace is not yet accepted; a PSF derivation can
fail when no photons reach the focal surface.

| Telescope family | `MODEL` | Native geometry represented |
| --- | --- | --- |
| LST | `LSTN-design` | Nominal curved finite primary panels, measured one-dimensional reflectivity, and circular focal bound. |
| MST | `MSTx-FlashCam` or `MSTx-NectarCam` | Nominal curved finite primary panels, measured one-dimensional reflectivity, circular focal bound, and model-supplied cylinders when present. |
| SST | `SSTS-design` | Rotationally symmetric M1, M2, and curved focal surface. |
| SCT | `SCTS-design` | Rotationally symmetric M1, M2, and curved focal surface; experimental, without an acceptance gate. |

The importer produces a provenance-recorded JSON intermediate representation
(`*.ir.json`). The compiler produces readable geometry (`*.geometry.json`) and
the small versioned surface table (`*.surfaces.csv`) consumed by the C++
tracer. These are generated from the selected model. `--scene-file` selects
that table.

## Nearby light sources

The illuminator samples rays from a finite position toward the dish with
inverse-square and projected-area weights. The laser samples a beam around its
given direction and divergence half-angle. Positions are metres in the
telescope frame.

```sh
obdeect-simtools-raytrace --scene-file "$MODEL.surfaces.csv" \
  --source illuminator --source-x-m 0 --source-y-m 0 --source-z-m 100 \
  --wavelength-nm 400 --photons 100000 --output illuminator.csv
obdeect-simtools-raytrace --scene-file "$MODEL.surfaces.csv" \
  --source laser --source-x-m 0 --source-y-m 0 --source-z-m 10 \
  --direction-x 0 --direction-y 0 --direction-z -1 \
  --divergence-deg 0.05 --wavelength-nm 355 \
  --photons 100000 --output laser.csv
```

These commands do not read model-defined flasher or laser source
configurations. They omit source visibility, spectrum, pulse shape,
attenuation, and atmospheric scattering. `--emission-time-ns` and
`--pulse-width-ns` add an ideal top-hat pulse for geometric timing studies.

## Commands

| Command | Output |
| --- | --- |
| `obdeect-import-simulation-models` | Optical geometry/response input selected from a production. |
| `obdeect-compile-scene` | Provenance-bound compiled geometry and native surface table. |
| `obdeect-simtools-raytrace` | Versioned photon-arrival CSV for a star, illuminator, or laser. |
| `obdeect-psf derive` | Weighted PSF, D80/R80, throughput, and loss summary. |
| `obdeect-photon-distributions` | Focal-plane and geometric arrival-time distributions. |
| `obdeect-plot-reference` | Compiled geometry, selected paths, or focal-plane image. |

Use `--help` on every command for its complete input/output contract.

## Results and limits

`obdeect-simtools-raytrace` writes `obdeect-arrival-v1` CSV. Each row records
source kind, wavelength, emission time, input weight, terminal status, optical
path length, interaction points, and reported incidence angles. Only detected
photons contribute to a PSF or arrival-time distribution.

Nominal LST/MST panels are traced as spherical facets using the catalogue panel
focal length. The focal bound does not model individual pixels. SST/SCT export
continuous aspheres without segment gaps or structural shadows. The native
path does not bind all incidence-dependent coatings, windows, noncylindrical
supports, alignment perturbations, or detector response. The production validation gate therefore
rejects these scenes. Use the standard `simtools`/`sim_telarray` backend for
CTAO production observables. See [status](docs/STATUS.md) and
[production validation](docs/PRODUCTION_VALIDATION.md) for the completion gate.
