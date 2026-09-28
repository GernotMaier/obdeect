# obdeect

`obdeect` traces optical photons to the focal surface. It is useful for model-derived nominal LST/MST panel studies, source studies, PSFs, and optical arrival-time distributions. It does not simulate camera electronics or sensor conversion.

For the supported workflows, see [docs/USE_CASES.md](docs/USE_CASES.md). That guide separates runnable native studies from production studies that must use `simtools`/`sim_telarray`.

## Install

```sh
python3 -m pip install obdeect-dev
```

From a source checkout:

```sh
cmake --preset debug
cmake --build --preset debug
python3 -m pip install -e .
```

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

## Results

`obdeect-simtools-raytrace` writes `obdeect-arrival-v1` CSV. Each row records source type, wavelength, emission time, input weight, terminal status, optical path length, interaction points, and primary/secondary/focal incidence angles. Only rows with `status=detected` contribute to a focal-plane PSF or arrival-time distribution.

## Limits

The native model-derived path currently exports nominal LST/MST panels and a finite focal boundary. It does not yet bind CTAO structural geometry, wavelength-dependent coatings, camera windows, mirror perturbations, or dual-mirror SST/SCT surfaces. Use the standard `simtools`/`sim_telarray` workflow for production studies until the matching scene and validation gate are available.
