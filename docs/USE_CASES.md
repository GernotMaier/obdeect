# Optical-study recipes

Example use cases for the `obdeect` ray tracing package.

Requires a simulation models as defined in the [CTAO simulation models repository](https://gitlab.cta-observatory.org/cta-science/simulations/simulation-model/simulation-models).
All examples use named arguments; run `--help` before adapting a recipe.

## Developer demo

This is a runnable kernel demonstration, not a CTAO telescope model. It traces
artificial photons through a simple MST-inspired spherical mirror, camera
shadow, and four mast supports. Use it to inspect the basic ray-tracing output
without downloading a production model.

```sh
obdeect-demo-mst --source star --field-x-deg 0.5 --field-y-deg -0.2 \
  --wavelength-nm 400 --photons 10000 --output demo-mst.csv
obdeect-plot-reference --input demo-mst.csv --view focal-plane \
  --output demo-mst-focal-plane.png
```

`obdeect-analytic-optics` is likewise a developer diagnostic. For imported
CTAO production geometry, use the compiled optical-model recipes below.

The following two environment variables might be useful:

```sh
export OBDEECT_SIMULATION_MODELS_PATH=../simulation-models
export OBDEECT_SIMULATION_MODELS_VERSION=7.0.0
```

## Define the optical model

The following command reads the simulation model and
compiles the optical model required by `obdeect` for the ray tracing.
This is for the collection of surfaces, materials, responses, and detector
geometry through which photons propagate.
The JSON output is the single canonical artifact: it contains the auditable
model and its directly traceable finite-surface representation.

Example:

```sh
obdeect-compile-optical-model \
  --source-root "$OBDEECT_SIMULATION_MODELS_PATH" --model LSTN-design \
  --version "$OBDEECT_SIMULATION_MODELS_VERSION" \
  --output lst.optical-model.json
```

## Runnable native studies

### Single-panel 2F test stand

Find the panel ID in `lst.geometry.json` under `primary.facets`. Place the source at the panel centre plus twice its `focal_length_m` along `nominal_normal`, and place a planar test screen at that point. Select the panel, then derive and plot the screen distribution.

```sh
obdeect-simtools-raytrace --optical-model lst.optical-model.json \
  --source illuminator --panel-id 42 --source-x-m -0.2 --source-y-m 0.1 --source-z-m 57.3 \
  --screen-x-m -0.2 --screen-y-m 0.1 --screen-z-m 57.3 --screen-radius-m 0.25 \
  --wavelength-nm 400 --photons 100000 --output panel-42.csv
obdeect-psf derive --input panel-42.csv --output panel-42-psf.json
obdeect-plot-reference --input panel-42.csv --view focal-plane --output panel-42-psf.png
```

Use coordinates computed from the selected panel; the numbers above are only an LST example. The custom screen replaces the compiled focal boundary for this command.

### LST/MST star PSF

Stars are plane waves. `field-x-deg` and `field-y-deg` define the two-dimensional offset. A comma-separated wavelength list is sampled evenly.

```sh
obdeect-simtools-raytrace --optical-model lst.optical-model.json --source star \
  --field-x-deg 0.5 --field-y-deg -0.2 --wavelength-nm 300,400,500 \
  --photons 100000 --output lst-star.csv
obdeect-psf derive --input lst-star.csv --output lst-star-psf.json
obdeect-plot-reference --input lst-star.csv --view focal-plane --output lst-star-psf.png
```

### Nearby illuminator, laser, and arrival time

An illuminator is a finite point source. A laser has an explicit origin, axis, and divergence. `--emission-time-ns` and `--pulse-width-ns` add a deterministic top-hat source pulse; the output time is source emission time plus vacuum geometric optical path time.

```sh
# Dish illuminator at 100 m.
obdeect-simtools-raytrace --optical-model lst.optical-model.json --source illuminator \
  --source-x-m 0 --source-y-m 0 --source-z-m 100 --wavelength-nm 400 \
  --emission-time-ns 20 --pulse-width-ns 4 --photons 100000 --output illuminator.csv

# Calibration/alignment laser.
obdeect-simtools-raytrace --optical-model lst.optical-model.json --source laser \
  --source-x-m 0 --source-y-m 0 --source-z-m 10 \
  --direction-x 0 --direction-y 0 --direction-z -1 --divergence-deg 0.05 \
  --wavelength-nm 355 --photons 100000 --output laser.csv

obdeect-photon-distributions --input illuminator.csv --output illuminator-distributions.svg
obdeect-plot-reference --input laser.csv --view rays --max-paths 300 --output laser-rays.png
```

The top-hat pulse is an optical timing input, not a replacement for a source-specific flasher or laser pulse model.

## Production CTAO studies

Use `simtools` with its default `sim_telarray` backend for full CTAO models. It is the validated path for LST, MST, SST, and SCT production configurations, including model-defined artificial sources. Select an actual telescope ID and site from the production; for example, `LSTN-01`/North, `MSTN-04`/North, `SSTS-04`/South, or `SCTS-01`/South.

### Dish PSF from stars

`simtools-validate-optics` runs a two-dimensional point-source offset scan and writes D80, effective area, effective focal length, and optional focal-plane images. Use a source distance of `10000` km for an effectively infinite-distance star.

```sh
simtools-validate-optics \
  --site South --telescope SSTS-04 --model_version 7.0.0 \
  --source_distance 10000 --max_offset 3 --offset_step 0.5 \
  --plot_images
```

Repeat for the selected LST, MST, or SCT telescope/site.

### Primary/secondary/focal incidence angles

`simtools-derive-incident-angle` traces point sources and writes the incident-angle tables and plots for the primary mirror, secondary mirror when present, and focal plane. `--debug_plots` adds two-dimensional hit maps.

```sh
simtools-derive-incident-angle \
  --site South --telescope SSTS-04 --model_version 7.0.0 \
  --source_distance 10000 --off_axis_angles 0 0.5 1 2 3 \
  --number_of_photons 1000000 --calculate_primary_secondary_angles --debug_plots
```

Run the same command with the selected LST, MST, or SCT telescope/site. Keep the production-generated tables; do not use native nominal LST/MST results as a substitute for this comparison.

### Flashers and model-defined illuminators

Use `simtools-simulate-flasher` with the light-source name defined by the selected model production. It resolves the source configuration from `simulation-models` rather than requiring a hand-written source geometry.

```sh
simtools-simulate-flasher \
  --run_mode full_simulation --site North --model_version 7.0.0 \
  --light_source LSFN-design --telescopes LSTN-01 --number_of_events 10
```

For model-defined lasers, use the `ls-beam` configuration generated by `simtools`/`sim_telarray`; do not treat the native geometric laser as an atmospheric-scattering calculation.

## Not yet runnable from this repository

These requests need geometry or response data that the current compiled optical model does not carry. They must not be approximated or labelled as CTAO results.

| Study | Missing compiled optical-model content |
| --- | --- |
| Camera-window incidence; flat or spherical windows | Window surface, curvature, thickness, and wavelength/incidence response. |
| Shadowing versus offset (1D/2D) | Camera housing, masts, baffles, and other obscurer solids. |
| Throughput versus offset and input spectrum | Mirror/coating response, alignment perturbations, obscurer geometry, and spectral weights. |
| Full structure-and-ray rendering | Compiled 3D component geometry and per-component interaction records. |
| Segmented SST/SCT structure, shadowing, and throughput | Segment placement/alignment, masts, camera housing, baffles, material binding, and validation fixtures. The M1/M2/focal aspheres are exported and traced. |

The existing native structure plot shows only compiled panel centres and the focal boundary. It intentionally does not draw unmodelled hardware.
