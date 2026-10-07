# Optical-study recipes

Example use cases for the `obdeect` ray tracing package.

Requires a simulation models as defined in the [CTAO simulation models repository](https://gitlab.cta-observatory.org/cta-science/simulations/simulation-model/simulation-models).
All examples use named arguments; run `--help` before adapting a recipe.

> [!WARNING]
> For current implementation limits and open tasks, read [STATUS.md](STATUS.md).

## Developer demo

Raytracing demonstration of a simple MST-inspired optical setup with a spherical mirror, camera shadow, and four mast supports.

Use it to inspect the basic ray-tracing output without downloading a production model.

```sh
obdeect-demo-mst --source star \
  --field-x-deg 0.5 --field-y-deg -0.2 \
  --wavelength-nm 400 --photons 10000 --output demo-mst.csv
obdeect-plot-reference --input demo-mst.csv --view focal-plane \
  --output demo-mst-focal-plane.png
```

`obdeect-analytic-optics` is likewise a developer diagnostic. For imported
CTAO production geometry, use the compiled optical-model recipes below.

## Realistic CTAO telescope studies

The following two environment variables might be useful:

```sh
export OBDEECT_SIMULATION_MODELS_PATH=../simulation-models
export OBDEECT_SIMULATION_MODELS_VERSION=7.0.0
```

### Define the optical model

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

### Plot the compiled optical model

Use the telescope plate for a model-faithful overview of the finite optical
geometry actually consumed by the tracer. It renders an orthographic assembly,
an axial section, the entrance pupil, and focal geometry. It does not add
camera supports, windows, or other structure absent from the selected model;
the footer declares the available geometry coverage.

```sh
obdeect-plot-reference --view telescope \
  --optical-model-json lst.optical-model.json \
  --output lst-telescope.png
```

For the panel layout as seen from the mirror, use the face-on view.  Each
finite panel face is drawn in its compiled tangent frame, so segmentation and
gaps remain visible.

```sh
obdeect-plot-reference --view pupil \
  --optical-model-json lst.optical-model.json \
  --output lst-primary-face.png
```

For a paper-style optical cross-section or a full-size CAD-like orthographic
view, select the corresponding model-derived renderer. `compiled-structure`,
`compiled-mirror`, and `compiled-3d` remain supported aliases.
`--view structure` also renders the compiled optical section when
`--optical-model-json` is supplied. Without an optical model, `structure` draws
the selected analytic reference outline; it is not an imported mechanical structure.

```sh
obdeect-plot-reference --view section --section-plane xz \
  --optical-model-json lst.optical-model.json \
  --output lst-optical-section.png

obdeect-plot-reference --view assembly-3d \
  --optical-model-json lst.optical-model.json \
  --output lst-optical-model-3d.png
```

### Recorded paths, terminal losses, and finite masks

Supplying a CSV to the telescope plate overlays its recorded flight segments and
interaction vertices. A cross marks the terminal vertex. Status controls colour;
source identity does not override that colour. Long source-flight segments are
cropped to a telescope-frame cube for display, without modifying the CSV or
solving another optical path. `--no-show-paths` retains the geometry-only plate.

```sh
obdeect-plot-reference --view telescope --input lst-star.csv \
  --optical-model-json lst.optical-model.json --path-selection stratified \
  --selection-seed 42 --max-paths 200 --context-radius-m 35 \
  --output lst-path-plate.svg

obdeect-plot-reference --view telescope --input lst-star.csv \
  --optical-model-json lst.optical-model.json --panel focal-plane-hits \
  --no-show-paths --bins 64 --output lst-geometry-and-psf.png
```

Use `--panel assembly`, `section`, `pupil`, or `focal-geometry` with the telescope
view to export one selected panel. In the pupil panel, supplied CSVs overlay
recorded entrance markers; in the focal geometry panel they overlay recorded
detected endpoints. Assembly and section panels overlay recorded flight
segments.

The default fourth panel contains detector geometry. Only the explicit
`--panel focal-plane-hits` option replaces it with weighted detected-hit density.
Parsed finite camera entrance footprints appear only in the focal panel; a
camera layout with centres and no footprint types does not create pixel polygons
or three-dimensional camera hardware. Segmented aspheric models use their
compiled hexagonal tangent-plane footprints or annular sectors, including
recorded radial and azimuthal gaps. Hexagonal curvature is not rendered; this
limitation is declared in the footer. Continuous aspheric contours use 256 radial
samples. Cylinder projection boundaries use 64 angular samples and capped
three-dimensional cylinders use 12 sides. The axial panel is an orthographic
projection of surface boundaries; cylinders entirely outside its selected plane
are omitted. It is labelled as a projection, rather than an exact mechanical cut.

A loss map uses each photon's recorded entrance vertex (point 0) and input weight.
Its numerator counts the input weight of terminal losses, and its denominator
counts all incident input weight in the same bin. Empty bins remain masked.
Detected photons with zero optical response remain detected outcomes; they do
not become geometric losses.

```sh
obdeect-plot-reference --view loss-map --input lst-star.csv \
  --optical-model-json lst.optical-model.json --bins 64 --output lst-losses.png

obdeect-plot-reference --view loss-map --input lst-star.csv \
  --optical-model-json lst.optical-model.json --status blocked_obscurer \
  --component-id 902 --output recorded-component-losses.svg
```

Status/component filters select the left-panel markers and the loss numerator;
the right-panel incident denominator retains all input photons. Component IDs
must be present in recorded diagnostics. Unknown or legacy IDs are never
assigned to nearby hardware by geometric inference. The legacy case is labelled
`component unknown`. Native versioned CSVs pass the shared optical-arrival
validator before plotting or PSF analysis; unversioned reference diagnostics
remain an explicit compatibility path.

`--path-selection first` retains the first-N behaviour. `stratified` distributes
slots across recorded status/source groups, ranked by a stable hash of photon
identity and the selection seed. `--colour-by wavelength`, `arrival-time`, or
`incidence-angle` requires the corresponding recorded finite values and shows a
numeric colour bar with units. `--view pupil --colour-by facet-z` colours aperture
centre height; `--label-panels` labels at most 300 segmented primary panels.
Plots always write to the output file without opening a window, regardless of
the configured Matplotlib backend. PNG and PDF use a file-only Matplotlib canvas.
Geometry SVGs have fixed artist
IDs and no creation timestamp; repeated exports are byte-stable for the same
inputs, options, and Matplotlib version. `--view focal-plane --output psf.svg`
uses the existing deterministic dependency-free weighted PSF summary.

### Single-panel 2F test stand

Find the panel ID in `lst.optical-model.json` under `primary.facets`. Place the source at the panel centre plus twice its `focal_length_m` along `nominal_normal`, and place a planar test screen at that point. Select the panel, then derive and plot the screen distribution.

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
obdeect-plot-reference --input lst-star.csv --view focal-plane \
  --optical-model-json lst.optical-model.json --output lst-star-psf.png
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

The compiled telescope plate shows finite panel boundaries, supplied detector surfaces, aspheric contours or finite masks, and supplied cylinders. Mechanical geometry absent from the compiled optical model is declared in the footer.
