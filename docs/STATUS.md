# obdeect status and implementation plan

## Scope and success metric

`obdeect` is a C++20 optical photon tracer. It must become a selectable
simtools backend for complete optical ray tracing, PSF analysis, effective
area/focal length, and focal/M1/M2 incidence angles. `sim_telarray` remains
the default backend and the reference implementation.

The scope ends at optical arrival at a physical detector surface. Camera
bodies, windows, filters, light guides, supports and obscurers are included
when they change photon transport, timing or shadowing.

Every release is judged by:

- **Accuracy:** matched simtools/sim_telarray observables and interaction
  identities within predeclared scientific tolerances.
- **Precision:** SI units, binary64 geometry, stable frame/time conventions,
  deterministic seeds and reproducible loss accounting.
- **Performance:** bounded memory, reusable compiled scenes, batch tracing and
  measured photons/s, compilation, I/O and peak RSS.

The importer accepts any compatible, explicitly selected simulation-models
production and schema. Production 7.0.0 is the first validation fixture, not a
version limit. ROBAST can supply an additional versioned geometry adapter when
it provides equivalent surfaces and provenance. SCT is a later validation
phase.

## Implemented

- Deterministic scalar tracing with normalized rays, finite facets, polynomial
  aspheres, nearest-hit ordering, bounded interactions, detector surfaces and
  terminal loss statuses.
- Optical material primitives: tabulated response, Fresnel/Snell refraction,
  total internal reflection, absorption and parallel windows.
- Star, finite-distance illuminator and laser source kernels are distinct:
  stars are plane waves with launch-plane phase timing, illuminators use a
  finite 3D position with inverse-square weighting, and lasers use an
  arbitrary launch-plane origin, axis and divergence cone. EventIO adapters
  decode CORSIKA 7 full, compact, split and 3D photon bunches.
- Model import provenance: selected production records, asset hashes,
  path-safety checks, current parameter-local and legacy shared asset layouts,
  nominal mirror footprints/normals, focal-plane layout geometry and
  unresolved-field reporting. ECSV mirror lists and focal-plane layouts from
  simulation-models compile for the LST and both MST reference variants.
  Imported scenes are not trace-ready production telescopes.
- Reference-manifest tooling and a 7.0.0 simtools/sim_telarray matrix with
  archived imaging lists, PSF containment radii and cumulative profiles.
- Reference structure, ray and focal-plane plotting plus machine-readable PSF
  summaries. Ray plots now distinguish stars, finite illuminators and lasers
  by source colour and annotate the mean focal incidence angle. Focal plots
  retain weighted x/y projections and containment summaries. These still use
  analytic structure outlines where compiled geometry is not available.
- The versioned `obdeect-arrival-v1` contract validates CSV records before
  analysis: finite SI scalars, non-negative weights, bounded throughput,
  terminal status, path length and interaction points. It exposes detected
  focal coordinates and optical weights.
- A native `obdeect-simtools-raytrace` reference CLI now emits that contract
  for star, finite-distance illuminator and laser inputs, with telescope,
  field, distance, wavelength, source-position and divergence options. It is a
  reference vertical slice; its records include focal, primary and secondary
  incidence angles, interaction points and an explicit `source_kind` field so
  downstream plots and analysis cannot confuse source classes. It does not yet
  compile simulation-model assets.
- A `SimulatorObdeect` runner and explicit `RayTracing` backend branch are now
  present in the sibling simtools checkout. They invoke the reference CLI and
  preserve sim_telarray as the default. The runner forwards telescope,
  offsets, source position, laser axis, distance, photon count, wavelength and
  divergence.
- `IncidentAnglesCalculator` now has an obdeect path using the shared arrival
  records, and `PSFImage` reads the same CSV contract for detected focal
  points. Both reject unsupported weighted PSF inputs instead of silently
  producing biased results.
- An optional simtools backend configuration keeps `sim_telarray` as the
  default and labels the present native path as a reference implementation.

## Remaining implementation, in order

### 1. Freeze contracts and references

1. Define versioned SI schemas for source batches, compiled scenes, arrivals,
   interactions and losses. Preserve run/event/telescope/bunch IDs, weights,
   wavelength, emission and arrival times, paths, surface ID and incidence
   angle.
2. Generate the same checksummed photon block for both engines. Record model
   checkout/version, site, telescope, atmosphere, seeds, enabled optics and
   tolerances. Use simtools wrappers wherever available.
3. Keep the 7.0.0 matrix as the first fixture, then add compatible model
   versions and schema adapters. Reject unknown required fields.

### 2. Compile realistic optical telescopes

1. Compile model-driven LST parabolic and MST classical/modified Davies--
   Cotton segmented reflectors: curvature, gaps, panel alignment/roughness,
   selectable panels, focal surface, camera body, windows, supports, masts and
   obscurers.
2. Compile SST dual-mirror M1/M2 aspheres, segmentation, holes, baffles,
   shadows and curved focal surface. Add SCT only with its own model and gate.
3. Bind optical coatings, filters, windows, light guides and active detector
   surfaces. Keep optical arrival, transmission and component losses separate.
   Do not model electronic response.
4. Keep the scene IR observatory-neutral. Add a ROBAST adapter only through
   the same surface/material/provenance contract.

### 3. Replace the simtools guards with a working backend

1. Implement `obdeect-simtools-raytrace` and `ObdeectSimulator` with selected
   model/site/version, zenith and two-dimensional offsets, source distance,
   wavelength, single-panel selection, alignment/roughness settings, seed and
   output paths. Preserve simtools force/test modes and concurrent offsets.
2. Complete the reference backend in `simtools.ray_tracing.RayTracing` for
   production scenes. Keep sim_telarray unchanged and default.
3. Emit a versioned arrival table and summary consumable without `rx`: focal
   x/y, throughput, optical time, path length, terminal loss, surface IDs and
   provenance.
4. Replace the `IncidentAnglesCalculator` guard and debug-trace column parsing
   with the same interaction records. Preserve focal, M1 and M2 positions,
   angles, optional mirror mode and offset scans.
5. Make `PSFImage`, effective-area, focal-length and mirror-panel PSF/RNDA
   analysis backend-neutral. Existing simtools ECSV results and plotting must
   work for both engines.

### 4. Implement CTAO light sources

1. **Stars:** match simtools `SimulatorRayTracing`: finite/infinite distance,
   zenith, two-dimensional offsets, launch sampling, wavelength and timing.
2. **Illuminators/flashers:** match `SimulatorLightEmission` and LightEmission
   `xyzls`/`ff-1m`: 3D position, visibility, pointing, spectrum, angular
   distribution, pulse shape, delay, photon levels and atmosphere.
3. **Lasers:** support centre-of-dish and arbitrary position, direction,
   divergence, spectrum and timing. Use LightEmission `ls-beam` as the
   qualified reference for beam attenuation and scattered light; keep
   atmospheric scattering a separate explicit model.
4. Validate source transforms, normalization, timing, wavelength sentinels and
   extinction exactly once before telescope optics.

### 5. Improve diagnostics and plotting

Plot the compiled 3D assembly and projections: mirror panels, M2, camera body
and window, light guide, focal surface, supports, masts, baffles and obscurers.
Overlay photon paths, component hits/losses and incidence angles; colour paths
by wavelength. Provide focal x/y, wavelength and time slices with centroid,
containment radii, throughput, timing and loss summaries. Never substitute
illustrative geometry for compiled geometry without saying so.

### 6. Validate accuracy, precision and performance

Run primitive/import tests first, then single-panel, full-telescope,
dual-mirror, spectral, source, EventIO and reproducibility fixtures. Compare
against simtools-driven sim_telarray with identical photon blocks and
conventions. Retain residual maps, first divergent interactions and loss
closures. Test photon order, block size and thread count. Publish tolerances,
Monte Carlo uncertainty, photons/s, peak RSS, compile time, I/O time and scene
cache time. Do not enable a telescope/source workflow until its fixture passes.

## Current blockers and completion gate

The simtools backend selector, container integration, reference CLI, incident
angle path and CSV PSF reader exist. The production model compiler still stops
before aligned geometry, structures, physical detector surfaces and material
bindings; those compiled scenes and their acceptance fixtures remain the next
blockers.

Claim production compatibility only after LST, MST and SST pass the complete
ray-tracing, PSF, effective-area/focal-length, source and incidence-angle
matrix. SCT requires a separate gate. Keep sim_telarray available as the
default and permanent comparison backend.
