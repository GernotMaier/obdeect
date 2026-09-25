# Project status and remaining work

**Status:** `obdeect` is a portable C++20 optical reference tracer.  It is not
yet a CTAO production simulator or a replacement for sim_telarray/ROBAST.
Current CTAO commands execute analytic optical prescriptions; they do not
construct the selected telescope from model assets.

The supported scope ends at optical arrival at a physical detector surface.
Sensor conversion may be added as a separate contract.  Electronics, triggers,
event writing, and array reconstruction are outside scope.

## Implemented baseline

- The core provides deterministic scalar tracing, normalized-ray validation,
  finite circular, square, and hexagonal facet intersections, simple detector
  surfaces, closed camera/mast obstructions for the toy scene, and
  component-specific terminal statuses.
- The analytic CTAO reference catalogue uses 6.3.0 values as a tested baseline.
  The production importer resolves an explicitly selected simulation-models
  production from the supplied checkout, records parameter records and asset
  SHA-256 hashes, retains metadata, constrains asset paths, and reports
  unresolved fields. It derives nominal single-reflector panel geometry but
  does not yet produce a trace-ready production scene.
- LST and MST have analytic reference prescriptions.  SST/SCT have tested
  two-mirror polynomial prescriptions, including corrected Schwarzschild--
  Couder reference-radius scaling and focal-plane vertices.  None is a
  production telescope implementation.
- Material kernels cover tabulated one-dimensional response, Fresnel/Snell
  refraction, total internal reflection, Beer--Lambert absorption, and a
  parallel window.  They are not connected to compiled telescope geometry.
- Source kernels provide deterministic star, finite-distance illuminator, and
  laser samples with wavefront timing, projected-pupil weighting, and
  beam-normal laser launch planes.  CSV input preserves photon and batch
  provenance.  An optional hessio probe checks the external library ABI; an
  EventIO reader is not implemented.
- Tests cover analytic primitives, material and source kernels, model
  prescriptions, importer validation, and toy-scene obstruction behavior.
- A reference-manifest tool freezes selected model revisions, input and
  executable hashes, commands, and comparison conventions. CSV and in-memory
  photon batches have an equivalence test; right-handed rigid frame transforms
  have a ray round-trip test. The full comparison result contract is pending.

## Plan progress

- Step 1: manifest creation and hash verification are implemented. Real
  reference commands and photon blocks are still needed to meet the exit gate.
- Step 2: input equivalence and frame round trips pass. The complete SI result
  and interaction schema remains open.
- Step 3: in progress. The compiler verifies imported records and assets and
  extracts mirror-list/segmentation footprints and camera pixel layouts.
  For single-reflector records it now derives unperturbed panel centres and
  normals using sim_telarray's `tel_setup_primary` prescription. These are
  nominal geometry only; run-specific random distance and alignment remain
  unresolved.
  Nested camera response files are hashed when present and missing references
  are reported. An explicit sim_telarray root resolves and hashes camera tables
  from its `cfg/CTA` search path. Normals, physical detector surfaces, structures, and materials
  remain unresolved, so these scenes are not trace-ready. The selected 6.3.0
  and 7.0.0 design records yield LST 198, MST 86, SST 18, and SCT 48/24
  primary/secondary mirror footprints; camera files yield 1855, 1764/1855,
  2048, and 11328 pixels respectively. The SCT camera references
  `Angular_response_MPPC_Prod3.dat`; it is present under the local
  `../sim_telarray/cfg/CTA` and resolves when that root is supplied. Patch
  productions 6.0.1, 6.0.2, 6.1.1, and 6.2.1 declare no primary geometry asset
  and cannot compile standalone scenes.
- Step 4: in progress. An immutable generic planar scene now traces the
  globally nearest finite mirror, detector, or opaque surface on each segment
  with a bounded interaction count. Analytic tests cover ordering before and
  after reflection, gaps, path/time, terminal IDs, and the interaction cap.
  Circular and annular even-polynomial aspheres can now participate in the
  same nearest-hit search, with tests for sag, normal, and a central hole.
  Segmented asphere masks, closed solids, active-area masks, and optical
  material behavior are not yet part of this scene.
- Step 10: analytic-reference CSVs can now produce a dependency-free weighted
  SVG focal-plane map annotated with centroid, D80, and throughput, plus a
  machine-readable PSF report. This does not validate production optics.

## Exact missing inputs and contracts

- Reference matrix: the selected scope is simulation-models 7.0.0, North and
  South. A sim_telarray **source checkout** is present locally, but a verified
  executable, hessio decoder revision, generated simtools config, resolved
  random seed, and checksummed shared photon blocks/commands are absent; no
  real comparison manifest can be frozen yet. This Python environment also
  lacks the dependencies needed to run simtools.
- Model data: the local 7.0.0 SCT camera response resolves from the explicit
  sim_telarray `cfg/CTA` root and is SHA-256 checked. Single-reflector nominal
  normals can be derived, but aligned
  panel normals require the selected run settings and random draws. The
  extracted camera pixels do not supply a compiled physical detector surface.
  Structure geometry and optical material semantics are also not yet compiled.

The audited 7.0.0 design matrix is LST North/South (198 panels, 1855 pixels),
MST FlashCam/NectarCam for both sites (86 panels, 1764/1855 pixels), SST South
(18 M1 segments, 2048 pixels), and SCT South (48 M1 plus 24 M2 segments,
11328 pixels). The checkout has no North SST/SCT design records. simtools'
`SimulatorRayTracing` writes a star file, generates a telescope/site config,
and runs sim_telarray with `IMAGING_LIST`, `random_state=none`, one telescope,
disabled camera filter and night-sky background, and 100000 photons per run
(5000 in test mode). It does not export the launched photon block; a matched
input comparison still needs a reproducible source block and recorded seed.
- Kernel: the generic scene represents finite planes and circular/annular
  polynomial aspheres. It still needs closed rods/caps and baffles, segmented
  asphere masks, pixel masks, material bindings, and bounded interaction
  records before it can replace the toy/segmented paths or trace a production
  model.

## External compatibility requirements

The local `../simulation-models` checkout contains the detailed, versioned
model data and is the primary import source.  A remote checkout is an allowed
acquisition route, but every production scene must record the resolved release,
parameter-record versions, asset paths, and SHA-256 hashes.  The compiler now
extracts mirror-list or segmentation footprints and camera pixel layouts and
derives nominal single-reflector panel normals. It stops before run-specific
alignment, physical detector surfaces, structure geometry, or material
semantics are established.

`../simtools` currently delegates these workflows to sim_telarray and its
`LightEmission` programs.  Their optical contracts define required obdeect
coverage:

| Workflow | Required obdeect result |
| --- | --- |
| `RayTracing`, `validate-optics`, and PSF tools | Stars at finite/infinite distance, two-dimensional off-axis scans, full-telescope and single-panel modes; focal-plane hits, D80/PSF, effective area, effective focal length, and optical-path timing. |
| Incident-angle derivation | Per-photon focal-surface, M1, and M2 hit positions and incidence angles for single- and dual-mirror models. |
| Mirror-panel PSF/RNDA optimisation | Selectable panel IDs and deterministic alignment/roughness perturbations, with per-panel PSF observables. |
| Flasher simulation | Flat-field source geometry, configured wavelength and photon-level sequences, pulse shape/offset, event mapping, and both full optical transport and direct-injection output. |
| Illuminator simulation | Positioned, oriented calibration sources; configured telescope visibility/layout; wavelength selection; finite-distance transport and timing. |
| CORSIKA transport | EventIO bunches and telescope metadata through the same compiled scene, preserving weights, wavelength sentinels, timing, and coordinate conventions. |
| Focal-plane diagnostics | Weighted focal-plane maps with x/y position as the base dimensions and optional wavelength and arrival-time axes; interactive or static projections must retain loss category, surface/pixel ID, field angle, and source provenance. |

The current LightEmission scope is limited to the programs simtools invokes:
`ff-1m` for flat-fielding and `xyzls` for illuminators.  The other
sim_telarray LightEmission executables (`fpls`, `xyls`, `ls-beam`, `pixled`,
`octo`, `nsbls`, and `fake-muon`) are outside the current implementation and
release scope.

The available production model records demonstrate that import must cover more
than mirror lists: camera configuration/filter and incidence response,
light-guide response, telescope transmission, camera-body shape, secondary
shadow geometry, focal-length/reflective perturbations, and family-specific
optics reference data.  Each field must be consumed, explicitly deferred with
a stated boundary, or rejected.

## Step-by-step implementation and validation plan

Each step has a deliverable and an exit criterion.  Complete the steps in
order; later comparisons are only meaningful when earlier geometry and frame
contracts are fixed.

1. **Freeze the reference matrix.** Record the sim_telarray release, hessio
   decoder, selected simulation-models checkout and production version, site,
   telescope variants, configuration overrides, seeds, atmosphere/extinction
   treatment, photon blocks, and coordinate frames.  The importer must support
   every compatible production version present in a supplied checkout; the
   recorded version makes a comparison reproducible without becoming a global
   version pin.  Create a separate `obdeect-tests` repository if it simplifies
   versioning large fixtures and running both engines.  Exit when a manifest
   can reproduce every reference command and verifies all hashes.
2. **Define the comparison contract.** Specify SI input/output schemas for
   photons, source descriptions, compiled scenes, arrivals, interactions, and
   terminal losses.  Include run/event/telescope/bunch provenance, weight,
   wavelength, emission time, geometric path, optical path, group time, and
   surface/component IDs.  Exit when CSV and in-memory inputs produce identical
   resolved batches and frame round trips preserve a physical ray.
3. **Complete model import.** Extend the provenance IR and compiler to resolve
   the selected production records and referenced assets for LST, both MST
   camera variants, SST, and SCT.  Version adapters by schema/field semantics,
   not by telescope name or a fixed production release.  Parse unit-bearing
   parameters, mirror/camera files, material tables, masks, and structure
   descriptions without inference.  Generate a consumed/deferred/rejected-field
   report and fail on unknown required fields.  Exit when fixtures from each
   supported schema/version have verified facet and detector counts, bounds,
   identities, and hashes.
4. **Build the generic scene kernel.** Replace the toy-specific compiled scene
   with immutable surfaces, transforms, finite facet masks, obscurers,
   material bindings, and detector surfaces.  Implement globally nearest
   intersections on every segment and bounded non-sequential tracing.  Exit
   when analytic primitive, ordering, gap, cap, and component-attribution
   tests pass before any telescope-level comparison.
5. **Deliver the single-reflector baseline.** Compile segmented parabolic LST
   and modified Davies--Cotton MST scenes including actual panel normals,
   alignment, camera choices, masts/CSS, focal surfaces, and gaps.  Add
   single-panel selection for the simtools mirror-PSF workflow.  Exit when
   on/off-axis position, D80, area, timing, shadow map, and panel-PSF fixtures
   agree with the recorded sim_telarray runs within predeclared tolerances.
6. **Deliver dual-mirror scenes.** Compile SST/SCT M1/M2 aspheres,
   segmentation, holes, masks, baffles, secondary shadowing, and curved
   detector/module geometry.  Exit when independent known rays and
   primary/secondary/focal-surface incidence-angle distributions agree with
   sim_telarray debug-trace output.
7. **Integrate transport through the detector boundary.** Bind measured
   coatings, windows/filters, Fresnel/Snell, absorption, light guides, and
   camera active-area masks to scene surfaces.  Keep optical arrival,
   concentrator loss, and sensor conversion distinct.  Exit when table-node,
   slab, Brewster, TIR, absorption, light-guide, and energy/loss-closure tests
   pass, followed by wavelength scans against reference fixtures.
8. **Implement source and input adapters.** Add production star, flasher,
   illuminator, laser, and EventIO sources to the same batch contract.
   Reproduce simtools `ff-1m` flat-field sequences and `xyzls` illuminator
   position/orientation/visibility behavior.  Exit when normalization, pulse,
   position, sample-count, launch-height, and EventIO-decoder golden tests
   pass.
9. **Expose simtools-facing tools.** Provide a stable native CLI/API for PSF,
   optics validation, incidence angles, single-panel PSF, flasher, and
   illuminator runs.  Adapt simtools only after the native outputs preserve the
   required metadata and analysis columns.  Exit when its existing workflow
   fixtures run against obdeect with no lossy format conversion.
10. **Deliver focal-plane analysis and plots.** Produce weighted two-dimensional
    focal-plane hit maps for every supported source and scene.  Support
    wavelength and arrival-time binning/slicing so the same result can be
    inspected as x/y, x/y/wavelength, x/y/time, or x/y/wavelength/time data.
    Include centroid, containment radii, throughput, timing moments,
    incidence angle, surface/pixel ID, and component-resolved losses in the
    underlying machine-readable output; plots must use the compiled detector
    geometry and expose their binning and weights.  Exit when static plot and
    data-product fixtures reproduce selected sim_telarray focal-plane,
    wavelength, and timing distributions.
11. **Make execution reproducible and scalable.** Add bounded SoA blocks,
    immutable compiled-scene reuse, validated acceleration, deterministic
    reductions, and counter-based stochastic streams.  Benchmark the scalar
    oracle before SIMD or GPU work.  Exit when photon order, block size, and
    thread count preserve results and benchmarks report photons/s, peak RSS,
    compilation, I/O, and kernel time.
12. **Run the release matrix in CI.** Compare both engines on fixed Cherenkov,
    star, flasher, illuminator, and laser fixtures across every supported
    telescope family, field angle, wavelength, and source distance.  Publish
    per-photon residuals plus D80/PSF, area, timing, incidence, camera-map,
    and loss-category residual distributions.  Freeze tolerances from the
    required science accuracy and Monte Carlo uncertainty before looking at
    candidate results.  Release only when all rows pass.

## sim_telarray verification design

Reference comparison is a feature of the implementation, not a final manual
check.  Create an `obdeect-tests` repository or equivalent versioned fixture
package that owns the reference manifests, input photon blocks, generated
sim_telarray configuration, obdeect scene IR, command lines, output hashes,
analysis code, expected results, and tolerances.  It must run both engines from
the same manifest and retain the exact software revision/container digest for
each result.

Where simtools provides the workflow, the harness must invoke its existing
implementation to run sim_telarray or its LightEmission package.  This applies
to `SimulatorRayTracing`, optics/incident-angle and mirror-panel tools,
`simulate-flasher`, and `simulate-illuminator`; it verifies the same production
entry points that users run.  Invoke sim_telarray or `ff-1m`/`xyzls` directly
only when simtools has no wrapper or does not expose the interaction-level
diagnostic needed by a fixture.  Record that exception and the direct command
in the fixture manifest.

### Matched inputs and conventions

For each fixture, the harness must make the following values identical or
record why they differ:

- simulation-models checkout, production version, selected telescope and site,
  parameter overrides, and every resolved input-file hash;
- telescope, camera, source, and array frames; position/direction handedness;
  mirror numbering; time origin; length, wavelength, and angle units;
- photon position, direction, wavelength, emission time, weight, and source
  identity; for CORSIKA, run/event/telescope/bunch IDs, reuse weight,
  emission data, and extinction convention;
- source geometry and sampling: star distance and field angle, flasher
  location/pulse/photon levels, and illuminator position/orientation/visibility
  and wavelengths;
- mirror alignment, focal-length/reflective perturbations, roughness, random
  seeds, and random-number policy; use deterministic settings wherever the
  reference supports them;
- enabled optical elements and their tables: mirrors, windows/filters,
  light guides, camera masks, pixel mapping, masts, baffles, M2, holes, and
  atmosphere.

The harness must reject a comparison when any required convention is absent or
when one engine applies a scalar transmission/extinction that the other already
applied explicitly.

### Test layers and required fixtures

| Layer | Fixture and comparison | Acceptance evidence |
| --- | --- | --- |
| Import | Every supported production/schema version and LST, both MST cameras, SST, and SCT | Exact model/parameter/asset hashes; field-coverage report; facet, detector, pixel, and obscurer counts; bounds and ID uniqueness. |
| Geometry oracle | Hand-authored rays for panels, gaps, camera caps, rods, M1/M2, holes, baffles, curved detector, and pixel boundaries | Exact terminal component, surface ID, intersection position, normal, reflected direction, and path length against analytic expectations. |
| Single-panel optics | simtools single-mirror runs at panel centre/edge and selected off-axis angles | Panel hit/miss, focal hit, D80, focal length, and reflection/roughness behavior compared with sim_telarray. |
| Full telescope star | On-axis and two-dimensional off-axis star grids at near and far source distances | Per-photon focal hits where deterministic pairing is available; otherwise weighted x/y residual maps, centroid, D50/D80, area, timing, and loss fractions. |
| Dual-mirror transport | SST/SCT rays that hit/miss M1 and M2, holes, baffles, and focal surface | Ordered interaction IDs, M1/M2/focal incidence angles, focal coordinates, geometric/optical path, and terminal losses. |
| Spectral transport | 300--700 nm nodes and interpolated wavelengths at normal, Brewster, and grazing incidence | Mirror/filter/window/light-guide transmission, detector arrival weight, wavelength-resolved area, and energy/loss closure. |
| `ff-1m` flat-fielding | simtools photon-level sequences, wavelengths, pulse settings, and full optical transport | Focal-plane illumination map, per-pixel/photoelectron-ready arrival statistics, pulse-time distribution, total throughput, and direct-injection contract where applicable. |
| `xyzls` illuminator | simtools layout/visibility fixtures at each supported illuminator position, orientation, and wavelength | Selected-telescope set, focal-plane x/y/wavelength/time maps, arrival time, weighted throughput, and shadow/loss distributions. |
| EventIO boundary | Standard, compact, and 3D CORSIKA bunch fixtures plus malformed/truncated inputs | Decoded metadata and resolved photons versus an independent decoder; identical obdeect trace after CSV and EventIO adaptation. |
| Reproducibility/performance | Photon order, block size, and thread-count permutations; representative shower-scale blocks | Bitwise identity for deterministic outputs or declared statistical equivalence for stochastic output; bounded memory and benchmark record. |

Use the existing simtools `SimulatorRayTracing`, incident-angle, mirror-panel
PSF, `simulate-flasher`, and `simulate-illuminator` fixtures as the initial
reference command set.  For focal-plane comparisons, obtain sim_telarray
debug-trace/imaging-list data when it exposes the needed interaction; otherwise
compare aggregate data products from identical input blocks and state the lost
per-photon observability in the manifest.

### Comparison procedure

1. Compile/import one scene and generate one resolved photon block from the
   fixture manifest.  Write both forms as immutable, checksummed artifacts.
2. Run simtools' existing wrapper for sim_telarray or LightEmission whenever it
   supports the fixture; otherwise run the recorded direct command.  Run
   obdeect with the same resolved block and compiled scene.  Capture stdout,
   configuration, environment revision, seed data, and terminal-loss summary.
3. Normalize outputs into a common SI comparison table with one row per
   resolved photon or aggregate bin.  Include input ID, terminal status,
   interaction sequence, focal x/y, wavelength, arrival time, path lengths,
   weight, surface/pixel ID, and loss component.
4. Compare deterministic geometry first: IDs/statuses, then positions,
   directions, normals, paths, and times.  Stop the fixture at the first
   divergent interaction and retain a ray report that identifies both scenes'
   preceding interactions.
5. Compare stochastic or unpaired samples as weighted distributions.  Require
   bin-integral/loss closure before comparing focal-plane x/y, wavelength, and
   time projections.  Report weighted residual maps, CDF differences for D80
   and timing, and per-loss-component residuals.
6. Write machine-readable results and static plots for every fixture.  A
   review must be able to inspect x/y, x/y/wavelength, x/y/time, and
   x/y/wavelength/time slices without rerunning the engines.

### Tolerances and failure policy

Use exact equality only for fixture identity, integer IDs, terminal categories,
and deterministic configuration fields.  Set geometry tolerances from the
conditioning of each primitive and the common binary64 implementation, before
examining the comparison result.  Set aggregate observable tolerances from the
required science accuracy plus quantified Monte Carlo uncertainty, using a
fixed confidence rule and effective weighted sample size.  Store both absolute
and relative limits, the sample size, and the uncertainty estimate in the
fixture manifest.

Do not loosen a tolerance, replace a golden result, or exclude a loss category
to make a failed comparison pass.  A failure must retain the model/asset hashes,
source block, normalized records, residual plots, terminal-loss table, and the
first divergent ray when available.  CI should run primitive/import fixtures
on every change; run a compact reference matrix on pull requests and the full
matrix on scheduled or release builds.

## Remaining work

### P0: production optical scenes

1. Build a versioned scene compiler that consumes selected model assets into an
   immutable, observatory-neutral scene IR.  It must compile transforms,
   facet geometry and boundaries, normals, focal lengths, alignment,
   focal-surface/camera geometry, structures, masks, material tables, and
   coordinate frames.  Preserve source parameters and fail on unsupported
   fields.
2. Implement real segmented LST and modified Davies--Cotton MST scenes,
   including facet curvature, gaps, missing facets, camera variants, alignment,
   CSS/structure shadows, and physical focal surfaces.
3. Implement SST/SCT M1/M2 segmentation, aspheres, holes, masks, baffles,
   secondary shadows, transforms, and curved aligned focal surfaces.
4. Generalize nearest-hit, bounded non-sequential transport to every compiled
   scene.  Test obstructions and material interactions on every path segment,
   and attribute all terminal losses.

### P1: transport, interfaces, and analysis

1. Integrate wavelength-, angle-, and polarization-dependent coatings and
   refractive geometry.  Add filters/lenses, dispersion and group delay, then
   separate detector-surface arrival, concentrator loss, and sensor conversion.
2. Replace the fixed four-point path record with bounded interaction diagnostics
   recording interaction ID, position, direction, wavelength, time, weight,
   and terminal status.
3. Complete configurable source spectra, pulses, pointing, yield, and
   beam-normal launch.  Specify a separate atmospheric propagation/scattering
   model before supporting laser side scatter.
4. Add an optional CORSIKA7 EventIO/hessio adapter behind the existing
   format-neutral photon contract.  Preserve all bunch metadata and unresolved
   wavelength sentinels; compare standard, compact, and 3D blocks with an
   independent decoder.
5. Stabilize the native scene/result contract, then expose C++ bindings and a
   Python API that execute the native tracer and preserve loss semantics.
6. Extend PSF analysis to two-dimensional field maps, material/wavelength
   scans, structure-shadow fractions, and primary/M2 incidence angles.  Add
   focal-plane x/y/wavelength/time data products and plots, with selectable
   projections and bins that preserve weights, pixel/surface IDs, timing, and
   loss categories.

### P2: reproducibility, performance, and presentation

1. Define weighted and stochastic tracing modes.  The toy parallel source
   currently uses a low-discrepancy golden-ratio sequence; replace it where
   stochastic sampling is required with a seeded, high-quality generator.
   Production stochastic streams should be counter based and keyed by
   run/event/telescope/photon/interaction/stream so results survive ordering,
   batching, and thread-count changes.
2. Compile immutable geometry and materials once; process bounded reusable
   SoA batches; add validated spatial acceleration; then add deterministic
   parallel execution.  Benchmark throughput, peak RSS, compilation, I/O, and
   kernel time separately.
3. Rename the `toy` executable/configuration/documentation as an
   artificial-light-source reference tracer if compatibility permits.  Render
   actual compiled obscuring geometry in diagnostics.
4. Strengthen tests with off-axis, grazing, edge, malformed-input, full-render,
   benchmark, and reference-observable coverage.  Run the project validation
   matrix before release.

## Production replacement gate

A sim_telarray/ROBAST-level claim requires each supported LST, MST, SST, and
SCT model to pass reproducible fixtures using the same photon blocks, source
wavefront, wavelength range, and atmospheric convention as the nominated
reference.  CI must retain toleranced comparisons for focal-plane positions,
PSF/D80, optical area, timing, material response, camera mapping, and each
terminal loss (facet gaps, M2, camera, masts, baffles, and holes).  Report
residual distributions and error budgets.  Do not state production
compatibility until these checks pass.

## Test identifiers

`T-OBS-004` is the closed-solid regression for a ray that reflects from M1 and
then enters the rear of the toy camera.  Related toy-obstruction identifiers
are `T-OBS-001` through `T-OBS-005` in `cpp/tests/test_core.cpp`; they are test
labels, not runtime requirements or external standards.
