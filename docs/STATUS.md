# obdeect status and remaining work

Reviewed 30 September 2026 against obdeect commit
`153b9e04bb0a98190b6be57093af3b145a90a733` and the local sibling simtools code.
This is the single current status and implementation backlog. Earlier plans,
checklists, and reviews are not part of the active documentation.

**Conclusion:** obdeect is a runnable nominal optical-geometry prototype. It
is not yet a validated CTAO production backend. LST/MST curved panels exist;
SST/SCT use continuous dual aspheres. Complete optical responses, structure,
source compatibility, shared-input comparisons, and performance parity remain
unfinished. The reviewed simtools integration is broken at its command and
configuration boundaries.

## Goal and boundaries

Trace optical photons through IACT optics, ending at a physical detector surface.
Support stars, nearby dish illuminators, calibration/alignment lasers, and
resolved CORSIKA photon bunches. Include mirrors, obscurers, camera bodies,
windows, filters, and light guides where they affect transport. Shower generation,
electronics, triggers, and event writing are outside this package.

The original requirements remain:

- Match sim_telarray accuracy for equivalent inputs, physics, and observables.
- Achieve at least sim_telarray CPU throughput for equivalent validated work.
- Derive CTAO optics from an explicitly selected simulation-models production.
  Version 7.0.0 is the first reference fixture, not an importer version limit.

Keep the C++20 core portable and independent of Python and external simulators.
Use binary64 geometry, explicit units/frames, immutable compiled model data,
bounded interactions, deterministic identities, and no per-photon allocation
or I/O in the transport kernel. CSV writing is a boundary cost to measure
separately. Keep telescope-specific import semantics outside the generic core.

sim_telarray remains the default simtools backend and compatibility reference.
ROBAST supplies independent detailed-geometry/window comparisons; IACTrace has
an existing simple spherical-mirror comparison script. Optional ROBAST import,
CAD/Embree and GPU work follow validated CPU functionality and use the same
generic model/result contracts.

## What exists today

“Runnable” means the code path exists. “Primitive only” means tested code is
not used by the imported-model command. Neither means production agreement.
Paths below are relative to obdeect unless explicitly labelled simtools.

| Capability | Reviewed state | Evidence |
| --- | --- | --- |
| C++ scalar kernels | Runnable; 14 configured Debug tests pass locally | `CMakeLists.txt`, `cpp/tests/` |
| Model import | Selected production, parameter/asset hashes, safe paths, nominal geometry, deferred-field reporting; schema coverage incomplete | `python/obdeect/model_import.py`, `optical_model_compiler.py` |
| LST/MST optics | Finite spherical panels with radius `2*f`, nominal centres/normals, circular detector envelope, supplied cylinders, 1D mirror reflectivity | `cpp/include/obdeect/optical_model_file.hpp`, `segmented_optical_model.hpp` |
| SST/SCT optics | Continuous M1/M2/focal aspheres, holes, 1D mirror responses; segment footprints do not constrain native transport | `optical_model_compiler.py`, `optical_model_file.hpp` |
| Sources | Plane-wave stars, finite isotropic illuminator with solid-angle weights, geometric laser, discrete wavelengths, top-hat timing | `cpp/include/obdeect/sources.hpp`, `cpp/src/simtools_raytrace_main.cpp` |
| Material/atmosphere physics | Primitive only: Snell/Fresnel, TIR, absorption, parallel slab, atmospheric helpers | `cpp/include/obdeect/materials.hpp`, `atmosphere.hpp` |
| Photon input | Memory/CSV and optional EventIO readers; public imported-model command generates its own photons | `photon_input.hpp`, `cpp/src/eventio_photon_input.cpp` |
| Results/analysis | Arrival-v1 CSV, weighted native PSF/D80, distributions, normalization; incomplete interaction/provenance contract | `python/obdeect/result_contract.py`, `analysis.py`, `arrival_normalizer.py` |
| Plotting | Actual recorded paths and compiled panel/focal views; hardware absent from the model cannot be shown | `python/obdeect/plotting.py`, `photon_distributions.py` |
| simtools backend | Selector, runner and consumers exist; runner cannot execute the current native command | Sibling `simtools/src/simtools/` |
| Reference runs | Derived sim_telarray PSF/CDF products and reference-run records; no completed matched production comparison established | `docs/reference/7.0.0/`, `python/obdeect/reference_run.py` |
| Packaging | Native executables in wheels, CMake install and CI workflows; PyPI publishing restricted to version tags | `pyproject.toml`, `.github/workflows/` |
| Python array API | Missing: Python launches executables; no nanobind bulk tracing module | `python/obdeect/cli.py`, `cpp/bindings/` |
| Performance parity | Unproven: exhaustive facet scans, scalar tracing, full-batch source allocation, no parity benchmark | `segmented_optical_model.hpp`, `sources.hpp` |

The generic `CompiledOpticalModel` supports framed planes, polynomial sag
surfaces and bounded nearest-hit tracing. It is not connected to the importer
or public model-derived command. Its `material_id` does not apply a material
response. This is not complete CTAO nonsequential transport.

## Verification and review limits

- Configured/built the Debug preset and ran CTest: **14/14 passed**. EventIO
  was disabled in this configuration.
- Traced the existing local `lst.optical-model.json` with 1,000 on-axis star
  photons: **768 detected**. The artifact was not regenerated or compared
  with sim_telarray in this review; this checks execution only.
- Replayed simtools' `--scene-file`: native command rejected it with exit 2.
  Its supported option is `--optical-model`.
- Changed a detector position in a temporary model without updating its hash:
  native loading and tracing succeeded. Hash-format validation does not verify
  compiled content integrity.
- A temporary C++ probe confirmed task 3: photon arrival weight was 0.8 but
  summary detected weight was 1.0 for a mirror with reflectivity 0.8.
- A temporary analytic probe confirmed an asphere limitation: a horizontal
  ray from `(0,0,1)` intersects `z=r²/4` at `(2,0,1)`, but
  `intersect_axisymmetric_mirror()` returned no hit. See task 8.
- Reviewed Python code/tests statically. Python tests, pre-commit, wheel builds,
  external comparisons and benchmarks were not rerun. Available Homebrew Python
  versions stop at 3.14; `AGENTS.md` asks for `>3.14`, while packaging/CI specify
  `>=3.14`/3.14. Resolve this inconsistency in task 12.

The archived review reports a previous five-model, 1,000-ray smoke run:
LST 767 hits, MST FlashCam 686, MST NectarCam 684, SST 123, SCT 0. These are
historical observations, not current regression targets. SCT reportedly reached
M1/M2 but missed the focal aperture. Its prescription/frame issue is unresolved;
do not tune geometry or widen the detector to reproduce a desired hit count.

Stored 7.0.0 sim_telarray products use 50,000 photons, a 10 km source,
20-degree zenith, zero field offset and seed 19780503. They disable camera
filter losses and set camera transmission to one, including focal crossings
that miss active pixels. These are optical PSF references, not full-camera
throughput or obdeect equivalence evidence. Incident-angle and field-scan
reference coverage remains incomplete.

## Ordered implementation tasks

All tasks are open. A handoff implements one task and its necessary tests,
preserves unrelated work, and updates this file with command/result evidence.
“Luna-sized” marks bounded starting work. Larger physics tasks need several
reviewed changes; their headings are not single-change assignments.

### 1. Repair the simtools command boundary — first, Luna-sized

**Defect:** sibling `simtools/simtel/simulator_obdeect.py` and
`simtools/ray_tracing/incident_angles.py` send `--scene-file`; native
`cpp/src/simtools_raytrace_main.cpp` accepts `--optical-model`.
`RayTracing._create_simulator()` also omits the model file required by
`SimulatorObdeect.__init__()`. The runner docstring names a removed compiler API.

**Implement:** pass the selected compiled model through both workflows, use
the current CLI option, and correct stale API descriptions. Resolve executables
through the installed package; introduce no path override or analytic fallback.
Keep the single-mirror guard until task 9 connects panel/screen options.

**Accept:** tests invoke the actual installed executable with a tiny model
through both workflows and validate arrivals. Missing input fails clearly;
force/test behavior and unique offset outputs remain correct; sim_telarray stays
default. Mocked command tests alone do not qualify the integration.

### 2. Make arrival parsing consistent — Luna-sized

**Defects:** `result_contract.read_arrivals()` allows duplicate IDs, a detected
row containing only its source vertex, and nonzero throughput on a lost row.
It accepts `escaped_optical model`, while C++ emits `escaped_optical_model`.
`analysis.py` and `plotting.py` have separate, partly different validation.

**Implement:** one version-aware validated reader, matching terminal enums,
unique identities within a declared batch, and status/vertex/weight invariants.
Define legitimate zero-response detector arrivals explicitly. Retain named
adapters for diagnostic formats where necessary; do not guess missing fields.

**Accept:** malformed-row tests cover these defects; native results, PSF,
distributions, normalization and plotting agree on acceptance/rejection.

### 3. Repair optical weight summaries — Luna-sized

**Defect:** segmented `cpp/include/obdeect/trace.hpp` writes
`input_weight * reflectivity` to the photon result but passes unattenuated input
to `TraceSummary.add()`. Its `detected_weight` is therefore wrong for lossy
mirrors. Terminal input-weight buckets do not explain partial absorption.

**Implement:** separate incident, arriving and lost optical weights, retaining
terminal counts and explicit partial-response loss accounting. Match C++ and
Python definitions.

**Accept:** unit incident weight and reflectivity 0.8 produce arrival 0.8 and
accounted loss 0.2; mixed detected/blocked/missed samples close the weight ledger.
Test count closure separately from weight closure.

### 4. Verify model integrity and readiness — Luna-sized first slice

**Defects:** native loading checks hash format without recomputing it and ignores
readiness/blockers. The compiler sets nominal `trace_ready`, never emits
`native_trace_ready`, and always adds production blockers. The validator requires
`native_trace_ready`, `primary.facets`, `detector` and `materials`; these do not
match the emitted artifact structure. The importer silently skips names outside
its allow-list, including any future optical parameter.

**Implement:** one schema and canonical hash policy for compiler/loader/validator;
separate nominal traceability from production completeness. Allow incomplete
models only in an explicitly labelled nominal workflow. Classify every source
parameter as applied, outside optical scope, or unsupported. Reject unknown
required optical semantics; record absent/not-applicable components explicitly.

**Accept:** stale-hash modifications are rejected, valid compiler output reloads,
nominal models cannot pass production checks, and schema/field errors name their
cause. First slice: hash round trip/tamper tests; do not remove physics blockers.

### 5. Strengthen the production comparator — Luna-sized initial checks

**Defects:** `production_validation.compare()` can pass an all-loss sample;
its reader accepts arbitrary nonempty statuses. Weight, wavelength, output
direction, surface identity and provenance are not compared. Matching losses
bypass numerical checks. “First divergence” names a photon/field, not an optical
interaction. CLI failure does not write the requested summary file.

**Implement:** fixture-declared detection/surface coverage, valid statuses/ranges,
source/model hash binding, weight/loss checks and persistent failure reports.
Allow legitimate all-loss fixtures only when declared. Add interaction comparison
after task 6; define comparable loss observables without requiring arbitrary
termination points to agree.

**Accept:** all-loss positive-acceptance fixtures, matching invalid statuses,
wrong weights and provenance fail; intentional negative fixtures pass their
declared criteria. Failed gates retain useful diagnostics.

### 6. Add shared photon replay and interaction records — foundational

**Gaps:** public CLI always generates sources despite existing readers. Output
holds at most four positions, without surface/component IDs, normals, directions,
pixel IDs or run/event/telescope/bunch identity. Incoming obscurer losses omit
travelled distance/end point; outgoing losses add distance without the hit vertex.
Secondary misses collapse into `missed_screen`. Normalization infers vacuum time.

**Implement:** versioned photon-batch/result contracts with explicit units/frames,
stable identities, weights, geometric/optical paths, actual arrival timing,
terminal component and ordered interactions. Connect memory/CSV input to the
model-derived CLI; sources become another input adapter. Share records between
CLI and library kernels before consolidating their duplicated trace loops.

**Accept:** frozen blocks replay unchanged, IDs/weights survive, every recorded
path ends at its actual interaction, and M1/M2/focal comparisons need no vertex
guessing. Reordering/chunking fixed blocks preserves results. Sampled diagnostics
have count/byte caps; default tracing does not retain full shower trajectories.

### 7. Complete single-reflector geometry and response

**Gaps:** spherical panels already exist. Missing are run-specific alignment,
distance/focal-length errors, roughness, degradation maps, camera housing,
quadrilateral/support/baffle geometry, physical pixel boundaries, filters/light
guides, and wavelength/incidence responses. A circular envelope is not an active
pixel map; centre-only layouts lack physical edges.

**Implement:** bind these model fields for LST and both MST cameras, with keyed
RNG/seed settings. Separate detector arrival, pixel gaps, concentrator losses
and later sensor conversion. Declare `simtel_compat` scalar transmission versus
`physical_3d` shadowing and reject double application. Windows use task 10.

**Accept:** isolated panel/boundary/perturbation/table/obscurer tests precede
full-dish tests; zero perturbations reproduce nominal geometry. Every optical
loss has a reason/component. Matched LST/MST fixtures pass task 11.

### 8. Complete dual-reflector optics; investigate SCT first

**Gaps:** continuous aspheres use fixed M1→M2→focal tracing. Segment boundaries,
gaps, transforms/errors, incoming secondary shadows, supports/baffles on every
leg and detector pixels are missing. SCT zero detections remain unexplained;
the production validator has no SCT family option. The asphere solver uses a
single Newton start at the vertex plane, rejects horizontal rays, and has no
bracketed nearest-root guarantee. Numerical failure is returned as an ordinary
miss, making physical and solver failures indistinguishable.

**Implement:** first isolate SCT prescription, reference-radius, sign/frame and
focal-boundary conventions using identical reference rays. Then bind finite M1/M2
segments and structures. Keep SST and SCT acceptance separate.
Make supported intersection domains explicit; add safeguarded root selection
and distinguish nonconvergence from a geometric miss where required.

**Accept:** report the first SCT divergent interaction without a manufactured
aperture fix. Independent monolithic asphere, grazing/multiple-root/near-parallel,
segment-edge, secondary-shadow and 2D field-map tests precede production
comparisons. Include source-distance/scale sweeps to check numerical precision.
SST passing does not qualify SCT.

### 9. Complete source semantics and simtools studies

**Gaps:** star `--distance-m` changes launch height, not wavefront curvature;
stars remain parallel. Laser launch radius equals the telescope pupil radius,
without independent beam size. Wavelengths cycle evenly rather than following
a model spectrum. Source profiles depend on total count; seeded/block-stable
production sampling, zenith/frame handling, model-defined sources and atmosphere
are missing. Native panel/test-screen options are not connected in simtools.

**Implement:** finite/infinite star modes, independent laser profile, model-defined
illuminator pointing/angular/spectral/time profile and emitted level; correct
frames, visibility and extinction exactly once. Atmospheric laser scattering
is a separate explicit model. Forward supported settings, panels and screens
through simtools.

**Accept:** finite-source conjugate/infinity-limit tests, absolute illuminator
normalization, beam-size and weighted spectral/time tests, seed/block invariance.
Compare qualified `xyzls`, `ff-1m`, `ls-beam` cases with matching semantics;
a geometric beam alone cannot establish atmospheric laser parity.

### 10. Bind physical windows and concentrators

**Gaps:** material primitives are disconnected from native imported-model tracing.
No production curved two-surface window, medium tracking, group delay or bounded
reflected concentrator path exists. Slab material transmission currently multiplies
measured response by Fresnel and bulk losses, risking double counting if the table
already describes the complete window.

**Implement:** model-derived interfaces/materials, Fresnel/absorption, medium
transitions and separate geometric/optical/group paths. Declare included losses
for measured tables. Use tabulated light-guide compatibility first; detailed
geometry and reflected branches require explicit selection and weight/depth bounds.

**Accept:** Snell/Fresnel/TIR, slab and real curved-interface tests, flat/zero/index-one
limits, energy closure and ROBAST comparisons. Filter-only sim_telarray output
cannot validate refractive spot shifts.

### 11. Establish scientific comparisons and backend-neutral analysis

**Gaps:** stored products are not a shared-input harness. Native weighted PSF
works, but simtools' `ray_tracing/psf_analysis.py` rejects nonuniform source
weights while ignoring varying throughput for accepted rows. Its count-based
PSF can therefore be biased even when source weights are equal; its default
area is a detection fraction with unit launch area. Production area, focal
length, panel PSF/RNDA and complete incidence scans are unqualified.

**Implement:** replay task 6 photons in both engines; provide a sim_telarray
adapter and recorded command executor. Freeze model/software/asset/source hashes,
frames, atmosphere, seeds, physics, observable definitions and tolerances before
acceptance runs. Use shared records for PSF, area, focal length, panels and
incidence. Area needs declared launch area/normalization, not just a hit fraction.

**Accept:** LST, both MST cameras and SST pass on/off-axis 2D, panel, spectral,
source/time and incidence matrices; SCT passes separately. Retain residuals,
component loss closure, weighted CDFs, uncertainties and first-divergent interactions.
Keep derived products in the repository; large raw inputs stay external with
hashes. Earlier tolerance examples are proposals, not accepted limits.

### 12. Complete test and installed-command checks — Luna-sized

**Gaps:** pytest functions coexist with `unittest.TestCase`; old `unittest discover`
instructions miss function tests. CI uses pytest. Python version policy conflicts.
Wheel smoke checks only analytic fallback. Release macOS metadata names obsolete
`obdeect_toy`/`obdeect_ctao`. EventIO is disabled in wheels and optional CMake
EventIO targets are not included in the current installation rules.

**Implement:** one full-suite command, consistent Python policy, accurate command
inventory/help, tiny compiled-model wheel smoke, explicit EventIO distribution
policy, source-independent CMake consumer and optional flag-on/off checks.

**Accept:** clean installs run every advertised `--help`, trace a small generic
compiled model and collect all Python tests. EventIO support/absence is reported
accurately. Separate local evidence from CI configuration. README edits remain
suggestions unless explicitly authorized.

### 13. Add Python bulk API, streaming and performance evidence

**Gaps:** no nanobind array API; CLI allocates full source vectors, traces scalarly,
scans all facets and writes every row. Candidate grid/BVH, thread/SIMD dispatch
and parity benchmarks are absent. Segmented model validation allocates an ID
set per block call.

**Implement:** stabilize contracts, then bind owned/borrowed arrays with lifetime/
contiguity tests, stream fixed blocks including EventIO, and validate immutable
models once. Measure exhaustive tracing before acceleration; retain it as the
nearest-hit oracle. Add threads/SIMD only with result reproducibility. Evaluate
direct simtools Python execution once the actual bulk API can be measured.

**Accept:** order/block/thread invariance, bounded RSS, no transport-loop allocations,
ownership checks and result checksums. Benchmark compile/cache, source/decode,
kernel, binding and I/O separately. Record host/compiler/flags, matched physics/
outputs/threads, warm-up and at least ten repeats with median/spread. Claim speed
parity only after optical parity for the same work.

### 14. Extend detailed geometry and diagnostics

**Gaps:** full assemblies, CAD import/BVH, component shadow maps and labelled
sampled interactions are not available in production transport.

**Implement:** start with a static four-panel diagnostic: orthographic assembly,
axial section, pupil polygons and focal geometry, sharing validated bounds and
provenance. Render existing cylinders and M1/M2/focal aspheres directly from
compiled data. Preserve current CLI views. This geometry-only first slice is
Luna-sized and can proceed independently of new physics; component-loss maps
depend on task 6. Use [the detailed plotting design](../TELESCOPE_PLOTTING_PLAN.md)
as a specification, with this document governing priorities and completion.
Then add validated analytic components and an audited transcription
of a small supplied ROBAST MST structure subset. Add optional meshes with unit/
frame/normal/degeneracy checks and hashes. Plot actual compiled geometry and
sampled interactions with wavelength/time/component filters and storage caps.
An interactive viewer is optional.

**Accept:** analytic/mesh nearest-hit agreement, ROBAST shadow maps, component
area attribution and clear missing-hardware labels. Display tessellation cannot
change physics. Preserve literature-based DC timing/off-axis, SC and SST-window
studies when exact inputs are reconstructable; otherwise label them exploratory.

## Work order and release gate

Start with **1–5** as separate Luna handoffs. Task 6 establishes contracts for
physics review; tasks 7–10 implement transport/sources. Task 11 grows alongside
each feature. Task 12 can proceed independently. Tasks 13–14 follow stable
contracts and measured needs. Analytic kernel/primitive tests precede scene,
cross-tool and performance tests; never relax tolerances or replace golden
products to hide a discrepancy.

A scientific release requires complete parameter mapping, interaction/loss
contracts, matched telescope/source/observable matrices, independent detailed
geometry checks, clean installed examples and reproducible performance evidence.
Unit-test success, plausible spots, model import, or cleared readiness flags
do not satisfy this gate.

## Documentation policy

Maintain this file for status, priorities, blockers and task completion.
`AGENTS.md` remains agent instructions; README is the entry point; CHANGELOG
is release history. Recipes, contracts and reference metadata retain technical
details rather than independent backlogs. The review covered the five root
Markdown files and all Markdown under obdeect, including the plotting design
that appeared during the review. That concurrent file was preserved as a detailed
specification; task 14 carries its current status. Archived claims are historical;
current code and fresh evidence govern status.

Suggested README correction: `--telescope MST` is an analytic installation smoke,
not a model-derived CTAO trace. Prefer a tiny compiled-model example and a link
here when README editing is authorized. The README was preserved in this review.
