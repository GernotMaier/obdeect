# obdeect current status

Updated 4 October 2026. This document describes the present capability, the evidence available in this checkout, and the work still required for scientific production use. It deliberately does not retain completed task history.

## Release position

obdeect can trace nominal, immutable optical models and replay a frozen CSV photon batch deterministically. It is **not a validated CTAO production backend**. sim_telarray remains the default simtools ray-tracing backend and the compatibility reference until the comparison matrix below is completed.

The present implementation is suitable for optical-kernel studies, model-import inspection, diagnostic plots, and controlled development comparisons. It must not be used to produce unqualified CTAO performance results.

## Available capability

| Area | Current capability | Evidence |
| --- | --- | --- |
| Model integrity | Canonical JSON SHA-256 is checked by native loading and by the production gate. Every imported parameter is classified as applied, outside optical scope, or unsupported. | `model_import.py`, `optical_model_compiler.py`, `optical_model_file.hpp` |
| Nominal LST/MST transport | Finite spherical facets, finite detector surfaces, cylinder obscurers, mirror response, shared scalar transport, and explicit optical loss accounting. | `segmented_path.hpp`, `trace.hpp` |
| Nominal SST/SCT transport | Bounded M1/M2/focal aspheres with the secondary frame transformed correctly, finite segment masks, incoming-M2 shadowing, angular mirror response, distinct secondary misses, and safeguarded nearest-root intersections. | `axisymmetric_optics.hpp`, `axisymmetric_segments.hpp`, `optical_model_file.hpp` |
| Sources | Replayable CSV photons; deterministic source IDs; plane-wave and finite-distance stars; point illuminators; independent laser beam radius/divergence; discrete spectra and emission profiles. | `photon_input.hpp`, `source_sampling.hpp` |
| Result contract | Stable identities, actual terminal position/direction/surface, arrival time, throughput, response and terminal losses, plus bounded interaction sidecars. | `result_contract.py`, `interaction_csv.hpp` |
| Analysis | Weighted focal-plane centroid, CDF and containment; declared-launch-area effective area; native and simtools consumers use arriving optical weights. | `analysis.py`, sibling `simtools` PSF analysis |
| Diagnostic plots | Rendered compiled geometry, path overlays, loss maps, focal-plane views, deterministic selection, SVG/PDF/PNG output, and explicit geometry-coverage labels. | `plotting.py`, `USE_CASES.md` |
| Python API | An immutable nanobind bulk API accepts contiguous NumPy arrays, releases the GIL while tracing, returns owned arrays, and preserves order/block invariance. | `cpp/bindings/core.cpp` |
| simtools integration | Both ray-tracing and incident-angle workflows pass a compiled optical model to the installed native command; sim_telarray remains the default. | sibling `simtools` runner tests and `test_simtools_installed_integration.py` |
| EventIO | Source builds can enable the EventIO reader and summary tool. Portable wheels intentionally exclude EventIO. | `CMakeLists.txt`, `FindEventio.cmake` |

`native_trace_ready` means that a nominal native trace can run. It is distinct from `production_trace_ready`; the compiler leaves the latter false whenever required geometry or optical response remains unbound. The native `--require-production-ready` option and the production comparator enforce that distinction.

## Verification currently recorded

- Debug CMake build and native tests: 20 configured tests pass, including installed CMake consumer, replay/chunk invariance, source sampling, optics, model hashes, and interaction diagnostics.
- EventIO-enabled source build and tests: 21 configured tests pass locally.
- Python source suite: 185 tests and 45 subtests pass; one optional installed-simtools test is skipped when its dependencies are absent. The installed-native integration itself has passed both simtools workflows in the simtools environment.
- Python formatting, linting, Towncrier fragment validation, and the repository pre-commit hooks pass after the current edits.
- An installed wheel smoke test invokes every public command and traces a tiny hashed generic optical model. It verifies an arriving mirror weight of 0.8.
- No matched sim_telarray photon-by-photon CTAO comparison has passed yet. Archived 7.0.0 PSF products remain diagnostic references, not acceptance data.

## Production requirements still open

The [five-model smoke record](review/nominal-optics.json) retains exact commands,
model-source revision and hashes, compiled optical-model hashes, and terminal
counts for 1,000 seeded photons per 7.0.0 model. Detected counts are LST 727,
MST FlashCam 336, MST NectarCam 691, SST 115, and SCT 599. These exercise the
nominal workflow; they are not sim_telarray acceptance results. Model compilation
retains measured reflectivity metadata; uncertainty columns are not applied as
random optical perturbations or interpreted as propagated errors.

1. Bind complete production optics for every telescope family: alignment and degradation, finite mechanical obscurers and baffles, physical pixel and concentrator boundaries, camera housing, windows, filters, and angle-dependent responses. Avoid applying a compatibility response together with the same physical loss.
2. Bind model-derived window and concentrator geometry to the generic material transport. Native and bulk tracing support bounded curved dielectric interfaces, phase/group indices, absorption and TIR; imported telescope models still lack the required finite geometry and material tables.
3. Establish the SCT prescription and frame convention with independent reference rays, then compare a full SCT field map separately from SST.
4. Extend the shared-photon sim_telarray diagnostic adapter to record complete terminal losses and downstream pixel/concentrator response. Its four real optical classes replay identically across repeated runs, but raw lost rows explicitly lack unavailable terminal fields. See [SIMTEL_REPLAY.md](SIMTEL_REPLAY.md) and its retained run record.
5. Run and retain the predeclared validation matrix for LST, both MST cameras, SST and SCT: on/off-axis two-dimensional fields, panel/gap/shadow scans, spectral and timing scans, finite sources, effective area, and loss closure.
6. Perform independent detailed-geometry checks with ROBAST for windows and structures. Report residuals, uncertainty studies, source/model hashes and first divergent interactions.
7. Execute the hash-bound [benchmark harness](BENCHMARK.md) only after optical equivalence: ten or more repetitions, matched physics and outputs, kernel/I/O separation, memory limits and block/thread reproducibility. No sim_telarray throughput-parity claim exists yet.

## Required comparison record

Use [REFERENCE_RUN.md](REFERENCE_RUN.md) to freeze and execute exact commands, input hashes, model selections, frames, atmosphere and seeds. Use [PRODUCTION_VALIDATION.md](PRODUCTION_VALIDATION.md) to define acceptance fixtures and tolerances before running either tracer. Retain the normalized tables, weighted residuals, loss ledger, interaction diagnostics and execution logs with the run record.

The production gate rejects an unready model, unbound source/model hashes, undeclared all-loss output, missing required surfaces, invalid optical weights, and missing comparison identities. Passing unit tests or a visually plausible focal image does not qualify a telescope model.
