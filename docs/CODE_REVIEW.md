# Code review notes

Review scope: C++ trace/source/scene interfaces, Python import/analysis/plot
commands, tests, examples, and documentation at this worktree's base commit.

## Fixed here

- The focal plot changed an explicit zero throughput into unit throughput.
- The mirror-list parser silently treated malformed optional height as zero.
- The path plot accepted non-finite vertices and invalid vertex counts.
- The README repeated its quick start, had an unclosed shell fence, and mixed
  build locations. It now points to short, runnable tutorials.
- Ray plots showed the long source flight at equal scale, hiding interactions
  near the mirror and camera. They now focus on the telescope region.

## Open design work

- The  reference, segmented, and general optical scenes each have separate trace
  loops. They cannot be deduplicated safely until they share one interaction
  record and terminal-status contract.
- Both native executables repeat CLI parsing and CSV writing. A common I/O
  layer would help when the trace output schema is stabilized.
- CTAO reference plots have no physical structure data. Full structure plots
  need a trace-ready compiled scene with supports, camera, and mirror panels.
- The general optical scene lacks curved surfaces, material bindings, and
  complete interaction diagnostics. See [STATUS.md](STATUS.md) for the
  validation plan and production blockers.

## September 2026 audit

This audit traced the installed commands through their C++ entry points, scene
compiler, arrival readers, PSF analysis, tests, and the available
`simulation-models` 7.0.0 production tables. The five-model smoke test below
uses 1,000 deterministic on-axis star samples per model. It checks that the
workflow runs; it is not a physical comparison against `sim_telarray`.

| Model | Import and native export | Detected / 1,000 | PSF derivation |
| --- | --- | ---: | --- |
| LSTN-design | Succeeded | 767 | Succeeded |
| MSTx-FlashCam | Succeeded | 686 | Succeeded |
| MSTx-NectarCam | Succeeded | 684 | Succeeded |
| SSTS-design | Succeeded | 123 | Succeeded |
| SCTS-design | Succeeded | 0 | Failed with no detected optical weight |

An additional SCT hit count found 794 M1 hits and 794 M2 hits, but no focal
hits. Sampled post-M2 rays crossed the focal vertex plane outside the exported
0.422 m focal radius. The next step is to verify the SCT prescription and
coordinate conventions against a frozen `sim_telarray` photon block. No SCT
acceptance claim follows from the successful export. All five compiled scenes have `native_trace_ready=false`
and nonempty `trace_blockers`; the production gate rejects them as designed.

### Execution and duplication map

| Path | Current use | Review conclusion |
| --- | --- | --- |
| `simtools_raytrace_main.cpp` with `scene_file.hpp` | Model-derived command used by the README | Active nominal path. It has a separate segmented and dual-asphere trace loop. |
| `reference_main.cpp` and `artificial_mst.hpp` | Developer demo and analytic tests | Retain as a diagnostic; keep it out of production claims. |
| `ctao_main.cpp`, `ctao_trace.hpp`, `ctao_models.hpp` | Analytic developer command; also the no-scene branch of the simtools tracer | Duplicate analytic route. Consolidate when the public command no longer needs that branch. |
| `optical_scene.hpp` and its `trace.hpp` overload | Generic nonsequential kernel and tests | No model importer or public command currently supplies it. Do not advertise its material/obscurer capabilities as active CTAO support. |
| `materials.hpp`, `atmosphere.hpp`, EventIO readers | Primitive tests and optional adapters | Not bound into the model-derived native CLI. Their presence does not imply coating, atmosphere, or EventIO support in that workflow. |
| `analysis.py`, `result_contract.py`, `plotting.py` | Three separate CSV consumers | Keep one versioned arrival contract and migrate consumers to it before changing the output schema. PSF analysis now rejects missing or impossible optical weights. |

The segmented primary and detector had duplicate aperture-containment code.
They now use one primitive; circular containment uses squared distance. The
star phase is covered by a common-wavefront test. For camera configuration
files with pixel entrance diameters, the circular native focal bound now
includes the outer edge of each known entrance. ECSV layouts that carry only
pixel centres still give a centre-only bound; their physical edge remains
unknown.

### Use-case gate

| Requested study | Native state | Missing before a CTAO production claim |
| --- | --- | --- |
| LST/MST star and PSF | Runnable with nominal planar panels | Panel curvature, optical response, structures, alignment, and matched reference observables. |
| SST/SCT star and PSF | Continuous dual-asphere export; SST smoke trace works, SCT smoke trace has no detections | Segment/structure geometry and a validated focal mapping; SCT requires its own first-hit and acceptance study. |
| Nearby illuminator | Generic finite point source runnable | Model-defined source position, pointing, angular/spectral/temporal distribution, and atmosphere. |
| Calibration laser | Generic geometric beam runnable | Model-defined laser configuration, attenuation, scattering, and response. |
| Throughput, shadows, incident-angle scans | Partial output fields, no qualified production workflow | Material bindings, obscurers, interaction identities, common photon blocks, and frozen comparison fixtures. |

### Next engineering steps

1. Introduce a shared interaction and terminal-result record before merging
   the three trace loops. Include surface IDs, normals, incoming/outgoing
   directions, path and time per interaction, and a reason for each loss.
2. Compile actual panel curvature and pixel/detector boundaries, then bind
   structures and wavelength/incidence response from the selected model.
   Keep unavailable fields as explicit trace blockers.
3. Feed an identical checksummed photon block to native and `sim_telarray`
   traces. Investigate the SCT zero-detection case at the first divergent
   surface. Record residuals and predeclared tolerances for each telescope and
   source class.
4. Only after those contracts are stable, consolidate the CLI parsers and CSV
   writers, remove the duplicate analytic branch from the public tracer, and
   benchmark photons/s and peak memory on the same scene and photon block.
