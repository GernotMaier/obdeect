# Review of optical transport changes

This review covers the native transport, compiled optical-model loading,
Python result consumers, diagnostic plots, and installed package boundary.
[STATUS.md](STATUS.md) records the remaining production requirements.

## Corrections and simplifications

| Finding | Result |
| --- | --- |
| Native command and library repeated segmented transport. | Both use `trace_segmented_path`, including actual terminal intersections and optical loss weights. |
| Root iteration could skip finite asphere intersections. | Bounded polynomial isolation selects the nearest accepted root and records numerical failures separately from physical misses. |
| A custom detector could reuse an obscurer ID or retain production qualification. | IDs include every component; panel and detector overrides invalidate production readiness. |
| Response and terminal losses were mixed with photon counts. | Independent arriving, response-loss, and terminal-loss weights close to incident weight. |
| Dual-reflector response export read the wrong primary object. | Compiler exports the primary mirror response and preserves complete wavelength/incidence grids. |
| The secondary polynomial used its local sag sign in the telescope frame. | The importer applies sim_telarray's secondary frame rotation; nominal SST and SCT rays now reach the detector. |
| Output paths could alias frozen inputs or each other. | Native commands reject path, symlink and hardlink aliases before opening output files. |
| Python consumers validated a complete arrival list before reparsing it. | Streaming arrival validation supports PSF analysis without retaining that list. |
| Direct comparison adapters bypassed some CSV invariants. | Both comparison paths share identity, optical scalar, direction, and loss checks. |
| Frozen comparison commands could modify their final input and still succeed. | Inputs are verified before and after each command; failures retain logs and a report. |
| Legacy timing and unavailable interactions could acquire production provenance. | Hash-bound normalization requires recorded timing and preserves missing interaction information. |
| Telescope and ray views repeated overlay drawing. | One renderer handles deterministic selection, interaction markers, status colours, and model provenance. |
| SVG images moved out-of-frame hits into detector edge bins. | Hits outside the displayed focal region are excluded from those bins. |
| Physical pixel entrances and dual-reflector segment masks were omitted. | Compiled finite pixel planes use an immutable stackless BVH; dual models retain explicit primary and secondary masks. |
| Camera and reflector offsets mixed coordinate origins. | Prescriptions, facet normals, focus offsets and physical camera planes share the optical reference frame. |
| Generic interfaces lacked material path and group timing. | One bounded transport records geometry, phase path, group delay, absorption, Fresnel transmission and TIR in native and bulk APIs. |
| Measured camera response and mirror scatter could be silently ignored. | Strict native schemas apply declared camera response once and keyed scatter independently of ordering and blocks. |
| Generated source batches were not reusable as reference input. | Resolved source export preserves photon context and every physical input at round-trip precision. |
| Docker cache entries accumulated in Actions storage. | Reusable layers now use release registry images and inline cache; 151 obsolete entries were removed. |
| Bulk API tests assumed a local build executable. | Installed wheel tests resolve the packaged native command and exercise native/bulk agreement. |

## Code retained deliberately

The analytic demonstrators remain independent test references. They do not
supply production telescope defaults. The generic nonsequential kernel and
slab primitive retain their analytic tests; their existence does not imply
that imported telescope windows or concentrators use them.

Model import, plotting, CSV diagnostics, and nanobind remain boundary code.
Transport uses the C++ standard library, immutable optical-model inputs, and
fixed interaction storage without per-photon allocation.

## Verification and limits

Focused tests exercise analytic intersections, angular response interpolation,
weight closure, model integrity, identity preservation, replay blocks, concurrent
bulk calls, frozen-input mutation, and plotted geometry. Full checks include
Debug and EventIO builds, the Python source suite, installed wheel smoke tests,
formatting, linting, and pre-commit hooks. Current measured results belong in
[STATUS.md](STATUS.md), rather than a second completion ledger here.

No matched production sim_telarray comparison or validated speed advantage is
established by these changes. See [PRODUCTION_VALIDATION.md](PRODUCTION_VALIDATION.md)
for the required common inputs, acceptance criteria, and evidence.

The [nominal five-model smoke record](review/nominal-optics.json) retains the
model hashes, commands and terminal counts for the corrected importer. Large
CSV outputs are reproducible from those commands and are not bundled here.
