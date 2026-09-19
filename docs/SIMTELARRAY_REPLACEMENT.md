# sim_telarray replacement gate

`obdeect` is not yet a sim_telarray replacement. A full replacement claim is
allowed only when every item below is implemented and passes its reference
test. This gate prevents the current toy and analytic kernels from being used
as an unvalidated CTAO production simulator.

| Capability | Required implementation | Required evidence | Current state |
| --- | --- | --- | --- |
| Photon input | CORSIKA7 EventIO reader preserving position, direction, wavelength, time, weight and event/telescope IDs | Golden CORSIKA block vs independent decoder | Not implemented |
| CTAO model import | Pinned simulation-models resolver with file hashes and explicit unsupported-field failures | LST/MST/SST/SCT IR goldens | Not implemented |
| Segmented reflector core | Finite transformed M1/M2 facets, nearest hit, per-facet normals and wavelength-dependent coating | Unit rays for aperture, ordering, reflection and coating interpolation | Finite circular facet kernel only; no model segment-list importer |
| LST | Parabolic dish, segmented facets, camera and structure; wavelength response | sim_telarray PSF/area/time fixtures | Analytic paraboloid kernel only |
| MST | Modified Davies-Cotton facet placement, 86-facet list semantics and CSS shadowing | MST reference D80, area, time and shadow maps | Central-facet relation only |
| SST/SCT | M1/M2 segmentation, holes, aspheres, baffles, curved focal surface | SSTS/SCTS observables against sim_telarray/ROBAST | Surface polynomial kernel only |
| Materials | Reflectivity, windows, Fresnel/Snell, bulk absorption and sensor response | Table, energy-closure and window tests | 1-D interpolation only |
| Camera | Focal surface, pixels, gaps, concentrator and detection separation | Pixel/PSF/throughput fixtures | Screen endpoint only |
| Structure | Analytic rods/plates/boxes and optional mesh backend with attribution | Component-resolved shadow maps | Toy camera/masts only |
| Performance | Chunked SoA, deterministic counter RNG, threads/SIMD and memory bounds | Reproducible benchmarks versus sim_telarray | Scalar reference only |

The current development rule is: implement one row end-to-end, add its lowest
level tests, compare a shared input array to a reference, then expose it to the
next layer. No compatibility summary may omit a terminal loss category or
apply a scalar transmission on top of explicit structure.

## What “at least sim_telarray/ROBAST level” means here

This project may use that description only after the following acceptance run
is reproducible for each supported LST, MST, SST and SCT production model:

1. Import the pinned model release into an immutable scene IR, record every
   input file SHA-256, and reject unknown or ignored fields.
2. Trace the same fixed photon blocks with this code and a nominated reference
   (sim_telarray and/or ROBAST), covering on-axis and off-axis 300--700 nm
   stars, flasher and laser configurations.
3. Publish area, D80/PSF, focal-plane coordinates, optical-path-time residuals
   and component-resolved losses (facet gaps, M2, camera, masts, baffles and
   holes), with agreed tolerances per observable.
4. Run camera mapping and photoelectron detection, preserving the distinction
   between optical survival, camera-window transmission, concentrator loss and
   sensor conversion.
5. Keep each reference fixture and tolerance in CI. A test must fail when a
   model, coordinate convention or material table changes.

Until this evidence exists, `obdeect_ctao` is a **reference-prescription
tracer**, not a production simulator. In particular, it must not be used to
generate CTAO instrument-response functions or claim agreement with either
reference package.
