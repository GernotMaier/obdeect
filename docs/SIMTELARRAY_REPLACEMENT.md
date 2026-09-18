# sim_telarray replacement gate

`obdeect` is not yet a sim_telarray replacement. A full replacement claim is
allowed only when every item below is implemented and passes its reference
test. This gate prevents the current toy and analytic kernels from being used
as an unvalidated CTAO production simulator.

| Capability | Required implementation | Required evidence | Current state |
| --- | --- | --- | --- |
| Photon input | CORSIKA7 EventIO reader preserving position, direction, wavelength, time, weight and event/telescope IDs | Golden CORSIKA block vs independent decoder | Not implemented |
| CTAO model import | Pinned simulation-models resolver with file hashes and explicit unsupported-field failures | LST/MST/SST/SCT IR goldens | Not implemented |
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
