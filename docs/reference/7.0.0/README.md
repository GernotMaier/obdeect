# Production 7.0.0 optical reference

These figures are **simtools/sim_telarray reference output**, not obdeect
production traces. simtools `validate-optics` generated the production 7.0.0
configs and star lists in test mode. The reference runs replay those configs
through sim_telarray with **50,000 photons** per telescope (10× the initial
test run), a 10 km source, 20° zenith, zero field offset, and seed 19780503.
The Podman launcher adds `-C random_seed=19780503` because simtools otherwise
leaves that seed on `auto`. The generated configs, star lists, imaging lists,
and sim_telarray logs remain in
`/private/tmp/obdeect-ref-50k-7.0.0/<site>_<model>/` on this workspace.
Direct replay was necessary because the updated local simtools checkout cannot
resolve this model checkout's legacy file paths; all seven configs retained
their original SHA-256 hashes. The sim_telarray logs confirm that it consumed
the 7.0.0 LST/MST/SST configs, selected site altitude, 50,000 photons, and
the seed. The configs use mirror classes 0/0/2 and focal lengths 28/16/2.15 m
for LST/MST/SST respectively.
The `validate-optics` settings disable the camera filter and night-sky
background and set camera transmission to one. The reported area is therefore
an optical ray-tracing area for this PSF fixture, not a full camera-throughput
prediction.
The [summary](summary.json) records revisions, hashes, focal-plane crossings,
effective area, and exact centroid-centred integration radii at 50%, 68%, 80%,
90%, 95%, and 99% containment. Each linked CSV gives the empirical cumulative
distribution every 0.1% in centimetres, including the enclosed photon count.
The five distinct compressed [imaging lists](imaging/) are archived here; the
North and South MST rows share identical photon lists. Recalculate and verify
all profiles and hashes with
`PYTHONPATH=python python -m obdeect.simtel_reference_psf docs/reference/7.0.0 --check`.
The summary keeps simtools' iterative D80 and also reports the exact empirical
`exact_d80_m`. The figures draw simtools' D80 circle.

| Telescope | North PSF / cumulative profile | South PSF / cumulative profile |
| --- | --- | --- |
| LST | [PSF](North_LSTN-design.png) / [CSV](profiles/North_LSTN-design.csv) | [PSF](South_LSTS-design.png) / [CSV](profiles/South_LSTS-design.csv) |
| MST FlashCam | [PSF](North_MSTx-FlashCam.png) / [CSV](profiles/North_MSTx-FlashCam.csv) | [PSF](South_MSTx-FlashCam.png) / [CSV](profiles/North_MSTx-FlashCam.csv) |
| MST NectarCam | [PSF](North_MSTx-NectarCam.png) / [CSV](profiles/North_MSTx-NectarCam.csv) | [PSF](South_MSTx-NectarCam.png) / [CSV](profiles/North_MSTx-NectarCam.csv) |
| SST | — | [PSF](South_SSTS-design.png) / [CSV](profiles/South_SSTS-design.csv) |

The imaging lists include focal-plane crossings that miss active pixels;
simtools computes the geometric PSF from all crossing positions. These
fixtures supply a production reference for implementing and checking the
obdeect scenes. They do not establish agreement between the two engines.
