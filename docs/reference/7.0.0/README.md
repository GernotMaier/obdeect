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
D80, and effective area. Each image shows the focal-plane x/y distribution
and the centroid-centred 80% containment circle in centimetres.

| Telescope | North | South |
| --- | --- | --- |
| LST | [PSF](North_LSTN-design.png) | [PSF](South_LSTS-design.png) |
| MST FlashCam | [PSF](North_MSTx-FlashCam.png) | [PSF](South_MSTx-FlashCam.png) |
| MST NectarCam | [PSF](North_MSTx-NectarCam.png) | [PSF](South_MSTx-NectarCam.png) |
| SST | — | [PSF](South_SSTS-design.png) |

The imaging lists include focal-plane crossings that miss active pixels;
simtools computes the geometric PSF from all crossing positions. These
fixtures supply a production reference for implementing and checking the
obdeect scenes. They do not establish agreement between the two engines.
