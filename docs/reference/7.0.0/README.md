# Production 7.0.0 optical reference

These figures are **simtools/sim_telarray reference output**, not obdeect
production traces. They use `validate-optics` in test mode: 5000 star photons,
10 km source distance, 20° zenith, zero field offset, and seed 19780503.
The Podman launcher adds `-C random_seed=19780503` because simtools otherwise
leaves that seed on `auto`. The generated configuration, stars file, compressed
imaging list, and sim_telarray log for each row remain in
`/private/tmp/obdeect-ref-fixed-7.0.0/<site>_<model>/` on this workspace.
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
