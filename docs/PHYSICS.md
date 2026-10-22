# Physics of the photon path

obdeect reproduces the optical part of sim_telarray: the same parameters give
the same geometry, random mirror errors with the same distributions, and the
same photon fate. This page lists what is computed, where it follows
sim_telarray (`sim_imaging.c`, `sim_config.c`), and where it deliberately
differs. The analytic and numerical checks are in `cpp/tests` and
`python/tests`; the comparisons with sim_telarray in [VALIDATION.md](VALIDATION.md).

## Compilation of the optical model

**Segmented dishes** (`mirror_class` 0): mirror centres from `mirror_list`, on
the dish shape (`parabolic_dish`, `dish_shape_length`) unless a z is given. Each
mirror is a sphere; its focal length is the catalogue value or, if zero, the
automatic prescription (distance to the focus or `mirror_focal_length`), with
`random_focal_length` errors (Gaussian clipped as in sim_telarray, or uniform).
Alignment follows `mirror_align_random_distance`, `_horizontal`, `_vertical`
with the zenith-dependent widths
`sqrt(r0² + d0²(cos θ - cos θ0)² + e0²(sin θ - sin θ0)²)` at the telescope
zenith θ. Mirror apertures: circle, hexagon (both orientations), square.

**Dual-mirror telescopes** (`mirror_class` 2): even polynomial surfaces for
primary, secondary and focal surface (`*_parameters`, optional `*_ref_radius`),
the secondary in sim_telarray's frame rotated by π about y, apertures from
diameters and holes, segmentation (`*_segmentation`) as gaps between segments.
sim_telarray does not misalign dual-mirror segments; neither does obdeect.

**Camera and focus.** Single-mirror telescopes: focal plane at
`focal_length + focus_offset - mirror_offset`. The zenith terms of
`focus_offset` apply only if both are non-zero, as in `sim_config.c`. Dual-mirror
telescopes: the focal surface polynomial, unbounded in radius.

**Shadows.**

- Camera body (`camera_body_diameter`, `_shape`, `camera_depth`): incoming
  photons are tested on the camera front and back planes (single mirror); the
  rays between primary and secondary are tested on the camera back plate and an
  open side wall (dual mirror).
- `secondary_mirror_shadow_diameter`/`_offset`: plate tested by incoming light,
  or the secondary itself if no diameter is given.
- `secondary_mirror_baffle`: open cone on the incoming and after-primary legs.
- Telescope structure (beyond sim_telarray, see the user guide): cylinders and
  convex polyhedra on every leg. With structure, `telescope_transmission` is
  not applied.

**Responses.** `mirror_reflectivity` (wavelength and optionally incidence
angle, resampled on a 1 nm grid), `mirror_degraded_reflection`, the secondary's
own tables, `telescope_transmission` (sim_telarray's angular parametrisation,
method and outer power), spatial degradation maps.

## Tracing a photon

1. The source supplies position, direction, wavelength, time, weight and id.
   Star light (`star_photons`) is sampled on sim_telarray's disk (radius 1.2
   times the dish radius, perpendicular to the optical axis); a source at finite
   distance emits from its position.
2. Incoming leg: nearest of primary mirror, secondary back side (if it blocks
   incoming light) and obscurers. Lost photons record where and on what.
3. Primary: nearest valid hit (facet grid or polynomial with Newton iterations
   from a paraboloid guess, with a bracketed fallback); specular reflection.
   Scattering by `mirror_reflection_random_angle` (two Gaussian components): for
   facets the reflected direction is rotated, for polynomial mirrors the normal
   is tilted by half the angle, as in sim_telarray.
4. Secondary (dual mirror): same, with the after-primary obscurers.
5. Focal surface: plane or polynomial; the photon's position, direction, time
   and incidence angles are recorded. obdeect stops at the camera entrance;
   pixels, light guides and the camera window are not part of the optical model.

Weights multiply the responses met (unless geometry only). Each random number
comes from a Philox-4x32-10 counter keyed by the seed, photon id and stream, so
results do not depend on threads, blocks or photon order. Times use the group
index stored in the model.

## Known differences from sim_telarray

- sim_telarray (`sim_imaging.c`, two places) rotates the scattered direction
  with `czn = -cyy*s3 + czz*c4`; it should use `s4`. obdeect uses the correct
  rotation; the difference is second order in the scattering angle.
- Telescope structure is traced only by obdeect.
- The camera front-face test on the dual-mirror leg is compiled out of
  sim_telarray; obdeect follows the compiled behaviour.
