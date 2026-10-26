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
- `secondary_mirror_baffle` (z1, z2, R1, wall thickness, R2): open cylinder or
  cone wall with annular end faces on the incoming and after-primary legs, as
  `baffle_intersection` in `sim_imaging.c`; a thickness not below R1 is read
  as the outer radius, as in `sim_config.c`.
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

## Calibration light

A point emitter sends N photons into 4 pi with the relative intensity I(theta)
per solid angle, theta measured from its axis. The photons crossing a target
disk (unit normal n, radius R, distance h of the emitter from its plane) are
sampled without waste: target points are uniform on the disk, and a photon to
a point at distance r carries the weight

```latex
w = \frac{I(\theta)}{I_\mathrm{max}} \left(\frac{h}{r}\right)^3 \le 1,
```

its fluence relative to the largest possible. The number of candidates is
Poisson with mean N I_max R^2 / (2 h^2 \int I \, d(1 - \cos\theta)), so
keeping each with the probability of its weight gives a Poisson number of
photons crossing the disk, as for the real source. Angular profiles are linear
in 1 - cos(theta), as in the sim_telarray LightEmission tables; emission times
follow the pulse table. simtools hands the photons reaching the camera to
sim_telarray, which simulates the camera with its optics bypassed
(`Bypass_Optics=1` for single-mirror telescopes, where the hand-off turns the
telescope frame by pi about x; `Bypass_Optics=2` for dual-mirror telescopes,
in sim_telarray's camera frame, the telescope frame rotated by the camera
rotation), with the photons on the plane through the focal-surface vertex.

## Known differences from sim_telarray

- sim_telarray (`sim_imaging.c`, two places) rotates the scattered direction
  with `czn = -cyy*s3 + czz*c4`; it should use `s4`. obdeect uses the correct
  rotation; the difference is second order in the scattering angle.
- Telescope structure is traced only by obdeect.
- The camera front-face test on the dual-mirror leg is compiled out of
  sim_telarray; obdeect follows the compiled behaviour.
- Flashers: the camera is where the optics put it, at `focal_length +
  focus_offset` (simtools places `ff-1m` there as well). `ff-1m` traces no
  structure.
- Illuminators: obdeect sends the light on straight lines; refraction over the
  few hundred metres to a telescope bends it by a few arcsec. `xyzls` applies
  CORSIKA's refraction and delay model, made for light from above; for nearly
  horizontal light (elevation below about 5 deg) its photon times are off by
  ns to tens of ns and its image can move by about 1 arcmin.
- Reflectivity tables over wavelength and incidence angle (2-D; the SST
  mirrors in simulation models 7.0.0): sim_telarray reads the reflectivity
  before ray tracing as R(lambda, 0 deg) (`sim_absorb.c`, `rpolate(.., 0., -1)`
  on a 2-D table interpolates at angle zero; the comment intends the upper
  envelope) and, after each reflection, multiplies by R(lambda, theta) /
  R_max(lambda) with the true envelope (`rpolate_1d(.., -1)`). Its
  reflectivity is therefore R(lambda, 0) R(lambda, theta) / R_max(lambda)
  instead of R(lambda, theta): for the SST at 400 nm 0.966 (primary) times
  0.945 (secondary), at 300 nm 0.944 times 0.607. obdeect applies
  R(lambda, theta).
