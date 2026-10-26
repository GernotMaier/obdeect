# Optical elements

An optical model describes the surfaces a photon can meet, their responses and
the focal surface where tracing ends. The elements below are grouped by their
effect on the light path.

## Reflecting elements

Mirrors redirect photons by specular reflection. Reflectivity and degradation
change their weights; scattering perturbs the reflected direction.

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

## Refracting elements

Refracting elements, such as lenses and camera windows, bend light at a material
boundary. obdeect currently does not trace them. The optical path ends at the
camera entrance; camera-window transmission belongs to the downstream camera
simulation. Refraction in air is also omitted.

The model's `propagation_group_index` sets photon flight times; it does not bend
rays. See [Scope and limitations](limitations.md) for the size of the omitted
effects and the assumptions behind them.

## Obscuring elements

Obscurers block light on specified legs of its path: `incoming`, `after_primary`
and `after_secondary`. They include camera bodies, secondary shadows, baffles
and telescope structure. A blocked photon records the obscurer's id as its
`terminal_surface`.

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

### Telescope structure

If a telescope model contains the parameter `telescope_structure`, the shadows
of its structure (masts, quadrupods, dish frame) are traced on all legs, and
`telescope_transmission` is not applied: the geometry replaces the
parameterised shadowing (as in sim_telarray's own mast option,
`RAYTRACING_INTERSECT_RODS`). The camera body shadow stays. The table is ECSV
in metres in the `primary_mirror_frame` (the frame of `mirror_list`, z = 0 at
the dish vertex), one row per element, by its `shape`:

- `cylinder`: axis end points `x1 y1 z1`, `x2 y2 z2` and `diameter`;
- `convex_polyhedron`: corner points `x1 y1 z1`, `x2 y2 z2`, ... (at least
  four); the element is the solid convex hull of its corners.

Columns `id`, `group` (optional) and `shape` are text. A row leaves the
columns it does not use as `nan`; unused points come last.

Which elements actually shade the light is shown by
`OpticalModel.obscurer_hits` and `obdeect obscurers` (see
[shadowing](shadowing.md)): per obscurer, the light whose path meets it, that
it absorbs, and that only it stops.

```{toctree}
:maxdepth: 1

shadowing
```

## Focal surface

The focal surface records arriving photons and ends the optical trace.

Single-mirror telescopes: focal plane at
`focal_length + focus_offset - mirror_offset`. The zenith terms of
`focus_offset` apply only if both are non-zero, as in `sim_config.c`. Dual-mirror
telescopes: the focal surface polynomial, unbounded in radius.

## Optical responses

`mirror_reflectivity` (wavelength and optionally incidence
angle, resampled on a 1 nm grid), `mirror_degraded_reflection`, the secondary's
own tables, `telescope_transmission` (sim_telarray's angular parametrisation,
method and outer power), spatial degradation maps.
