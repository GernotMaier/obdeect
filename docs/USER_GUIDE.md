# obdeect user guide

obdeect traces photons through the optics of an imaging atmospheric Cherenkov
telescope: from the light source, over the primary (and secondary) mirror, to
the camera entrance (the focal surface). It is a C++20 kernel with a Python
package and the command `obdeect`.

```sh
python -m pip install obdeect-dev
```

## Command line

```sh
# Compile an optical model from a simulation-models checkout.
obdeect compile --simulation-models ../simulation-models --model LSTN-01 \
  --version 7.0.0 --seed 1 --zenith-deg 20 --output lst.json

# Trace a star 0.5 deg off axis; print the image summary, keep every photon.
obdeect trace --model lst.json --seed 2 --photons 100000 \
  --field-deg 0.5 0 --output lst-arrivals.npz

# Draw the dish seen from the sky and the focal-plane image.
obdeect plot --model lst.json --arrivals lst-arrivals.npz --output lst.png
```

`compile` reads full and patch productions of the simulation models; a patch
production is laid over the productions in its `model_version_history`, as in
simtools. `--seed fixes the random mirror errors (focal
lengths, alignment); `--zenith-deg` is the telescope zenith angle used by
zenith-dependent focus offsets and alignment errors. `trace --seed` fixes source
sampling and mirror scattering. `--geometry-only` leaves the photon weights
unchanged (no reflectivity, telescope transmission or degradation). Results do
not depend on `--threads`, nor on photon order or block size.

## Python

```python
from obdeect import DETECTED, OpticalModel

model = OpticalModel.from_simulation_models(
    "../simulation-models", "LSTN-01", "7.0.0", configuration_seed=1, zenith_deg=20.0
)
photons = model.star_photons(100_000, seed=2, towards_source=(0.0087, 0.0, 0.99996))
arrivals = model.trace(photons, seed=2)
image = arrivals["position_m"][arrivals["status"] == DETECTED]
```

`OpticalModel.from_records(records, read_file, ...)` compiles parameter records
held in memory (simtools uses it); `save`/`load` store the optical model as JSON.

`trace` takes photons as arrays: `position_m` (n, 3), `direction` (n, 3, unit
vectors), `wavelength_nm`, `time_ns`, `weight` and `photon_id` (uint64; the
random numbers of a photon depend only on the seed and its id). It returns, per
photon:

| Array | Meaning |
| --- | --- |
| `status` | index into `STATUS_NAMES`: detected, missed_primary, missed_secondary, missed_focal_surface, blocked, invalid_input, numerical_failure |
| `terminal_surface` | id of the surface the photon ended on (focal surface, obscurer, mirror) |
| `position_m`, `direction` | at the focal surface, or where the photon was lost |
| `path_length_m`, `time_ns` | geometric path; time with the model's group index |
| `weight` | input weight times the applied responses; zero when lost |
| `primary_point_m`, `primary_surface`, `primary_incidence_deg` | reflection on the primary |
| `secondary_point_m`, `secondary_surface`, `secondary_incidence_deg` | reflection on the secondary (NaN for single-mirror telescopes) |
| `focal_incidence_deg` | angle to the focal-surface normal |

By default `trace` returns every array above for every photon. Writing 14
arrays costs memory bandwidth, so ask for what you use:

```python
detected = model.trace(
    photons, seed=2,
    outputs=("position_m", "direction", "time_ns", "weight", "index"),
    detected_only=True,
)
```

`outputs` lists the arrays to compute and return (an unknown name raises
`ValueError` and lists the valid ones). With `detected_only=True` only photons
with status `detected` are returned, in input order; the additional array
`index` (valid only with `detected_only`) gives the row of each in the input.
The values are those of the full trace, for any `threads`.

`classify_losses=False` (default `True`) reports a photon that misses the
primary mirror as `missed_primary`, even where an obscurer on the incoming leg
would have blocked it first; the search for such obscurers is then skipped.
Detected photons are the same either way.

### Many small traces: `ThreadPool` and `trace_many`

`trace(threads=...)` starts its threads on every call, which is costly when an
event is traced as many small blocks. A `ThreadPool` keeps its threads between
calls: `ThreadPool(n)` (0 uses every core) runs one job at a time, and
`model.trace(photons, seed=2, pool=pool)` uses it instead of `threads`.
`trace_many` traces the photons of several optical models, as one job, in a
single call:

```python
from obdeect import ThreadPool, trace_many

pool = ThreadPool(0)
results = trace_many(
    [
        (model_a, photons_a, {"outputs": ("position_m", "weight", "index"), "detected_only": True}),
        (model_b, photons_b, {"apply_responses": False}),
    ],
    seed=2,
    pool=pool,
)
```

Each request is `(model, photons, options)` with the options `apply_responses`,
`direct`, `outputs`, `detected_only` and `classify_losses` of `trace`; the result is the list of
the dicts that one `trace` call per request returns, with the same values for
any pool size. A call must not use the same pool from inside another call that
runs on it.

## Shower photon bunches: `trace_bunches`

Air-shower photons come as bunches: `size` photons (may be fractional) with one
position, direction, time, emission height and wavelength. `trace_bunches`
thins and traces them without making every photon a full trace:

```python
out = model.trace_bunches(
    {"x": x, "y": y, "cx": cx, "cy": cy, "time_ns": t, "zem": zem,
     "size": size, "wavelength_nm": wavelength, "bunch_id": bunch_id},
    seed=run_seed_of_this_telescope,
    frame_rotation=rotation,          # 3 x 3, input frame to telescope frame
    frame_offset_m=offset_m,
    length_unit_m=0.01,               # metres per input length unit (cm here)
    plane_z=0.0,                      # height of the input plane, input units
    acceptance=(wavelength_nm_axis, zem_m_axis, values),   # optional, in [0, 1]
    wavelength_density=(nm, density),                      # optional
    pool=pool,                        # or threads=...
)
```

The inputs are float32 arrays (`bunch_id` is uint64). Photons move towards -z of
the input frame: `cz = -sqrt(1 - cx^2 - cy^2)`. Position is
`R (u x, u y, u plane_z) + offset` and direction `R (cx, cy, cz)`, with
`u = length_unit_m`; `R` must be orthonormal. A bunch of `wavelength_nm = 0`
takes a wavelength per photon from `wavelength_density` (linear between the
nodes, over nm).

For each photon of a bunch (`ceil(size)` of them; the last one counts with the
fraction `size - floor(size)`) obdeect

1. looks up the external `acceptance` (atmosphere, quantum efficiency, ...; a
   table you supply, with wavelength in nm and emission height in metres; none
   means 1) and the upper bound of the model's responses at that wavelength
   (`response_bound`: reflectivities, degradation maps, telescope
   transmission);
2. keeps the photon with probability `fraction * acceptance * bound`, so that
   most lost photons cost no ray tracing;
3. traces the kept photon with the responses, and with `roulette=True`
   (default) keeps it again with probability `weight / bound`. Survivors have
   weight 1 and their expected number equals the sum of weights of a full
   trace. With `roulette=False` all of them stay and carry `weight / bound`.

The result is a dict of columns in input order, bunch by bunch and photon by
photon: `position_m`, `direction`, `time_ns`, `weight`, `wavelength_nm`,
`bunch_index` (row of the bunch in the input) and `sub_photon` (number within
the bunch). It does not depend on `threads`, on `pool`, or on how the bunches
are split over calls. Invalid input (a frame that is not orthonormal, an
acceptance outside [0, 1], a bunch of wavelength 0 without a density, a bad
size or direction) raises `ValueError` before any tracing.

The random numbers of a photon depend on `(seed, bunch_id, sub_photon)`.
**Bunch ids must be unique per seed.** The recommended id is
`(event_number << 32) | bunch_index_in_event`, and a different `seed` per
telescope, derived by the caller from the run seed and the telescope number
(otherwise two telescopes would draw the same numbers for equal ids).
`classify_losses=False` is as for `trace`.

## Calibration light

`emitter_photons` makes the light of a point emitter, such as a flasher or an
illuminator, that crosses a target disk. The emitter has an angular profile
(relative intensity per solid angle against the angle to its axis) and a pulse
shape (relative rate against time), both tables:

```python
camera = model.focal_surface_z_m
photons = model.emitter_photons(
    1e6,                                   # photons emitted into 4 pi
    seed=3,
    position_m=(0.2, 0.2, camera - 16.5),   # telescope frame
    intensity=([0.0, 11.0, 11.0001], [1.0, 1.0, 0.0]),  # 11 deg filled cone
    pulse=([-1.25, 1.25], [1.0, 1.0]),     # 2.5 ns flat pulse
    wavelength_nm=385.0,
    disk=((0.0, 0.0, camera), (0.0, 0.0, 1.0), 1.2),  # the camera
)
arrivals = model.trace(photons, seed=3, direct=True)
```

The number of photons is Poisson distributed and their weights are in [0, 1]:
keeping each photon with the probability of its weight gives the photons that
cross the disk (see [PHYSICS.md](PHYSICS.md)). Without `disk`, the target is the
star-light disk in front of the telescope. `direct=True` traces light that
starts inside the telescope straight to the focal surface: no mirror and no
telescope transmission act, the obscurers of that leg do.

A bright source (an illuminator emits about 10^10 photons per pulse) gives
more photons than fit in memory at once. `trace_emitter` takes the same
arguments as `emitter_photons` plus those of `trace`, makes and traces the
photons `block_size` at a time, and keeps only those reaching the focal surface
with a positive weight (`position_m`, `direction`, `time_ns`, `weight`,
`wavelength_nm`). `trace_emitter_blocks` yields them block by block, so that
memory stays bounded by one block. The photons do not depend on `block_size`:

```python
arrivals = model.trace_emitter(1e10, seed=3, position_m=(0.0, 0.0, 600.0),
                               intensity=([0.0, 180.0], [1.0, 1.0]), axis=(0.0, 0.0, -1.0))
```

## Frames and units

All lengths are metres, angles in degrees unless named `_rad`. The telescope
frame is sim_telarray's: +z along the optical axis towards the sky, the dish
vertex at z = -`mirror_offset`, and x, y as in the mirror lists and the camera
photon lists of sim_telarray (with the telescope pointing above the horizon, x
points down). `star_photons` samples a disk of radius 1.2 times the dish radius
around the optical axis, perpendicular to it and centred on the dish vertex, as
sim_telarray does for star light; `sampling_area_m2` is its area. The distance
of a point source is measured from the disk centre.

## Optical model format

An optical model (`obdeect.optical-model.v2`) is one JSON document:

- `format`, `sha256` (of the document without it), `provenance` (parameter
  versions, file hashes, seed, zenith angle, structure flag), `notes`;
- `sampling`: the star-light disk;
- `optics`: `primary`, optional `secondary` (each either `facets` or an
  `asphere` with `segments`, a `reflectivity` table, `scatter`, optional
  `degradation` map), `focal_surface` (plane or asphere), `obscurers`
  (cylinder, polyhedron, plate, frustum; each acting on the legs `incoming`,
  `after_primary`, `after_secondary`), `telescope_transmission`,
  `propagation_group_index`.

The C++ loader accepts only known keys, checks every value and names the JSON
path of the first error. Surface and obscurer ids are unique.

## Telescope structure

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

## Build options

- `-DOBDEECT_USE_FMA=ON` (CMake, default OFF) evaluates the aspheric polynomials
  with `std::fma`. Use it only on hardware with fused multiply-add: without it
  `std::fma` is emulated and very slow. Results differ from builds without it
  at the rounding level. Default builds, wheels and CI do not use it.
- Results also remain platform-dependent in the last bits, because the math
  library (`std::pow`, `std::hypot`, trigonometric functions) differs between
  platforms; see `VALIDATION.md`, section "Rounding-level changes".
