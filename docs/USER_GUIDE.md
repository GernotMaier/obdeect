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
