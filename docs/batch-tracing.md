# Batch tracing

## Reuse worker threads

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

## Shower photon bunches

Air-shower photons come as bunches: `size` photons (may be fractional) with one
position, direction, time, emission height and wavelength. `trace_bunches`
thins and traces them without making every photon a full trace:

```python
out = model.trace_bunches(
    {
        "x": x,
        "y": y,
        "cx": cx,
        "cy": cy,
        "time_ns": t,
        "zem": zem,
        "size": size,
        "wavelength_nm": wavelength,
        "bunch_id": bunch_id,
    },
    seed=run_seed_of_this_telescope,
    frame_rotation=rotation,  # 3 x 3, input frame to telescope frame
    frame_offset_m=offset_m,
    length_unit_m=0.01,  # metres per input length unit (cm here)
    plane_z=0.0,  # height of the input plane, input units
    acceptance=(wavelength_nm_axis, zem_m_axis, values),  # optional, in [0, 1]
    wavelength_density=(nm, density),  # optional
    pool=pool,  # or threads=...
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
