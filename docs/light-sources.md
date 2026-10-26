# Light sources

Choose a source by the light path you want to study. obdeect provides star and
point-source sampling, calibration emitters, and an interface for externally
generated shower photons.

| Source | Typical use | Interface |
| --- | --- | --- |
| Star at infinity | Point-spread function, effective area, off-axis response | `star_photons` |
| Point source at finite distance | Focus and distance-dependent imaging | `star_photons(distance_m=...)` |
| Flasher | Direct illumination of the camera | `emitter_photons`, then `trace(direct=True)` |
| Illuminator | External calibration light passing through the mirrors | `emitter_photons`, then `trace` |
| Air-shower photon bunches | Cherenkov light from an external shower simulation | `trace_bunches` |

## Stars and finite-distance point sources

`star_photons` samples the optical model's entrance disk. At infinity, rays
are parallel; at a finite distance, their directions follow the point source.
`towards_source` points towards the source in the telescope frame, whose +z
axis points towards the sky. `distance_m` is measured from the dish vertex.

```python
# On-axis star at infinity (the default).
photons = model.star_photons(100_000, seed=2, wavelength_nm=400.0)

# On-axis point source 10 km from the dish vertex.
photons = model.star_photons(100_000, seed=2, distance_m=10_000.0)
arrivals = model.trace(photons, seed=2)
```

Use `towards_source` for an off-axis source. On the command line, use
`obdeect trace --field-deg X Y` and optionally `--distance-m`.
This sampler fixes the number of incident photons; it does not describe an
emitter's luminosity or pulse shape.

## Flashers and illuminators

Both use a point emitter with a position, wavelength, angular intensity profile
and optional pulse shape. `emitted` specifies the mean number of photons emitted
into all directions. Sampling targets a disk and uses photon weights to represent
the illumination; it does not generate every emitted photon.

A **flasher** illuminates the camera directly from inside the telescope. Set
the target disk at the camera and trace with `direct=True`: obscurers are tested,
but mirrors and telescope transmission are bypassed.

An **illuminator** shines onto the telescope from outside. Target the entrance
disk (the default) and use the normal mirror-tracing path. For bright pulses,
`trace_emitter` or `trace_emitter_blocks` generates and traces photons in blocks.

For example, a flasher with an 11-degree cone and a 2.5 ns flat pulse:

```python
camera = model.focal_surface_z_m
photons = model.emitter_photons(
    1e6,  # photons emitted into 4 pi
    seed=3,
    position_m=(0.2, 0.2, camera - 16.5),  # telescope frame
    intensity=([0.0, 11.0, 11.0001], [1.0, 1.0, 0.0]),  # 11 deg filled cone
    pulse=([-1.25, 1.25], [1.0, 1.0]),  # 2.5 ns flat pulse
    wavelength_nm=385.0,
    disk=((0.0, 0.0, camera), (0.0, 0.0, 1.0), 1.2),  # the camera
)
arrivals = model.trace(photons, seed=3, direct=True)
```

The candidate count is Poisson distributed and weights lie in [0, 1]. Keeping
each candidate with the probability of its weight gives the photons crossing
the target disk. The sum of weights estimates their mean number.

## Air-shower photons and custom inputs

`trace_bunches` accepts photon bunches produced by an external air-shower
simulation. Supply their positions, directions, times, wavelengths, sizes and
unique bunch ids, together with an explicit frame transformation and length
unit. obdeect samples and traces these photons; it does not generate the shower.
See [Shower photon bunches](batch-tracing.md#shower-photon-bunches).

For other sources, supply photon arrays directly to `trace`. The required fields
and output arrays are described in [Python usage](python-guide.md).
