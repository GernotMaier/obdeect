# Optics

obdeect traces photons from an entrance surface to the camera entrance.
The optical model describes mirrors, obscurers, the focal surface and optical
responses. The CTAO importer follows the sim_telarray conventions described below;
comparisons and their limits are recorded in [Validation](VALIDATION.md).

```{toctree}
:maxdepth: 1

optical-elements
telescopes
shadowing
sim-telarray-differences
limitations
```

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
