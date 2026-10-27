# Scope and limitations

obdeect traces light through the optical model to the camera entrance. The
following effects are outside the current trace. Their impact depends on the
instrument and its optical-model parameters.

## Misalignment of aspheric (dual-mirror) mirrors

**What:** the primary and secondary mirrors of the Schwarzschild-Couder
telescopes (SST, SCT) are even polynomial surfaces (`mirror_class` 2) built
from segments. In the real telescope each segment, and the secondary as a
whole, can be tilted or shifted with respect to the design.

The alignment errors (`mirror_align_random_distance`, `_horizontal`,
`_vertical`) and `random_focal_length` act only on the facets of segmented
dishes (`mirror_class` 0). Dual-mirror surfaces are traced at their design
position; `primary_segmentation` and `secondary_segmentation` define gaps
between segments. The surface-scattering parameter
`mirror_reflection_random_angle` tilts the surface normal at each reflection.
See [validation](VALIDATION.md) for the resulting dual-mirror PSF comparisons.


## Curvature and refraction of the camera window (LST and others)

**What:** a camera window is a transparent plate (flat or curved) in front of
the pixels. Refraction in it shifts each photon sideways, moves the focus
backwards, and delays the photon. A curved window additionally acts as a weak
lens.

obdeect stops at the camera entrance (the focal surface). The window and its
transmission belong to the downstream camera simulation; window geometry and
refraction are not traced.

**Size of the effects** for a flat plate of thickness t and refractive index n
(PMMA n = 1.49 as an example; the thickness of each camera window has to be
taken from the camera documentation):

| t | lateral shift at 10 / 16 / 25 deg incidence | focus shift t(1 - 1/n) | extra delay (n_g - 1) t / c |
| --- | --- | --- | --- |
| 3 mm | 0.17 / 0.28 / 0.46 mm | 1.0 mm | 5 ps |
| 5 mm | 0.29 / 0.47 / 0.77 mm | 1.6 mm | 8 ps |
| 10 mm | 0.58 / 0.95 / 1.5 mm | 3.3 mm | 17 ps |

The lateral shift is t sin θ (1 - cos θ / sqrt(n² - sin² θ)); the delay uses a
group index n_g ≈ 1.5. The mean incidence angle on the LST focal plane is
16 deg ([validation](VALIDATION.md)). For comparison: the LST D80 is about
3 cm, its pixels about 5 cm wide, and a focus error δ blurs the image by about
δ D / f = δ / 1.2 for the LST. A shell of uniform thickness t and radius of
curvature R has an optical power of order (n - 1) t / R²: a focal length of
kilometres for millimetres and metres, negligible against 28 m.

The focus shift is the only one of these that is not negligible. It belongs
with the camera position: where `focus_offset` (6.55 cm for the LST) was tuned
to a measured PSF, the window is already in it, and adding the window to the
ray tracing would count it twice unless `focus_offset` were re-derived at the
same time. The origin of each telescope's `focus_offset` has to be audited
before any window refraction is added.


## Other effects not traced

- **Telescope structure of LST, SST and SCT.** Only the MST has a
  `telescope_structure` table in the models examined. The other telescopes use
  `telescope_transmission`, a parametrisation of shadowing.
- **Pixels, light guides, photodetectors:** handled by downstream camera simulation.
- **Refraction in air** between an illuminator and the telescope: a few arcsec
  over hundreds of metres.
- **Diffraction and polarisation:** irrelevant at the size of IACT mirrors and
  for unpolarised Cherenkov light.
