# Not implemented, and why

obdeect replaces the ray tracing of sim_telarray. The effects below are
deliberately left out: sim_telarray does not trace them either, so leaving
them out keeps the two backends equivalent, and their size is small compared
with the uncertainties of the parameters that describe them. Each section says
what would change if they were added.

## Misalignment of aspheric (dual-mirror) mirrors

**What:** the primary and secondary mirrors of the Schwarzschild-Couder
telescopes (SST, SCT) are even polynomial surfaces (`mirror_class` 2) built
from segments. In the real telescope each segment, and the secondary as a
whole, can be tilted or shifted with respect to the design.

**sim_telarray:** the alignment errors (`mirror_align_random_distance`,
`_horizontal`, `_vertical`) and `random_focal_length` act only on the facets of
segmented dishes (`mirror_class` 0; `sim_config.c`). Dual-mirror surfaces are
traced at their design position; `primary_segmentation` and
`secondary_segmentation` only define the gaps between segments. The measured
PSF of the dual-mirror telescopes is reproduced by the surface-scattering
parameter `mirror_reflection_random_angle`, which tilts the surface normal at
every reflection.

**obdeect:** the same. The dual-mirror PSF agrees with sim_telarray within
0.15% (SST) and 0.10% (SCT) in D80 ([validation](VALIDATION.md)).

**If it were added:** segment tilts and secondary decentring need per-segment
frames and two random streams (segment, then secondary) keyed by the
configuration seed, as for facets. Without measured misalignment
distributions, the parameters would have to be tuned to the measured PSF,
which the random reflection angle already does. Moving PSF width from
scattering to misalignment changes the PSF shape (misalignment broadens
differently with off-axis angle), so it should come with a PSF measurement
that distinguishes the two. Not needed for sim_telarray equivalence.

## Curvature and refraction of the camera window (LST and others)

**What:** a camera window is a transparent plate (flat or curved) in front of
the pixels. Refraction in it shifts each photon sideways, moves the focus
backwards, and delays the photon. A curved window additionally acts as a weak
lens.

**sim_telarray:** the window is part of the camera simulation, as a
transmission only (`camera_filter`, over wavelength and, where given, incidence
angle; `camera_transmission`). It has no geometry and no refraction.

**obdeect:** stops at the camera entrance (the focal surface); everything
after it, the window included, belongs to the camera simulation. simtools hands
the photons to sim_telarray with the optics bypassed, where `camera_filter`
acts as before. Nothing is lost relative to sim_telarray.

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

**If it were added:** a window surface (plane or sphere) with thickness and an
index table on the leg after the last mirror, applied after the focal-surface
crossing has been redefined as the window entrance; the camera position of
every telescope audited at the same time. Not needed for sim_telarray
equivalence.

## Other effects not traced

- **Telescope structure of LST, SST and SCT.** Only the MST has a
  `telescope_structure` table (simulation-models branch
  `simulation-model-masts`). The other telescopes use `telescope_transmission`,
  sim_telarray's parametrised shadowing, in both backends. Tracing their masts
  needs their geometry from the instrument teams (the sources of the shadowing
  studies that the parametrisations were fitted to).
- **Pixels, light guides, photodetectors:** camera simulation (sim_telarray).
- **Refraction in air** between an illuminator and the telescope: a few arcsec
  over hundreds of metres ([physics](PHYSICS.md#known-differences-from-sim_telarray)).
- **Diffraction and polarisation:** irrelevant at the size of IACT mirrors and
  for unpolarised Cherenkov light.
