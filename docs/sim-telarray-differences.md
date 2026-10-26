# Differences from sim_telarray

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
