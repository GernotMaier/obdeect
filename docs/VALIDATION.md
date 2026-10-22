# Validation against sim_telarray

obdeect is validated in three layers:

1. **Kernels and primitives** (`cpp/tests`): random numbers (Philox known
   answers), tables, apertures, facet grid against brute force, polynomial
   surfaces against brute force, obscurers, sources, transmission; traced
   photons in analytic optics; model loading and its error messages.
2. **Compilation** (`python/tests`): parameter records to optical models —
   Davies-Cotton geometry, seeded mirror errors, focus offsets, camera body,
   dual-mirror surfaces, telescope structure, unit handling, store and reload.
3. **Equivalence with sim_telarray** (below): the same simtools applications
   run with both backends on the same simulation models.

## How to repeat the comparison

In a simtools environment with sim_telarray and obdeect installed, for each
telescope and backend:

```sh
simtools-validate-optics --site North --telescope LSTN-01 --model_version 6.0.2 \
  --zenith_angle 20 --source_distance 10 --max_offset 4 --offset_step 1 \
  --offset_directions N,E --ray_tracing_backend obdeect --obdeect_seed 11 \
  --output_path out/LSTN-01/obdeect
```

and the same with `--ray_tracing_backend sim_telarray`, then compare the
`results/ray_tracing_*.ecsv` tables (D80, effective area, effective focal
length per offset). Incident angles: `simtools-derive-incident-angle` with both
backends. Both backends draw 20 mirror-error realisations of 100,000 photons per
offset (2 million photons), so random mirror errors and photon statistics
limit the agreement; repeating with other seeds shows the spread.

## Results (simulation models 6.0.2, zenith 20 deg, star at 10 km, 0-4 deg)

Largest deviation of obdeect from sim_telarray over nine offsets:

| Telescope | D80 | Effective area | Effective focal length |
| --- | --- | --- | --- |
| LSTN-01 | 0.3% | 0.25% | 0.02% |
| MSTN-01 (NectarCam) | 0.9% | 0.25% | 0.1% (realisation spread) |
| MSTS-01 (FlashCam) | 0.5% | 0.3% | 0.06% (realisation spread) |
| SSTS-01 | 0.15% | 0.3% | < 0.01% |
| SCTS-01 (one-pixel camera) | 0.15% | 0.13% | < 0.01% |

- Star at 10,000 km (LSTN-01, 0-2 deg): D80 within 0.15%, area within 0.13%.
- Per mirror (MSTN-01, on axis, 86 mirrors): the fractions of photons
  reflected by each mirror agree within statistics (pulls: mean 0.05,
  standard deviation 0.96), including the mirrors partly shadowed by the camera.
- Single-mirror mode (29 MSTN-01 mirrors): mean D80 1.503 cm for both, spread
  0.183 cm (sim_telarray) and 0.160 cm (obdeect).
- Incident angles (SSTS-01, 1 deg, 10^6 photons, per photon): mean angles on
  the focal surface 46.656/46.658 deg, primary 9.899/9.898 deg, secondary
  33.184/33.183 deg; in bins of the primary radius they agree within 0.05 deg.
  LSTN-01 focal surface: 15.81/15.80 deg.
- Frames: image centroids of off-axis stars agree in sign and size (LSTN-01,
  2 deg: -102.42 and -102.40 cm) without any frame conversion.

The PSF applications, LSTN-01, measured PSF curves of the simtools test
resources:

| Application | sim_telarray | obdeect |
| --- | --- | --- |
| `validate-cumulative-psf` (6.0.2, 20 deg) | D80 3.051 cm, area 365.8 m², 31 s | D80 3.059 cm, area 365.4 m², 3 s |
| `derive-psf-parameters` (6.0.2, 20 deg) | RMSD 0.0131, D80 3.053 cm | RMSD 0.0127, D80 3.059 cm |
| `derive-mirror-rnda` (7.0.0, test mode, 5 mirrors) | fitted D80 within 0.9-3.8% of measured, 117 s | within 0.8-4.9%, 19 s |

Both backends start `derive-psf-parameters` below its RMSD threshold, so the
optimisation stops at once. Its final D80-versus-offset plot (46 offsets of 2
million photons) needs more than 2 GB of memory with either backend, because
simtools keeps every photon image; in `--test` mode the whole application
takes 11 s (sim_telarray) and 6 s (obdeect). With 5000 photons per evaluation the
`derive-mirror-rnda` fits differ by their noise; with a fixed `--obdeect_seed`
every evaluation reuses the same photons, so the obdeect objective is smooth.

The comparison found and fixed one convention: sim_telarray measures the
distance of a point source from the dish vertex, not from the telescope frame
origin (0.02% plate-scale difference for the SCT at 10 km).

## Telescope structure (beyond sim_telarray)

With the MSTx-FlashCam structure prototype (simulation-models branch
`simulation-model-masts`), the traced structure transmits 0.895 of the light on
axis and 0.882 at 4 deg, against 0.898 and 0.885 from `telescope_transmission`.
This needs one prototype error fixed: `dish-structure-prism-046` has its fourth
vertex at x = -5.48 m instead of +5.48 m, which spans a slab across the dish and
lowers the transmission to 0.856.
