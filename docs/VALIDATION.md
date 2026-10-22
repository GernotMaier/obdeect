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
| LSTS-01 | 0.4% | 0.11% | 0.01% |
| LSTN-design | 0.5% | 0.15% | 0.01% |
| LSTS-design | 0.4% | 0.10% | 0.01% |
| MSTx-FlashCam | 1.1% | 0.24% | 0.1% (realisation spread) |
| MSTx-NectarCam | 0.5% | 0.23% | 0.04% |
| SSTS-design | 0.12% | 0.45% | 0.01% |
| SCTS-design (one-pixel camera) | 0.10% | 0.22% | < 0.01% |

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
| `derive-mirror-rnda` (7.0.0, full statistics) | first 20 mirrors: within 4.2% (median 1.7%), 848 s | all 198 mirrors: within 5% (median 2.4%), 241 s |

Both backends start `derive-psf-parameters` below its RMSD threshold, so the
optimisation stops at once. Its final D80-versus-offset plot (46 offsets in four
directions, 2 million photons each) needed more than 2 GB of memory with either
backend, because simtools kept every photon image (fixed in simtools#2664): at
full statistics it takes 2120 s with sim_telarray (RMSD 0.0132 to 0.0101 in
one step, peak 0.3 GB) and 115 s with obdeect (0.0127 to 0.0107, peak 1.45
GB); in `--test` mode the whole application takes 11 s
(sim_telarray) and 6 s (obdeect). With 5000 photons per evaluation the
`derive-mirror-rnda` fits differ by their noise; with a fixed `--obdeect_seed`
every evaluation reuses the same photons, so the obdeect objective is smooth.
On the 20 mirrors fitted by both, the fitted random reflection angle of obdeect
is 2.8% below sim_telarray's on average, with a spread of 6.6% between mirrors
(fit noise). A full sim_telarray fit needs `--max_workers 1`: parallel mirror
jobs share one model directory and rewrite its files while others read them.

The comparison found and fixed one convention: sim_telarray measures the
distance of a point source from the dish vertex, not from the telescope frame
origin (0.02% plate-scale difference for the SCT at 10 km).

## Telescope structure (beyond sim_telarray)

`telescope_transmission` is sim_telarray's fit of the ratio of the effective
mirror area of a detailed ray tracing (with masts) to its own (camera body
only), versus the off-axis angle. The same ratio from obdeect, the same photons
traced with and without the structure (geometry only, mean of four azimuths),
for MSTx-FlashCam (`telescope_structure` 1.0.0 on simulation-models branch
`simulation-model-masts`):

| Off axis | 0 deg | 1 deg | 2 deg | 3 deg | 4 deg |
| --- | --- | --- | --- | --- | --- |
| obdeect structure | 0.8945 | 0.8917 | 0.8897 | 0.8866 | 0.8855 |
| `telescope_transmission` 2.0.0 | 0.8980 | 0.8967 | 0.8939 | 0.8898 | 0.8846 |

`simtools-validate-optics` for MSTS-01 (production 7.0.0, both backends)
agrees in effective area within 0.4% on axis and 1.2% up to 4 deg; obdeect
shows the azimuthal dependence (up to 2% between N and E) that the radial
parameterisation cannot. The structure also widens D80 by 1.4% on axis (same
photons with and without structure), since it shades the dish unevenly.

## Calibration light (Phase 2)

`simtools-simulate-flasher --run_mode full_simulation` and
`simtools-simulate-illuminator` with `--ray_tracing_backend obdeect`: obdeect
makes and traces the light, simtools writes the photons at the camera entrance
as photon bunches, and sim_telarray simulates the camera with
`Bypass_Optics=1`. The reference is `ff-1m` or `xyzls` followed by sim_telarray
(simulation models 7.0.0).

- Photon bunch file: the bunches of an `ff-1m` file, written again by simtools,
  give an identical sim_telarray output (photoelectrons of every pixel).
- Flasher (MSFx-FlashCam on MSTS-04, 200 pulses of 10^6 photons):
  photoelectrons 0.24% and charge 0.46% below `ff-1m`, as expected from the
  camera position (`focus_offset` 2.55 cm, 0.31%); time centroid of the camera
  sum 9.80 and 9.79 samples. Two `ff-1m` runs differ per pixel by 1.2 times
  Poisson, as obdeect does from `ff-1m`.
- Illuminator, mirrors without random errors (10^10 photons):

| Case | Photoelectrons sim_telarray / obdeect | Image centre [m] sim_telarray / obdeect |
| --- | --- | --- |
| ILLN-01 to MSTN-02 (axes offset 1.6 m) | 58969 / 58843 | (0.0122, -0.0025) / (0.0119, -0.0017) |
| ILLN-02 to LSTN-01 | 56902 / 56802 | (-0.0084, 0.0024) / (-0.0106, 0.0023) |
| LSTN-01, source 2 deg off axis (camera rotated 10.9 deg) | 62343 / 61844 | (-0.8679, -0.3747) / (-0.8701, -0.3726) |

  Image widths agree within 2%. Single illuminator pulses fluctuate more than
  Poisson: two sim_telarray runs of the off-axis case differ by 1.9% in total
  and per pixel by 2.1 times Poisson, obdeect from sim_telarray by 0.8% and 2.3.
- With random mirror errors, each realisation shapes the defocused image of a
  source at 600 m: two runs of either code differ per pixel by 7 to 9 times
  Poisson; the totals agree within 0.6%.
