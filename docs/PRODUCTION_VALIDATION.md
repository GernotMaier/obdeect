# Scientific validation of optical transport

[STATUS.md](STATUS.md) records the production requirements that remain open.
The comparison gate is a measurement tool, not a way to qualify incomplete
physics. Imported models currently remain nominal.

## Prepare the same experiment in both tracers

Select the simulation-models version and telescope variant explicitly. Freeze
one resolved photon input with positions in metres, unit directions, wavelengths
in nanometres, emission times in nanoseconds, input weights, and stable
run/event/array/telescope/bunch/photon identities. Set the coordinate transform,
atmospheric attenuation, detector acceptance, mirror losses, and seeds once.

Native replay uses the `CsvPhotonReader` header shown in
`cpp/include/obdeect/photon_input.hpp`. Pass it to the installed command:

```sh
obdeect-compile-optical-model --source-root /path/to/simulation-models \
  --model LSTN-design --version 7.0.0 --output lst.optical-model.json
obdeect-simtools-raytrace --optical-model lst.optical-model.json \
  --photon-input common-photons.csv --input-block-size 4096 \
  --output obdeect-arrivals.csv
```

This is a nominal replay. Adding `--require-trace-ready` to compilation or
`--require-production-ready` to tracing rejects unresolved production inputs.
Keep those gates for qualification runs; do not clear blockers to obtain a pass.

sim_telarray's star-light imaging list is a detected-photon list, with positions
in centimetres and no complete shared-photon terminal record. The current
`simtel_reference_psf.py` computes diagnostic PSFs from these lists; it cannot
recover undetected IDs or interactions. A production comparison requires a
sim_telarray adapter that actually records those quantities from the same
resolved input. Never assign sequential IDs to detections to manufacture
photon-by-photon correspondence. The local reference implementation is
`../sim_telarray/common/sim_imaging.c`; configuration/source generation is in
`../sim_telarray/common/sim_config.c`. The simtools runners are in
`../simtools/src/simtools/simtel/`.

The available container exposes
`/workdir/simulation_software/sim_telarray/bin/sim_telarray` from release
`2025.246.0`. Probe it with:

```sh
podman run --rm ghcr.io/gammasim/simtools-dev:latest \
  /workdir/simulation_software/sim_telarray/bin/sim_telarray -h
```

It currently contains no resolved `.lis` or `.iact` input, so it cannot produce
a matched validation record by itself. Supply a frozen common input and the
adapter before claiming a comparison result.

Record each executable and script, full argv, working directory, environment
settings, model/source hashes, and all input files with
[REFERENCE_RUN.md](REFERENCE_RUN.md). The reference executor runs those recorded
commands and preserves stdout/stderr and failures. The [shared-photon reference adapter](SIMTEL_REPLAY.md) exercises all four real optical classes, with stable input identities and explicit raw lost rows. It still needs complete terminal-loss instrumentation and downstream response before production normalization; it does not provide qualified comparison data.

## Declare acceptance before measuring residuals

Normalize both outputs to one CSV row per resolved input photon. Required base
columns are:

```text
photon_id,status,focal_x_m,focal_y_m,path_length_m,arrival_time_ns,incidence_primary_deg,incidence_secondary_deg,incidence_focal_deg
```

Production rows additionally require `wavelength_nm`, `source_weight`,
`throughput`, `terminal_surface_id`, `final_dx`, `final_dy`, `final_dz`,
`optical_model_sha256`, and `source_sha256`. Preserve full identity-context
columns whenever supplied. Ordered `interaction_surface_ids` support the first
recorded-surface divergence; full interaction positions/normals remain a separate
requirement. For detected photons, throughput is arriving optical weight; a
zero-response detector crossing is legitimate and must not be changed to a miss.

A fixture declares the experiment's coverage. For example, replace these hashes
and surface IDs with the actual frozen values:

```json
{
  "minimum_detected": 100,
  "required_surfaces": ["198"],
  "optical_model_sha256": "<64 lowercase hexadecimal characters>",
  "source_sha256": "<64 lowercase hexadecimal characters>"
}
```

A negative experiment can declare `minimum_detected: 0`,
`required_surfaces: []`, and `allow_all_loss: true`. Positive experiments must
include useful detector coverage. An empty/all-loss comparison cannot establish
optical image agreement.

Define tolerances from analytic precision and reference uncertainty studies.
All four base tolerance keys are required: `focal_position_m`, `path_length_m`,
`arrival_time_ns`, and `incidence_deg`. Numerical values are scientific choices,
not defaults to adjust until a run passes. Freeze them alongside the fixture.

```sh
obdeect-normalize-arrivals --input obdeect-arrivals.csv \
  --optical-model-sha256 "$MODEL_SHA256" --source-sha256 "$SOURCE_SHA256" \
  --output obdeect.normalized.csv
obdeect-validate-production --optical-model lst.optical-model.json \
  --telescope-family LST --fixture fixture.json --tolerances tolerances.json \
  --obdeect-arrivals obdeect.normalized.csv --simtel-arrivals simtel.normalized.csv \
  --output validation-summary.json
```

Declare `MODEL_SHA256` from the compiled content hash, not the hash of its
pretty-printed file. `SOURCE_SHA256` binds the frozen common photon bytes. The
normalizer preserves recorded arrival time and requires it when binding source
or optical-model hashes. Vacuum timing from older records is only a diagnostic
compatibility conversion without those bindings. Missing interaction sequences
remain unavailable rather than becoming invented empty sequences. Lost rays do not acquire invented
surface IDs. Failure summaries are written to the requested output path.

## Evidence to show the astronomy community

| Study | Demonstration |
| --- | --- |
| Analytic reflector/interface limits | Image position, timing, Snell/Fresnel and energy closure agree with independent calculations. |
| Shared photons and interaction sequences | The same photon reaches the same surfaces; identify the first divergent interaction. |
| On-axis and two-dimensional field maps | Weighted centroids, D50/D68/D80/D95 and radial CDF residuals over the stated field. |
| Panel, gap and shadow boundaries | Edge scans distinguish mirror gaps, obscurers, detector gaps and concentrator losses. |
| Spectral and angular response | Throughput residuals versus wavelength and angle, with partial optical losses included. |
| Finite sources and timing | Source-distance scans, absolute illuminator normalization, arrival-time distributions and path residuals. |
| Effective collecting area | Declared launch area times arriving/input optical weight; account for every lost component. |
| Windows and guides | Spot shifts, optical/group delay, transmission and independent ROBAST comparisons. |
| Reproducibility | Reordered photons, different blocks and thread counts give the same identity-keyed results. |
| Performance | At least ten repeats, median/spread and memory use for matched validated physics; separate kernel and I/O. |

Run these studies separately for LST, MST FlashCam, MST NectarCam, SST and SCT,
with exact production variants recorded. SST evidence does not qualify SCT.
Store residual tables, weighted CDFs, loss closure, uncertainties, first divergent
interactions, compiler/host details and hashes. Retain large raw files externally
with retrievable locations and hashes. Current passing unit tests and historical
PSF figures do not constitute this completed matrix.
