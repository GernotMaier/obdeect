# Production validation

An obdeect trace is not a production CTAO result until this gate passes for
each supported telescope family. The required 7.0.0 matrix is LSTN-design,
MSTx-FlashCam, MSTx-NectarCam, and SSTS-design; run it for the selected site
and model variant rather than substituting an analytic prescription.

## Inputs

Use the same resolved photon block for both engines. Before comparison, adapt
each engine output to an ASCII CSV with exactly one row per resolved photon and
these columns:

```text
photon_id,status,focal_x_m,focal_y_m,path_length_m,arrival_time_ns,incidence_primary_deg,incidence_secondary_deg,incidence_focal_deg
```

`photon_id` is the stable identity from that common source block. Positions
are metres, paths are optical/geometric path metres as declared by the source
adapter, and incidence angles are degrees. Do not compare independently
sampled Monte Carlo runs: the gate rejects missing or added photon identities.

Declare tolerances before running, for example:

```json
{
  "focal_position_m": 0.001,
  "path_length_m": 0.001,
  "arrival_time_ns": 0.01,
  "incidence_deg": 0.01
}
```

Those values are an example only. The accepted values belong in the frozen
reference manifest after the uncertainty study; they must not be relaxed to
make a failed comparison pass.

## Reproducible command sequence

```sh
export OBDEECT_SIMULATION_MODELS_PATH=/path/to/simulation-models
SIMTEL_ROOT=/path/to/sim_telarray
MODEL=LSTN-design

obdeect-import-simulation-models \
  --source-root "$OBDEECT_SIMULATION_MODELS_PATH" --model "$MODEL" --version 7.0.0 \
  --output "$MODEL.ir.json"
obdeect-compile-optical-model \
  --input "$MODEL.ir.json" --source-root "$OBDEECT_SIMULATION_MODELS_PATH" \
  --simtel-root "$SIMTEL_ROOT" --require-trace-ready \
  --output "$MODEL.optical-model.json" --native-output "$MODEL.optical-model.csv"

# Generate the two normalised comparison tables from the frozen common source block.
obdeect-normalize-arrivals --input obdeect-arrivals.csv --output obdeect.normalized.csv
obdeect-validate-production \
  --optical-model "$MODEL.optical-model.json" --telescope-family LST \
  --obdeect-arrivals obdeect.normalized.csv \
  --simtel-arrivals simtel.normalized.csv \
  --tolerances tolerances.json --output validation-summary.json
```

The compiler gate fails before tracing if it cannot bind finite primary and
secondary surfaces where needed, detector geometry, material responses, or
other required optical transport. The comparison gate then fails on the first
different photon identity, terminal status, focal position, path length, or
incidence angle, and reports maximum residuals and the first divergence.

Run the matrix for on-axis and two-dimensional field offsets, selected panels,
spectral samples, stars, illuminators, lasers, EventIO input, deterministic
repeats, and the effective-area/focal-length/PSF/incidence observables. Store
the source block hashes, model and software revisions, atmosphere, seed,
commands, tolerances, and generated summaries in an
`obdeect-reference-manifest` before interpreting any result as validation.
