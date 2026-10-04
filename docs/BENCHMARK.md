# Matched optical performance measurements

`python -m obdeect.benchmark` executes explicit adapters without a shell. Each
adapter must produce a normalized comparison CSV accepted by
`production_validation.read_comparison_table`. The benchmark applies the same
predeclared fixture hashes, detection/surface coverage and optical tolerances as
[production validation](PRODUCTION_VALIDATION.md). It cannot qualify an unready
optical model, repair simulator output or infer matched physics from two command
names.

The configuration requires these fields:

```json
{
  "optical_model": "/absolute/compiled.optical-model.json",
  "source": "/absolute/resolved-photons.csv",
  "fixture": {
    "minimum_detected": 1,
    "required_surfaces": ["901"],
    "optical_model_sha256": "<canonical compiled optical-model SHA256>",
    "source_sha256": "<resolved photon file SHA256>"
  },
  "tolerances": {
    "focal_position_m": 1e-9,
    "path_length_m": 1e-9,
    "arrival_time_ns": 1e-9,
    "incidence_deg": 1e-9
  },
  "provenance": {
    "reference_run": "<frozen reference-run identifier>",
    "software_revisions": "<exact baseline and candidate revisions>",
    "physics": "<matched geometry, response, atmosphere and seed conventions>"
  },
  "repeats": 10,
  "timeout_s": 3600,
  "variants": [
    {"threads": 1, "block_size": 4096},
    {"threads": 4, "block_size": 16384}
  ],
  "inputs": ["/absolute/baseline-adapter.py", "/absolute/candidate-adapter.py"],
  "engines": {
    "baseline": {
      "argv": ["/absolute/python3.14", "/absolute/baseline-adapter.py",
               "--output", "{output}", "--threads", "{threads}",
               "--block-size", "{block_size}"],
      "cwd": "/absolute/run-directory",
      "environment": {},
      "timing_metrics": false
    },
    "candidate": {
      "argv": ["/absolute/python3.14", "/absolute/candidate-adapter.py",
               "--output", "{output}", "--threads", "{threads}",
               "--block-size", "{block_size}"],
      "cwd": "/absolute/run-directory",
      "environment": {},
      "timing_metrics": false
    }
  }
}
```

These illustrative tolerance values are not production acceptance criteria.
Choose tolerances and coverage from the scientific validation matrix before
running either engine. List every adapter, linked simulator binary,
configuration, data table and other command input in `inputs`; the command
executables, optical model and resolved source are included automatically.
Freeze the reference run separately when it needs model-checkout revision and
production selection verification.

```sh
PYTHONPATH=python python -m obdeect.benchmark --action freeze \
  --input benchmark-config.json --output frozen-benchmark.json
PYTHONPATH=python python -m obdeect.benchmark --action verify \
  --input frozen-benchmark.json
PYTHONPATH=python python -m obdeect.benchmark --action execute \
  --input frozen-benchmark.json --output /absolute/new-benchmark-directory
```

At least ten paired repetitions and two different thread counts and block sizes
are required. The candidate command must actually contain both variant
placeholders. Execution order alternates baseline/candidate to reduce a
consistent ordering bias. Frozen hashes are verified before and after every
command. Every result is compared photon by photon to its matched engine, and
every repeated/variant result must be exactly identical as normalized records
within that engine, including recorded terminal-loss diagnostics. Row order may
change. The first optical mismatch or reproducibility failure stops the run,
retains logs and a failed `benchmark.json`, and publishes no speedup summary.

A successful report records all invocations, output hashes, matched residuals,
wall time, user/system CPU and peak measured-process RSS for each run, then
median, minimum, maximum and standard deviation for each variant. Wall-time
speedup includes process startup, input and output. It is not kernel speedup.

For instrumented adapters, set `timing_metrics` to `true` and include a
`{metrics}` argument in `argv`. Write a JSON file with exactly `kernel_s` and
`io_s`, both finite nonnegative seconds from separate measured phases; their
sum must fit inside measured wall time. The harness reports their distributions
and kernel speedup only when both adapters supply these measurements and the
candidate median kernel time is positive. It never subtracts wall times to
invent an I/O estimate.

Resource accounting currently requires macOS or Linux `wait4`. Peak RSS is the
maximum resident set of the measured process, not simultaneous summed RSS of a
multiprocess tree. Inherited environment variables are not recorded; freeze all
physics-relevant environment overrides explicitly. Use sufficiently long runs
that process startup and the 1 ms child-completion polling interval are small
compared with the timing being assessed.

## Independent ROBAST checks

The local `../ROBAST` checkout was inspected at revision
`9ba75dfadcde6da25c4875e34b5669aabd266c2e`. ROOT (`root` and `root-config`) is absent
from PATH and the standard Homebrew installation locations, and the checkout
contains no built ROBAST library. No independent detailed-geometry comparison
was executed. A future ROBAST adapter must export the same normalized identity,
frame, wavelength, weights, directions, timing and interaction diagnostics and
bind the exact source/model hashes before this harness can compare or time it.
The available source checkout is not independent validation evidence.
