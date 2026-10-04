# Reference run record (step 1)

For current implementation limits and open tasks, read [STATUS.md](STATUS.md).
This file retains technical recipes or contracts, not the current completion ledger.

This is a reproducibility record for a matched optical comparison. It fixes
the sim_telarray run configuration, telescope model, site and atmosphere,
random seeds, common input photon lists, software revisions, and commands.
Freezing the record does not trace photons or produce a PSF result. The
`execute` action runs the recorded commands and retains their logs.

A comparison begins with a JSON configuration containing these required keys:

- `sim_telarray_release`, `hessio_decoder`, `site`, `production_version`: exact selected identifiers.
- `telescope_variants`: production table names present in that version of the supplied checkout.
- `configuration_overrides`, `seeds`, `atmosphere_extinction`, `coordinate_frames`, `software_revisions`: explicit, nonempty objects recording the values used by both engines. `software_revisions` records the exact simtools/obdeect revisions or container digests used by the commands.
- `photon_blocks`: absolute paths to resolved input blocks.
- `commands`: ordered objects with `name`, `argv`, `cwd`, `environment`, and `inputs`. The executable, working directory, and input paths must be absolute. Put referenced configuration files and scripts in `inputs`; record every override in `configuration_overrides`.

The `obdeect-reference-run` tool records the exact sim_telarray comparison setup: selected telescope model records, the model checkout revision, and SHA-256 hashes for all input files and command executables. Verification resolves the selected production again and fails if a file or revision changed. The caller supplies version identifiers; the tool does not claim to detect an installed sim_telarray or hessio release.

```sh
obdeect-reference-run --action freeze --input reference.json \
  --model-root /path/to/simulation-models --output frozen.json
obdeect-reference-run --action verify --input frozen.json \
  --model-root /path/to/simulation-models
```

The reference-run record is a prerequisite for comparisons, not a reference result. A comparison harness must execute the recorded commands and retain their outputs, hashes, tolerances, and environment details. No production comparison is available until the reference executables and photon blocks are supplied.

## Execute the frozen commands

```sh
obdeect-reference-run --action execute --input frozen.json \
  --model-root /path/to/simulation-models --output comparison-run-001
```

The output directory must be new. Each command runs with its recorded argument
list and working directory, without a shell. Its declared environment settings
replace the corresponding inherited settings. Record all settings relevant to
optical physics in the configuration; inherited credentials are not written to
reports. Inputs and executables are checked against the frozen SHA-256 hashes
before and after each command. Commands must not modify a frozen input in place;
even a successful exit is rejected when an input hash changes.

The runner writes `000.stdout.txt`, `000.stderr.txt`, and subsequent numbered
logs, plus `execution.json` containing commands, return codes, elapsed time, and
log hashes. A nonzero exit or timeout stops the sequence and retains the failure
report. Success means the recorded programs executed successfully; the numerical
comparison gate must also pass before agreement can be claimed.
