# Reference manifest (step 1)

A comparison begins with a JSON configuration containing these required keys:

- `sim_telarray_release`, `hessio_decoder`, `site`, `production_version`: exact selected identifiers.
- `telescope_variants`: production table names present in that version of the supplied checkout.
- `configuration_overrides`, `seeds`, `atmosphere_extinction`, `coordinate_frames`, `software_revisions`: explicit, nonempty objects recording the values used by both engines. `software_revisions` records the exact simtools/obdeect revisions or container digests used by the commands.
- `photon_blocks`: absolute paths to resolved input blocks.
- `commands`: ordered objects with `name`, `argv`, `cwd`, `environment`, and `inputs`. The executable, working directory, and input paths must be absolute. Put referenced configuration files and scripts in `inputs`; record every override in `configuration_overrides`.

The `obdeect-reference-manifest` tool freezes the selected model records and assets, the model checkout revision, and SHA-256 hashes for all listed files and command executables. Verification resolves the selected production again and fails if a file or revision changed. The caller supplies version identifiers; the tool does not claim to detect an installed sim_telarray or hessio release.

```sh
obdeect-reference-manifest freeze reference.json \
  --model-root /path/to/simulation-models --output frozen.json
obdeect-reference-manifest verify frozen.json \
  --model-root /path/to/simulation-models
```

The manifest is a prerequisite for comparisons, not a reference result. A comparison harness must execute the recorded commands and retain their outputs, hashes, tolerances, and environment details. No production comparison is available until the reference executables and photon blocks are supplied.
