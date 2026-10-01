# Project-owned model fixtures

Only tiny, redistributable small fixtures belong here. CTAO simulation-models,
CORSIKA files and external telescope data are resolved by explicit importers
and are never copied into this repository without provenance and licence review.

`obdeect-compile-optical-model` reads a separately obtained, versioned
`simulation-models` checkout and records relative source paths and SHA-256
hashes for the production manifest, every parameter record, and each selected
asset in the compiled optical-model output.
