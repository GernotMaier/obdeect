# Project-owned model fixtures

Only tiny, redistributable small fixtures belong here. CTAO simulation-models,
CORSIKA files and external telescope data are resolved by explicit importers
and are never copied into this repository without provenance and licence review.

`tools/import_simulation_models.py` reads a separately obtained,
versioned `simulation-models` checkout and emits a JSON IR outside this
directory by default. The IR contains relative source paths and SHA-256 hashes
for the production manifest, every parameter record and each selected asset.
