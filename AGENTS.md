# obdeect Agent Guide

## Scope

`obdeect` is a portable C++20 optical ray tracer with a thin Python interface.
It traces photons from a defined entrance surface to a focal/sensor surface.
Air-shower generation, electronics, triggers, and event writing are out of
scope. Read the relevant code and tests first; use the workspace plans and
validation matrix for architectural or scientific changes.

## Design rules

- Never ever modify the main README.md. Suggest instead changes, don't apply them.
- Terminology: this project is for gamma-ray astronomy. Always call the
  compiled telescope geometry and optical-response artifact an "optical model";
  never call it a "scene" in user-facing text, documentation,
  diagnostics, or new APIs.
- Keep it simple, efficient, deterministic, and portable. Prefer small,
  explicit data structures and testable physical kernels over frameworks,
  global state, hidden defaults, or abstraction in the hot path.
- Keep the core C++ standard-library-only. Python, plotting, file import, and
  external-simulator dependencies belong at the boundary, never in a trace
  kernel.
- Use immutable compiled optical models and contiguous SoA photon buffers. After
  output buffers are sized, the per-photon path must not allocate, perform I/O,
  call Python, or dispatch virtually.
- Make units, frames, random seeds, statuses, and loss accounting explicit.
  Reject invalid or unsupported input; never silently discard fields.

## Data, not telescope-specific code

- Do not hardwire CTAO names, telescope geometry, model versions, file paths,
  or optical/physics parameters into the core, public API, or defaults.
- Geometry, response tables, trace modes, and conventions are validated input
  data. Model import produces a generic, provenance-recorded IR.
- CTAO support is an explicit importer/adapter or versioned fixture only. It
  must not be required to build, test, or run the generic engine, and it must
  never be an implicit fallback.
- Use python >3.14. Do not fall back to earlier versions.

## Changes and validation

- Make the narrowest complete change and preserve unrelated work.
- Add or update a focused test for every behavioural change. Establish analytic
  kernel and primitive tests before optical-model, cross-tool, or performance tests.
- Preserve results across photon ordering, blocks, and thread counts. Record
  provenance for fixtures and cross-tool comparisons.
- Do not weaken tolerances or replace golden data to hide a discrepancy. Keep
  examples and fixtures clearly labelled and configuration-driven.
- When renaming a public concept, update implementation symbols, imports,
  exception types, serialized keys, fixture variables, and test call sites
  together. For scene-to-optical-model renames, do not leave stale references
  such as `write_native_scene`, `SceneCompileError`, `scene_sha256`, or local
  variables named `scene` in optical-model tests unless they intentionally test
  backward compatibility.

## Build, test, and lint

Run from this directory. Configure before static analysis so the compilation
database is current.

```bash
cmake --preset debug
cmake --build --preset debug
ctest --test-dir build/debug --output-on-failure
python -m unittest discover -s python/tests
ruff format --check python
ruff check python
find cpp -type f \( -name '*.cpp' -o -name '*.hpp' \) -print0 | \
xargs -0 clang-format --dry-run --Werror
```

Ruff is configured with a 100-character line limit. Run `ruff format` before
`ruff check`, and manually wrap long strings or boolean expressions that the
formatter leaves unchanged.

After renames, search the affected Python files for removed symbols and run
the focused tests before the full suite. A clean handoff requires both
`ruff format --check python` and `ruff check python`; do not rely on a passing
test run alone because stale test references can fail collection or linting.

Run `clang-tidy -p build/debug cpp/src/<changed-file>.cpp` for changed C++
translation units when available, and `pre-commit run --all-files` before
handoff. Start with the smallest relevant test; run the full checks for shared
kernels, interfaces, configuration, or tooling changes.
