# obdeect Agent Guide

## Scope

`obdeect` is a portable C++20 optical ray tracer with a thin Python interface.
It traces photons from a defined entrance surface to a focal/sensor surface.
Air-shower generation, electronics, triggers, and event writing are out of
scope. Read the relevant code and tests first; use `docs/VALIDATION.md` and
`docs/optical-models.md` for architectural or scientific changes.

## Layout

- `cpp/include/obdeect/`: C++ core headers (kernels, optical models, public API).
- `cpp/bindings/core.cpp`: nanobind extension; the only C++ file that may
  depend on Python.
- `cpp/tests/`: native tests (`test_*.cpp`), shared helpers (`check.hpp`),
  and fixtures in `cpp/tests/data/`.
- `cpp/bench/`: benchmarks and the replay comparison tool.
- `python/obdeect/`: Python boundary (compile, trace, formats, plotting, CLI).
- `python/tests/`: pytest suite; `telescopes.py` holds shared test helpers.
- `docs/`: Sphinx sources. `docs/changes/` is towncrier input (see below).

## Design rules

- Never ever modify the main README.md. Suggest instead changes, don't apply them.
- Terminology: this project is for gamma-ray astronomy. Always call the
  compiled telescope geometry and optical-response artifact an "optical model".
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
- Use Python >=3.14 (as in `pyproject.toml`). Do not fall back to earlier versions.

## Changes and validation

- Make the narrowest complete change and preserve unrelated work.
- Do not add anything under `docs/changes/` for current work. Leave that
  directory untouched unless the user explicitly asks for a changelog or
  release-note change.
- Add or update a focused test for every behavioural change. Establish analytic
  kernel and primitive tests before optical-model, cross-tool, or performance tests.
- Preserve results across photon ordering, blocks, and thread counts. Record
  provenance for fixtures and cross-tool comparisons.
- Do not weaken tolerances or replace golden data to hide a discrepancy. Keep
  examples and fixtures clearly labelled and configuration-driven.
- When renaming a public concept, update implementation symbols, imports,
  exception types, serialized keys, fixture variables, and test call sites
  together. Do not leave obsolete geometry terminology in optical-model tests.

## Build, test, and lint

Run from this directory. Configure before static analysis so the compilation
database is current. The `debug` preset enables `OBDEECT_BUILD_TESTING` and
writes to `build/debug/`; wheels build with `OBDEECT_BUILD_TESTING=OFF`.

```bash
python -m pip install -e '.[dev]'   # editable install; rebuilds the extension on C++ changes
cmake --preset debug
cmake --build --preset debug
ctest --test-dir build/debug --output-on-failure
python -m pytest python/tests
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

Run `clang-tidy -p build/debug <changed-file>` for changed C++ translation units
(`cpp/tests/`, `cpp/bench/`, `cpp/bindings/`, or `cpp/include/obdeect/` headers)
when available, and `pre-commit run --all-files` before handoff. Pre-commit runs
ruff, codespell, markdownlint, yamllint, actionlint, pyproject-fmt, shellcheck,
and the towncrier check. Start with the smallest relevant test; run the full
checks for shared kernels, interfaces, configuration, or tooling changes.

CI also builds and tests a `-DOBDEECT_USE_FMA=ON` Release configuration and an
AddressSanitizer/UBSan Debug configuration. Run those locally for changes to
floating-point kernels or memory-handling code:

```bash
cmake -S . -B build-sanitize -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_FLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer" \
  -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=address,undefined"
cmake --build build-sanitize --parallel
ctest --test-dir build-sanitize --output-on-failure
```

Docs build with warnings as errors: `sphinx-build -W --keep-going -b html docs docs/_build/html`.
Regenerate figures only from new validation results, via `docs/figures/make_figures.py`.
