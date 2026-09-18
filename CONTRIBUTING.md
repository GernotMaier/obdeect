# Contributing to obdeect

`obdeect` is deliberately IACT-only: changes must concern optical transport of
Cherenkov photons through imaging atmospheric Cherenkov telescope optics,
camera, or structure. Do not add generic telescope or detector simulation
features without an explicit design decision.

## Local checks

Install CMake, Ninja, a C++20 compiler, Python 3.10+, and the development
tools:

```bash
python -m pip install --upgrade pip ruff
python -m pip install .
cmake --preset debug
cmake --build --preset debug
ctest --test-dir build/debug --output-on-failure
python -m unittest python/tests/test_plot_toy_mst.py
ruff format --check python
ruff check python
clang-format --dry-run --Werror cpp/include/obdeect/toy_mst.hpp cpp/src/toy_main.cpp cpp/tests/test_core.cpp
```

## Change rules

1. Add a focused test before or alongside every observable behaviour change.
2. State coordinate frames and units at every public boundary; use metres,
   nanometres and nanoseconds.
3. Keep the core free of Python, ROOT, EventIO and plotting dependencies.
4. Preserve deterministic inputs and output identifiers. Do not silently drop
   a photon or a loss category.
5. Treat the current continuous spherical dish and four-mast structure as a
   demonstrator, not a CTAO production model. Compatibility work needs a
   pinned reference fixture and an acceptance tolerance.

## Pull requests

Keep pull requests small. Explain the physical convention, list the test IDs
or tests changed, and attach a generated diagnostic plot only when it helps
review an optical effect. The CI suite must pass before merge.
