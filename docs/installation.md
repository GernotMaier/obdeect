# Installation

Install the Python package and the `obdeect` command in a virtual environment:

```sh
python3.14 -m venv .venv
source .venv/bin/activate
python -m pip install obdeect-dev
obdeect --help
```

The package requires Python 3.14 or later. If no wheel is available for your
platform, installation builds the extension and requires a C++20 compiler.

## Install from source

From an obdeect checkout, install the package in editable mode:

```sh
python -m pip install -e ".[dev]"
```

To build the documentation, also install the documentation dependencies:

```sh
python -m pip install -e ".[docs]"
sphinx-build -W --keep-going -b html docs docs/_build/html
```

## Build options

- `-DOBDEECT_USE_FMA=ON` (CMake, default OFF) evaluates the aspheric polynomials
  with `std::fma`. Use it only on hardware with fused multiply-add: without it
  `std::fma` is emulated and very slow. Results differ from builds without it
  at the rounding level. It is disabled by default.
- Results also remain platform-dependent in the last bits, because the math
  library (`std::pow`, `std::hypot`, trigonometric functions) differs between
  platforms; see [Validation](VALIDATION.md#rounding-level-changes).
