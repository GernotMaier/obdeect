# obdeect

`obdeect` is a C++20 optical ray tracing prototype for imaging atmospheric
Cherenkov telescopes. The C++ core uses only the standard library. Python
provides plotting, PSF analysis, and model import tools.

It currently supports **nominal geometry studies**, not full CTAO production simulations.

For example use cases, see [optical-study recipes](docs/USE_CASES.md).

## User Installation

Install from PyPI with:

```sh
python -m venv .venv
source .venv/bin/activate
python -m pip install obdeect-dev
```

(requires python 3.14 or higher)

Run a simple test to verify the installation.
The following command simulates a simple ray tracing scenario of a mid-size telescope:

```sh
obdeect-simtools-raytrace --telescope MST --photons 10000 --output trace.csv
```

## Developer Installation

```sh
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'
cmake --preset debug   # add '--target clean' if you want to clean the build directory before building
cmake --build --preset debug
```

Executables are in `build/debug/` and tests in `build/debug/tests/`.

Test your installation with:

```sh
ctest --test-dir build/debug --output-on-failure
python -m unittest discover -s python/tests
```

## License and Citation

BSD-3-Clause license. Citation metadata is in [CITATION.cff](CITATION.cff).

## Generative AI disclosure

Generative AI tools (mostly ChatGPT 5.6) were used to write the entire code of this project. All AI-assisted outputs were reviewed, validated, and, where necessary, modified by the authors to ensure accuracy and reliability.
