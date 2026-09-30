# obdeect

`obdeect` is a C++20 optical ray tracing prototype for imaging atmospheric
Cherenkov telescopes. The C++ core uses only the standard library. Python
provides plotting, PSF analysis, and model import tools.

It currently supports **nominal geometry studies**, not full CTAO production simulations.

For panel tests and more plots, see [optical-study recipes](docs/USE_CASES.md).

## User Installation

Install from PyPI with:

```sh
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

Run a simple test to verify the installation:

```sh
obdeect-simtools-raytrace --telescope MST --photons 10000 --output trace.csv
```

## Developer Installation

```sh
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'
cmake --preset debug
cmake --build --preset debug
```

Executables are in `build/debug/` and tests in `build/debug/tests/`.

Test your installation with:

```sh
ctest --test-dir build/debug --output-on-failure
python -m unittest discover -s python/tests
```

## Recipes

See [optical-study recipes](docs/USE_CASES.md).

## License and Citation

BSD-3-Clause license. Citation metadata is in [CITATION.cff](CITATION.cff).

## Generative AI disclosure

Generative AI tools were used to write much of this project; outputs were
reviewed and validated by the authors.

Generative AI tools (mostly ChatGPT 5.6) were used to write the entire code of this project. All AI-assisted outputs were reviewed, validated, and, where necessary, modified by the authors to ensure accuracy and reliability.
