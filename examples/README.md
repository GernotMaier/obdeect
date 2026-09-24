# Runnable optical and source examples

From the repository root:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
python -m pip install .
python examples/run_examples.py --build build --output out/examples --photons 1000
```

The script runs a toy telescope with a star, point flasher and divergent laser;
an ideal parabolic LST on/off axis; and the continuous MST sphere reference.
Each case writes CSV, a PNG of the actual paths, and status counts in
`summary.json`. The script fails on subprocess errors, malformed CSV,
missing photons, or missing focal-plane hits for cases that are expected to detect light.

These examples exercise analytic geometry and source sampling. Their focal
screen counts are not photoelectron counts or CTAO effective areas. The point
source currently carries relative inverse-square weights, not an absolutely
normalized flasher yield; the laser is a sampled beam, not atmospheric laser
scattering. See the production review for the missing calibration physics.
SST/SCT full-model examples need curved focal surfaces, segment gaps and
obscurations before their outputs can be used as reference results.
