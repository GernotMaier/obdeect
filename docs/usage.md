# Command-line usage

Compile an optical model once, then reuse it for tracing and plotting.
The example below uses an explicit CTAO simulation-models checkout, telescope
and version; choose values available in your checkout. Plotting requires the
`obdeect[plot]` extra.

```sh
# Compile an optical model from a simulation-models checkout.
obdeect compile --simulation-models ../simulation-models --model LSTN-01 \
  --version 7.0.0 --seed 1 --zenith-deg 20 --output lst.json

# Trace a star 0.5 deg off axis; print the image summary, keep every photon.
obdeect trace --model lst.json --seed 2 --photons 100000 \
  --field-deg 0.5 0 --output lst-arrivals.npz

# Draw the dish seen from the sky and the focal-plane image.
obdeect plot --model lst.json --arrivals lst-arrivals.npz --output lst.png

# Which obscurers (camera, baffle, masts) shade a star 2 deg off axis.
obdeect obscurers --model lst.json --seed 2 --photons 100000 --field-deg 2 0
```

`compile` reads full and patch productions of the simulation models; a patch
production is laid over the productions in its `model_version_history`. `--seed` fixes the random mirror errors (focal
lengths, alignment); `--zenith-deg` is the telescope zenith angle used by
zenith-dependent focus offsets and alignment errors. `trace --seed` fixes source
sampling and mirror scattering. `--geometry-only` leaves the photon weights
unchanged (no reflectivity, telescope transmission or degradation). Results do
not depend on `--threads`, nor on photon order or block size.

Use `obdeect --help` or `obdeect trace --help` for the available options.
For structural shadow diagnostics, see [Shadowing](shadowing.md).
