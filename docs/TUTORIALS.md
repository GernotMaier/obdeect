# Tutorials

Run these commands from the repository root after [building the project](../README.md).

## 1. Explore the reference telescope

```sh
obdeect-plot-reference --view structure --telescope reference-mst --output  reference_structure.png
./build/debug/obdeect_reference --photons 10000 --output  reference.csv
obdeect-plot-reference  reference.csv --view rays --max-paths 150 --output  reference_rays.png
obdeect-plot-reference  reference.csv --view focal-plane --bins 80 --output  reference_focal.png
```

The structure view has x-z and y-z projections. Rays are coloured by terminal
status; blue rays reach the focal screen. The focal image contains only
`detected` rows, weighted by `source_weight × throughput`.

## 2. Change the source

```sh
./build/debug/obdeect_reference --source star --field-x-deg 0.5 --output off_axis.csv
./build/debug/obdeect_reference --source illuminator --distance-m 50 --output flasher.csv
./build/debug/obdeect_reference --source laser --distance-m 50 --divergence-deg 0.1 --output laser.csv
obdeect-plot-reference off_axis.csv --view rays --output off_axis.png
```

The star is a plane wave. The illuminator is a finite distance point source
with relative inverse square weights. The laser samples a beam with the given
angular divergence. All three are artificial calibration sources.

## 3. Compare analytic optical references

```sh
./build/debug/obdeect_ctao --telescope LST --photons 10000 --output lst.csv
obdeect-plot-reference --telescope LST --view structure --output lst_outline.png
obdeect-plot-reference lst.csv --telescope LST --view rays --output lst_rays.png
obdeect-plot-reference lst.csv --telescope LST --view focal-plane --output lst_focal.png
```

`MST` is another tested continuous-mirror baseline. LST/MST structure views
are optical surface outlines only. SST/SCT outlines are marked unavailable
until their two mirror geometry is validated.

## 4. Measure a PSF and scan field angle

```sh
obdeect-psf derive lst.csv --field-x-deg 0 --output lst_psf.json
obdeect-psf scan --executable ./build/debug/obdeect_reference --photons 10000 \
  --field-x-deg 0 0.5 1.0 --output-dir out/scan --plot out/scan.png
```

`derive` reports weighted centroid, R80/D80, throughput, and terminal losses.
The scan writes one trace per offset and a summary table.

## 5. Inspect a production model

```sh
obdeect-import-simulation-models /path/to/simulation-models LSTN-design \
  --version 6.3.0 --output lst.ir.json
obdeect-compile-scene lst.ir.json --source-root /path/to/simulation-models \
  --output lst.scene.json
```

The IR hashes selected records and assets. The compiler extracts supported
mirror and camera geometry and lists unresolved trace requirements. It does
not produce a trace executable or a validated production model.
