# Telescopes

Nominal optical layouts of the CTAO telescope models, drawn with
`simtools-plot-telescope-structure`. Each drawing shows the view from the sky
and the side, with distances in metres. Mirrors are blue, camera housings dark
grey, and other structures light grey.

These are optical-model drawings, not complete engineering models. Only
components present in the input data are shown; detailed support structures
are currently available for the MST. Shadowing is shown separately in
[Shadowing](shadowing.md).

## Large-Sized Telescope (LST)

```{figure} _static/figures/telescope_lst.png
:alt: LST mirror facets and camera housing, viewed from the sky and the side.

LSTN-01: a single segmented primary mirror and the camera at its focus.
```

## Medium-Sized Telescope (MST)

::::{tab-set}

:::{tab-item} NectarCam

```{figure} _static/figures/telescope_mst_nectarcam.png
:alt: MST with NectarCam, viewed from the sky and the side.

MSTN-01: the segmented primary mirror and NectarCam housing.
```

:::

:::{tab-item} FlashCam

```{figure} _static/figures/telescope_mst_flashcam.png
:alt: MST with FlashCam and its support structure, viewed from the sky and the side.

MSTS-01: the segmented primary mirror, FlashCam housing and telescope structure.
```

:::

::::

## Small-Sized Telescope (SST)

```{figure} _static/figures/telescope_sst.png
:alt: SST primary and secondary mirrors, baffle and camera, viewed from the sky and the side.

SSTS-01: a dual-mirror optical layout with a curved focal surface.
```

## Schwarzschild-Couder Telescope (SCT)

```{figure} _static/figures/telescope_sct.png
:alt: SCT primary and secondary mirrors, baffle and camera, viewed from the sky and the side.

SCTS-design: the dual-mirror SCT design model.
```

All drawings use simulation models 7.0.0. The source revisions, input and image
hashes are recorded in the [figure provenance](_static/figures/telescopes-provenance.json).
See the [reproduction workflow](https://github.com/GernotMaier/obdeect/blob/main/docs/figures/telescopes.yml).

:::{dropdown} Reproduce the drawings

Run the checked-in workflow in a simtools environment with obdeect installed.
Use the simulation-models Git revision recorded in the provenance file; working-tree
edits are ignored.

```sh
simtools-run-application --config_file docs/figures/telescopes.yml \
  --simulation_models_git_path /path/to/simulation-models \
  --simulation_models_git_revision 2ed8d06347c8552367288eacab28ec0b75aaf01c \
  --replace __SIMULATION_MODELS__=/path/to/simulation-models \
  --replace __MODEL_REVISION__=2ed8d06347c8552367288eacab28ec0b75aaf01c \
  --replace __OUTPUT__=/path/to/output
```

Each output directory contains `plot_telescope_structure_<telescope>.png`.
Only the final workflow step enables shadowing, for the FlashCam MST.
:::
