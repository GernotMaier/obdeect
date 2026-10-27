# Shadowing

Camera bodies, baffles and telescope structures can block photons on their way
to the focal surface. These **obscurers** act on the `incoming`, `after_primary`
or `after_secondary` legs. A blocked photon records the obscurer's id in
`terminal_surface`.

## MST structure

```{figure} _static/figures/mst_shadowing.png
:alt: MST structure coloured by the fraction of light each element alone blocks.

MSTS-01 (FlashCam), simulation models 7.0.0: an on-axis star at infinity,
1,000,000 photons, zenith angle 20 degrees and source seed 11. [1](#figure-provenance)
```

Dotted elements are never crossed in this sample. Blue shows the fraction of unobscured light that each element alone blocks;
grey elements stop no light on their own in this sample. The nominal layout is
drawn by simtools, with shadowing calculated by obdeect.

## Which elements matter?

`OpticalModel.obscurer_hits` traces photons without obscurers, then checks every
obscurer along paths that reach the focal surface. It reports summed input
weights for each element:

| Quantity | Meaning |
| --- | --- |
| `crossed` | Light whose path meets the element. |
| `first` | Light absorbed first by this element in the full optical model. |
| `sole` | Light stopped only by this element; removing it would recover this light. |

```python
from obdeect import OpticalModel

model = OpticalModel.load("mst.json")
photons = model.star_photons(1_000_000, seed=11)
hits = model.obscurer_hits(photons, seed=12)
share = hits["sole"] / hits["unobscured"]
```

The command-line equivalent is:

```sh
obdeect obscurers --model mst.json --seed 11 --photons 1000000 --field-deg 0 0
```

This diagnostic checks every obscurer on every leg and is intended for studies.
An element that blocks nothing in one sample may matter at another source
position. Check the full field of view and internal calibration sources before
removing elements. The [optical elements](optical-elements.md#telescope-structure)
section describes the structure table.

## Figure provenance

[1] [MST shadowing figure](_static/figures/telescopes-provenance.json).
