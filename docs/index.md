# obdeect

obdeect is an optical ray tracer for imaging atmospheric Cherenkov
telescopes and one of the ray-tracing backends of [simtools](https://github.com/gammasim/simtools).

Photons might originate from the sky, a light source or a shower photon bunch, over the
mirrors, past the telescope structure, to the camera entrance. It is a C++20
library with a Python interface and a command line. The telescope is described
by an optical model compiled from the CTAO simulation
models.

::::{grid} 1 2 2 4
:gutter: 3

:::{grid-item-card} {fas}`bolt;pst-color-primary` User Guide
:link: user-guide
:link-type: doc
:class-card: sd-shadow-sm

Learn how to get started as a user.
This guid will cover how to install obdeect, and introduce the main functionalities.
:::

:::{grid-item-card} {fas}`eye;pst-color-primary` Optics and Validation
:link: physics-validation
:link-type: doc
:class-card: sd-shadow-sm

Learn how obdeect builds the optical model and traces photons through it.
Find validation results and comparisons with other ray-tracing backends.
:::

:::{grid-item-card} {fas}`book;pst-color-primary`  Python API
:link: api
:link-type: doc
:class-card: sd-shadow-sm
Learn about the `obdeect` Python API.
:::

:::{grid-item-card} {fas}`hammer;pst-color-primary` Changelog
:link: changelog
:link-type: doc
:class-card: sd-shadow-sm
Learn what has changed in each release of obdeect.
:::

::::

```{toctree}
:hidden:
:maxdepth: 1

user-guide
physics-validation
api
changelog
```
