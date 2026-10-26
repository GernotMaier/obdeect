# obdeect

```{image} _static/images/obdeect_21ratio_small.jpg
:alt: obdeect — optical ray tracing for gamma-ray astronomy
:align: center
:class: homepage-logo
```

Optical ray tracing for gamma-ray astronomy, with a C++20 core, a Python
interface and a command line. Trace photons through mirrors and telescope
structure to the camera entrance, using a compiled optical model.

:::{container} homepage-actions

```{button-link} https://github.com/GernotMaier/obdeect
:color: primary
:shadow:

{fab}`github` GitHub Repository
```

```{button-link} https://github.com/GernotMaier/obdeect/blob/main/LICENSE
:color: secondary
:outline:

{fas}`scale-balanced` License BSD-3
```

:::

---

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} {fas}`bolt;pst-color-primary` User Guide
:link: user-guide
:link-type: doc
:class-card: sd-shadow-sm

Install obdeect and run your first trace from the command line or Python.
:::

:::{grid-item-card} {fas}`eye;pst-color-primary` Optics
:link: optics
:link-type: doc
:class-card: sd-shadow-sm

Understand mirrors, photon paths, optical responses and shadowing.
:::

:::{grid-item-card} {fas}`check-circle;pst-color-primary` Validation
:link: VALIDATION
:link-type: doc
:class-card: sd-shadow-sm

Explore accuracy checks, sim_telarray comparisons and measured performance.
:::

:::{grid-item-card} {fas}`book;pst-color-primary` Python API
:link: api
:link-type: doc
:class-card: sd-shadow-sm

Look up optical-model methods, compilation and analysis helpers.
:::

::::

See the [changelog](changelog.md) for release history.

```{toctree}
:hidden:
:maxdepth: 2

user-guide
optics
VALIDATION
api
changelog
```
