# Validation status

This file reports executed checks, their exact scope, and the boundary between
verified kernels and deferred CTAO MST observables. It does not treat a passing
toy-model test as validation of a different telescope model.

## Executed 2026-09-18

### IACTrace analytic cross-check — pass

`validation/compare_iactrace_spherical.py` ran the same 256 deterministic,
parallel rays in obdeect and local IACTrace `0.10.3`. Both traces used binary64
and the shared geometry was deliberately limited to one continuous spherical
primary (`R=9.75 m`, aperture radius `6 m`) and a plane at `R/2=4.875 m`.
Masts and camera shadowing were disabled because the reference scene contains
no matching obstruction definition.

| Quantity | Result | Acceptance |
| --- | ---: | ---: |
| maximum primary hit-point residual | `2.89e-15 m` | `<= 1e-10 m` |
| maximum focal-plane hit-point residual | `3.16e-14 m` | `<= 1e-10 m` |

The result validates the shared ray/sphere-root/reflection/plane-intersection
path. It does not validate throughput, wavelength response, timing, facets,
or an MST configuration.

Repeat after building `obdeect_toy`:

```bash
python validation/compare_iactrace_spherical.py \
  --executable build/debug/obdeect_toy --photons 256 --tolerance-m 1e-10
```

The script needs an environment containing IACTrace and JAX; it intentionally
is not a normal project dependency.

### Local MST literature check — pass at the supported scope

Garczarczyk, *MST MC parameters* (27 September 2017), p. 2, specifies a
modified Davies-Cotton MST with a 13.77-m dish, 16.00-m
central-mirror-to-focal-plane distance, spherical 1.2-m hexagonal facets, and
up to 90 mirrors. Page 19 lists the central facet radius of curvature
`32.14 m` and focal length `16.07 m`, satisfying the single-spherical-facet
relation `f=R/2`. The PDF was examined in the parent research workspace and
is not redistributed in this repository.

`T-LIT-001` in `cpp/tests/test_literature_mst.cpp` traces the corresponding
on-axis, one-facet analytic case and verifies a focal-plane-origin landing.
This test passes.

## Full MST literature validation — not yet eligible

The same paper specifies effective mirror area above `88 m²`, optical
`theta_80 < 0.18°`, and RMS optical time spread below `0.8 ns` for the
modified Davies-Cotton telescope. Garczarczyk, *MST Optics Parameters*
(20 March 2020), pp. 1–2, also reports the 86-facet CSS design, approximately
`90.68 m²` effective mirror area, `15.37 m²` total shadowing, and off-axis
PSF values. This slide deck likewise is not redistributed here.

The current executable is a continuous spherical dish with an illustrative
camera disk and four masts. It has no hexagonal facet layout, Davies-Cotton
dish placement, CTAO camera geometry, cable/rope/frame model, reflectivity
tables, camera pixels, or optical path-time distribution. Therefore it cannot
be compared to those published area/PSF/time values, and no “good MST result”
is claimed for them.

The next eligible validation gate is to import a pinned 86-facet MST geometry
and matched CSS data, then compare equal input rays at 0° and the published
off-axis angles for effective area, `theta_80`, component shadow maps and RMS
time spread. That gate is already the Step-2/Step-4/Step-5 requirement in the
project plan; it must remain failing/unsupported until those features exist.
