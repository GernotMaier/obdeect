# Photon transport through an optical model

This guide describes the implemented optical path from a source ray to its
terminal state. Geometry and response functions come from the compiled optical
model. A response multiplies accumulated optical throughput, while source
weight remains a separate quantity in the output. A response does not by itself
decide whether a ray intersects a later surface. A geometrical obstruction
terminates the path.

The reflector trace begins at the supplied photon position and direction; it
does not generate an air shower or bend rays through an atmospheric profile.
For reflector models, an imported propagation group index affects arrival time,
not the ray trajectory, unless atmospheric media are explicitly part of a
non-sequential optical model.

The stages below describe the segmented and axisymmetric telescope tracers.
An optical model may omit an optional response, in which case that response is
unity. The non-sequential material tracer follows the same input and accounting
conventions but advances through its declared surfaces and media instead of a
fixed primary/secondary/detector sequence.

## 1. Source ray and input frame

The source adapter supplies a position, unit direction, wavelength, emission
time, weight, and photon identity. A plane-wave star has parallel directions;
its field angles set the direction and phase across the sampled entrance disk.
A finite-distance star or illuminator has position-dependent directions. A
laser samples directions within its configured divergence. Supplied photon
records are replayed without resampling. Ground-frame replay is transformed to
the telescope frame at the input boundary.

**Verification.** `cpp/tests/test_sources.cpp`,
`cpp/tests/test_source_sampling.cpp`, and `cpp/tests/test_photon_input.cpp`
check source geometry, keyed reproducibility, and record preservation.

**Visual check.** Save resolved source samples with `--photon-output` when
reproducibility needs the exact entrance rays. Overlay the resulting arrival
records on the model to inspect their trajectories:

```sh
obdeect-plot-reference --view telescope --input arrivals.csv \
  --optical-model-json telescope.optical-model.json --path-selection stratified \
  --max-paths 200 --output paths.png
```

## 2. Entrance transmission and incoming obscuration

For reflector models, an optional direction-dependent telescope transmission
multiplies the initial throughput. The incident ray is then tested against
incoming cylinders, opaque surfaces, and configured shadow planes. The nearest
blocking surface before the primary mirror terminates the photon as
`blocked_obscurer`. A ray that reaches no primary aperture is `missed_primary`.
Obscurers are absorbers in this transport; they do not produce reflected rays.

**Verification.** `cpp/tests/test_core.cpp` exercises support and camera
obscuration; `cpp/tests/test_optical_model.cpp` checks nearest-hit ordering and
loss attribution.

**Visual check.** Inspect blocked paths and their recorded component IDs on the
telescope plate, or map blocked input weight by entrance position:

```sh
obdeect-plot-reference --view loss-map --input arrivals.csv \
  --optical-model-json telescope.optical-model.json --status blocked_obscurer \
  --output incoming-losses.png
```

## 3. Primary-mirror intersection and reflection

The tracer finds the nearest valid hit on a finite primary facet or the
axisymmetric primary profile. Aperture shape, gaps, holes, surface sag, and the
local normal determine whether and where the ray intersects. Reflection is
specular about that normal, with configured stochastic mirror scatter applied
to the outgoing direction. Wavelength- and incidence-dependent reflectivity
and any spatial degradation map multiply throughput at the hit. The incidence
angle is measured from the local normal.

**Verification.** `cpp/tests/test_facets.cpp` and
`cpp/tests/test_segmented_optical_model.cpp` check aperture boundaries, hit
selection, and reflection; `cpp/tests/test_optical_response.cpp` checks response
interpolation. Scatter and degradation import have focused tests in
`python/tests/test_optical_model_compiler.py`.

**Visual check.** Compare the compiled pupil with the sampled entrance points,
then inspect recorded primary interactions and the focal-plane distribution:

```sh
obdeect-plot-reference --view pupil --optical-model-json telescope.optical-model.json \
  --output pupil.png
obdeect-plot-reference --view telescope --input arrivals.csv \
  --optical-model-json telescope.optical-model.json --path-selection stratified \
  --max-paths 200 --output primary-paths.png
obdeect-plot-reference --view telescope --panel focal-plane-hits --input arrivals.csv \
  --optical-model-json telescope.optical-model.json --output focal-hits.png
```

## 4. Secondary mirror and inter-mirror transport

For a dual-reflector model, the reflected ray propagates to the secondary
A configured incoming secondary shadow is tested on the incident ray before the primary. After primary reflection, inter-mirror obscurers are tested along the segment to the secondary; the nearest obstruction terminates the path.
miss is recorded separately. At a secondary hit, its local normal sets the
reflection, configured scatter perturbs the outgoing direction, and the
wavelength/incidence reflectivity and spatial degradation multiply throughput.
Single-mirror models skip this stage.

**Verification.** `cpp/tests/test_axisymmetric_optics.cpp` checks analytic
primary/secondary intersections and aperture failures.
`cpp/tests/test_optical_model.cpp` checks finite masks and ordered transport;
`python/tests/test_optical_model_compiler.py` checks imported dual-reflector
geometry, scatter, and obscuration.

**Visual check.** Use the telescope assembly or axial-section panel with
recorded paths. Confirm that rays intercept the declared secondary surface,
that shadow losses stop on the obstruction, and that surviving rays continue
toward the focal surface. A `missed_secondary` status identifies rays that
leave the primary without a secondary hit.

```sh
obdeect-plot-reference --view telescope --panel section --section-plane xz \
  --input arrivals.csv --optical-model-json telescope.optical-model.json \
  --path-selection stratified --max-paths 200 --output dual-section.png
```

## 5. Focal surface, pixels, and camera response

The final ray is intersected with the model's focal surface or finite detector
surfaces. Pixel assignment follows the compiled aperture geometry and declared
assignment policy. A detector hit is recorded even when its optical throughput
has fallen to zero. At that hit, configured camera degradation, pixel
acceptance/response, and camera angular response multiply throughput. The
`--focal-surface-image` diagnostic instead stops at the continuous focal
prescription before pixel acceptance and camera response; it is not the normal
detector result.

**Verification.** `cpp/tests/test_detector_planes.cpp` and
`cpp/tests/test_segmented_optical_model.cpp` check finite detector apertures and
assignment. `python/tests/test_camera_surfaces.py` checks pixel placement, and
`python/tests/test_bulk_tracing.py` checks measured camera response and keyed
scatter. `python/tests/test_result_contract.py` checks the distinction between
zero response and geometric loss.

**Visual check.** Draw the compiled focal geometry and compare it with weighted
detected-hit density using the `focal-plane-hits` panel. Check detector
boundaries and enabled pixel footprints in the focal-geometry panel. The
continuous-focal diagnostic should be compared separately from a normal trace.

```sh
obdeect-plot-reference --view telescope --panel focal-geometry \
  --optical-model-json telescope.optical-model.json --output focal-geometry.png
obdeect-plot-reference --view telescope --panel focal-plane-hits --input arrivals.csv \
  --optical-model-json telescope.optical-model.json --output focal-hits.png
```

## 6. Refractive surfaces, timing, and terminal accounting

The non-sequential tracer advances to the nearest declared surface. In a
material it accumulates geometric path length, phase optical path, group delay,
and bulk absorption. At a refractive interface it computes transmitted or
reflected directions from the media on either side. Except under total internal
reflection, transport follows the Snell-transmitted direction. With
`additional_coating` semantics, throughput includes the unpolarized Fresnel
power transmission (the mean of the s- and p-polarized terms) and the configured
coating response; `complete_interface` response tables replace the Fresnel
factor. Polarization is not propagated and Fresnel reflection is not sampled as
a stochastic branch. Total internal reflection follows the reflected direction
and retains the photon in its incident medium. Mirror reflectivity and opaque
obscurers act at their corresponding surfaces. The declared interaction limit
prevents an unbounded path.

For reflector-only traces, arrival time is emission time plus geometric path
length times the configured propagation group index divided by *c*. For
material transport, arrival time uses the accumulated group delay. The output
retains response loss separately from terminal geometric loss, together with
path vertices, surface IDs, and the terminal status. `missed_primary`,
`missed_primary`, `missed_secondary`, `missed_screen`, `no_detector`, and `blocked_obscurer` identify geometric termination; non-sequential paths can also escape the model or material, be absorbed, or reach the interaction limit.
remains a detection with zero response, not a geometric loss.

**Verification.** `cpp/tests/test_refractive_transport.cpp` checks interface
directions, transmission, and material timing; `cpp/tests/test_optical_model.cpp`
checks interaction limits and path accounting. Python result-contract tests
check that input weight closes into surviving throughput, response loss, and
terminal loss.

**Visual check.** Use `--view rays` to inspect recorded vertices and colour
paths by the stored arrival time, wavelength, or incidence angle. The weighted
distribution command is limited to focal-plane position and vacuum-geometric
arrival-time summaries; it is not a material group-delay diagnostic. For
modelled refractive stacks, check that paths cross the declared interfaces in
order and compare the output's `optical_path_m` and `arrival_time_ns` with
analytic layer-by-layer expectations.

```sh
obdeect-plot-reference --view rays --input arrivals.csv --colour-by arrival-time \
  --path-selection stratified --max-paths 300 --output timed-rays.png
obdeect-photon-distributions --input arrivals.csv --output detected-distributions.svg
```

## What a plot establishes

A plot checks spatial consistency and makes missed or blocked paths visible; it
does not establish response normalization or equivalence to another simulator.
Use analytic tests for individual surfaces and response kernels, then compare
the same photon inputs and compiled model against the reference implementation
under predeclared tolerances. See [PRODUCTION_VALIDATION.md](PRODUCTION_VALIDATION.md)
and [REFERENCE_RUN.md](REFERENCE_RUN.md) for that numerical evidence.
