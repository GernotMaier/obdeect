# Telescope plotting implementation plan

## Purpose and decision

Implement a **model-faithful telescope diagnostic plate** in
`python/obdeect/plotting.py`.  Its primary output is a static, headless PNG or
SVG for a compiled optical model, with optional overlays from one trace CSV.
It must answer, at a glance:

1. What physically defined surfaces, panels, obscurers, and detector geometry
   did this trace use?
2. How are those components arranged in the telescope frame?
3. Where did the selected photons interact or terminate, and what hardware
   caused recorded losses?

This is not an illustrative telescope drawing.  A component is rendered only
when it is represented in the selected compiled model.  Missing mechanical
geometry is declared in the figure footer and is never filled in with a
generic LST/MST/SST/SCT outline.  That rule is already established by
`docs/USE_CASES.md`, `TODO.md`, and the importer's provenance-first design.

The target is more useful than the current independent
`compiled-structure`, `compiled-mirror`, and `compiled-3d` images: one
consistent telescope plate, with drill-down views for the pupil and the
focal plane.  Retain those existing views and their command-line behaviour as
compatibility aliases until the new tests and documentation are in place.

## Design evidence from `../literature`

The following visual language should be adopted, not copied mechanically.

| Reference | Visual lesson to apply |
| --- | --- |
| `app_82_36_2016.pdf`, pp. 2-3, Figs. 1-3 | A perspective computational model and a labelled side view make the primary, M2, camera, support masts, trusses, and focal surface legible.  The labels identify physical roles, not a guessed telescope family. |
| `app_28_10_2007.pdf`, p. 6, Fig. 2; p. 14, Figs. 8-9 | A clean axial cross-section establishes P1/M1, M2, focal surface, aperture, and optical distances.  Ray/spot plots are separate diagnostic panels, rather than decoration over a CAD view. |
| `SFegan_optics_2006.pdf`, p. 21, Figs. 14-15 | The face-on pupil must draw each finite mirror polygon exactly; segmentation, gaps, central holes, and panel identity are scientifically visible.  Do not replace panels with centre dots. |
| `2309.09560.pdf`, pp. 3-6, Figs. 1-6 | A full 3-D assembly and a simplified optical model must be visibly distinguishable.  Full structure produces azimuthally asymmetric shadowing, so visualisation needs an explicit aperture/loss map rather than an axisymmetric surrogate. |
| `Appendix-LST.pdf`, p. 6; `Appendix-MST-Structure.pdf`, pp. 6 and 12-13; `Appendix-MST-SC.pdf`, p. 6; `Appendix-SST-Structure.pdf`, pp. 11 and 16-17 | Plot model-derived panel placement, focal-plane geometry, PSF, and transmission as distinct observables.  Use metre coordinates and title/annotation provenance, not CTAO-specific hard-coded dimensions. |

The visual hierarchy therefore is: exact optics first, model-defined
obscuration second, optional path/loss evidence third, and metrics/metadata
last.  Decorative ground, tower, mount, or CAD surfaces are out of scope
unless their finite geometry becomes a compiled input.

## Scope, non-goals, and invariants

### In scope

- Static Matplotlib/Agg plots, using the existing optional dependency and
  headless execution path.
- Segmented models: exact primary facet polygons, detector polygons, and
  finite cylinder obscurers already consumed by the tracer.
- Axisymmetric models: model-derived primary, secondary, and focal-surface
  contours sampled from their compiled polynomial coefficients.  They are
  labelled as rotationally symmetric surfaces, not as segmented hardware.
- Orthographic 3-D, axial x-z/y-z projections, a face-on pupil, focal-plane
  geometry, and optional recorded trace paths/losses.
- A deterministic SVG summary for focal-plane data when the output suffix is
  `.svg`; retain the current dependency-free focal PSF path.

### Not in scope for the first implementation

- Inventing camera bodies, M2 supports, baffles, windows, light guides,
  quadrilateral obscurers, towers, or counterweights from a telescope name or
  nominal prescription.
- An interactive viewer, CAD/mesh importer, PyVista/VTK dependency, or WebGL
  export.
- Re-tracing photons inside plotting, deriving a component hit from a line
  intersection, or presenting a physically unvalidated `reference-mst` sketch
  as a compiled CTAO telescope.
- Changing C++ trace decisions solely to support cosmetics.  A new recorded
  interaction field is justified only when it records an actual already-made
  trace decision.

### Invariants

- All coordinates are telescope-frame metres.  The optical-axis direction and
  origin are named in every figure; do not assume a global astronomical frame.
- Equal data ranges and equal aspect ratios are used within a projection; the
  camera is never visually enlarged to make it readable.
- Geometry comes from one immutable `obdeect.compiled-optical-model.v1`
  document.  The model SHA-256 and model/version provenance appear in the
  footer.
- The figure is deterministic for identical model, CSV, selection options,
  Matplotlib version, and output size.  No global random sampling.
- A rendered path is a recorded path vertex sequence from the CSV.  It is not
  an interpolated optical solution.

## Intended output

### Default telescope plate

Add `--view telescope` as the primary compiled-model view.  With a compiled
model and no CSV, it produces this 2 x 2 layout:

```text
+---------------------------+---------------------------+
| A. orthographic assembly  | B. axial optical section  |
| exact panels, detector,   | x-z (or chosen plane),    |
| supplied obscurers        | surfaces and scale bar    |
+---------------------------+---------------------------+
| C. entrance pupil         | D. focal surface / camera |
| exact facet polygons,     | detector outline/pixels;  |
| holes and obscurer shadow | optical context only      |
+---------------------------+---------------------------+
footer: frame, units, model + SHA prefix, geometry coverage, availability
```

- **A, assembly:** Orthographic, azimuth -55 degrees, elevation 24 degrees,
  with `set_box_aspect` from the true global bounds.  Primary facets have a
  muted blue-grey face, thin charcoal edges, and no height colour bar by
  default.  M2 (when compiled) is a darker blue-grey.  Detector/focal surface
  is amber; opaque structures are muted red-brown at 60 percent opacity;
  selected paths use high-contrast status colours.  The legend contains roles
  and counts, never an unqualified component name.
- **B, axial section:** Show x-z by default and offer `--section-plane xz|yz`.
  Draw surface cross-sections rather than centre scatter.  Add a thin dashed
  optical axis at x=0.  Include compact labels `M1`, `M2`, `detector`, and
  `opaque cylinder #N` only for actual components.  If a section cuts a
  cylinder, draw the exact projected capsule; if it is out of plane, do not
  pretend it is in the cut.
- **C, pupil:** View towards the incoming aperture, with x right and y up.
  Use the exact `compiled_facet_polygons` shapes and true polygon orientation.
  Draw the detector/M2 aperture projected into the pupil only if the compiled
  geometry can supply it; otherwise omit it.  For cylinder obscurers, show
  their projected physical footprint only when a finite projection is
  well-defined; otherwise list them in the coverage footer.  Provide a
  separate loss overlay described below instead of manufacturing a shadow.
- **D, focal geometry:** Plot all actual detector surfaces and, when present,
  the parsed camera-pixel polygons.  The default is a clean geometry view.  A
  trace CSV replaces it with the weighted focal-plane PSF only with the
  explicit `--panel focal-plane-hits`; never silently mix geometry and density
  colour scales.

Use a white background, dark grey axes/labels, restrained semantic colours,
and readable 9-11 pt type.  Do not use a rainbow colour map for component
identity.  A colour bar is only allowed for a numeric quantity (for example,
facet z, wavelength, time, incidence angle, or weighted density) and must
state its units.

### Drill-down views

Keep current choices but normalise their responsibilities:

| CLI view | Deliverable |
| --- | --- |
| `telescope` | The four-panel diagnostic plate above. |
| `pupil` (successor to `compiled-mirror`) | Full-size exact primary layout, panel ID tooltip-free labels only when `--label-panels` is given, optional scalar colouring. |
| `assembly-3d` (successor to `compiled-3d`) | Full-size orthographic geometry, no synthetic elements. |
| `section` | Full-size x-z or y-z optical section for paper-style schematic comparison. |
| `rays` | Existing two axial projections, upgraded to use compiled geometry when supplied. |
| `loss-map` | Entrance-pupil positions of actual terminal outcomes from a trace CSV. |
| `focal-plane` | Existing weighted hit density/projections and PSF statistics. |

Deprecation can be documented, but do not remove `structure`,
`compiled-structure`, `compiled-mirror`, or `compiled-3d` in this work.
Implement each legacy name as an explicit dispatch to the corresponding new
renderer rather than maintaining two drawing code paths.

## Data contract and required compiler work

### 1. Introduce a renderer-facing optical-model adapter

Create a small, pure-Python `PlotScene` dataclass and loader in
`python/obdeect/plotting.py` (or `plot_scene.py` if the module becomes too
large).  It must be built strictly from the compiled JSON and validate before
rendering.  Suggested records:

```python
@dataclass(frozen=True)
class PlotSurface:
    id: int
    role: Literal["primary", "secondary", "detector"]
    kind: Literal["facet", "axisymmetric"]
    vertices_m: tuple[tuple[float, float, float], ...] | None
    axial_profile_m: tuple[tuple[float, float], ...] | None

@dataclass(frozen=True)
class PlotObscurer:
    id: int
    first_endpoint_m: tuple[float, float, float]
    second_endpoint_m: tuple[float, float, float]
    diameter_m: float

@dataclass(frozen=True)
class PlotScene:
    provenance: Mapping[str, object]
    model_sha256: str
    kind: Literal["segmented", "axisymmetric"]
    surfaces: tuple[PlotSurface, ...]
    obscurers: tuple[PlotObscurer, ...]
    camera_pixels: tuple[PlotPolygon, ...]
    unavailable_roles: tuple[str, ...]
```

The loader owns all JSON-schema checks, finite-value checks, normal/tangent
normalisation, and component-bound computation.  Render functions receive
only `PlotScene`, never raw dictionaries.  Reuse `_vector3`, `_normalised`,
`_cross`, `_aperture_local_vertices`, `compiled_facet_polygons`, and
`_detector_polygons` after tightening their return types, so every polygon has
one source of truth.

### 2. Extend compiled JSON only where actual information exists

`trace_model` already contains enough geometry for segmented M1 facets,
detector surfaces, and cylinder obscurers.  Promote their renderer roles in a
backward-compatible optional `plot_geometry` section rather than relying on
undocumented field names.  It is a denormalised rendering index, not a second
physical model:

```json
"plot_geometry": {
  "schema_version": 1,
  "frame": {"origin": "primary_vertex", "axes": "+x,+y,+z", "unit": "m"},
  "components": [
    {"id": 0, "role": "primary", "source": "trace_model.primary_facets"},
    {"id": 901, "role": "detector", "source": "trace_model.detector_surfaces"},
    {"id": 902, "role": "opaque_cylinder", "source": "trace_model.cylinder_obscurers"}
  ],
  "unavailable_roles": ["camera_body", "window", "m2_support", "baffle"]
}
```

Rules for the compiler:

- Put a role in `components` only if the complete finite geometry consumed by
  the tracer is available.
- Add an unavailable role only from the compiler's existing deferred or
  `trace_blockers` evidence.  Do not list a component merely because typical
  telescopes have one.
- Preserve v1 input compatibility: a plot loader falls back to established
  `trace_model` paths and reports `plot_geometry unavailable (legacy model)`.
- For axisymmetric models, provide component roles for M1, M2, and focal
  surface and reuse the exact compiled `vertex_z_m`, radii, scale, and
  coefficient arrays.  Do not tessellate a new approximate shape in the
  compiler; the renderer samples the defined profile at a documented fixed
  radial count (e.g. 256).
- Do not claim camera housing/window/pixel footprints are 3-D surfaces from
  the existing 2-D camera layout.  Pixel centres are eligible only for panel D
  and only when their finite footprint type is parsed.

Future mechanical import can add `finite_mesh`, `cylinder`, `frustum`, and
`polygonal_prism` component records to this section.  It must include a
component ID, role, closed finite geometry, material/opacity semantics, and
provenance.  The generic renderer then adds support/mast/baffle views without
telescope-family branches.

### 3. Expand trace diagnostics before adding component-loss plots

The present CSV has ragged vertices and terminal status but not a terminal
component ID or an explicit interaction sequence.  Add an optional,
versioned diagnostic extension (for example, `interaction_surface_ids` and
`terminal_surface_id`) only after the tracer records the actual winning
intersection.  A loss map may group legacy CSVs by terminal `status`, but it
must state `component unknown` rather than infer a mast/camera from geometry.

For new CSVs, record per photon:

- `entrance_x_m`, `entrance_y_m` (or define point 0 as the entrance point in
  the contract, and validate that interpretation),
- an ordered interaction surface-ID list aligned with the stored vertices,
- `terminal_surface_id` when a finite compiled component caused termination,
- existing source, wavelength, time, weight, and terminal status.

This gives the loss-map panel a factual basis and enables deterministic filters
such as `--status blocked_obscurer` and `--component-id 902`.

## Rendering implementation

### Geometry builders

1. **Segmented facets:** retain every finite aperture polygon exactly.  For
   curved panels, the view may use boundary vertices in the panel tangent
   plane for a face-on outline; in 3-D, sample the already compiled spherical
   cap only if all necessary curvature data are present.  Otherwise show the
   true planar aperture boundary and annotate `surface curvature not rendered`.
2. **Detector surfaces:** use the same aperture factory and normal/tangent
   frame as facets.  Do not use camera-pixel-centre extent as a replacement
   for a detector polygon where `trace_model.detector_surfaces` exists.
3. **Axisymmetric surfaces:** generate a radial profile from the coefficients
   and revolve it only for the 3-D artist.  The section uses the exact sampled
   profile; the pupil uses concentric finite annuli.  State sample count in
   code/comments and test it independently of visual pixels.
4. **Cylinders:** render a finite capsule in 2-D projections and a 12-sided
   or configurable low-poly cylinder in 3-D.  End caps are shown; an infinite
   line is forbidden.  Use the same endpoints and diameter as the trace model.
5. **Bounds:** calculate the union of all finite vertices/endpoints plus radii.
   Add a fixed 5 percent or 0.25 m minimum margin.  Set all three 3-D limits
   from the largest span so aspect is physically honest.  Raise a clear error
   for an empty or non-finite optical model.

### Paths, selection, and status semantics

- `--input TRACE.csv` activates an overlay; no input means geometry-only.
- `--max-paths` is a deterministic first-N compatibility default.  Add
  `--path-selection first|stratified` and `--selection-seed`; `stratified`
  reserves slots by terminal status and source kind using a deterministic hash
  of `photon_id`, never global RNG state.
- Draw status with a colour-blind-safe palette and a dedicated legend:
  detected blue; blocked/opaque red; missed-primary grey; missed-detector
  purple; invalid input black.  Source kind may control line style only when
  status colour is already needed, preventing a misleading double colour
  encoding.
- Put a marker at each recorded interaction, optionally labelled with a
  verified surface ID for up to a small cap.  A terminal marker is distinct
  from an intermediate reflection.
- Clip long source-flight segments at an explicitly rendered `--context-radius-m`
  telescope box and add `incoming path clipped` to the legend.  This retains
  the existing protection against distant sources hiding the telescope.

### Pupil loss map

Implement `--view loss-map --input TRACE.csv` after the diagnostic extension.
It has two panels:

- left: physical pupil geometry plus entrance points, coloured by terminal
  status/component;
- right: hexbin or fixed 2-D weighted loss fraction, with a labelled
  `lost input weight / incident-bin weight` scale and bins with no incident
  photons masked rather than zero-filled.

Use the raw entrance points and source weights.  Never estimate a shadow map
from outlines or use detected-only density as an effective-area plot.  Plot
azimuthal cuts only when a scan file explicitly supplies field coordinates;
an individual trace does not imply a field-angle scan.

### Titles, legends, and footer

Use this title format:

```text
<model>/<version> compiled telescope geometry [segmented|axisymmetric]
```

Use a single concise footer, for example:

```text
Telescope frame (x, y, z), metres | model SHA256: 0123456789ab |
rendered: 198 M1 facets, 1 detector, 4 opaque cylinders |
not modelled: camera body, M2 support
```

The `not modelled` clause appears only when the compiler can establish it;
otherwise use `mechanical geometry coverage not declared`.  A single large
annotation replaces a panel only when nothing renderable is present.

## CLI and API changes

Extend `main()` with these validated options:

```text
--view telescope|pupil|assembly-3d|section|loss-map|rays|focal-plane|...
--optical-model-json PATH                     # required for all geometry views
--input PATH                                  # required for paths, loss map, focal plane
--section-plane xz|yz                         # default xz
--panel assembly|section|pupil|focal-geometry|focal-plane-hits
--show-paths / --no-show-paths                # telescope default: show if CSV supplied
--path-selection first|stratified
--selection-seed INTEGER
--component-id INTEGER                        # filter validated recorded component IDs
--status NAME [NAME ...]
--colour-by status|wavelength|arrival-time|incidence-angle|facet-z
--label-panels                                # pupil only; disabled over a safe panel-count limit
--context-radius-m POSITIVE
```

Validation rules:

- reject an input-dependent option without `--input`;
- reject a view/model-kind pair that cannot be represented (for example,
  `--label-panels` on an axisymmetric model);
- reject `--component-id` when the CSV has no recorded component identifiers;
- reject an unsupported output suffix; document PNG, PDF, and SVG semantics;
- preserve the existing `--focal-plane` shorthand as a deprecated alias only
  when it is unambiguous;
- no CLI option changes trace coordinates, units, or the input model.

Keep `TELESCOPE_NAMES` only for the legacy analytic outline mode.  Compiled
views title themselves from JSON provenance; a user-supplied telescope name
cannot relabel or override a compiled optical model.

## Implementation sequence for Luna

### Phase 0 - establish fixtures and compatibility tests

1. Do not touch the user-owned modified files in the working tree without
   reconciling them.  Start by recording their diff and isolate this work.
2. Add tiny deterministic JSON fixtures under `python/tests/fixtures/`:
   - segmented single facet + detector;
   - segmented hexagonal mini-dish with a finite cylinder obscurer;
   - axisymmetric M1/M2/focal model with an annular central hole;
   - legacy v1 model lacking `plot_geometry`;
   - trace CSV with detected, missed, and blocked outcomes plus known entrance
     positions/component IDs for the new contract.
3. Preserve `test_plot_reference_mst.py` coverage while extracting geometry
   helpers into unit-testable functions.

### Phase 1 - pure optical-model loader and exact geometry

1. Implement `load_plot_scene()` and component dataclasses with strict error
   messages that name the JSON path and invalid field.
2. Refactor current facet/detector polygon helpers into the loader and add
   axisymmetric profile sampling with unit tests for vertex/radius/coefficient
   transforms.
3. Implement shared bounds, projection, semantic style, legend, provenance
   footer, and output helper.  Ensure no plot function reads raw JSON after
   the loader returns.
4. Add `plot_geometry` emission to `optical_model_compiler.py` only after the
   loader can consume both legacy and new models.

### Phase 2 - geometry views

1. Implement `pupil`, then `section`, then `assembly-3d`; validate each with
   fixture geometry before composing the telescope plate.
2. Make `compiled-mirror` and `compiled-3d` dispatch to the new renderers.
   Make `compiled-structure` dispatch to `section` until a purpose-specific
   compatibility layout is needed.
3. Compose `telescope` with `GridSpec`, deliberately sharing only compatible
   limits; the pupil and focal plane retain their own equal-scale bounds.
4. Test PNG and SVG/PDF export from an Agg backend, including a model with no
   camera pixels and a model with no renderable mechanical geometry.

### Phase 3 - trace-aware plots

1. Make rays render against `PlotScene` whenever `--optical-model-json` is
   supplied, while retaining the current analytic reference outline for its
   explicit legacy mode.
2. Implement deterministic path selection and status/source legend logic.
3. Extend the result contract and C++/CSV writer for actual interaction and
   terminal component IDs.  Version the contract, provide a migration reader,
   and update `analysis.py` without changing mathematical PSF results.
4. Implement `loss-map` only after those fields are genuinely emitted.
5. Add the `focal-plane-hits` panel mode and keep its density colour bar and
   optical metrics separate from geometry colours.

### Phase 4 - documentation and validation

1. Update `docs/USE_CASES.md` with one command and one expected artefact per
   renderer.  Do not modify `README.md` (per `AGENTS.md`).
2. State geometry coverage and model limits prominently in the docs and every
   screenshot caption.
3. Render documented fixture outputs in CI or a dedicated test target.
4. Run the focused tests first, then the repository checks required by
   `AGENTS.md`.

## Test matrix and acceptance criteria

| ID | Level | Assertion |
| --- | --- | --- |
| T-PLOT-001 | unit | A square, circle, and both hex orientation polygon have correct vertices in their supplied tangent frame. |
| T-PLOT-002 | unit | Optical-model loader rejects non-finite coordinates, duplicate IDs, zero cylinders, invalid normals/tangents, and unsupported schema versions. |
| T-PLOT-003 | unit | Axisymmetric profile sampling preserves vertex z, inner/outer radii, and polynomial endpoint values. |
| T-PLOT-004 | unit | Bounds include panel vertices, detector vertices, and the complete cylinder radius; all 3-D axes use equal span. |
| T-PLOT-005 | unit | Legacy compiled JSON loads with an explicit coverage warning and no invented components. |
| T-PLOT-006 | integration | The segmented pupil SVG/PNG contains the known exact number of polygon paths, the central gap, and correct x/y units. |
| T-PLOT-007 | integration | The axisymmetric section contains M1, M2, and focal profiles derived from its model, not an LST/MST reference curve. |
| T-PLOT-008 | integration | Telescope plate has four named panels, one provenance footer, a physical equal-aspect pupil, and no `Unvalidated` analytic outline. |
| T-PLOT-009 | integration | A cylinder obscurer appears at correct projected endpoints/diameter in 2-D and 3-D. |
| T-PLOT-010 | integration | Identical CSV/model/seed produce byte-stable SVG and stable raster pixels within a documented tolerance. |
| T-PLOT-011 | integration | `--status` and `--component-id` filter only recorded paths and reject unavailable terminal component data. |
| T-PLOT-012 | contract | A terminal component ID equals the trace kernel's selected first blocker; paths and interaction lists agree in length/order. |
| T-PLOT-013 | scientific | Weighted loss-map numerator/denominator closes to the trace summary per status, including zero-throughput detected photons. |
| T-PLOT-014 | regression | Existing focal-plane centroid, D80, throughput, and zero-throughput tests still pass unchanged. |
| T-PLOT-015 | visual QA | Render the mini-dish, segmented, and axisymmetric fixtures; inspect panel spacing, labels, clipping, legends, colour bars, and small-component legibility at 1600 px width. |

Acceptance requires all of the following:

- a model with only M1 and detector visibly says that other mechanical
  geometry is undeclared or unavailable, rather than showing representative
  supports;
- the pupil uses panel polygons, not panel-centre scatter;
- a dual-mirror fixture visibly contains distinct M1/M2/focal surfaces in the
  section and assembly;
- every path overlay is trace-recorded and every component-specific loss is
  backed by a recorded component ID;
- no existing `focal-plane` analysis result changes merely because plotting
  was refactored;
- all modified Python code satisfies the Python 3.14 and Ruff requirements in
  `AGENTS.md`.

## Suggested first implementation command

After Phase 1 fixtures exist, the first reviewable image should be produced
without any telescope-name fallback:

```sh
PYTHONPATH=python python -m obdeect.plotting \
  --view telescope \
  --optical-model-json python/tests/fixtures/segmented-mini.optical-model.json \
  --output /tmp/segmented-mini-telescope.png
```

The first trace-overlay review adds the same model and a fixture CSV with an
explicit path selection seed.  Do not use a production CTAO telescope image
as the first golden baseline: a small fixture is the only practical way to
verify exact geometry and status semantics.
