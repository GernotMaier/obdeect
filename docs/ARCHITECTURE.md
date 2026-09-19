# Core architecture status

The project follows the `Plan_v0.md` target layout. The C++20 core is split by
responsibility so a photon kernel never needs Python, YAML, EventIO or a
plotting dependency.

| Component | Current implementation |
| --- | --- |
| `math.hpp` | `Vec3`, dot/cross/norm and checked direction normalisation. |
| `photon_buffer.hpp`, `abi.hpp` | SoA input contract, result buffers, status values and validation. |
| `source.hpp` | Deterministic 400-nm star, point illuminator and laser sources. |
| `tables.hpp` | Immutable no-extrapolation 1-D response interpolation. |
| `intersections.hpp`, `geometry.hpp` | Plane, sphere/cap, disk, finite cylinder and axisymmetric-surface dispatch. |
| `interactions.hpp` | Checked specular reflection. |
| `axisymmetric_optics.hpp` | Paraboloid/even-polynomial SC surfaces and forward Newton intersection. |
| `facets.hpp` | Finite circular facet intersection, nearest-hit selection and tabulated coating response. |
| `scene.hpp`, `trace.hpp` | Immutable directed toy-scene compilation and scalar SoA block tracing. |
| `diagnostics.hpp` | Status/weight closure summary. |
| `model_import.hpp` | Canonical CTAO model provenance/import target. |

The existing `toy_mst.hpp` is a minimal concrete scene used by the executable.
It remains intentionally separate from LST/MST/SST/SCT production model
import. A production scene requires the model’s segment list, transforms,
apertures, structure, material tables and reference fixture before it can be
declared compatible.

## Deliberately deferred interfaces

`cpp/bindings/`, `python/obdeect/api.py`, CORSIKA7 EventIO input and CTAO JSON
parsing are integration layers, not hot-path physics. They remain deferred
until the C++ scene contract stabilises. Their absence is explicit rather than
silently offering a Python API that cannot execute the native trace.
