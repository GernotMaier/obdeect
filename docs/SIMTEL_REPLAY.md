# Replay resolved photons through sim_telarray

`tools/build_simtel_replay.py` builds a diagnostic executable against a supplied
sim_telarray source checkout and its matching object files. It copies the main
and imaging translation units into a separate build directory. The simulator
checkout is read-only. Source anchors must match the optical-backend boundary;
unknown source layouts fail explicitly.

Supply a JSON build recipe with four lists: `compile_c`, `compile_cxx`, `link`,
and `inputs`. Copy the exact compiler flags, include paths, linked objects and
libraries from that simulator build. Include `-DDEBUG_TRACE_99` when compiling
the copied C units. Use these argv placeholders:

| Placeholder | Meaning |
| --- | --- |
| `{source}`, `{object}` | Translation unit and output object for each compilation |
| `{include}` | obdeect's generic C++ headers, for the C++ replay driver |
| `{main_object}`, `{imaging_object}`, `{driver}` | Replacement objects in the final link |
| `{executable}` | Output executable |

Use a C++20 compiler for `compile_cxx` and the final link. Do not link the
original `sim_telarray.o` or `sim_imaging.o` alongside their replacements.
`inputs` must list all remaining objects, libraries, configuration headers and
toolchain inputs needed to reproduce the build. The builder retains commands,
input/executable SHA-256 values and compiler output in `build-record.json` and
`build.log`. Compilation failures retain the log.

```sh
python tools/build_simtel_replay.py --simtel-root /path/to/sim_telarray \
  --output-directory build/simtel-reference --build-recipe build-recipe.json
```

Freeze the resolved source before comparing tracers. Native generation can
export the complete pre-trace input, including undetected photons:

```sh
obdeect-simtools-raytrace --optical-model telescope.json --source star \
  --star-mode finite --distance-m 10000 --entrance-z-m 50 --ray-tracing-seed 21 \
  --photons 10000 --photon-output common-photons.csv --output candidate.csv
```

The optical model stores `random_seeds.detector_configuration_seed`. Set it with
`obdeect-compile-optical-model --detector-configuration-seed N`. This seed is
used once while compiling random panel positions, focal lengths, and alignment;
the resulting geometry is saved in the optical model and stays fixed for every
trace. `--ray-tracing-seed` is a separate run-level seed used for source
sampling and per-photon mirror-surface scatter. The arrival CSV records both
`detector_configuration_seed` and `ray_tracing_seed` on every row. Reusing the
same optical model, run seed, and photon IDs reproduces the same randomized
trace; changing the run seed does not alter the compiled mirror geometry.
Both command-line seeds default to zero, so production runs should set and
retain their chosen values explicitly. `--source-seed` remains a deprecated
compatibility alias for `--ray-tracing-seed`; new commands should use the
ray-tracing name because the value also controls per-photon optical scatter.

| Effect | Seed | When it is used | What changing it changes |
| --- | --- | --- | --- |
| Random panel positions, focal lengths, and panel alignment | `detector_configuration_seed` (`--detector-configuration-seed`) | Once, during optical-model compilation | The compiled mirror geometry; the result is stored in the optical model. |
| Randomized source sampling | `ray_tracing_seed` (`--ray-tracing-seed`) | While generating source photons for a run | The sampled source-photon inputs. Preloaded photon inputs are not resampled. |
| Per-photon mirror-surface scatter | `ray_tracing_seed` (`--ray-tracing-seed`) | During tracing, at each mirror encounter when scatter is configured | The stochastic surface deflections; it does not change the compiled geometry. |

The configuration seed fixes telescope geometry; the run seed controls random
draws made for that trace. To reproduce a run, retain the optical model, run
seed, and photon identities (or the exported photon inputs).

Run the reference from its normal configuration/data working directory. Supply
the simulator preprocessor's absolute command because the replay executable
lives outside the simulator installation:

```sh
SIMTEL_CONFIG_PREPROCESSOR='/path/to/sim_telarray/bin/pfp -v -I.' \
OBDEECT_REPLAY_INPUT=/absolute/common-photons.csv \
OBDEECT_REPLAY_OUTPUT=/absolute/reference.csv \
OBDEECT_REPLAY_TELESCOPE_INDEX=0 \
/absolute/build/simtel-reference/simtel-replay --quiet -c telescope.cfg /dev/null
```

The explicit telescope index selects the configured optics; input telescope
identities remain unchanged. One replay contains one telescope identity but can
span events, array reuses and bunches. Duplicate full identities, unresolved
wavelengths and output/input aliases are rejected. Input positions use the
same telescope-local metres as native replay. The boundary converts to
centimetres and applies the inverse of sim_telarray's actual ground-to-telescope
transform. Finite emission distance is preserved; unavailable distance uses a
declared effectively infinite `1e30 cm` bound. No source photons are resampled.

The adapter calls the actual legacy optical backend. It records camera-frame
hits, direction, measured travel/arrival time, relative efficiency, mirror
index, actual telescope-frame primary/secondary reflection points and incidence
cosines. Lost photons retain their input identities and actual imaging-source
loss branch, with unavailable terminal measurements left empty. Resetting
reflection diagnostics before each photon prevents stale data from appearing
as a new interaction.

This raw table is **not production-normalized**. The legacy handoff does not
expose classified terminal surfaces/positions for losses; its relative efficiency
excludes nominal spectral response applied by the caller, and pixel/concentrator
transport lies downstream. Do not relabel raw camera hits as qualified detector
arrivals or derive missing lost-ray endpoints. Complete those reference boundary
measurements before using [PRODUCTION_VALIDATION.md](PRODUCTION_VALIDATION.md)
or publishing the speedup measured by [BENCHMARK.md](BENCHMARK.md).

[The local four-class record](review/simtel-replay.json) retains hashes, exact
commands and terminal counts for 1,000 common photons per class, each replayed
twice with byte-identical seeded output. It exercises real segmented,
paraboloid, Fresnel and secondary optics using the simulator's generic regression
configuration. It does not qualify CTAO production variants.
