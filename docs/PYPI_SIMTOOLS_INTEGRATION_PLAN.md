# PyPI distribution and simtools integration plan

## Goal and container impact

Publish `obdeect-dev` as a pip-installable distribution whose supported wheels
contain the compiled C++ ray-tracing engine, its runnable executables, Python
helpers, and the C++ headers needed by downstream developers. A user must be
able to run `pip install obdeect-dev` and then a real trace through an
`obdeect` command without a source checkout or local CMake build. Make
simtools invoke that same packaged engine. An obdeect installation should then
need no obdeect Docker image on supported platforms.

This removes the `ghcr.io/gammasim/obdeect` image dependency from simtools dev
and production images. It does **not** remove the need for simtools container
images when running the other simulation software: those images still supply
CORSIKA, sim_telarray, hessio, and their runtime dependencies. Local simtools
users who already have those dependencies may install obdeect with pip without
using an obdeect image.

## Implemented in this repository state

- The distribution metadata is now `obdeect-dev`; its wheels are built with
  scikit-build-core and include the C++ executables, headers, and an exported
  CMake target.
- `obdeect-simtools-raytrace`, `obdeect-reference`, and `obdeect-ctao` launch the
  packaged native binaries. `obdeect.executable_path()` is the discovery API
  used by simtools.
- simtools has a `gammasimtools[obdeect]` extra, resolves the installed wheel,
  and records the pinned wheel version in dependency-manifest schema 0.2.0.
  Its Docker files no longer copy an obdeect runtime image; they install
  `obdeect-dev`.
- simtools has a `SimulatorObdeect` reference runner, can execute the packaged
  CLI from `RayTracing`, reads `obdeect-arrival-v1` files in `PSFImage`, and
  computes reference incident-angle tables from the same records. Weighted
  arrival PSF inputs are rejected explicitly.
- A model adapter can export nominal LST/MST single-reflector scenes as the
  strict, dependency-free `obdeect-scene-v1` surface table. The C++ executable
  parses that table, traces the imported finite facets and detector boundary,
  and emits the same arrival contract. simtools accepts the file through
  `obdeect_scene_file`; provenance remains attached to the scene file.
- The installed wheel's exported CMake target has been verified from an
  external consumer project, and simtools includes a published-wheel smoke
  integration test that is skipped when `obdeect-dev` is not installed.
- A GitHub Actions workflow builds platform wheels, smoke-tests the installed
  command, builds an sdist, publishes stable GitHub Releases through PyPI
  trusted publishing, and routes prereleases/manual candidates to TestPyPI.
  Tag pushes remain available for build-only checks.

The full production scene/arrival adapter remains a separate scientific
integration gate. The packaged C++ reference tracer and nominal imported scene
path are runnable and usable independently, but they do not by themselves
make every simtools production workflow equivalent to sim_telarray.

## Review outcome

The plan remains because the production telescope-scene work is not complete.
The current workspace has no validated aligned LST/MST/SST scene adapter,
material/obscurer model, or sim_telarray comparison fixtures. Removing this
file would incorrectly imply production equivalence. The package and reference
backend work below is implemented; the remaining gates are intentionally kept
as the acceptance checklist.

## PyPI registration and release procedure

The publishing workflow uses PyPI Trusted Publishing. It requests a short-lived
OIDC credential through the `pypa/gh-action-pypi-publish` action, so no PyPI API
token is stored in GitHub. Trusted Publisher setup is a one-time PyPI account
operation.

For the current repository (`github.com/GernotMaier/obdeect`):

1. Confirm that the project name `obdeect-dev` is available on PyPI. If it does
   not exist, open **Account settings → Publishing → Add a new pending
   publisher**. If it already exists, open that project's **Publishing** page
   instead.
2. Select **GitHub Actions** and enter repository owner `GernotMaier`,
   repository name `obdeect`, workflow filename `.github/workflows/pypi.yml`,
   and environment `pypi`. For a new project, enter `obdeect-dev` as the
   project name and add the pending publisher.
3. In GitHub, create an environment named `pypi` under repository settings.
   Required reviewers can be added there if releases should require a human
   approval. The workflow already grants `id-token: write` only to its publish
   job.
4. Create a GitHub Release whose tag is `v<version>` (for example `v0.1.0`).
   The `release: published` event builds Linux and macOS wheels, builds the
   source distribution, runs the wheel smoke test, and publishes stable
   artifacts after the build jobs pass. A tag push alone only exercises the
   build jobs. GitHub prereleases are sent to TestPyPI instead of production
   PyPI.
5. Verify the result with `python -m pip install obdeect-dev` and a small
   `obdeect-simtools-raytrace` run. PyPI versions are immutable; use a new
   version if an upload partially succeeds.

For prerelease testing, create the same Trusted Publisher on TestPyPI with
environment `testpypi`, then publish a GitHub prerelease. A manual workflow run
can also set the `publish_testpypi` input. TestPyPI uses its own project and
account configuration; its packages are not production PyPI releases.

When the repository moves to `github.com/gammasim/obdeect`, add a second
Trusted Publisher for owner `gammasim`, repository `obdeect`, and the same
workflow/environment before publishing from the new location. Keep the PyPI
distribution name `obdeect-dev` unless a separately planned package rename is
approved; transferring a GitHub repository does not rename a PyPI project.

## Current state

- `obdeect/pyproject.toml` uses scikit-build-core to build the Python helpers
  and native wheel payload. `obdeect/CMakeLists.txt` builds `obdeect_reference`,
  `obdeect_ctao`, and `obdeect-simtools-raytrace`.
- The EventIO reader is an optional CMake target and requires a separate
  hessio/EventIO C library. The current Dockerfile does not enable it.
- The Python native binding directory is reserved but has no binding yet. A
  binding is not needed to distribute and invoke the C++ executable.
- simtools still accepts an explicit legacy `obdeect_path`, but its packaged
  path resolution is implemented. Both simtools Dockerfiles install a pinned
  `obdeect-dev` wheel instead of copying `/opt/obdeect` from an obdeect image.
- The simtools obdeect path runs the packaged reference CLI and accepts a
  provenance-bound nominal scene file. It still does not compile run-specific
  alignment, obscurers, materials, camera surfaces, or dual-mirror geometry;
  production scene and arrival validation remain scientific integration gates.

## Target architecture

1. The initial PyPI distribution is named `obdeect-dev`, subject to confirming
   name availability and project ownership before publishing. This temporary
   distribution name reflects that the repository will move to the gammasim
   organization later. Keep the import package name `obdeect` and existing
   `obdeect-*` commands, so the install command is `pip install obdeect-dev`
   while Python code continues to use `import obdeect`. When the repository
   moves, plan any PyPI distribution-name transition explicitly: a repository
   transfer does not automatically rename an existing PyPI project.
2. Build the C++ code with CMake during wheel creation, using a CMake-aware
   PEP 517 backend such as scikit-build-core. Install native executables into
   `obdeect` package data (for example, `obdeect/_native/`), expose installed
   ray-tracing commands that forward to them, and provide a small Python
   function that returns the installed tracer path for simtools. The Python
   command is a launcher; tracing remains in C++. Do not invoke a compiler at
   import time or download binaries during installation.
3. Produce correctly tagged platform wheels for the supported Linux and macOS
   architectures, including the Linux environment used by simtools. Add
   Windows only after the C++ build and runtime tests pass there. Because the
   first release has no Python native extension, evaluate a `py3-none-<platform>`
   wheel tag; never label a wheel containing executables as `*-any`.
4. Publish a source distribution that contains all Python and C++ build
   inputs. Keep the standalone CMake build working. Ship the C++ headers and a
   relocatable CMake package/export in the wheel after testing an out-of-tree
   C++ consumer. These development artifacts are distinct from the compiled
   runtime needed to execute ray tracing, but are part of the requested full
   C++ package.
5. Keep telescope/model files as versioned inputs. The wheel contains engine
   code, not a silent snapshot of mutable production models.

## Implementation order

### 1. Define the release contract in obdeect

- Record the exact binaries and headers to ship, supported platforms, minimum
  OS/CPU targets, and whether the EventIO reader is part of the first release.
  The wheel must include the complete standalone tracer. A release described
  as a simtools alternative must also include a production-ready
  simtools-compatible executable; do not use the current reference-only CLI
  to satisfy that requirement.
- Decide how EventIO is supplied. Its current CMake target depends on an
  external hessio library. For an all-in-one EventIO wheel, audit the library
  license, bundle or statically link a compatible build, repair native-library
  dependencies, and test reading a real fixture on every target platform.
  Otherwise document EventIO as an optional adapter and have simtools prepare
  the tracer's photon input using its existing EventIO handling. Do not claim
  that the optional reader is bundled if it is absent.
- Define a versioned CLI/input/output contract and expose the obdeect version
  and contract version. This allows simtools to reject an incompatible wheel.

### 2. Build and verify obdeect wheels

- Migrate `pyproject.toml` to the selected CMake-aware backend and adjust
  `CMakeLists.txt` so wheel installs include native artifacts while ordinary
  `cmake --install` remains useful. Ensure the sdist is self-contained.
- Add a package-level executable locator and narrow command wrappers. The
  locator must report a clear error for an unsupported platform or a damaged
  installation. Check executable permissions and native-library resolution.
- Use CI to build wheels for the supported operating systems/architectures,
  test each *built wheel in a fresh environment*, and run a native smoke trace
  through the installed command plus a Python import/CLI check. Test sdist
  rebuild separately. On Linux, check
  the wheel against its declared manylinux/musllinux policy and repair bundled
  native dependencies where required. Ensure no accidental local CPU flags.
- Publish release candidates to TestPyPI, then publish a published GitHub
  Release to PyPI using a narrowly scoped trusted publishing workflow. Confirm
  PyPI name ownership before configuring the publisher. Record the released
  wheel version, artifact hash, compiler/build options, and source revision.

### 3. Make simtools consume the installed package

- Add a bounded `obdeect-dev` distribution requirement, preferably as a
  `gammasimtools[obdeect]` extra while `sim_telarray` remains the default.
  Production container builds explicitly install this extra with a pinned
  compatible version. If obdeect becomes mandatory for all simtools users,
  move the requirement to core dependencies only after that decision.
- Change `settings.obdeect_exe` to use the installed package's executable
  locator. Retain an explicit executable override for local development and a
  deprecation period for `obdeect_path`; remove the path argument, environment
  variable, and old tests after the transition. Report package/version and
  supported contract when obdeect is selected.
- Finish the remaining production scene and arrival contract already identified in
  `obdeect/docs/STATUS.md`: imported/aligned scenes, photon source data,
  detector arrivals, loss accounting, and versioned output. The packaged
  reference runner and CSV PSF/incident-angle readers are implemented; keep
  sim_telarray as the default until the model-driven adapter and comparative
  validation pass.
- The packaged reference and nominal scene paths are covered by the wheel
  smoke trace, scene compiler tests, the `SimulatorObdeect` command test, and
  the `PSFImage`/arrival-contract tests.
  Add a published-wheel simtools integration fixture and a cross-tool fixture
  for each production telescope class before claiming production equivalence.

### 4. Remove the obdeect image plumbing

- Remove the obdeect `FROM`, `COPY`, image args, path variables, labels, and
  image inputs from both simtools Dockerfiles and the corresponding build
  workflows. Install the pinned wheel in the simtools Python environment.
- Replace `obdeect_image` in simtools' dependency manifest with the installed
  obdeect distribution version, wheel hash or artifact identity, source
  revision, and CLI contract version. Version the manifest schema if its
  serialized shape changes; keep readers of older manifests working.
- Retire obdeect's runtime `Dockerfile` and its image-publishing workflow only
  after wheel-consuming simtools dev and production images pass integration
  checks. Keeping containers as CI builders for manylinux wheels is compatible
  with removing the separately published obdeect runtime image.

## Release gates

1. `pip install obdeect-dev` installs a wheel on every declared platform without
   requiring a local C++ compiler. An installed obdeect command runs the C++
   tracer end to end, reads a small scene and photon source, and writes a
   validated result in a clean environment. Unsupported platforms receive a
   clear installation/build outcome.
2. A fresh install of `gammasimtools[obdeect]` resolves the executable and
   identifies its version with no obdeect path or image. Existing sim_telarray
   workflows continue to pass.
3. A simtools dev and production image build succeeds without the obdeect
   image stage and records obdeect wheel provenance. CORSIKA and sim_telarray
   stages remain as required by their own workflows.
4. Only after scene, arrival, and cross-tool validation may simtools expose
   obdeect as a working production backend. Release packaging can happen
   earlier, but must label the CLI as reference-only until that gate passes.
5. An external CMake project can find the packaged C++ target and compile
   against installed headers (and any compiled libraries the target requires)
   without accessing the obdeect source checkout.

## References

- [Python wheel platform tags](https://packaging.python.org/en/latest/specifications/platform-compatibility-tags/)
- [scikit-build-core CMake install directories](https://scikit-build-core.readthedocs.io/en/stable/guide/cmakelists.html#install-directories)
- [cibuildwheel](https://cibuildwheel.pypa.io/en/stable/)
- [PyPI trusted publishing](https://docs.pypi.org/trusted-publishers/)
