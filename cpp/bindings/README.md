# Native binding boundary

`core.cpp` builds the optional nanobind bulk API, included in portable wheels.
`obdeect._core.OpticalModel(path)` loads and validates an immutable compiled
optical model. Its `trace` method uses the same C++ transport as the native
command; Python does not reimplement optical physics.

Inputs are contiguous CPU NumPy arrays: positions and unit directions have shape
`(N, 3)`, wavelength, emission time and source weight have shape `(N,)`, all in
binary64. Photon IDs are unique unsigned 64-bit integers with shape `(N,)`.
Coordinates are telescope-local metres, wavelengths nanometres and times
nanoseconds. Conversion and copying are never implicit.

The API borrows input arrays and releases the GIL during transport. Callers must
keep those arrays unchanged until the call returns. Independent calls may share
the immutable optical model. Returned arrays own their storage and remain valid
after the model and result dictionary are deleted. `STATUS_NAMES` maps numeric
statuses to the shared terminal names.

This API traces nominal optical models. It does not qualify production optics or
provide material group delay. The current arrival time uses vacuum flight time.
