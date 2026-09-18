# Native binding boundary

This directory is reserved for the nanobind module described in the project
plan. It is intentionally not compiled yet: the C++ `PhotonBlockView` and
`TraceResult` contracts must remain stable before a Python binary API is
introduced. Python must not reimplement or alter core optical physics.
