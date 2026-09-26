# Code review notes

Review scope: C++ trace/source/scene interfaces, Python import/analysis/plot
commands, tests, examples, and documentation at this worktree's base commit.

## Fixed here

- The focal plot changed an explicit zero throughput into unit throughput.
- The mirror-list parser silently treated malformed optional height as zero.
- The path plot accepted non-finite vertices and invalid vertex counts.
- The README repeated its quick start, had an unclosed shell fence, and mixed
  build locations. It now points to short, runnable tutorials.
- Ray plots showed the long source flight at equal scale, hiding interactions
  near the mirror and camera. They now focus on the telescope region.

## Open design work

- The toy, segmented, and general optical scenes each have separate trace
  loops. They cannot be deduplicated safely until they share one interaction
  record and terminal-status contract.
- Both native executables repeat CLI parsing and CSV writing. A common I/O
  layer would help when the trace output schema is stabilized.
- CTAO reference plots have no physical structure data. Full structure plots
  need a trace-ready compiled scene with supports, camera, and mirror panels.
- The general optical scene lacks curved surfaces, material bindings, and
  complete interaction diagnostics. See [STATUS.md](STATUS.md) for the
  validation plan and production blockers.
