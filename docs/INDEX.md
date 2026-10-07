# Reading the obdeect documentation

Start with [USE_CASES.md](USE_CASES.md) for commands that answer an optical
question. Read [STATUS.md](STATUS.md) to find what has been demonstrated and
what still prevents production use.

| Document | What you use it for |
| --- | --- |
| [USE_CASES.md](USE_CASES.md) | Run a trace, measure its image, or draw the actual compiled telescope geometry. |
| [STATUS.md](STATUS.md) | Check implementation limits and the remaining production requirements. |
| [COMPARISON_WITH_SIM_TELARRAY.md](COMPARISON_WITH_SIM_TELARRAY.md) | Understand what constitutes an equivalent scientific comparison. |
| [PRODUCTION_VALIDATION.md](PRODUCTION_VALIDATION.md) | Prepare comparison tables, declare acceptance criteria, and run the numerical gate. |
| [REFERENCE_RUN.md](REFERENCE_RUN.md) | Freeze input hashes and exact commands, execute them, and retain evidence. |
| [SIMTEL_REPLAY.md](SIMTEL_REPLAY.md) | Build and run the actual sim_telarray shared-photon diagnostic adapter. |
| [BENCHMARK.md](BENCHMARK.md) | Measure matched performance after every optical comparison passes. |
| [reference/7.0.0/README.md](reference/7.0.0/README.md) | Inspect previously archived sim_telarray PSF products and their limitations. |
| [CODE_REVIEW.md](CODE_REVIEW.md) | Check the correctness fixes and simplifications reviewed in the optical transport changes. |
| [TELESCOPE_PLOTTING_PLAN.md](../TELESCOPE_PLOTTING_PLAN.md) | Review the detailed telescope-diagnostic design and its acceptance tests. |

An optical model is a file containing the telescope surfaces, their optical
responses, and the provenance of the input data. A nominal model can demonstrate
geometry and numerical transport while still lacking production physics. A
reference-run record describes the exact experiment; its existence is not evidence
that two tracers agree. A comparison summary supplies measured residuals and states
whether the predeclared acceptance criteria were met.
