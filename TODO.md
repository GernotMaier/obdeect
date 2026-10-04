
obdeect needs to get production ready to be used as replacement for sim_telarray in simtools for ray tracing.

There the absolute top priorities are:

- implement all features required to do this and be able run efficient and precise ray traceing for all CTAO telescope types.
- implement (if needed in a separate obdeect-tests repository in the directory ../obdeect-tests) convincing validation scripts with detailed comparison against sim_telarray. We need to convince the community, meany a lot of astronomers.

Follow up the following task - if possible in parallel:

- review the code and update wording, program names, etc to use gamma-ray astronomer language (not competing geeks)
- in docs, I understand docs/STATUS.md and docs/USE_CASES.md - but all the other documents are hard to understand. What is it?
- implement all items in docs/STATUS.md. All - if you stop before your done, you failed.
- ensure that TELESCOPE_PLOTTING_PLAN.md is implemented
- work hard to develop a good description how we validate against sim_telarray. You have the sim_telarray code in ../sim_telarray and also ../simtools. Ensure that you document exact scripts, commands, input. We need to be convincing! Develop a list of items we show to demonstrate that obdeects work.
- review every single line of code and try to reduce it. There are too many lines!
- ensure that we have meaningful and signal rich tests for both C++ and python code
- add towncrier and enforce changelogs for new merge requests.

Don't stop until you are done. Don't.  Work and do the job.

Feel free to open as many branches you need (but don't and never use worktree).
