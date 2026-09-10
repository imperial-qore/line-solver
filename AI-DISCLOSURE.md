# AI Disclosure

## Statement

During the development of LINE the authors used Anthropic's Claude models,
namely Opus 4.0, Opus 4.5, Opus 4.8, Opus 5 and Fable 5, accessed through Claude
Code, to assist with implementation, porting across the MATLAB, Java, Python and
C++ codebases, test authoring, refactoring and documentation drafting. The tool
is provided as-is, without warranty of any kind, under the terms of its license.

## Scope

Claude assistance covers source code, tests, examples and documentation across
all four codebases (`matlab/`, `jar/`, `python/`, `cpp/`), together with the
knowledge base under `_kb/` and the build and test harnesses. No Claude model is
an author or a copyright holder of any part of LINE; authorship is recorded in
`AUTHORS`, and licensing in `LICENSE`.

The analytical methods, algorithms and their published references are the work
of their respective authors and are cited in `BIBLIOGRAPHY.md`, in
`doc/latex/biblio.bib` and through the `citations()` registry exposed by each
solver. Model assistance concerns the realization of those methods in code, not
their derivation.

## Provenance

The last stable release developed without any AI assistance is **2.0.39**.
Every release after 2.0.39 was developed with the assistance described above.
Users needing a model-assistance-free baseline should use 2.0.39.

## Verification

LINE is validated by cross-language consistency testing, by the JUnit, pytest
and MATLAB test suites, and by comparison against independent tools (JMT, LQNS,
lqsim, external stochastic Petri net tools) where an external oracle exists. Algorithms taken from the
scientific literature are checked against the numerical examples reported in the
corresponding papers, in addition to being cross-checked against simulation.
Numerical results are checked against the MATLAB reference implementation, which
remains the ground truth when the codebases disagree.
