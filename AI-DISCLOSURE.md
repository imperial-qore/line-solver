# AI Disclosure

## Statement

During the development of LINE the authors used Anthropic's Claude models,
namely Opus 4.x, Sonnet 4.x, Haiku 4.x, Opus 5.x, Sonnet 5.x, and Fable 5.x, 
accessed through Claude Code, to assist with implementation, porting across 
the MATLAB, Java, Python and C++ codebases, test authoring, refactoring and 
documentation drafting. Much more rarely, OpenAI's Codex has been used to
fix small bugs. The tool is provided as-is, without warranty of any kind, 
under the terms of its license.

## Scope

Claude assistance covers source code, tests, examples and documentation across
all four codebases (`matlab/`, `jar/`, `python/`, `cpp/`), together with the
documentation and the build and test harnesses. No LLM is an author or a 
copyright holder of any part of LINE; authorship is recorded in `AUTHORS`, 
and licensing in `LICENSE`.

The analytical methods, algorithms and their published references are the work
of their respective authors and are cited in `BIBLIOGRAPHY.md`, in
`doc/latex/biblio.bib` and through the `citations()` registry exposed by each
solver. Model assistance concerns the realization of those methods in code, not
their derivation.

## Provenance

The 3.x codebase makes a heavy-duty usage of AI assistance. The last stable 
release developed without any AI assistance is **2.0.39**. Users needing a 
AI-assistance free baseline should use 2.0.39, available on Sourceforge.

Every release after 2.0.39 was developed with the AI assistance described above.
Code added or changed in versions after 2.0.39 have little to none human-written 
code. Debugging is mostly automated. New features have been conceptualized and
orchestrated by a human and validated against examples given in original papers 
and/or simulation and/or external tools.

LDES and the C++ codebase are 100% AI-generated. The majority of the python 
codebase is also an automated AI-generated port. Most of MATLAB and Java code is
instead originally implemented by humans, but heavily debugged, refactored, and 
further extended via AI-assistances after version 2.0.39.

## Verification

LINE is validated by cross-language consistency testing, by the JUnit, pytest
and MATLAB test suites, and by comparison against independent tools (JMT, LQNS,
lqsim, external stochastic Petri net tools) where an external oracle exists. 
Algorithms taken from the scientific literature are checked against the numerical 
examples reported in the corresponding papers, in addition to being cross-checked 
against simulation. Numerical results across the codebases are checked against 
the MATLAB reference implementation, which remains the ground truth when 
the codebases disagree.
