---
name: port-audit
description: Audit a LINE API family, module, or symbol across MATLAB, JAR, Python-native, and C++, then port what is genuinely missing and verify leaf-first against MATLAB. Use when the user invokes $port-audit, names a parity family to audit, or requests an end-to-end cross-codebase port audit in line-dev.git.
---

# Port Audit

Treat text following the explicit `$port-audit` invocation as the target API
family, module, or symbol. Audit it across MATLAB, JAR, Python-native and C++,
then port what is genuinely missing. If no target is supplied and none is clear
from the request, ask which family to audit before doing anything else.

The default scope of new work is `m,j,p,c` (CLAUDE.md). An audit that ends with
"C++ still lacks it" is not done, it is a report.

Each phase below exists because skipping it has already produced a wrong result
in this repository. Do not compress them.

## Phase 0: Establish a fresh baseline

**A parity number is only as fresh as the binary that produced it.** A probe
measured against a stale `common/line-cli` routed five defects, and four were
already fixed in the working tree.

Before measuring anything:

```bash
cd /home/gcasale/Dropbox/code/line-dev.git   # ALWAYS absolute; cwd persists between Bash calls
find cpp/include -name '*.h' -newer common/line-cli | wc -l   # must be 0
```

If nonzero, rebuild into a private dir before quoting any number:

```bash
PORT_AUDIT_TMP=$(mktemp -d)
cmake -S cpp -B "$PORT_AUDIT_TMP/build" -DCMAKE_BUILD_TYPE=Release
tsp -f cmake --build "$PORT_AUDIT_TMP/build" -j4 --target line_mp_tests
```

Every C++ compile and every run of a built binary goes through `tsp`, no
exceptions. `tsp -f` prints only the job id on stdout; read the output back with
`tsp -c <id>`. `tsp -f` exit 0 is NOT proof the job ran, so check the artifact
exists and its mtime is fresh. The same rule holds for the JAR: a `lang='java'`
only failure is a stale `common/jline.jar` until proven otherwise.

Quote the artifact mtime beside any number you report.

## Phase 1: Screen by name, with word boundaries

Build the candidate list by diffing symbol names across the four trees. The
naming conventions are MATLAB/Python `snake_case`, JAR `Camel_case` class files,
C++ `snake_case` in namespaced headers.

**A substring match hides the shorter name.** A present `mapqn_bnd_lr_pf` made an
absent `mapqn_bnd_lr` (34K of Java) read as present. Always search with explicit
boundaries:

```bash
grep -rnaE '(^|[^a-z0-9_])NAME([^a-z0-9_]|$)' matlab/src jar/src/main/java python/line_solver cpp/include
```

Use `-a`: grep treats some `.m` files as binary and silently reports nothing.

**Acronyms break the case converter both ways.** `SnHasDPSPRIO` is
`sn_has_dps_prio`, which neither `sn_has_dpsprio` nor `sn_has_d_p_s_p_r_i_o`
matches. Generate several candidate spellings per symbol and expect false hits.

The screen is a candidate list, not an answer.

## Phase 2: Classify every candidate before porting one

Roughly half of a name-diff's hits are not gaps. Verified categories, each seen
in this tree:

| Category | Example | Action |
|---|---|---|
| Verbatim duplicate | `cache_rayint` == `cache_spm`; `Cache_xi_bvh` == `Cache_xi_iter` | none |
| Alias / renamed | `Map_piq` is `map_prob`; `Mapqn_lpmodel` is `lp::LpModel` | none |
| Nested function | `cache_ttl_tree` lives inside `cache_ttl_lrua.m` | none |
| Folded into a solver header | `FJConvert`/`FJValidation` are in `solver_mam_fj.h`; `SnToAG` in `solver_mam_ag.h` | none |
| Language idiom | `IdcFunction` is a Java functional interface, i.e. `std::function` | none |
| Real gap | the algorithm exists nowhere on the target side | port |

**Grep each candidate for its ALGORITHM, not its name, then read both sides
before porting.** A missing symbol in one file is not missing behaviour: work
funnels through shared dispatchers here, and an audit that read a symbol's
absence as a behavioural absence produced two false "silent wrong number"
findings out of four. Trace the call graph to where the value is produced.

Report a refutation as a valid result. Four fixes to already-correct code is a
worse outcome than the original bug.

## Phase 3: Choose the reference, and distrust the convenient one

MATLAB is ground truth. The JAR is usually the cleanest thing to port FROM
because it is already imperative and typed, which is exactly why it is dangerous:
some JAR api classes are ad-hoc SUBSTITUTES for the MATLAB algorithm, not ports.
Porting one propagates a wrong answer into a new codebase.

Known substitutes (re-verify, do not assume this list is closed):

- `jline.api.wf.*`: all four `Wf_pattern_updater` convolutions return
  `params.get(0)` unchanged; `removeMatrixRows` returns inside its first
  iteration so it drops at most one row; `Wf_sequence_detector.validateSequence`
  returns false for every chain of length >= 2 because `jline.util.Pair` has no
  `equals`/`hashCode`. MATLAB has no `api/wf` at all, so **Python is the
  reference for `api/wf`**.
- `jline.api.mam.Map_block`: MATLAB carries the Maple closed form for E1, E2, E3
  and G2; the JAR discards E3 and G2 for an ad-hoc hyperexponential. Transcribe
  the MATLAB expressions MECHANICALLY (a `^` -> `pow` script), never by hand and
  never from the JAR.

**The tell is always a signature carrying parameters the body never mentions.**
Read the MATLAB and the candidate reference side by side and diff their
signatures before writing a line.

If MATLAB and the chosen reference disagree and you cannot rule which is right,
stop and ask. Do not pick the one that is easier to port.

## Phase 4: Port

- Complete implementations only. No stubs, no partial paths with a placeholder
  comment, no defensive bounds/try-catch/null checks that mask the real problem.
- Java in `jar/` must compile on Java 8: no `var`, no `List.of`, no records, no
  `Stream.toList()`, no switch arrows, no pattern-matching `instanceof`.
- Python-native must not touch JPype or a JVM.
- In-method comments are at most one line; longer rationale goes to `_kb/`.
- ASCII only, no `---` and no non-keyboard punctuation, in code, comments and
  commit messages alike.
- Porting a call site into an explicit dispatch is not a behaviour fix: it must
  keep the SAME path, or every golden re-baselines.

## Phase 5: Verify leaf-first

**Probe the leaf, then walk outwards.** The symptom of a bad port is a plausible
number, not a crash: a corrupted ME representation reported QLen 0.5 on
M/Pareto/1, which is exactly the utilization (the infinite-server answer),
against the exact Pollaczek-Khinchine 0.95. Nothing threw. Feeding the inner
algorithms directly (`mmapph1fcfs_ncmean`, `qbd_mapmap1`, `ph_equilibrium`) found
the corruption in one step; reading the station solver top-down would not have.

So:

1. Call the ported function directly on a case with a closed form or an
   independently derived value. Prefer an exact canary (a Geo/Geo tandem is exact
   by discrete Burke; M/G/1 has P-K).
2. Only then compare through the solver.
3. Compare against MATLAB with full precision. The readable AvgTable prints
   `%12.6g`, which alone yields ~1e-6 apparent deviation; use `-o json`.
4. Key any comparison table on (model, solver, **method**) and pass `--method`.
   Keying on (model, solver) lets the last method row silently win, so you
   compare a default against, say, `amva`.

Focused tests only. Do NOT start a full suite without explicit authorization;
that includes `make.sh -t` (it builds and runs all 164 test TUs), `mvn test`
without `-Dtest=`, `run-tests.sh`, a whole-directory `pytest`, `allTests*`, and
the `parity/` sweep. For focused C++ work compile just the TUs you need:

```bash
SP=<the session scratchpad directory>
rm -f $SP/t
tsp -f g++ -std=c++17 -O0 -ffp-contract=off \
  -I cpp/include -I cpp/third_party -DLINE_MP_HAVE_LAPACK \
  -DLINE_MP_CLI_BINARY='"/dev/null"' \
  -DLINE_MP_REPO_ROOT='"/home/gcasale/Dropbox/code/line-dev.git"' \
  cpp/tests/test_main.cpp cpp/tests/test_<yours>.cpp -llapack -o $SP/t
test -x $SP/t && tsp -f $SP/t -tc="*<filter>*"
```

`test_main.cpp` holds the doctest main; the library is header-only, so two TUs
is ~11 s against ~1 h for the Release suite. Linking `cpp/src/cli/line_cli.cpp`
by hand additionally needs `src/reg/api_dispatch.cpp` and
`src/api/mam/iltcme_table.cpp`, or you get undefined references after a ~30 s
compile.

Never chain a run onto a piped compile: `... | tail && ./binary` always runs the
binary, because the pipe's exit status is `tail`'s, so a failed compile is
invisible and you read the PREVIOUS binary's output. Quote `-tc` filters, and
read the case count: a doctest run exiting 0 in 0.03 s matched NOTHING.

At the end of the cycle, ASK whether the broader suites should run.

## Phase 6: Definition of done

The port is not done until all of these are true:

- [ ] All four codebases carry it, or the exception is named explicitly with what
      is missing. Do not report `m,j,p` work as complete.
- [ ] `_kb/` updated in the SAME change: the relevant page, and `_kb/index.md`
      if a page was added or removed. No changelog entry -- `_kb/log.md` was
      removed 2026-08-18 and must not be recreated.
- [ ] `.citations()` registry updated if a method, analyzer, transformation or
      approximation was added or renamed, with the `biblio.bib` key.
- [ ] Any new invariant, silent-failure mode or gotcha written down, in `_kb/`
      if it is about the code, or as a memory file if it is about how to work.
- [ ] No test assertion, expected value or tolerance was changed. If one looks
      wrong, raise it; do not edit it.
- [ ] Commit message proposed and approved, prefix in fixed order `c,m,j,p`
      (e.g. `c,m,j,p feat: ...`). Never commit without explicit consent.

## Report format

Close with a table, one row per candidate:

| Symbol | m | j | p | c | Verdict | Evidence |
|---|---|---|---|---|---|---|

`Verdict` is one of: ported, already present, alias of X, duplicate of X, nested
in X, refuted (behaviour lives in dispatcher Y). `Evidence` is a file:line, never
a recollection. State the binary mtime you measured against, and list explicitly
anything left unported and why.
