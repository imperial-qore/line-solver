---
name: tests
description: Run all eleven LINE test phases sequentially across Java, C++, Python native and wrappers, MATLAB native and dispatch, plus numeric and JSON parity, then report every failure without modifying code. Use only when the user invokes $tests or explicitly requests the complete multi-hour LINE test matrix.
---

# Full Test Suite: Java + C++ + Python + MATLAB + Parity

Run all LINE solver tests in sequence, observe and analyze results on-the-go, and produce a final summary of all failures and errors. Each phase is driven by its codebase's own `run-tests.sh` entry point, so the scripts stay the single source of truth for how that codebase is tested.

Order rationale: the four NATIVE suites run first, in the order j, p, c, m -- Java (1), Python native (2), C++ (3), MATLAB (4) -- so each codebase is characterized on its own implementation before anything cross-backend runs. The two build phases sit inside that block: Java (1) writes `common/jline.jar` and C++ (3) writes `common/line-cli` / `common/ldes`, so both still precede every consumer except Phase 2, whose `SolverLDES` rows read whatever `common/ldes` the tree already holds. The two parity runs follow (5, 6), since parity is the cross-codebase check and is most informative once every codebase has been characterized. The five WRAPPER (dispatch) passes run LAST (7-11), after Parity (JSON): each still has its own native baseline established earlier in the run (Phase 2 for the Python files, Phase 4 for the MATLAB suite), so a dispatch failure remains attributable, and the parity verdicts land before the longest block of the run.

Cost note: the full MATLAB `allTests` suite executes four times per invocation (Phase 4 under `lang=matlab`; Phases 9-11 under `lang=java`, `lang=python` and `lang=cpp`), and Phase 3 recompiles ~280 C++ test translation units before running them. Expect a multi-hour run.

## Phase 1: Java Tests

Run the JAR test suite via its script:
```bash
/home/gcasale/Dropbox/code/line-dev.git/jar/run-tests.sh 2>&1
```

Use a **7200000ms (2 hour) timeout**.

This runs `mvn test -DskipTests=false -Dstyle.color=always` from `jar/` and tees to `jar/logs/test_run_<timestamp>.txt`. It accepts an optional `-Dtest` filter as `$1`; do not pass one here.

Note on coverage: the script does NOT pass `-DexcludedGroups=no-ci`, so the pom's default `slow` exclusion stays in force and the **slow test group does not run** in this phase. To include it, invoke maven directly with `-Dtmp -DskipTests=false -DexcludedGroups=no-ci`. The script also omits `-Dtmp`, so build output goes to `jar/target/`, not `/tmp/target-jar/`; avoid running this phase concurrently with another build that shares `jar/target/`.

Analyze the output:
- Look for Surefire result lines: `Tests run:`, `FAILED`, `ERROR`
- Track each test class/method and its status
- Note the final Surefire summary (e.g., `Tests run: X, Failures: Y, Errors: Z, Skipped: W`)

## Phase 2: Python Native Tests

Run the Python native test suite:
```bash
/home/gcasale/Dropbox/code/line-dev.git/python/run-tests.sh 2>&1
```

Use a **7200000ms (2 hour) timeout**.

This runs pytest with `LINE_SOLVER_LANG=python` and `LINE_NATIVE_STRICT=true`, and tees to `python/logs/test_run_<timestamp>.txt`. Phases 7 and 8 run these same test files over the JAR and the C++ binary. Analyze the output:
- Look for pytest result lines: `PASSED`, `FAILED`, `ERROR`, `SKIPPED`
- Track each test function name and its status
- Note the final pytest summary line (e.g., `X passed, Y failed, Z errors`)

## Phase 3: C++ Tests (doctest)

Run the C++ doctest suite, which also rebuilds and installs `common/line-cli` and `common/ldes` for Phases 4, 5, 6, 8 and 11:
```bash
/home/gcasale/Dropbox/code/line-dev.git/cpp/run-tests.sh 2>&1
```

Use a **21600000ms (6 hour) timeout**: the script reconfigures from scratch and recompiles every file in `cpp/tests/` before running them.

It drives `cpp/make.sh -f -t`, the only sanctioned build entry point (every compiler invocation goes through `tsp` behind the shared cmake flock, at `-O0`). The fresh configure is REQUIRED, not cautious: `cpp/CMakeLists` globs `tests/*.cpp` at configure time with no `CONFIGURE_DEPENDS`, so a reused directory silently omits every test file added since. The build directory is its own (`/tmp/line-mp-suite-build`, override `LINE_CPP_BUILD_DIR`) so the deletion cannot land under another session's build. It tees to `cpp/logs/test_run_<timestamp>.txt`. Analyze the output:
- doctest reports per-case failures inline, with the file and line of the failed assertion
- The tail carries `[doctest] test cases:` / `assertions:` / `Status:`; quote those counts rather than recomputing them, since they are what distinguishes a green suite from a binary that was never built
- A compile error here blocks Phases 8 and 11: report it as a build failure, not as a test failure
- This phase now runs AFTER Phase 2, so Phase 2's `SolverLDES` rows used whatever `common/ldes` the tree already held

## Phase 4: MATLAB Tests

Run the MATLAB test suite:
```bash
/home/gcasale/Dropbox/code/line-dev.git/matlab/run-tests.sh 2>&1
```

The MATLAB backend is chosen automatically: `matlab/matlab_docker.sh` runs this
host's MATLAB inside a container when the image and the tools the phases shell
out to are available, and the native MATLAB otherwise. The line
`matlab_guard: MATLAB backend = ...` records which one ran. Containerized
instances take no host-wide MATLAB lock, so this phase no longer queues behind a
parity run; `LINE_MATLAB_DOCKER=0` forces the native path. Phases 10 and 11
always run native: the image carries no CPython, and it cannot load `line-cli`
at all (missing `liblapack.so.3`, and a glibc older than the host's). See
`_kb/08-build-and-test.md`.

Use a **7200000ms (2 hour) timeout**.

This runs MATLAB's `allTests` from `line-test.git` under the default `lang=matlab`. Analyze the output as it streams:
- Look for lines indicating test failures (e.g., `FAILED`, `AssertionError`, `Error using`, non-zero exit)
- Look for lines indicating test passes (e.g., `PASSED`, `Verification passed`)
- Track each test case name and its pass/fail/error status
- If the script exits with a non-zero code, note overall failure

## Phase 5: Parity Tests (numeric)

Run the assertion-based numeric sub-suite of `parity-static/`: each example has a single checked-in golden table and every codebase row is asserted against it within an explicit tolerance (`parity-static/tolerances.py`). Output is pure pass/fail, no manual inspection required.

```bash
cd /home/gcasale/Dropbox/code/line-dev.git/parity-static
./run_parity_static.sh 2>&1
```

Use a **600000ms (10 minute) timeout**. A per-row subprocess timeout can be raised with `./run_parity_static.sh --timeout 600`. Output is tee'd to `parity-static/logs/parity_static_<timestamp>.txt`.

Analyze the output:
- Each test id is `<example>[<ROW>]` (e.g. `cqn_oneline[MATLAB-WRAPPER]`), so a failure already pinpoints the example and row
- For each failure, relay the offending cells from the assertion message
- Track the pytest tally: total / passed / failed / skipped

Notes:
- `*-WRAPPER` rows use `common/jline.jar` and are skipped if that JAR is absent.
- To regenerate goldens (agreement-gated), use `./run_parity_static.sh --gen [example ...]`. Not part of this run.
- Never widen tolerances to force a pass; fix the divergence.

## Phase 6: Parity Tests (JSON round-trip)

Run the JSON round-trip sub-suite of `parity-static/`, which covers model serialization and the `CROSS:*` loads:
```bash
cd /home/gcasale/Dropbox/code/line-dev.git/parity-static
./run_parity_json.sh 2>&1
```

Use a **600000ms (10 minute) timeout**. A per-row subprocess timeout can be raised with `./run_parity_json.sh --timeout 600`. This runs `pytest tests/test_parity_json.py -v`, tee'd to `parity-static/logs/parity_json_<timestamp>.txt`. Never reconstruct a tally from `.pytest_cache/v/cache/lastfailed`: pytest prunes it only for tests the last run selected, so it mixes ids across runs and records no pass count.

Analyze the output with the same pytest parsing as Phase 5.

## Phase 7: Python Wrapper Tests (lang=java)

Run the same test files as Phase 2, with each solver delegating to the canonical `common/jline.jar` via JSON subprocess instead of the native implementation:
```bash
/home/gcasale/Dropbox/code/line-dev.git/python/run_tests_wrapper.sh 2>&1
```

Use a **7200000ms (2 hour) timeout**.

This sets `LINE_SOLVER_LANG=java` and `LINE_SEED=23000`, and requires both a JVM (`$LINE_JAVA`/`$JAVA_HOME`) and `common/jline.jar` (`$LINE_JLINE_JAR`); it exits non-zero with a clear message if either is missing. Note this is the lang=java dispatch, NOT the retired JPype wrapper. The script copies the JAR to a temp file and pins `LINE_JLINE_JAR` to it, so a concurrent JAR rebuild cannot invalidate the run mid-flight.

It runs `tests/` with the notebook suites ignored, then each notebook suite in its own pytest process, prints a per-session plus cumulative tally, and tees to `python/logs/test_run_wrapper_java_<timestamp>.txt`. Analyze the output:
- Look for pytest result lines: `PASSED`, `FAILED`, `ERROR`, `SKIPPED`
- Track each test function name and its status
- Read the final `TOTAL (N selected)` line, which aggregates all sessions

## Phase 8: Python Wrapper Tests (lang=cpp)

Run the same test files again, with each solver delegating to the C++ binary `common/line-cli` via JSON subprocess:
```bash
/home/gcasale/Dropbox/code/line-dev.git/python/run_tests_wrapper_cpp.sh 2>&1
```

Use a **7200000ms (2 hour) timeout**.

This sets `LINE_SOLVER_LANG=cpp` and `LINE_SEED=23000`, preflights that `line-cli` exists AND runs on this host, and pins `LINE_CLI_BINARY` to a frozen private copy so a concurrent `cpp/make.sh` cannot invalidate the run mid-flight. `cpp_dispatch` degrades to native Python for ONE reason only -- the binary is absent or unrunnable -- so a construct the C++ analyzer refuses surfaces here as a reported failure instead of being silently answered by a different engine. Solvers the port does not carry (JMT, LQNS, QNS) are unaffected by the variable. Same session split and tally as Phase 7; log at `python/logs/test_run_wrapper_cpp_<timestamp>.txt`.

## Phases 9-11: MATLAB Dispatch Tests (lang=java, lang=python, lang=cpp)

Run the Phase 4 suite three times more, with solvers delegating to the JAR, to native Python and to the C++ binary instead of running the MATLAB implementation. These are three separate scripts; run all three, in this order:
```bash
/home/gcasale/Dropbox/code/line-dev.git/matlab/run_tests_wrapper_java.sh 2>&1
```
```bash
/home/gcasale/Dropbox/code/line-dev.git/matlab/run_tests_wrapper_python.sh 2>&1
```
```bash
/home/gcasale/Dropbox/code/line-dev.git/matlab/run_tests_wrapper_cpp.sh 2>&1
```

Use a **7200000ms (2 hour) timeout** for each: each is one full `allTests` pass in its own MATLAB process.

The backend is selected by the `LINEDefaultLang` global read in `SolverOptions.m`, set in the base workspace immediately before `allTests`: `lang=java` dispatches via JLINE.m to `common/jline.jar`, `lang=python` via PYLINE.m to the native `line_solver`, `lang=cpp` via CPPLINE.m to `common/line-cli`. All three fail fast if the JAR is absent (`lineStart` errors without it, whatever the lang; override with `$JLINE_JAR`); the python one additionally pins `OMP/OPENBLAS/MKL/NUMEXPR_NUM_THREADS=1` to stop numpy from oversubscribing MATLAB's worker threads, and the cpp one preflights and freezes `line-cli` (`$LINE_CLI_BINARY`). Logs land in `matlab/logs/test_run_wrapper_<lang>_<timestamp>.txt`, one per script.

Expected failures, not regressions: JMT rejects `lang='python'`, and RL routing plus transient LQN `getTranAvg` are intentionally not bridged to native Python. `SolverLDES` pins `options.lang='java'` regardless of the global, so its rows are identical in both passes.

Analyze the output:
- Each script tees its own log and exits non-zero on failure; report the three passes as separate rows
- Attribute every failure to its pass (`lang=java` vs `lang=python` vs `lang=cpp`)
- A test that passes in Phase 4 but fails here is a dispatch/marshalling defect, not a solver defect

## Final Summary

After all eleven phases complete, produce a unified report:

### 1. Phase Summary Table

```
| Phase                | Total | Passed | Failed | Errors | Skipped |
|----------------------|-------|--------|--------|--------|---------|
| Java (mvn)           |   ... |    ... |    ... |    ... |     ... |
| Python Native        |   ... |    ... |    ... |    ... |     ... |
| C++ (doctest)        |   ... |    ... |    ... |    ... |     ... |
| MATLAB Tests         |   ... |    ... |    ... |    ... |     ... |
| Parity (numeric)     |   ... |    ... |    ... |    ... |     ... |
| Parity (JSON)        |   ... |    ... |    ... |    ... |     ... |
| Python lang=java     |   ... |    ... |    ... |    ... |     ... |
| Python lang=cpp      |   ... |    ... |    ... |    ... |     ... |
| MATLAB lang=java     |   ... |    ... |    ... |    ... |     ... |
| MATLAB lang=python   |   ... |    ... |    ... |    ... |     ... |
| MATLAB lang=cpp      |   ... |    ... |    ... |    ... |     ... |
```

### 2. Failed/Errored Tests Table

List every test that failed or errored across all phases:

```
| Phase           | Test Name                        | Status | Details                              |
|-----------------|----------------------------------|--------|--------------------------------------|
| Java            | <test_class#method>              | FAIL   | <assertion or error detail>          |
| C++             | <doctest case>                   | FAIL   | <file:line and failed assertion>     |
| Python Native   | <test_function>                  | FAIL   | <error message or assertion detail>  |
| Python lang=java| <test_function>                  | FAIL   | <exception type and message>         |
| Python lang=cpp | <test_function>                  | FAIL   | <exception type and message>         |
| MATLAB          | <test_name>                      | FAIL   | <error message or assertion detail>  |
| MATLAB lang=java| <test_name>                      | FAIL   | <error message or assertion detail>  |
| MATLAB lang=py  | <test_name>                      | FAIL   | <error message or assertion detail>  |
| MATLAB lang=cpp | <test_name>                      | FAIL   | <error message or assertion detail>  |
| Parity (numeric)| <example>[<ROW>]                 | FAIL   | <which metric/codebase diverged>     |
| Parity (JSON)   | <example>[<ROW>]                 | FAIL   | <which field diverged on round-trip> |
```

Report backend divergences as such, not as test bugs: a test failing in Phase 7/8 but passing in Phase 2 (or vice versa) is a Python native-vs-JAR or native-vs-C++ divergence; a test failing in Phases 9-11 but passing in Phase 4 is a MATLAB dispatch/marshalling defect. Phases 7-8 are the Python files under `lang=java` and `lang=cpp`, and Phases 9-11 the MATLAB suite under `lang=java`, `lang=python` and `lang=cpp`; report every pass as its own row.

### 3. Proposed Corrective Actions

For each failure/error, propose what should be investigated or fixed. Do NOT execute any fixes — only propose them. Consider:
- Whether the failure is a known issue (check MEMORY.md and CLAUDE.md)
- Whether the failure indicates a regression vs. a pre-existing problem
- Which codebase likely has the bug (MATLAB is ground truth for parity)
- Whether the fix is in test code vs. implementation code

## Critical Rules

- **Do NOT attempt to shorten or optimize test execution time.** Let every test run as defined.
- **Do NOT skip any phase.** All eleven phases must run to completion (or timeout).
- **Do NOT modify any test files or source code.** This is an observation-only run.
- **Do NOT batch or parallelize tests within a phase.** Run each phase sequentially as its script defines.
- **Analyze output on-the-go** to avoid losing context from long outputs. Summarize intermediate results after each phase before moving to the next.
- **If a phase times out**, report what was completed and move on to the next phase.
