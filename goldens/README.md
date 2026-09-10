# `goldens/` -- shared cross-codebase golden results

One checked-in expected-result table per example, agreed by every codebase and
asserted by all of them. This is **data**, not a harness: nothing here runs, and
nothing here is specific to MATLAB, the JAR, Python or C++.

## Why it lives in line-dev.git

It used to live in `line-test.git/parity-static/`, beside the harness that
generated it. That was fine while ONE harness read it. It stopped being fine as
soon as the native suites started asserting against it: `python/tests/`,
`jar/src/test/` and `cpp/tests/` are all in **this** repo, so each had to locate
a sibling `line-test.git` checkout and skip when it was absent --
and a parity test that skips reads exactly like a parity test that passed.
`python/tests/test_lqnx_cache_roundtrip.py` carried that skip for months.
Moved here on 2026-08-19; that test is now unconditional.

The MATLAB suite is the one consumer still outside this repo (it lives in
`line-test.git/test/testsParity/`), and it reaches in through `$LINE_DEV` -- the
mirror of how `matlab/run-tests.sh` already reaches out to `line-test.git`.
`ParityGoldens.root` resolves `$LINE_GOLDENS`, then `$LINE_DEV/goldens`, then
`lineRootFolder()/../goldens`, and RAISES rather than skipping when none of them
exists -- for the reason in the paragraph above.

The four consumers, and the rows each asserts:

| Suite | File | Rows | Runs by default? |
|---|---|---|---|
| python | `python/tests/parity/test_parity_static.py` | PYTHON, P2J, P2C | yes |
| jar | `jar/src/test/java/jline/parity/ParityStaticTest.java` | JAVA | yes |
| cpp | `cpp/tests/test_parity_static.cpp` | CPP | yes -- opt out with `LINE_CPP_PARITY=0` |
| MATLAB | `line-test.git/test/testsParity/` | MATLAB, M2J, M2P, M2C | yes |

Each carries its own thin reader for `tolerance_policy.json` and its own port of
the comparator, and every one of them RAISES rather than skipping when the
goldens are unreachable.

A `FAMILY:method` key names a method no example constructs, so a row that only
ran the example would report it missing. What each consumer does with one:

| Suite | `FAMILY:method` |
|---|---|
| python | **replays** it (`rows.method_keys` -> `record_example._sweep`) |
| MATLAB | **replays** it (`ParityRows.methodKeys` -> `ParityRows.sweepModel`) |
| jar | **excuses** it by name, printed as `PARITY UNSWEPT` |
| cpp | **excuses** it by name, printed as `UNSWEPT` |

The golden itself is the list, so an example with no qualified key behaves
exactly as before and costs nothing.

The jar and cpp twins excuse rather than replay because a twin's body has no
model handle -- `ParityExample.Body` is a `void run()`, and `cpp/examples` is the
same shape -- so there is nothing in those suites to solve the method against.
Closing that means giving the twin a model handle, across every `*Twins` file,
and is a separate change. It is a bounded gap, not a silent one: the numbers are
still checked, by the generator, whose JAVA and CPP rows drive the same pairs
through `jline.cli.LineCLI` and `line-cli`, and whose gate refuses to write a
key either engine disagreed on. What those two suites lose is the RE-assertion
of that agreement on every later test run.

The generating harness
(`line-test.git/parity-static/tests/test_parity_static.py`) excuses them too,
for the same structural reason: it runs each example once, to produce the
goldens; it is not their measurer.

## Contents

| Path | What it is |
|---|---|
| `baselines/<example>.json` | the golden result table for one example, keyed by solver then by `(Station, JobClass)` then by metric |
| `examples.txt` | the 210 examples the numeric parity suite covers |
| `json_examples.txt` | the 115 examples the JSON round-trip suite covers |
| `test_tolerances.json` | per-example absolute tolerance; `0.0` means the codebases agree exactly |
| `tolerance_policy.json` | the tolerance **policy**: class defaults, the near-zero rule, and every justified per-`(example, solver)` override |
| `corpus.json` | facts about the corpus that more than one suite must agree on: which examples build SEVERAL models, and why each one is on that list |
| `method_coverage.json` | which examples are swept over **every runnable solver method**, and why each one was picked |

A golden carries a `class` map beside `solvers`, marking each solver
`deterministic` or `stochastic`. That is what selects the tolerance class, so a
fixed-seed simulator is never held to the 1e-6 a closed-form solver is.

## The solver key: `FAMILY`, and `FAMILY:method`

A key in `solvers` is one of two things.

**`MVA`** -- a bare FAMILY -- means *that family's default path, as the example
itself drove it*. This is what every key has always meant and every one of the
209 baselines is keyed this way, unchanged. Note the "as the example drove it":
`cqn_oneline` writes `MVA(model, method='exact')` and keys the result `MVA`, so
a bare key is not a claim about the family default in the abstract, it is the
answer the example's own solve returned.

**`MVA:schmidt`** -- `FAMILY:method` -- names one runnable method of that family.
It is an *additional* key beside the bare one, never a replacement, and it comes
from the METHOD SWEEP below. The method half is the method name `model.findSolver()`
prints, minus its family prefix: `mva.schmidt` becomes `MVA:schmidt`,
`nc.lc.ue` becomes `NC:lc.ue`.

Reading one is nearly free, because every consumer already iterates
`golden['solvers'].items()` and keys its tolerance lookup on that string. Two
rules make a colon key behave:

* **The tolerance CLASS resolves by the FAMILY.** `SSA:nrm` is a simulator
  because `SSA` is; reading the whole key against `stochasticSolvers` would hold
  every swept simulator to 1e-6 and quarantine it for being what it is. A
  golden's own `class` entry still wins where it has one, and it usually does:
  the generator records the class per pair from findSolver's own `Class` column,
  which is finer than the family -- `NC:sampling` and `NC:mcmc` are simulation
  arms of an otherwise deterministic family.
* **A tolerance OVERRIDE matches the EXACT key.** An override on `(cqn_oneline,
  MVA)` was justified by measuring MVA's default arm; spreading it over the
  family's other sixteen method names would widen sixteen gates nobody measured. To
  grant slack to one method, name it in full:

  ```json
  { "example": "cqn_oneline", "solver": "MVA:schmidt",
    "rel": 0.02,
    "note": "...", "why": "the measurement that bought it" }
  ```

  Same table, same `note`/`why` requirement, same review. Prefer quarantining
  the pair over loosening it, exactly as for a family.

The ensemble spelling `LN(NC)` is untouched and means something different: a
member SOLVER, where the method half of a swept key is a LAYERING STRATEGY.
`findSolver()` on a layered model reports `ln.srvn`, `ln.flat`, `ln.moment3` and
so on, so `LN:srvn` names how the layers were cut and `LN(NC)` names what solved
them. The layer engine is not part of the key.

## The method sweep

`MVA` appears in 115 of the 209 baselines, `JMT` in 68, `CTMC` in 33 -- and until
this extension **no non-default method of any family was compared across
codebases at all**. That is not academic: native Python's 16-name closed-population
AMVA family returned `[2, 0]` where the answer is `[1.4118, 0.5882]`, a 100%
error on any class-switching model, and every parity run stayed green for as
long as the defect existed because the goldens only ever exercised MVA's default
arm. Meanwhile `model.help()` advertises about a hundred runnable method names
per model and tells the caller each one can be acted on.

**Coverage is every runnable method, on a curated subset of examples.** Sweeping
the whole corpus would be roughly 25,000 solves per codebase per run and is not
worth it: the 209 baselines already cover breadth of MODELS, and what was
missing is breadth of METHODS. A method that is wrong is wrong on a SHAPE, so
the subset spans the shapes that make methods behave differently. The list, with
each pick's justification, is `method_coverage.json`; in brief:

| Example | Shape | Why this one |
|---|---|---|
| `test_gallery_mm1` | open, single class, M/M/1 | the exact answer is closed form, so a wrong method disagrees with arithmetic rather than with another approximation |
| `oqn_oneline` | open, two classes, PS + INF delay | everything that closes an open chain, with no closed population in the way |
| `cqn_repairmen` | closed, single class | every population recursion collapses to one chain, so a one-chain mistake has nowhere to hide |
| `cqn_oneline` | closed, two classes | the shape the AMVA family is written for, and the one whose defect motivated this |
| `mqn_basic` | mixed open + closed | the only shape where a method must combine an open and a closed chain in one fixed point |
| `cqn_multiserver` | multiserver FCFS, four classes, a `Disabled` pair | every multiserver correction is a distinct code path |
| `cqn_multichain_cs` | class switching into three chains | the construction the AMVA defect corrupted, with too many chains to get right by symmetry |
| `cs_single_diamond` | class switching by ROUTING | a different chain construction from the ClassSwitch-node one |
| `cache_replc_lru` | cache node, LRU | the cache arms of NC and the RMF/fluid methods are reached by nothing else |
| `fj_basic_closed` | fork-join, closed | the fj transformation runs before the solver, so every method sees a model the example did not build |
| `fcr_mm1kdrop` | finite buffer with drop | a method that ignores the buffer returns a plausible wrong answer |
| `ld_multiserver_ps` | load-dependent, multiserver PS | the LD arms of MVA and NC are separate implementations |
| `prio_hol_open` | HOL priority | only some methods carry priority, and those disagree the most |
| `spn_basic_closed` | stochastic Petri net | an SPN reaches the solvers through a different translation |
| `lqn_twotasks` | layered | an `ln.*` method is a LAYERING STRATEGY (`srvn`, `flat`, `moment3`), compared nowhere else |

**How a pair is asked for.** `SolverAUTO(model, 'mva.schmidt')` -- literally the
method name in the `Method` column of `model.findSolver()`. It exists with that
spelling in all four codebases, so ONE list of method names can be asked of every row;
going through each family's own constructor would need a per-family translation
table per codebase, and a table that drifted would silently sweep the wrong
method. It also means the sweep tests exactly what `findSolver` advertises.

**The runnable set is enumerated once**, from native Python's `findSolver()` on
the model the example builds, and handed to every row. A per-codebase
enumeration would make the rows structurally incomparable: a pair one codebase
never ran is an absence, not a disagreement, and the gate cannot tell those
apart.

**Simulation-class pairs get a shared seed and run length**
(`simulationDefaults` in `method_coverage.json`), given identically to every
row. Left to each engine's default the four codebases walk four different sample
paths, which is not a disagreement about the model.

**How much this is.** Measured with `example_methods.py --runnable` on
2026-09-04: the fifteen examples report **1121** runnable pairs, **1022** of them
after `excludeFamilies` drops `qns` and `lqns` -- 828 deterministic and 194
simulation-class. That is at most 1022 new `FAMILY:method` tables, and about
9,200 solves for a generation run across the nine rows. Widest is
`cqn_repairmen` at 140 pairs over nine families; narrowest is `lqn_twotasks` at
10, where only `ln` and `ldes` accept a layered model. Fewer than 1022 will
actually be written: the gate quarantines every pair the rows disagree on, and a
run that writes fewer is the gate doing its job, not a defect.

Two kinds of quarantine are expected on the first run and are findings, not
faults. A pair a codebase cannot ASK FOR at all: `line-cli` accepts only
`default`, `moment3`, `mwba.upper` and `mwba.lower` as `--method` on the layered
path (the layer engine is `--layer-solver`), so the `ln.srvn` / `ln.flat`
layering strategies native Python reports as runnable are unreachable from the
C++ CLI and that row reports `SWEEP-NA` by name. And a pair a codebase REFUSES
on this model: `line-cli` declines `mva.mva`, `mva.sum`, `mva.esum` and
`mva.amva` on a one-station one-class model that `findSolver()` marks runnable.
Both print the engine's own sentence, and both leave the key unwritten rather
than written from the rows that managed it.

### The one rule applies unchanged

A swept pair is goldened only where the rows actually agreed within tolerance,
with the values taken from the MATLAB run. A pair that disagreed, or that only
the MATLAB row could produce, is QUARANTINED and printed -- never written with
inflated slack, and never invented. Merging is additive: `gen_baselines.py
--methods` copies every existing key through untouched and only adds or replaces
`FAMILY:method` ones, because a bare key is the example's own solve and this
mode never runs the example.

### Asking for the coverage

The selection is DATA. To sweep more (or less), edit `method_coverage.json` --
add an entry with its `shape` and its `why`, or add an `excludeFamilies` /
`excludeMethods` entry with the reason. No generator change is involved.
`python/tests/parity/test_method_coverage.py` asserts the file stays
well-formed: every named example real, goldened and justified, every exclusion
naming a family that exists.

## The one rule

**A golden is never per-codebase.** `gen_baselines.py` writes a solver's table
only when the rows actually agreed within tolerance, taking the values from the
MATLAB run (the project's ground truth per `CLAUDE.md`); a solver whose rows
disagreed is quarantined and printed, never written with inflated slack. That is
the whole reason these files can be asserted by four codebases at once, and it is
why an expected value here must never be edited to make a run go green.

## Regenerating

From a checkout where MATLAB, `common/jline.jar` and `common/line-cli` are all
available:

```bash
../line-test.git/parity-static/run_parity_static.sh --gen cqn_oneline   # one
../line-test.git/parity-static/run_parity_static.sh --gen               # all
```

The METHOD SWEEP is its own mode, because it answers a different question and
because it MERGES rather than rewrites:

```bash
../line-test.git/parity-static/run_parity_static.sh --gen-methods              # the curated subset
../line-test.git/parity-static/run_parity_static.sh --gen-methods cqn_oneline  # one of them
```

Naming an example that is not in `method_coverage.json` is reported and skipped,
not swept: an example outside the curated set has no justification recorded, and
that record is the point of keeping the selection as data.

The generator resolves this directory through
`parity-static/goldens.py`: `$LINE_GOLDENS`, else `$LINE_DEV/goldens`, else the
sibling checkout. It raises rather than degrading to a skip when it finds
nothing.

## The tolerance policy

`tolerance_policy.json` was `parity-static/tolerances.py` until 2026-08-19. It
became data for the same reason the baselines are here: four codebases now
assert their own parity rows, and four hand-kept copies of a tolerance table are
four tables that drift apart. Each language carries a thin reader --
`python/tests/parity/tolerance.py` is the Python one,
`line-test.git/test/testsParity/ParityTolerance.m` the MATLAB one, and
`parity-static/tolerances.py` is now a re-export of the Python reader, so the
harness in line-test.git keeps its API.

Every override carries `note` (the one line a reviewer scans) and `why` (the
measurement that bought it). The `why` is the audit trail, not decoration: an
override without one is how a gate gets widened and nobody can later say what
for. `python/tests/parity/test_tolerance_policy.py` asserts that every override
still names an example that has a golden and a solver that golden actually
produces -- a stale override reads as a justified allowance while granting
nothing.
