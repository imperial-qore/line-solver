# MATLAB-to-C++ Port Audit

Audit date: 2026-09-23  
Repository revision: `be82992e8905`  
Scope: MATLAB ground truth compared with JAR, Python-native, and C++  
Status refreshed: 2026-10-04 at `2b44c234e8` (fixes in `82ab24f433` and the 2026-10-03/04 merges up to `2b44c234e8`)

## Verdict

The C++ port is substantial, not a stub, but it is not fully MATLAB-complete.
The remaining gaps generally fail explicitly with `UnsupportedError`; this audit
found no silent do-nothing stub in the registered core C++ API.

The live gaps found on 2026-09-23 were MAM priority MMAP/PH/1 analysis and
priority CDFs, CTMC breakdowns, CTMC `rate_sched`, advanced `mmap_compress`
methods, SSA NRM sub-engines, the LQN `refpath` interlock method, LDES
probability queries, and secondary MVA/SMC branches. **As of 2026-10-03 all of
these are fixed.** What remains open:

- Native C++ layered LDES: FIXED 2026-10-04. Processor Util under replication
  is the total over the copies in both engines (JAR task-Util host multiplicity
  also fixed), and the C++ budget counts activity completions with the JAR's
  20% warmup reset. Both engines count per-observation measures by root cohort
  (a cycle straddling the reset counts nowhere, flow balance exact; no delayed
  hit for a read parked before the reset) and report every QLen as a time
  integral; the JAR's `.lqnx` reader keeps a task's own replication.
  Forwarding, partial AND-join quorums (both
  engines), entry R, task Util and expected-wait balking were fixed 2026-10-03.
  Server pools on a task and setup on an INF task stay refused, as in Java.
- Flat LDES stream seeding (C++ and Java): FIXED 2026-10-03. Equally spaced
  MRG32k3a seed offsets made the streams linearly dependent and biased every flat
  3+-way fork-join by 0.2-0.4% in X (|t| 9-17); both flat engines now hash
  (seed, offset) with the layered engine's splitmix64 (`rng::SeedHash`, Java
  `LdesStreamSeeds`) and still share one sample path.
- LN interlock baselines: FIXED 2026-10-04. `none`, `ilrate` and `refpath` match
  MATLAB in j/p/c on the lqns interlock models 10..19 (Markov overtaking, phase-2
  call residence and think-time tail, `post` update order, Python's 0-based
  interlock tables). C++ is at machine precision on 14, 15-split and 18 (raw
  iterates in the moving average; `sn_has_multi_class_heter_fcfs` now relative
  FineTol in m/j/p/c, a one-ulp rate spread had flipped exact MVA vs AMVA), and
  replication > 1 closes the think time per replica and reports totals in j/p/c
  (`lqn_repl_thinktime`). Left: JAR 14/15-split/18 at 3e-8 or less, Python 14 at
  1.3e-7, all under the test tolerance.
- PAS: the C++ reader now treats an absent `oiCutoffs` as unclamped, MATLAB
  LDES refuses an uncapped PAS buffer by name, and the MATLAB, JAR and C++
  writers emit Python's table extent instead of clamping open classes at 10
  (all FIXED 2026-10-03). The follow-ups are FIXED too (2026-10-03): the class-
  and joint-dependence tables are written over the reachable extent
  (`sn.classcap`/`sn.cap`) in all four writers, and the C++ reader reads an
  absent dependence `cutoffs` as unclamped; save -> load -> save is a fixed
  point in all four (search restricted to the buffer, table miss 0 everywhere);
  MATLAB `getProb*`, `sample*`, reward and `initFromSolver` paths reach the
  named refusal through `SolverLDES.getStruct`. Only the uncapped
  non-saturating case (open class, no buffer) is still saved over the 10-box.

These are explicit omissions rather than placeholder implementations that return
plausible but incorrect values.

## Baseline and verification

The production C++ CLI was fresh relative to the C++ headers:

- `common/line-cli` mtime: `2026-09-23 04:44:22.850652271 +0100`
- C++ headers newer than `common/line-cli`: `0`
- `common/jline.jar` mtime: `2026-09-21 11:58:10.300882047 +0100`

The C++ API registry reported 777 registered functions. A focused registry test
was compiled during the audit and passed:

- Test cases: 7 passed, 0 failed, 0 skipped.
- Assertions: 6,463 passed, 0 failed.

This verifies registry population, uniqueness, reference-file correspondence,
and the direct same-name MATLAB/C++ entries covered by the registry. It does not
prove behavioral parity for every option or solver branch.

The complete multi-hour LINE test matrix was not run because the audit was
report-only and no authorization was given for the full suite.

**Re-verification 2026-10-03** (Dramatiq run `20261003_160550`, master at
`bdbab7eee3`): phase 1 JAR build PASS (144 s), phase 2 C++ build PASS (1047 s,
picard01), phase 3 Java tests 4996 run / 2 failed (picard06), phase 4 C++ tests
3876 cases / 2 failed, 933534 assertions / 15 failed (picard04). Every failure
was a seeded LDES reference that the stream-seeding fix moved: the shared
goldens of `tut02_mg1_multiclass_solvers` and `oqn_cs_routing` (Java and C++
produce identical new values) and `test_cli_ldes_prob.cpp`. Those references
were regenerated (goldens through `run_parity_static.sh --gen`), with the
user's approval. Phases 5-8 were not run.

**Re-verification 2026-10-04** (Dramatiq run `20261004_070441`, master at
`2b44c234e8`, line-test.git at `87d2aace`): all eight phases PASS. Phase 1 JAR
build (199 s), phase 2 C++ build (2126 s, picard01), phase 3 Java tests 5015 run
/ 0 failed / 0 skipped (picard06), phase 4 C++ tests 3891 cases / 0 failed
(picard04), phase 6 MATLAB 1060 tests / 0 failed (picard07, getting-started
soft checks now counted), phases 7 and 8 parity PASS (picard01, picard02, with
SSA now in the `oqn_cs_routing` golden at 1e6 samples). Phase 5 Python 5801
passed / 1 failed: `ld_class_dependence.ipynb` lost its Jupyter kernel to a ZMQ
port collision ("Address already in use") and passed on two isolated reruns.
Seeded references moved by the 2026-10-04 fixes were regenerated with the
user's approval: four MATLAB example `.mat` files and four getting-started
tutorials (LDES entries only), `oqn_cs_routing.mat` and its golden (SSA only),
and `SolverLDESLayeredReplicationTest` E1 throughput (7.63111 -> 7.630918, exact
MVA 7.626340).

## Exact-name API screen

Every MATLAB `matlab/src/api/**/*.m` basename was screened against the C++
include and source trees using explicit word boundaries. There were 41 MATLAB
names without an exact C++ name match:

| Classification | Count | Meaning |
|---|---:|---|
| `mexify_*` build scripts | 15 | Build infrastructure, not runtime API gaps |
| Tests, benchmarks, demos, profiles, and validation scripts | 13 | Test/development infrastructure |
| Private helpers folded into callers | 4 | Behavior is present without a public C++ symbol |
| Environment helper | 1 | `lineGetAvailableMemory`; not a solver algorithm |
| Algorithms present under another name or folded into a dispatcher | 6 | Name-only differences |
| Genuine public API or feature omissions | 2 | `infer_quick_model` and `lqn_ref_routes` |
| **Total** | **41** | |

The four private helpers folded into callers are:

- `bicgstab_iterate`
- `foxglynn_weights`
- `gmres_iterate`
- `explicit_signedlogsumexp`

The six apparent algorithm misses whose behavior is present are:

- `cache_build_item_graphs`: folded into the RMF implementation as
  `rmf_detail::build_item_graphs`.
- `fj_simplex`: provided by the generic `lp::simplex_solve` implementation.
- `laplace_invert_cme`: represented by `LaplaceMethod::Cme` in the common
  Laplace inversion dispatcher.
- `ctmc_ssg_reachability`: folded into `reachable_space_generator` and CTMC
  state-space aggregation.
- `sn_has_immfeed`: represented by `NetworkStruct::has_immediate_feedback`.
- `sn_has_quorum_join`: represented by partial-join feature detection and the
  `NetworkStruct::quorum_joins` data.

The two genuine exact-name omissions are:

1. `infer_quick_model`: MATLAB and Python provide this convenience model
   factory. C++ can build the same networks programmatically but lacks the
   convenience API.
   **Status 2026-09-24: fixed** (`cpp/include/line/api/infer/infer_quick_model.h`,
   plus a JAR twin `jline.api.infer.InferQuickModel`).
2. `lqn_ref_routes`: the route analysis supporting
   `config.interlock_method='refpath'` exists only in MATLAB.
   **Status 2026-09-25: fixed** in all three ports
   (`cpp/include/line/api/lqn/lqn_ref_routes.h`, `jline.api.lqn.LqnRefRoutes`,
   `python/line_solver/api/lqn/lqn_ref_routes.py`). `refpath` matches MATLAB to
   1e-7 or better on the lqns interlock models 10, 12, 14, 15-split and 18. On
   the phase-2 models it falls back to `ilrate`, whose baseline diverged from
   MATLAB in every port until 2026-10-03, when it was fixed in j/p/c (see
   `_kb/07`).

## Live solver and algorithm gaps

### MAM priority queues and response-time distributions

**Status 2026-09-24: fixed.** `cpp/include/line/api/mam/mmapph1prio.h` ports
both analyzers; `solver_mam_basic` and `solver_mam_passage_time` serve distinct
priorities under HOL and FCFSPRPRIO, matching MATLAB SolverMAM to 12 digits
(`line-test.git/cpp/tests/test_mmapph1prio.cpp`, `test_mam_priority.cpp`).

MATLAB, Java, and Python provide the `MMAPPH1NPPR` and `MMAPPH1PRPR` analyzers.
C++ did not. A MAM model with distinct priorities under HOL or preemptive
priority reaches an explicit refusal rather than being silently treated as
ordinary FCFS/HOL.

Evidence:

- `cpp/include/line/solvers/mam/solver_mam_basic.h:424`
- `cpp/include/line/solvers/mam/solver_mam_basic.h:655`
- `cpp/include/line/solvers/mam/solver_mam_passage_time.h:34`
- `matlab/src/solvers/MAM/solver_mam_basic.m:366`
- `matlab/src/solvers/MAM/solver_mam_passage_time.m:83`

### CTMC server breakdowns

**Status 2026-09-24: fixed.** C++ SolverCTMC carries the up/down marker and the
FAILURE/REPAIR transitions; the default initial state now starts UP in all four
codebases (`line-test.git/cpp/tests/test_ctmc_breakdown.cpp`).

The C++ language and model layers carry breakdown metadata, and native LDES can
simulate breakdowns. SolverCTMC nevertheless refuses a network with a server
breakdown because its CTMC state representation does not yet carry the
up/down marker and FAILURE/REPAIR transitions used by the reference.

Evidence:

- `cpp/include/line/solvers/ctmc/solver_ctmc_analyzer.h:301`
- `cpp/include/line/lang/qn/network_struct.h:439`
- `cpp/include/line/lang/qn/network_builder.h:778`
- `python/line_solver/api/state/after_event_station.py:69`

### CTMC time-varying rate schedules

**Status 2026-09-24: fixed.** `CtmcOptions::rate_sched` rebuilds the transient
generator per segment (`line-test.git/cpp/tests/test_ctmc_rate_sched.cpp`).

MATLAB, Java, and Python rebuild the transient generator from
`options.config.rate_sched`. C++ solves only the constant-rate transient
generator. `CtmcOptions` does not expose the field, so the unsupported option
cannot currently be requested and silently ignored.

Evidence:

- `cpp/include/line/solvers/ctmc/solver_ctmc_transient.h:118`
- `matlab/src/solvers/CTMC/solver_ctmc_transient_analyzer.m:48`
- `jar/src/main/java/jline/solvers/ctmc/SolverCTMC.java:3988`
- `python/line_solver/solvers/solver_ctmc/solver_ctmc.py:4905`

### Advanced MMAP compression and fitting

**Status 2026-09-24: fixed** (all seven methods, `mmap_super_safe` order-2 and the npfqn Mixture merge; JAR/python gaps in `_kb/07`).

C++ implements the order-1 mixture path of `mmap_compress`. It explicitly
refuses:

- `mixture.order2`
- `mamap2`
- `mamap2.fb`
- `m3pp.approx_cov`
- `m3pp.approx_ag`
- `m3pp.exact_delta`
- `m3pp.approx_delta`

The unported dependency families include `mmap_mixture_fit_mmap`,
`mamap2m_fit_mmap`, `mamap2m_fit_gamma_fb_mmap`, and
`m3pp2m_fitc_theoretical`.

Evidence:

- `cpp/include/line/api/mam/mmap_compress.h:188`
- `cpp/include/line/api/mam/mmap_compress.h:210`
- `cpp/include/line/api/npfqn/npfqn_traffic_merge.h:211`

### SSA NRM sub-engine

**Status 2026-09-24: fixed.** Polling, Cache, DROP/WAITQ regions, PAS/OI and
cd/jd scaling are ported to `solver_ssa_nrm.h`. The NRM gate now refuses only
the reference's own set (Fork/Join, BAS/BBS/RSRD regions, global dependence).
Tests are in `line-test.git/cpp/tests/test_ssa_nrm_subengines.cpp`. What was
originally reported:

The C++ next-reaction-method engine is present, but not every MATLAB sub-engine
is ported. Explicitly refused cases include:

- Polling scheduling and its controller state.
- Cache access and replacement behavior.
- Finite-capacity-region admission and release machinery.
- Several scheduling policies and specialized rate laws.
- Some fork-join cases, which are also outside the reference NRM's eligibility
  envelope.

Evidence:

- `cpp/include/line/solvers/ssa/ssa_dispatch.h:176`
- `cpp/include/line/solvers/ssa/ssa_dispatch.h:193`
- `cpp/include/line/solvers/ssa/ssa_dispatch.h:347`

### LDES probability queries

**Status 2026-09-24: fixed** (`-s ldes -a prob` served off the state histogram at `sn_declared_marginal`; matches MATLAB exactly).

The C++ CLI does not port `-s ldes -a prob`. The reference methods weigh a
trajectory against the model's current state, while the C++ `NetworkStruct`
does not carry the required current-state row. CTMC and SSA probability queries
are suggested instead.

Evidence:

- `cpp/src/cli/line_cli.cpp:6465`

### Native LDES coverage

The native C++ LDES engine refuses unsupported node kinds and recommends the
subprocess client where that client can cover the model. Remaining native-only
gaps include expected-wait balking and layered-network constructs such as cache
tasks, admission constraints, heterogeneous server pools, and non-sequential
precedence.

**Status 2026-09-25: layered part fixed.** The native layered engine
(`ldes_ln_engine.h`) now covers activity graphs (OR, LOOP, AND fork/join,
phase-2 early reply), cache tasks (RR/FIFO/SFIFO/LRU/HLRU/CLIMB/QLRU, levels,
delayed-hit retrieval), admission constraints (host and task rows, FIFO drain,
deadlock reported) and heterogeneous server pools. Each agrees with
`Solver_ssj_ln` under Welch |t|<3 over 6 seeds. Async calls, activity think
time, open arrivals and PS hosts landed 2026-09-26; forwarding and partial
AND-join quorums (implemented in both engines) on 2026-10-03, when entry R
(thread acquisition, cold start charged to the waiter) and task Util (derived
X*D over the host multiplicity) were also aligned to Java and the LN streams
were decorrelated (equally spaced MRG seeds had biased every AND join by
0.4-2%). Server pools on a task and setup on an INF task are refused by name,
as `Solver_ssj_ln` refuses them, and so is a quorum join reachable after an
earlier join (in both engines). Processor Util under replication (a total over
the copies) and the 20% warmup over activity completions were aligned
2026-10-04; no definitional divergence remains.

**Status 2026-10-03: expected-wait balking fixed.** The flat engine runs
EXPECTED_WAIT and COMBINED balking natively, with the Java engine's rule (one
job-count table for all three strategies, first match, one routing draw per
match) and its sample path: six seeds agree with `common/ldes.jar` to 1e-14 with
identical balk counts. See `_kb/09-ldes-and-cache.md`.

Evidence:

- `cpp/include/line/solvers/ldes/ldes_engine.h:254`
- `cpp/include/line/solvers/ldes/ldes_engine.h:1253`
- `cpp/include/line/solvers/ldes/ldes_ln_engine.h:284`

### Secondary MVA, QNA, SMC, and inference branches

Additional explicit limitations include:

- MVA/AMVA multiserver approximations, scheduling combinations, and high-
  variability methods that are not implemented in the C++ dispatcher.
  **Status 2026-09-24: fixed.** The multiserver (default, softmin, seidmann,
  suri, conway, krzesinski), highvar (hvmva), HOL (cl, shadow) and scheduling
  arms already matched MATLAB; the C++ refusals are the ones MATLAB raises too.
  Two linearizer-arm gaps were real and are now ported: `erlang` is aliased to
  `default` (C++ refused it in `ms_term`), and an unconverged amvald iterate that
  breaks the closed populations is re-solved with the Conway rule.
- The closed-chain state-encoding layer needed by the reference QNA branch.
  **Not a live gap:** MATLAB's QNA featset refuses ClosedClass and
  SelfLoopingClass (`SolverMVA.m:360`), as the C++ one does, so the closed arms
  of `solver_qna.m` run only under `setChecks(false)`. The refusal stays.
- SMC/MG1 alternatives `MG1_NI`, `MG1_RR`, and `MG1_IS`; the FI and CR methods
  used by the active ETAQA path are present. **Status 2026-09-24: fixed** (NI, RR, IS and every
  MG1_Shifts branch, all matching MATLAB; see `_kb/07`).
- MATLAB inference's MCI Gibbs branch is not ported. This is not treated as a
  live gap because MATLAB hard-codes the TE path and the dead MCI call has an
  argument-count defect.

Evidence:

- `cpp/include/line/solvers/mva/solver_mva.h:379`
- `cpp/include/line/solvers/mva/solver_qna.h:90`
- `cpp/include/line/lib/smc/mg1.h:878`
- `cpp/include/line/api/infer/infer_gibbs.h:34`

## Stub assessment

### No silent core C++ stub found

The source screen covered `TODO`, `FIXME`, `stub`, `not ported`, and
`not implemented` markers, followed by inspection of the associated dispatch
paths. No registered core C++ implementation was found that simply returned a
dummy success value or silently skipped a requested algorithm.

The common pattern is an explicit `UnsupportedError`, often with the missing
reference behavior named in the message.

### `solver_ag_autocat` is a deliberate tombstone

`solver_ag_autocat` always throws, which superficially looks like a C++ stub.
It is not a live parity regression: the reference's autocat route is dead,
`exact` was moved out to the legacy repository, and the live reference path
falls back to the ported INAP variants.

Evidence:

- `cpp/include/line/solvers/ag/solver_ag_autocat.h:60`

### MATLAB contains an actual SSA stub

`solver_ssa_nrm_space_analyzer.m` declares nine outputs but its body contains
only a scheduling gate and a debug line, then ends without assigning the
outputs. Calling it raises an output-not-assigned error. C++ deliberately ports
the working state-space analyzer branch instead of reproducing this stub.

Evidence:

- `matlab/src/solvers/SSA/solver_ssa_nrm_space_analyzer.m:1`
- `cpp/include/line/solvers/ssa/solver_ssa_nrm_space.h:908`

## Examples and documentation

Every MATLAB basic, advanced, and getting-started example name had a C++
word-boundary match. The gallery model coverage is therefore broad. However,
the C++ examples contain 54 `TODO(cpp)` markers across 13 translation units. **Status 2026-09-25: 3 remain in 2 files.**
**Status 2026-10-03: none remain.** `passAndSwap.cpp` now catches and prints
the reference's LDES refusal of an unset PAS buffer, which `ldes::solver_ldes`
raises through `io::require_finite_pas_buffers`; `workflowModels.cpp` samples
through the new `Workflow::sample` on MATLAB's settings, matching in
distribution rather than in RNG stream. LDES, LQNS,
JMT, CTMC transient, initial-state, LDES warm start, TBI cells and
`SolverENV.getGenerator` arms were ported.
They mostly mark omitted solver invocations or presentation helpers, including
LDES, LQNS, JMT, initial-state, environment, workflow, and sampling arms. These
are partial example ports, not evidence that the corresponding core model
builder is absent.

Some comments and help text are stale. In particular, CLI help states that
QN2LQN is not ported, while a real C++ implementation exists and is used by the
QNS wrapper.

**Status 2026-09-24: the stale QN2LQN help is fixed.**

Evidence:

- Stale help: `cpp/src/cli/line_cli.cpp:10112`
- Implementation: `cpp/include/line/io/qn2lqn.h:64`
- Consumer: `cpp/include/line/solvers/wrappers/lqns/lqns_qnsolver.h:66` (QNS was merged into SolverLQNS on 2026-09-28)

## Conclusion

C++ should be described as a broad beta-quality port with a healthy API
registry and substantial solver coverage. It is not a thin facade and is not
generally stubbed. Every solver and API gap listed on 2026-09-23 has since been
closed, as have the native LDES forwarding, quorum, balking, stream-seeding and
example gaps and the LN interlock baselines. What remains is small: two layered
LDES measure conventions (replicated processor Util, warmup), the C++ LN
think-time divisor under replication.

The Java and C++ suites were run on 2026-10-03 (see "Baseline and
verification"); the MATLAB, Python and cross-codebase parity phases (5-8) were
not.

## Cross-codebase report

Legend: `yes` = implemented; `partial` = some methods or submodes; `no` = absent;
`dead` = present only in an unreachable or obsolete reference path; `mixed` =
the grouped family varies by subfeature.

| Symbol or family | m | j | p | c | Verdict | Evidence |
|---|---|---|---|---|---|---|
| `cache_build_item_graphs` | yes | mixed | mixed | yes | Nested/folded in `rmf_detail::build_item_graphs` | `cpp/include/line/api/cache/cache_miss_rmf.h:552` |
| `fj_simplex` | yes | mixed | mixed | yes | Alias of generic `lp::simplex_solve` | `cpp/include/line/util/simplex.h:286`; `cpp/include/line/api/fj/fj_tsm_capacity.h:136` |
| `laplace_invert_cme` | yes | yes | yes | yes | Folded into the Laplace dispatcher | `cpp/include/line/api/lti/laplace_invert.h:418` |
| `ctmc_ssg_reachability` | yes | yes | yes | yes | Refuted; behavior lives in CTMC state-space generation and aggregation | `cpp/include/line/solvers/ctmc/solver_ctmc.h:1121`; `cpp/include/line/solvers/ctmc/solver_ctmc.h:1275` |
| `sn_has_immfeed` | yes | mixed | yes | yes | Alias of `NetworkStruct::has_immediate_feedback` | `cpp/include/line/lang/qn/network_struct.h:1632` |
| `sn_has_quorum_join` | yes | mixed | yes | yes | Refuted; represented by partial-join feature detection and `quorum_joins` | `cpp/include/line/lang/qn/feature_set.h:943`; `cpp/include/line/lang/qn/network_struct.h:2686` |
| `infer_quick_model` | yes | yes | yes | yes | Fixed 2026-09-24, with a JAR twin `jline.api.infer.InferQuickModel` | `cpp/include/line/api/infer/infer_quick_model.h:1`; `matlab/src/api/infer/infer_quick_model.m:1`; `python/line_solver/inference/api/infer_quick_model.py:6` |
| `lqn_ref_routes` / `refpath` | yes | yes | yes | yes | Fixed 2026-09-25; `none`/`ilrate` baselines fixed 2026-10-03, so the phase-2 fallback matches too | `cpp/include/line/api/lqn/lqn_ref_routes.h:1`; `jar/src/main/java/jline/api/lqn/LqnRefRoutes.java:1`; `python/line_solver/api/lqn/lqn_ref_routes.py:1`; `matlab/src/api/lqn/lqn_ref_routes.m:8`; `matlab/src/solvers/LN/@SolverLN/buildLayersRecursive.m:1409` |
| MAM priority MMAP/PH/1 and CDF | yes | yes | yes | yes | Fixed 2026-09-24 (`mmapph1prio.h`), matches MATLAB to 12 digits | `cpp/include/line/api/mam/mmapph1prio.h:1`; `cpp/include/line/solvers/mam/solver_mam_basic.h:655`; `cpp/include/line/solvers/mam/solver_mam_passage_time.h:34` |
| CTMC server breakdown/repair | yes | yes | yes | yes | Fixed 2026-09-24; default initial state starts UP in all four | `cpp/include/line/solvers/ctmc/solver_ctmc_analyzer.h:301` |
| CTMC transient `rate_sched` | yes | yes | yes | yes | Fixed 2026-09-24 (`CtmcOptions::rate_sched`) | `cpp/include/line/solvers/ctmc/solver_ctmc_transient.h:118` |
| Advanced `mmap_compress` methods | yes | yes | yes | yes | Fixed 2026-09-24: all seven methods; JAR/python gaps in `_kb/07` | `cpp/include/line/api/mam/mmap_compress.h:188` |
| SSA NRM polling/cache/regions/schedules | yes | yes | partial | yes | Fixed 2026-09-24: all five sub-engines ported | `cpp/include/line/solvers/ssa/ssa_dispatch.h:176`; `cpp/include/line/solvers/ssa/ssa_dispatch.h:347` |
| LDES current-state probabilities | yes | yes | yes | yes | Fixed 2026-09-24: served off the state histogram at `sn_declared_marginal` | `cpp/src/cli/line_cli.cpp:6465` |
| Native LDES specialized constructs | yes | yes | yes | yes | Fixed 2026-10-03: forwarding, partial AND quorum (both engines), balking; task pools / INF-task setup refused as in Java; flat stream seeding fixed 2026-10-03 (hashed seeds, both engines); LN replicated processor Util, 20% warmup, root-cohort counting and integral QLen aligned 2026-10-04 | `cpp/include/line/solvers/ldes/ldes_engine.h:254`; `cpp/include/line/solvers/ldes/ldes_ln_engine.h:284` |
| MVA/QNA/SMC alternative branches | yes | mixed | mixed | yes | Fixed 2026-09-24: MG1 methods and the two linearizer multiserver gaps ported; QNA closed arm is dead in MATLAB too | `cpp/include/line/solvers/mva/solver_qna.h:90`; `cpp/include/line/lib/smc/mg1.h:878` |
| `solver_ssa_nrm_space_analyzer` | stub | mixed | mixed | yes | MATLAB is stubbed; C++ ports the working analyzer branch | `matlab/src/solvers/SSA/solver_ssa_nrm_space_analyzer.m:1`; `cpp/include/line/solvers/ssa/solver_ssa_nrm_space.h:908` |
| `solver_ag_autocat` | dead | no | no | dead | Deliberate refusal-only tombstone for a dead reference route | `cpp/include/line/solvers/ag/solver_ag_autocat.h:60` |
| Example gallery solver arms | yes | mixed | mixed | yes | Fixed 2026-10-03: no `TODO(cpp)` markers remain | `cpp/examples/` |
