# `cpp/examples`: the example gallery in C++

The port of `matlab/examples/` and `python/examples/`. Python is the closer
reference of the two (same file names, same helper spellings); MATLAB is ground
truth whenever the two disagree.

## Layout and naming

- One translation unit per **reference directory**, not per example: a TU here
  instantiates the solver templates it calls, and 300 of them cost more to
  compile than the whole library. `examples/basic/openQN.cpp` carries every
  example of `python/examples/basic/openQN/`.
- One function per example, named **exactly** as the reference file
  (`cqn_bas_blocking`, `gallery_mm1`, `tut01_mm1_basics`).
- Each is registered with `LINE_EXAMPLE("<reference dir>", <fn>)` so
  `line-examples <name>` runs it and `line-examples --group basic/openQN` runs
  the directory. Group strings are the reference path (`basic/openQN`,
  `gallery`, `gettingstarted`).
- Model factories that other files reuse (the `gallery`) are declared in a
  header beside them; a self-contained example needs no header.

## What an example does

It reproduces the reference script's `__main__` block: build the model, run the
solvers it runs, print what it prints. Helpers in `examples_common.h`:

| helper | use |
|---|---|
| `section("MVA")` | the `print('\nSOLVER: MVA')` banner |
| `print_avg(sn, r)` | `getAvgTable` for anything returning `mva::AvgResult` (MVA, NC, MAM, BA, CTMC) |
| `print_avg_sim(sn, r)` | the same for `FluidSolution` / `SsaSolution`, which carry no ResidT or ArvR matrix |
| `kv("System response time", v)` | one labelled scalar |
| `na("JMT", "...")` | a solver this port does not carry, refused BY NAME |
| `serial(P, {a, b, c})` / `serial(P, r, {...})` | `Network.serialRouting` |
| `cyclic(P, {...})` | the same closed into a cycle, which is what `serial_routing` does when the last node is not a Sink |
| `read_trace(path)` | the samples behind `Replayer(file)` |

`LINE_EXAMPLES_REPO_ROOT` is defined by the build, so a data file is read from
the repository rather than copied.

## The library facade, and why this header still exists

Since 2026-08-31 the library itself carries the user-facing API --
`line/lang/qn/nodes.h` (`Queue`, `Delay`, `Source`, `ClosedClass`, ...),
`line/lang/distributions.h` (`Exp`, `Erlang`, `HyperExp`, ...) and
`line/solvers/solver.h` (`SolverMVA(model).avg_table()`) -- so a new example
should be written against those and reads like its Python twin line for line:

```cpp
Network model("model");
Queue  queue1(model, "Queue1", SchedStrategy::PS);
ClosedClass closed(model, "ClosedClass", 2, delay, 0);
queue1.set_service(closed, Exp(1.0));
AvgTable t = SolverMVA(model).avg_table();
```

`example_util.h` keeps the printing helpers and the analyses the facade does not
carry. Its own `solve_avg`/`AvgTable` predate the facade and are the same
computation; prefer the library spelling in new code.

## The wrapper solvers

**JMT AND LDES RUN HERE.** `example_util` exposes `jmt_avg(model, opts)` and
`ldes_avg(model, opts)`, which forward to `jmt::solver_jmt_run_analyzer` and
`ldes::solver_ldes` -- the same engines the reference drives, through the same
`common/JMT.jar` and `common/ldes`. An example whose reference calls one of them
CALLS IT, and a refusal saying "no C++ simulation engine" is stale text: it was
true before the wrappers landed and nothing re-reads a refusal when the arm
arrives. LQNS and QNS stay out; they analyse a LayeredNetwork, not the `Net`
this facade takes.

When something genuinely is missing, an example keeps the model and does **both**
of these:

1. keeps the reference's own call **as a comment**, marked `TODO(cpp)`, so the
   source records what the reference runs and the port's gap is greppable;
2. prints a runtime refusal through `na("SOLVER", "<why>")`, so a reader of the
   OUTPUT sees exactly what was not answered and a sweep cannot mistake silence
   for coverage.

```cpp
// TODO(cpp): pr = JMT(model, 'seed', 23000).getProbAggr(station)
na("JMT", "SolverJMT getProbAggr weighs the trajectory against the model's CURRENT state, "
          "which this port's NetworkStruct does not carry");
```

**NAME THE THING THAT IS ACTUALLY MISSING, and name it narrowly.** "No C++
simulation engine" was written of a getter, an option and a constraint, and in
each case the engine was there -- what was missing was the getter, the option or
the constraint. A refusal is a claim about another component and it DECAYS; check
the sibling function in the same file before believing one.

**Never** substitute another solver for a refused one: the number would be
plausible, wrong and unattributed, which is exactly what
`line-test.git/parity/test_example.sh:126-140` refuses to do. The same applies
to a method C++ does not implement -- refuse it by name rather than falling back
to the default.

## Randomness

`random.seed(n)` / `random.random()` is **reproduced exactly**, not frozen:
`py_random.h` carries CPython's MT19937 with its `init_by_array` seeding and its
`genrand_res53` double, plus `py_round` for Python's banker's rounding. Verified
digit for digit against CPython on several seeds. Use it wherever the reference
draws parameters:

```cpp
PyRandom rng(seed);                       // random.seed(seed)
const double n = py_round(rng.random() * M + 2);   // round(random.random()*M+2)
```

Freezing the drawn values as literals is wrong for a **parameterized** factory:
`gallery_cqn(M)` draws `M + 1` of them, so literals would be correct only at the
default `M` and silently wrong elsewhere.

Two streams this does NOT cover, because they are different generators:
`numpy.random` (scalar-seeded with `init_genrand`) and MATLAB's `rand`. A model
drawn from `numpy` -- `MAP.rand(seed=...)` is one -- carries its realized values
as literals, with a comment naming the call and the seed. Extract them by
running the reference, never by drawing a new sample.

MATLAB and Python already disagree on every one of these models, since their
streams differ; this port follows Python, and says so where it matters.

## Distribution fitters

`Erlang.fitMeanAndOrder`, `HyperExp.fitMeanAndSCV[Balanced]`,
`Coxian.fitMeanAndSCV`, `Cox2/Coxian.fitCentral`, `APH.fitMeanAndSCV`,
`APH.fitCentral`, `Gamma.fitMeanAndSCV` and `Pareto.fitMeanAndSCV` live in
`include/line/lang/dist_fitters.h` (namespace `line::lang`). Everything else is
a `Distrib<double>` factory in `lang_types.h`: `exp_rate`, `exp_mean`, `det`,
`erlang`, `erlang_fit`, `hyperexp`, `coxian`, `cox2`, `phase_type`, `map_dist`,
`uniform`, `pareto`, `gamma_dist`, `weibull`, `lognormal`, `replayer`,
`immediate`, `disabled_dist`, `nhpp`, `mapt`, `pht`.

## Building and checking one file

Compiles go through `tsp`, and `tsp` prints a job id rather than the compiler
output, so read the output back with `tsp -c <id>`:

```bash
cd cpp
export TSP_CLASS=cpp TSP_SLOTS_CPP=4
ID=$(tsp -f g++ -std=c++17 -fsyntax-only -O0 -I include -I third_party -I examples \
     -DLINE_EXAMPLES_REPO_ROOT='"'$PWD/..'"' examples/basic/openQN.cpp 2>/dev/null)
tsp -c $ID            # empty output means it compiled
```

Building the whole binary is `cmake --build <dir> --target line-examples`, again
through `tsp`. `line-examples --list` enumerates what is registered.
