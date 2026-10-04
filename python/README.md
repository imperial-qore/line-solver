# LINE Solver for Python 
This folder includes the Python version of the [LINE solver](http://line-solver.sf.net).

LDES, the discrete-event simulation engine of LINE, is built to span the entire feature space of the package, so that every approximation LINE publishes can be validated against a sample path of the same model, expressed once. Several analytical solvers are included alongside it to accelerate the evaluation of specific classes of models, and their solutions can in turn warm-start an LDES run, placing it near steady state already at the start of the simulation.

## Installation
Requirements: Python 3.11 or later.

Install from PyPI with:
```
pip install line-solver
```
Alternatively, install from source by running `pip install .` in this folder.

After installing, run `line-install` (or `python3 -c "import line_solver; line_solver.line_install()"`) to check optional dependencies. It warns when Java or the symbolic backend is missing without failing.

## Solving with the Java or C++ engine (`lang=`)

**The same Python model can be solved by any of LINE's engines.** Every solver constructor
takes a `lang` keyword that decides which codebase actually does the work; nothing else in
the script changes:

```python
from line_solver import *

MVA(model).avg_table()                  # native Python (default)
MVA(model, lang="java").avg_table()     # solved by the Java JAR (jline.jar)
MVA(model, lang="cpp").avg_table()      # solved by the C++ engine (line-cli)
```

Both delegated routes use the same transport: the model is serialized to `model.json`, the
other codebase's command-line front end solves it, and its JSON answer is read back into the
usual result object, so `avg_table()` and the other getters behave as before. The banner a
run prints names the engine that produced the numbers (`lang: cpp`). Export
`LINE_SOLVER_LANG=java` (or `cpp`) to change the default for a whole session without editing
any script.

|  | `lang="python"` (default) | `lang="java"` | `lang="cpp"` |
|---|---|---|---|
| Engine | native Python | `jline.jar` | `line-cli` |
| Needs | nothing | a Java runtime | the `line-cli` binary |
| Solvers | all | all | `MVA`, `NC`, `CTMC`, `MAM`, `FLD`, `SSA`, `BA`, `AG`, `AUTO`, `JMT`, plus `LN` (steady state) and `ENV` |

`lang="java"` looks for the JAR in `$LINE_JLINE_JAR`, inside the installed package, and in a
checkout's `common/`, and downloads it from SourceForge when none is found, so in practice it
works out of the box wherever `java` is on the path. It serves every solver, including the
`LDES` and `LQNS` wrappers.

`lang="cpp"` locates its binary through `$LINE_CLI_BINARY`, then a checkout's `common/`, then
`cpp/build/`, then `PATH`. Build it from the source tree with:

```
cmake -S cpp -B cpp/build -DCMAKE_BUILD_TYPE=Release && cmake --build cpp/build --target line-cli
```

`LDES` is deliberately absent from its list, because the Python `LDES` wrapper is already a
subprocess client of that same C++ engine, and `LQNS` wraps an external tool the C++ front end
does not carry. `lang="cpp"` also accepts `arith="exact"` or `arith="real:<digits>"`, which
runs the solve in exact or extended-precision arithmetic instead of IEEE double; an analyzer
needing transcendental functions (`NC`, which forms throughputs as `exp(lG(N-1) - lG(N))`)
refuses `exact` by name and takes `real:<digits>`.

A missing or unrunnable binary is the one failure that falls back to native Python, with a
warning. Every other failure is raised rather than silently answered by a different engine,
since `lang` is a statement about what produced the numbers.

## Symbolic backend (optional)
The symbolic methods of `SolverCTMC`/`SolverFLD` use a SageMath service packaged as the Docker image [`imperialqore/line-sage-rest`](https://hub.docker.com/r/imperialqore/line-sage-rest). It is not pulled automatically on first use, so obtain it once with `docker pull imperialqore/line-sage-rest:latest`; thereafter LINE starts and stops a container on its own. Without it, symbolic analysis falls back on `sympy`.

## Documentation
The Python syntax is nearly identical to the MATLAB one, see for example the scripts in the Python `examples/gettingstarted/` folder compared to the ones in the corresponding MATLAB `examples/gettingstarted/` folder.

A Python version of the [manual](https://line-solver.sourceforge.net/doc/LINE-user-python.pdf) is also available.

## Example
From a source checkout, solve a simple M/M/1 model with 50% utilization by running ```python3 mm1.py``` in this folder. The script uses the deterministic MVA solver, so after an analysis banner (whose elapsed time varies) it prints the following pandas DataFrame. After a `pip install`, the tutorials ship inside the package and run from any directory as modules, e.g. ```python3 -m line_solver.examples.gettingstarted.tut01_mm1_basics```.
```
Station  JobClass  QLen  Util  RespT  ResidT  ArvR  Tput
Source   Class1       0     0      0       0     0     1
Queue    Class1       1   0.5      1       1     1     1
```
Alternatively, you can open and run mm1.ipynb in Jupyter.

## Getting Started Examples
The `examples/gettingstarted/` folder contains tutorial examples demonstrating key LINE features. They are installed with the package as `line_solver.examples.gettingstarted`, so an installed copy runs them without a checkout:
```
python3 -m line_solver.examples.gettingstarted.tut03_repairmen
```

## License
This package is released as open source under the [BSD-3 license](http://opensource.org/licenses/BSD-3-Clause).
