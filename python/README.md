# LINE Solver for Python 
This folder includes the Python version of the [LINE solver](http://line-solver.sf.net).

## Installation
Requirements: Python 3.11 or later.

Install from PyPI with:
```
pip install line-solver
```
Alternatively, install from source by running `pip install .` in this folder.

After installing, run `line-install` (or `python3 -c "import line_solver; line_solver.line_install()"`) to check optional dependencies. It warns when Java or the symbolic backend is missing without failing.

## Symbolic backend (optional)
The symbolic methods of `SolverCTMC`/`SolverFluid` use a SageMath service packaged as the Docker image [`imperialqore/line-sage-rest`](https://hub.docker.com/r/imperialqore/line-sage-rest). It is not pulled automatically on first use, so obtain it once with `docker pull imperialqore/line-sage-rest:latest`; thereafter LINE starts and stops a container on its own. Without it, symbolic analysis falls back on `sympy`.

## Documentation
The Python syntax is nearly identical to the MATLAB one, see for example the scripts in the Python `examples/gettingstarted/` folder compared to the ones in the corresponding MATLAB `examples/gettingstarted/` folder.

A Python version of the [manual](https://line-solver.sourceforge.net/doc/LINE-python.pdf) is also available.

## Example
Solve a simple M/M/1 model with 50% utilization running: ```python3 mm1.py```. You should then get as output the following pandas DataFrame
```
Station  JobClass     QLen     Util   RespT  ResidT    ArvR    Tput
Source   Class1          0        0       0       0       0     1.0
Queue    Class1        1.00     0.50    1.00    1.00     1.0     1.0
```
Alternatively, you can open and run mm1.ipynb in Jupyter.

## Getting Started Examples
The `examples/gettingstarted/` folder contains tutorial examples demonstrating key LINE features.

## License
This package is released as open source under the [BSD-3 license](http://opensource.org/licenses/BSD-3-Clause).
