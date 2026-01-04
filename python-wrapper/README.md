# LINE Solver for Python 
This folder includes the Python version of the [LINE solver](https://sourceforge.net/p/line-solver/code/ci/master).

## Installation
Requirements: Python 3.10 or later versions; a Java Virtual Machine, either Sun/Oracle JDK/JRE Variant or OpenJDK. 

On Windows, make sure that JAVA_HOME is correctly configured and java.exe is available on the system path.

## Documentation
The Python syntax is nearly identical to the MATLAB one, see for example the scripts in the Python [gettingstarted/](https://sourceforge.net/p/line-solver/code/ci/master/tree/python/gettingstarted) folder compared to the ones in the corresponding MATLAB [gettingstarted/](https://sourceforge.net/p/line-solver/code/ci/master/tree/matlab/examples/gettingstarted) folder.

A Python version of the [manual](https://sourceforge.net/p/line-solver/code/ci/master/blob/main/doc/LINE-python.pdf) is also available.

## Example
Solve a simple M/M/1 model with 50% utilization running: ```python3 mm1.py```. You should then get as output the following pandas DataFrame
```
    Station   JobClass   QLen  Util  RespT  ResidT  Tput
0  mySource     Class1    0.0   0.0    0.0     0.0   0.5
1   myQueue     Class1    1.0   0.5    2.0     2.0   0.5
```
Alternatively, you can open and run mm1.ipynb in Jupyter.

## Getting Started Examples
The `examples/gettingstarted/` folder contains tutorial examples demonstrating key LINE features. Notable examples include:

- **Tutorial 9 (tut09_dep_process_analysis.ipynb)**: Departure process analysis comparing simulated and theoretical squared coefficient of variation (SCV) using Marshall's formula. This example demonstrates validation methodology by computing relative error between simulation and analytical results.

## License
This package is released as open source under the [BSD-3 license](https://raw.githubusercontent.com/imperial-qore/line-solver/main/python/LICENSE).

## Version
This version is an early alpha release with support for basic models with open and closed classes. MVA, Fluid, MAM, and JMT solvers are mostly functional. Examples that have still some missing functionalities are marked as incomplete by a warning in the header.

