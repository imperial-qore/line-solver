# LINE Solver

[![Download LINE](https://img.shields.io/sourceforge/dt/line-solver.svg)](https://sourceforge.net/projects/line-solver/files/latest/download)
[![Download LINE](https://a.fsdn.com/con/app/sf-download-button)](https://sourceforge.net/projects/line-solver/files/latest/download)

LINE is an open-source queueing network solver for performance and reliability analysis. Visit the [LINE website](http://line-solver.sf.net) for more information. 

## Available Versions

| Version | Requirements | Maturity | Documentation |
|---------|--------------|----------|---------------|
| [MATLAB](matlab/) | MATLAB | Stable | [Manual](https://line-solver.sourceforge.net/index.html#papers) |
| [Java/Kotlin](jar/) | Java SE 8+ | Stable | [Manual](https://line-solver.sourceforge.net/index.html#papers) |
| [Python Wrapper](python-wrapper/) | Python 3.10+, Java SE 8+ | Beta | [Manual](https://line-solver.sourceforge.net/index.html#papers) |
| [Python Native](python/) | Python 3.10+ | Alpha | [Manual](https://line-solver.sourceforge.net/index.html#papers) |

The Python Wrapper interfaces with the Java/Kotlin JAR via JPype, allowing Python users to leverage the faster, more mature JAR-based solvers while using familiar Python syntax. The JAR implementation offers better performance than the native Python version for large-scale models. 

## License

LINE is released under the [BSD-3 license](LICENSE).

## Acknowledgement

The development of LINE has been partially funded by the European Commission grants FP7-318484 (MODAClouds), H2020-644869 (DICE), H2020-825040 (RADON), and by the EPSRC grant EP/M009211/1 (OptiMAM).
