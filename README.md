# LINE Solver

[![Download LINE](https://img.shields.io/sourceforge/dt/line-solver.svg)](https://sourceforge.net/projects/line-solver/files/latest/download)
[![Download LINE](https://a.fsdn.com/con/app/sf-download-button)](https://sourceforge.net/projects/line-solver/files/latest/download)

LINE is an open-source queueing network solver for performance and reliability analysis. Visit the [LINE website](http://line-solver.sf.net) for more information. 

## Available Versions

| Version | Requirements | Maturity | Manual | API Reference |
|---------|--------------|----------|--------|---------------|
| [MATLAB](matlab/) | MATLAB | Stable | [PDF](https://line-solver.sourceforge.net/doc/LINE-matlab.pdf) | [Doxygen](https://line-solver.sourceforge.net/doxygen/index.html) |
| [Java/Kotlin](jar/) | Java SE 8+ | Stable | [PDF](https://line-solver.sourceforge.net/doc/LINE-java.pdf) | [Javadoc](https://line-solver.sourceforge.net/javadoc/index.html) |
| [Python Wrapper](python-wrapper/) | Python 3.10+, Java SE 8+ | Beta | [PDF](https://line-solver.sourceforge.net/doc/LINE-python.pdf) | [Sphinx](https://line-solver.sourceforge.net/sphinx/index.html) |
| [Python Native](python/) | Python 3.10+ | Alpha | [PDF](https://line-solver.sourceforge.net/doc/LINE-python.pdf) | [Sphinx](https://line-solver.sourceforge.net/sphinx/index.html) |

The Python Wrapper interfaces with the Java/Kotlin JAR via JPype, allowing Python users to leverage the faster, more mature JAR-based solvers while using familiar Python syntax. The JAR implementation offers better performance than the native Python version for large-scale models.

## Command-Line Interface

The `line-cli.py` script provides a standalone command-line interface for solving queueing network models without writing code. It wraps the Java/Kotlin JAR and supports multiple solvers, input formats (JSIM, JSIMG, LQNX), and output formats (table, JSON, CSV). Run `python line-cli.py solve model.jsimg` to solve a model, `python line-cli.py list solvers` to see available solvers, or `python line-cli.py info` for command line options and features. The script can also start WebSocket or REST API servers for integration with other tools.

## License

LINE is released under the [BSD-3 license](LICENSE).

## Acknowledgement

The development of LINE has been partially funded by the European Commission grants FP7-318484 (MODAClouds), H2020-644869 (DICE), H2020-825040 (RADON), and by the EPSRC grant EP/M009211/1 (OptiMAM).
