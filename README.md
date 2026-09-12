# LINE Solver

[![Download LINE](https://img.shields.io/sourceforge/dt/line-solver.svg)](https://sourceforge.net/projects/line-solver/files/latest/download)
[![Download LINE](https://a.fsdn.com/con/app/sf-download-button)](https://sourceforge.net/projects/line-solver/files/latest/download)

LINE is an open-source software package to analyze queueing models via analytical methods or simulation. The package is developed by the [QORE lab](https://qore.doc.ic.ac.uk/) at Imperial College London and distributed under the BSD-3 license.

The package offers solution algorithms for queueing systems (e.g., M/M/1, M/M/k, M/G/1, ...), queueing networks, layered queueing networks, and queueing models in random environments. Models are solved in LINE either natively or via external solvers, such as [JMT](http://jmt.sourceforge.net/), [LQNS](http://www.sce.carleton.ca/rads/lqns/), [MAMSolver](https://www.cs.wm.edu/MAMSolver/), [Q-MAM](https://win.uantwerpen.be/~vanhoudt/), [SMCSolver](https://win.uantwerpen.be/~vanhoudt/), and [BuTools](http://webspn.hit.bme.hu/~telek/tools/butools/). See the [LINE website](http://line-solver.sf.net).

## Available Versions

| Version | Folder                       | Requirements             | Maturity | Manual | API Reference |
|---------|------------------------------|--------------------------|----------|--------|---------------|
| [MATLAB](matlab/) | matlab/                      | MATLAB                   | Stable       | [PDF](https://line-solver.sourceforge.net/doc/LINE-user-matlab.pdf), [Primer](https://line-solver.sourceforge.net/doc/LINE-primer-matlab.pdf) | [Sphinx](https://line-solver.sourceforge.net/sphinx-matlab/index.html) |
| [Java](jar/) | jar/                        | Java SE 8+               | Stable       | [PDF](https://line-solver.sourceforge.net/doc/LINE-user-java.pdf), [Primer](https://line-solver.sourceforge.net/doc/LINE-primer-java.pdf) | [Javadoc](https://line-solver.sourceforge.net/javadoc/index.html) |
| [Python Native](python/) | python/ | Python 3.11+             | Stable                    | [PDF](https://line-solver.sourceforge.net/doc/LINE-user-python.pdf), [Primer](https://line-solver.sourceforge.net/doc/LINE-primer-python.pdf) | [Sphinx](https://line-solver.sourceforge.net/sphinx/index.html) |
| [C++](cpp/) | cpp/ | C++17 compiler | Beta | [PDF](https://line-solver.sourceforge.net/doc/LINE-user-cpp.pdf), [Primer](https://line-solver.sourceforge.net/doc/LINE-primer-cpp.pdf) | [Doxygen](https://line-solver.sourceforge.net/doxygen-cpp/index.html) |

The `jar/` folder contains the canonical Java implementation, building `common/jline.jar`. The JAR is callable from any JVM language and is the portable option: one self-contained file, no build step, any platform with a JVM. Native Python users should use the `python/` folder.

The `cpp/` folder holds a header-only C++ port (`cpp/include/line/`), the `line-cli` binary and the native LDES simulation engine. It is the fastest implementation and ships as source; build it with `cpp/make.sh -O` for an optimized build.

## Command-Line Interface

The `line-cli.py` script solves models from the command line. It wraps the Java JAR and supports multiple solvers, LINE's native JSON format, and table, JSON and CSV output. Run
```
python line-cli.py solve example.json --solver mva
```
to solve a model,
```
python line-cli.py list solvers
```
to see available solvers, or
```
python line-cli.py info
```
for command line options and features. The script can also start a WebSocket server; an HTTP REST API is available in `io/rest-api/`.

LINE's native model format is a JSON shared across the MATLAB, Java and Python codebases. See [`example.json`](example.json) (a queueing network) and [`example_lqn.json`](example_lqn.json) (a layered network) in the repository root. It is specified by [`doc/line-model.schema.json`](doc/line-model.schema.json) (JSON Schema) and documented in the "JSON model format" appendix of each manual. `line-cli.py` also accepts JMT's [JSIMG](https://jmt.sourceforge.net/Papers/JMT_system_Manual.pdf#page=7) and LQNS's [LQNX](https://github.com/layeredqueuing/V6/blob/master/xml/lqn.xsd) formats.

## MCP Integration (for LLM-based Analysis)

LINE ships a [Model Context Protocol (MCP)](https://modelcontextprotocol.io/) server, for building and solving models from LLM clients such as [Claude Code](https://claude.ai/code) and [Claude Desktop](https://claude.ai/download). `pip install line-solver` installs the server and the MCP SDK and puts a `line-mcp` command on your path; point your MCP client at that command. See the [MCP Getting Started Guide](https://line-solver.sourceforge.net/doc/LINE-mcp.pdf) for setup and examples.

## License

LINE is released under the [BSD-3 license](LICENSE). It also invokes and
embeds third-party software; see
[THIRD-PARTY-NOTICES.md](THIRD-PARTY-NOTICES.md) for attributions, per-component
license terms, and which external solvers are redistributed rather than
installed by the user.

## Acknowledgement

LINE has been partially funded by the European Commission grants FP7-318484 (MODAClouds), H2020-644869 (DICE), H2020-825040 (RADON), and by the EPSRC grant EP/M009211/1 (OptiMAM).
