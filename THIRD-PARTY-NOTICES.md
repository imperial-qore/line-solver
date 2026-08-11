# Third-Party Notices

LINE is distributed under the BSD 3-Clause license in `LICENSE`. It also
invokes, embeds, or depends on software written by other groups. This file
records that software: what it is, who wrote it, under which terms it reaches
the user, and where it sits in this tree. It is informational and does not
modify any upstream license.

Provenance in this file was read from the artifacts in this repository (license
files, POM metadata, package metadata, source headers) at the time of writing.
Entries marked `ACTION` lack an in-tree license statement and must be resolved
against the upstream project before the next release. Do not reword an author
list, a URL, or a license name in this file from memory: re-read the source.

## 0. Invariant: LINE stays BSD-3-Clause

No component listed here may impose copyleft obligations on LINE's own source.
This is a hard constraint on what may enter the tree, not a preference. A
candidate dependency or embedded library is admissible only if its terms are
permissive (BSD, MIT, Apache-2.0, or equivalent), or if it is invoked strictly
as a separate process.

The current position, component class by component class:

| Class | Terms | Effect on LINE's license |
|---|---|---|
| Embedded numerical libraries, section 3 | Permissive where stated; none is copyleft | None |
| Java dependencies, section 4 | Apache-2.0, MIT, BSD; EPL-2.0 in test scope only | None |
| Python dependencies, section 5 | BSD-3-Clause, Matplotlib/PSF | None |
| JMT, section 2.1 | GPL, redistributed as a binary in `common/` | None while it remains a separate program, see below |

**JMT is the only copyleft item, and it does not propagate.** LINE never links
JMT, never imports its classes, and never derives source from it: `SolverJMT`
writes a model file, runs JMT in a separate process, and reads the results
back. That is aggregation of two independent programs, so the GPL governs
`common/JMT.jar` alone and LINE's own source stays BSD-3-Clause. What the GPL
does attach is a distribution obligation on the archive we ship, which is the
open item in section 2.1. Fetching JMT from upstream on first use discharges
that obligation and removes the only copyleft artifact from the repository.

Two rules preserve the invariant:

1. **Never link, import, or copy source from a copyleft component.** Copyleft
   tools are reachable only across a process boundary, through files and
   command lines.
2. **Never take a copyleft library as a dependency of the MATLAB, Java, or
   Python packages**, not even an unused declared one.

## 1. Distribution policy for external solvers

LINE's wrapper solvers call external tools **as separate processes**, over
files and command lines. No external solver is linked into `jline.jar` or
imported as a library, so each tool remains a separate work under its own
terms.

Two rules follow, and both are deliberate:

1. **Nothing is redistributed without an explicit redistribution grant.** A
   tool whose upstream states no license, or states terms that do not permit
   redistribution, is never shipped in this repository. The user obtains it
   from the upstream site.
2. **Every wrapper acknowledges its upstream at run time.** The first time a
   wrapper solver runs in a session it prints the upstream authors and the
   official project URL (`line_ack`, mirrored in
   `matlab/src/io/line_ack.m`, `jline.io.InputOutput.line_ack`, and
   `line_solver.api.io.logging.line_ack`).

Wrapper solvers maintained outside this tree supply their own acknowledgement
text through the same interface; their tools are intentionally not named in
this codebase.

### 1.1 External solvers invoked by LINE

| Tool | Authors / group | Upstream | Redistributed here? |
|---|---|---|---|
| JMT (Java Modelling Tools) | M. Bertoli, G. Casale, G. Serazzi | http://jmt.sourceforge.net/ | Yes, see 2.1 |
| LQNS / LQSIM | G. Franks, M. Woodside et al., Real-Time and Distributed Systems Group, Carleton University | http://www.layeredqueues.org/ | No, user installs |
| qnsolver (QNS) | Part of the LQNS distribution, same authors | http://www.layeredqueues.org/ | No, user installs |
| Wrapper solvers distributed outside this tree | See each tool's own site | Reported by the wrapper at run time | No, user installs |

LQNS is not redistributed for a substantive reason: its published sources
carry a bare copyright notice with no license grant, so we have no permission
to redistribute the binaries. Users must obtain LQNS from the upstream site.

## 2. Binaries redistributed in `common/`

| Artifact | Contents | Terms |
|---|---|---|
| `common/JMT.jar` | JMT, repackaged build | GPL, see 2.1 `ACTION` |
| `common/jline.jar` | LINE itself, plus the shaded Java dependencies of section 4 | BSD-3-Clause (LINE) plus section 4 |
| `common/ldes.jar`, `common/ldes` | LDES engine, a shade repackaging of `jline.jar` | Same as `jline.jar` |
| `common/line-viewer.jar` | LINE model viewer, a fat JAR bundling its Kotlin/Compose and Java dependencies | BSD-3-Clause (LINE) plus bundled dependency terms |

### 2.1 JMT: open compliance item

`common/JMT.jar` is not a stock upstream release. Its Maven metadata reads
`jmt:jmt-light`, built 2021-12-20, and it additionally bundles EJML and
`org.qore:kpc-toolbox-java`. The upstream project states that "Java Modelling
Tools (JMT) is a suite of applications released under GPL license"; the
homepage does not state the GPL version.

`ACTION`: the redistributed archive contains no copy of the GPL and no written
offer of corresponding source. Before the next release, either

- ship the upstream license text alongside the JAR, record the exact JMT
  version and revision the build derives from, and provide the corresponding
  source or a written offer for it, together with a description of the
  modifications that produced `jmt-light`; or
- stop bundling and fetch the official release from
  http://jmt.sourceforge.net/ on first use. This is the preferred option: it
  removes the compliance burden and returns download counts to the upstream
  project, which is itself a form of credit.

## 3. Numerical libraries embedded in the source tree

These are libraries written by other groups whose sources (or ports of them)
live inside this repository. Where a port exists, the MATLAB, Java, and Python
copies are listed together: the attribution is the same in all three.

| Library | Authors / group | License as shipped | In-tree paths |
|---|---|---|---|
| KPC-Toolbox | G. Casale, E. Smirni (College of William and Mary; Imperial College London) | BSD-3-Clause, `matlab/lib/kpctoolbox/LICENSE.TXT` | `matlab/lib/kpctoolbox/`, `jar/src/main/java/jline/lib/kpctoolbox/`, `python/line_solver/lib/kpctoolbox/` |
| M3A | A. Sansottera, G. Casale, P. Cremonesi | BSD-3-Clause, `matlab/lib/m3a/LICENSE` | `matlab/lib/m3a/`, `jar/src/main/java/jline/lib/m3a/`, `python/line_solver/lib/m3a/` |
| BuTools 2.0 | BuTools developers (webspn.hit.bme.hu) | Not stated in tree `ACTION` | `matlab/lib/thirdparty/BUTools/`, `jar/src/main/java/jline/lib/butools/`, `python/line_solver/lib/thirdparty/butools/` |
| SMCSolver (QBD, M/G/1, GI/M/1 solvers) | D. A. Bini, B. Meini, S. Steffe, B. Van Houdt | Not stated in tree `ACTION` | `matlab/lib/thirdparty/smcsolver/`, `matlab/lib/thirdparty/QBDfiles/`, `matlab/lib/thirdparty/MG1files/`, `jar/src/main/java/jline/lib/smc/`, `python/line_solver/lib/thirdparty/smc/` |
| Q-MAM | B. Van Houdt (University of Antwerp) | Inconsistent in tree `ACTION` | `matlab/lib/thirdparty/QMAM/`, `jar/src/main/java/jline/lib/qmam/`, `python/line_solver/lib/thirdparty/qmam/` |
| MAMSolver (ETAQA algorithms) | A. Riska, E. Smirni (College of William and Mary) | Not stated in tree `ACTION` | `matlab/lib/thirdparty/MAMSolver/` |
| FJ_codes | Z. Qiu, J. F. Perez, P. Harrison | BSD-3-Clause, Copyright 2015 Imperial College London, `python/line_solver/lib/thirdparty/fj/LICENSE.txt` | `matlab/lib/thirdparty/FJ_codes/`, `jar/src/main/java/jline/lib/fjcodes/`, `python/line_solver/lib/thirdparty/fj/` |
| aoi-fluid toolbox | O. Dogan, N. Akar, E. U. Atay | BSD-2-Clause, stated in the source headers | `matlab/lib/aoi/`, `matlab/lib/thirdparty/aoi/`, `python/line_solver/lib/thirdparty/aoi/` |
| rmf_tool | N. Gast, E. Rodriguez (Inria) | MIT, `python/line_solver/lib/rmf_tool/LICENSE` | `python/line_solver/lib/rmf_tool/`, `jar/src/main/java/jline/lib/rmf/` |
| ILT-CME inverse Laplace coefficients | CME inverse-Laplace project | Not stated in tree `ACTION` | `matlab/lib/thirdparty/iltcme/`, `python/line_solver/lib/thirdparty/iltcme/` |
| libQBD | S. Astaf'ev, IAMR Karelian Research Centre RAS | BSD-3-Clause, `matlab/lib/thirdparty/libQBD/LICENSE` | `matlab/lib/thirdparty/libQBD/` |
| LSODA (ODE integrator port) | S. Frost | MIT, `matlab/lib/thirdparty/lsoda/LICENSE` | `matlab/lib/thirdparty/lsoda/`, `python/line_solver/lib/lsoda.py`; see also `lsoda-java` in section 4 |
| freeLYAP and related Lyapunov solvers | A. Townsend (freeLYAP); Q. Fang, Nedialko, The MathWorks Inc., F. Glineur (contributed files) | MIT, `matlab/lib/thirdparty/lyap/lyap.LICENSE.txt`; BSD-3-Clause, `matlab/lib/thirdparty/lyap/license.txt` | `matlab/lib/thirdparty/lyap/` |
| DERIVEST suite | J. D'Errico | Not stated in tree `ACTION` | `matlab/lib/m3a/lib/thirdparty/derivestsuite/` |
| QP (quadratic programming) | A. Barraud | BSD-3-Clause, `matlab/lib/m3a/lib/thirdparty/QP/license.txt` | `matlab/lib/m3a/lib/thirdparty/QP/` |
| MAPMsG | O. Gursoy | BSD-3-Clause, `matlab/lib/thirdparty/MAPMsG/license.txt` | `matlab/lib/thirdparty/MAPMsG/`, `python/line_solver/lib/thirdparty/mapmsg/` |
| Hurst estimators | C. Chen | BSD-3-Clause, `matlab/lib/thirdparty/hurst_estimators/license.txt` | `matlab/lib/thirdparty/hurst_estimators/`, `python/line_solver/lib/thirdparty/hurst_estimators/` |
| graph_cc (connected components) | M. Vedenyov | BSD-3-Clause, `matlab/lib/thirdparty/graph_cc/license.txt` | `matlab/lib/thirdparty/graph_cc/` |
| uniqueperms | J. D'Errico | Author notice only, `matlab/lib/thirdparty/uniqueperms/LICENSE.txt` `ACTION` | `matlab/lib/thirdparty/uniqueperms/`, `python/line_solver/lib/thirdparty/uniqueperms/` |
| xml_read / xml_write | J. Tuszynski | BSD-3-Clause, `matlab/lib/thirdparty/xml_read/LICENSE.txt` | `matlab/lib/thirdparty/xml_read/` |

Notes on the `ACTION` rows:

- **BuTools, SMCSolver, MAMSolver, ILT-CME, DERIVEST**: the copies in this tree
  carry no license file. Obtain a written statement of terms from each upstream
  project, add the license file next to the code, and record it here.
- **Q-MAM**: the shipped license files disagree between copies.
  `matlab/lib/thirdparty/QMAM/license.txt` is a BSD notice naming an unrelated
  author (a MathWorks File Exchange utility that travelled with the directory),
  while `python/line_solver/lib/thirdparty/qmam/LICENSE.txt` is the Apache
  License 2.0. Neither states terms for Q-MAM itself. Confirm the terms with
  the upstream author and make the two copies agree.
- **uniqueperms**: the file contains an author and release notice, not a
  license grant.

## 4. Java dependencies

Declared in `jar/pom.xml`. Under the bundled-deps profile (`-P b`) these are
shaded into `common/jline.jar` and are therefore redistributed in binary form.
License names below were read from the dependency POM metadata (resolving
parent POMs) or from the `META-INF` license files inside the artifacts.

| Artifact | Version | License |
|---|---|---|
| org.apache.commons:commons-math3 | 3.6.1 | Apache-2.0 |
| org.apache.commons:commons-lang3 | 3.14.0 | Apache-2.0 |
| commons-io:commons-io | 2.11.0 | Apache-2.0 |
| commons-cli:commons-cli | 1.8.0 | Apache-2.0 |
| us.hebi.matlab.mat:mfl-core, mfl-ejml | 0.5.15 | Apache-2.0 |
| org.ejml:ejml-all | 0.41 | Apache-2.0 |
| xerces:xercesImpl | 2.12.2 | Apache-2.0 |
| xalan:xalan, xalan:serializer | 2.7.3 | Apache-2.0 |
| com.google.code.gson:gson | 2.11.0 | Apache-2.0 |
| ca.umontreal.iro.simul:ssj | 3.3.2 | Apache-2.0 |
| io.grpc:grpc-netty-shaded, grpc-protobuf, grpc-stub | 1.60.0 | Apache-2.0 |
| io.opentelemetry.proto:opentelemetry-proto | 1.0.0-alpha | Apache-2.0 |
| com.formdev:flatlaf | 3.2.5 | Apache-2.0 |
| de.xypron.jcobyla:jcobyla | 1.3 | MIT |
| org.java-websocket:Java-WebSocket | 1.5.1 | MIT |
| com.quantego:josqp | 0.6.5 | MIT |
| org.apfloat:apfloat | 1.10.1 | MIT |
| org.slf4j:slf4j-nop | 1.7.36 | MIT |
| net.sf.jung:jung-api, jung-graph-impl, jung-algorithms, jung-visualization | 2.1.1 | BSD |
| lsoda-java:lsoda | 1.0 | MIT, see section 3; served from https://github.com/imperial-qore/lsoda-java |
| org.junit.jupiter:junit-jupiter-* | 5.10.1 | EPL-2.0, test scope, not redistributed |

## 5. Python dependencies

Declared in `python/pyproject.toml` and installed by the user's package
manager; not redistributed in this repository. Licenses read from the
installed package metadata.

| Package | Constraint | License | Scope |
|---|---|---|---|
| numpy | ^2.3.3 | BSD-3-Clause | required |
| scipy | ^1.16.2 | BSD-3-Clause | required |
| pandas | ^2.3.2 | BSD-3-Clause | required |
| matplotlib | ^3.10.6 | Matplotlib License, PSF-based | required |
| websockets | ^12.0 | BSD-3-Clause | required |
| torch | >=2.0.0 | BSD-3-Clause | optional, `native` / `gpu` extras |
| sympy | >=1.10 | BSD-3-Clause | optional, `symbolic` extra |
| pytest, nbformat, nbconvert, ipykernel | see `pyproject.toml` | see each project | development only |

Every entry is permissive; no copyleft package is present in the Python
dependency set. Keep it that way: see section 0.

## 6. Citing the tools LINE relies on

Acknowledgement in a run log is the weakest form of credit. Where a result in a
paper was produced through a wrapper solver or an embedded library, cite the
upstream work as well as LINE.

The run-time acknowledgement makes this mechanical. On its first call in a
session each wrapper prints the upstream authors, the official project page,
and the canonical paper to cite:

```
SolverJMT delegates to Java Modelling Tools (JMT), by M. Bertoli, G. Casale,
G. Serazzi (Politecnico di Milano, Imperial College London). Please acknowledge
the JMT authors: http://jmt.sourceforge.net/
  Cite: M. Bertoli, G. Casale, G. Serazzi. "The JMT Simulator for Performance
  Evaluation of Non-Product-Form Queueing Networks". Proc. of the 40th Annual
  Simulation Symposium (ANSS), pp. 3-10, 2007.
```

The same reference is available in machine-readable form, so it can go into a
`.bib` file without being retyped:

| Codebase | Call |
|---|---|
| MATLAB | `line_citation('JMT')` |
| Java | `jline.io.InputOutput.line_citation("JMT")` |
| Python | `from line_solver.api.io import line_citation; line_citation('JMT')` |

`'JMT'`, `'LQNS'` and `'QNS'` are recognised; `'QNS'` returns the LQNS
reference, qnsolver being part of that distribution. The keys match
`doc/latex/biblio.bib` and `BIBLIOGRAPHY.md`. References for the embedded
numerical libraries of section 3 are collected in `BIBLIOGRAPHY.md`.

## 7. Reporting an attribution error

If you are an author of any software listed here and the attribution, the
license, or the redistribution status is wrong, this is a defect and will be
corrected. Open an issue on the LINE repository or contact the maintainers
through the LINE website, http://line-solver.sf.net.
