# Third-Party Notices

LINE is distributed under the BSD 3-Clause license (`LICENSE`). It also invokes,
embeds, or depends on software written by other groups; this file records what
that software is, who wrote it, under which terms it reaches the user, and where
it sits in this tree. It is informational and does not modify any upstream
license. Provenance was read from in-tree artifacts (license files, POM/package
metadata, source headers); rows marked `ACTION` lack an in-tree license
statement and must be resolved upstream before release. Do not reword an author
list, URL, or license name from memory — re-read the source.

## 0. Invariant: LINE stays BSD-3-Clause

No component here may impose copyleft on LINE's own source. A dependency or
embedded library is admissible only if its terms are permissive (BSD, MIT,
Apache-2.0, or equivalent) or if it is invoked strictly as a separate process.

**JMT is the only copyleft item, and it does not propagate.** LINE never links,
imports, or derives source from JMT: `SolverJMT` writes a model file, runs JMT
as a separate process, and reads results back. The GPL therefore governs
`common/JMT.jar` alone (as aggregation), not LINE's source. Two rules preserve
this: never link/import/copy source from a copyleft component (reach it only
across a process boundary), and never take a copyleft library as a declared
dependency of the MATLAB, Java, or Python packages.

## 1. External solvers invoked by LINE

Wrapper solvers call external tools **as separate processes** over files and
command lines; none is linked or imported. Nothing is redistributed without an
explicit redistribution grant — a tool whose upstream permits no redistribution
is obtained by the user from upstream. Each wrapper prints its upstream authors
and project URL on first use (`line_ack`, mirrored in `matlab/src/io/line_ack.m`,
`jline.io.InputOutput.line_ack`, `line_solver.api.io.logging.line_ack`).

| Tool | Authors / group | Upstream | Redistributed here? |
|---|---|---|---|
| JMT (Java Modelling Tools) | M. Bertoli, G. Casale, G. Serazzi | http://jmt.sourceforge.net/ | Yes, see 2.1 |
| LQNS / LQSIM | G. Franks, M. Woodside et al., Carleton University | http://www.layeredqueues.org/ | No, user installs (upstream states no license grant) |
| qnsolver (QNS) | Part of the LQNS distribution, same authors | http://www.layeredqueues.org/ | No, user installs |
| Wrapper solvers distributed outside this tree | See each tool's own site | Reported by the wrapper at run time | No, user installs |

## 2. Binaries redistributed in `common/`

| Artifact | Contents | Terms |
|---|---|---|
| `common/JMT.jar` | JMT, repackaged build (`jmt:jmt-light`, bundles EJML and `org.qore:kpc-toolbox-java`) | GPLv2-or-later, see 2.1 |
| `common/jline.jar` | LINE plus the shaded Java dependencies of section 4 | BSD-3-Clause plus section 4 |
| `common/ldes.jar`, `common/ldes` | LDES engine, a shade repackaging of `jline.jar` | Same as `jline.jar` |
| `common/line-viewer.jar` | LINE model viewer, fat JAR bundling its Kotlin/Compose and Java deps | BSD-3-Clause plus bundled dependency terms |

### 2.1 JMT source publication

Per-file GPLv2-or-later headers let the modified build coexist with Apache-2.0
EJML (as GPLv3). The corresponding source of the published JAR is shipped beside
it as `JMT-custom-sources.zip` (GPLv2 §3(a)); both are built and uploaded by
`upload-jmt.sh`, which refuses to upload unless the zip carries `COPYING` and the
`patches/` record, uploads source before binary, and pins the git revision in
`BUILD-INFO`. `ACTION`: the in-tree `common/JMT.jar` and the copy fetched by
`jmtGetPath`/`SysUtils` still travel without the license text next to them — ship
`COPYING` in `common/` or have the download fetch the sources zip.

## 3. Numerical libraries embedded in the source tree

Where a port exists, the MATLAB/Java/Python copies share one attribution.

Two entries were RETIRED on 2026-08-19 rather than resolved, and the rule they
set is the one to follow for the remaining `ACTION` rows. **DERIVEST** (its one
live caller was `map_count_moment`) and **uniqueperms** (state-space
enumeration) both lacked a license grant — uniqueperms shipped an author byline
and nothing else, DERIVEST no statement at all — so each was replaced by LINE-owned
code written from the published method — the block-Toeplitz derivative of a
truncated polynomial, and multiset permutation enumeration (Knuth's Algorithm L
in python; in MATLAB a preallocated expansion over leading values, which is what
reproduces the row order the state builders depend on) — not by transliterating
the vendored source. That distinction is the whole point: a line-by-line rewrite
would still be a derivative work and would still owe attribution, so removal
only clears an obligation when the replacement is written from the algorithm.
A related lesson for this page: uniqueperms had a FOURTH copy that no row here
listed, `jline.util.Maths.uniquePerms`, a line-by-line transliteration into
Java. When adding a row, grep for the algorithm under every naming convention
the codebases use, not only for the upstream file name.

| Library | Authors / group | License as shipped | In-tree paths |
|---|---|---|---|
| KPC-Toolbox | G. Casale, E. Smirni | BSD-3-Clause, `matlab/lib/kpctoolbox/LICENSE.TXT` | `matlab/lib/kpctoolbox/`, `jar/.../lib/kpctoolbox/`, `python/.../lib/kpctoolbox/` |
| M3A | A. Sansottera, G. Casale, P. Cremonesi | BSD-3-Clause, `matlab/lib/m3a/LICENSE` | `matlab/lib/m3a/`, `jar/.../lib/m3a/`, `python/.../lib/m3a/` |
| BuTools 2.0 | BuTools developers (webspn.hit.bme.hu) | Not stated in tree or upstream `ACTION` | `matlab/lib/thirdparty/BUTools/`, `jar/.../lib/butools/`, `python/.../lib/thirdparty/butools/` |
| SMCSolver (QBD, M/G/1, GI/M/1) | D. A. Bini, B. Meini, S. Steffe, B. Van Houdt | Not stated in tree `ACTION` | `matlab/lib/thirdparty/{smcsolver,QBDfiles,MG1files}/`, `jar/.../lib/smc/`, `python/.../lib/thirdparty/smc/` |
| Q-MAM | B. Van Houdt, J. F. Perez, J. Van Velthoven | Apache-2.0, `matlab/lib/thirdparty/QMAM/README.md` and upstream | `matlab/lib/thirdparty/QMAM/`, `jar/.../lib/qmam/`, `python/.../lib/thirdparty/qmam/` |
| MAMSolver (ETAQA) | A. Riska, E. Smirni | Not stated in tree `ACTION` | `matlab/lib/thirdparty/MAMSolver/` |
| FJ_codes | Z. Qiu, J. F. Perez, P. Harrison | BSD-3-Clause, ©2015 Imperial College London, `python/.../lib/thirdparty/fj/LICENSE.txt` | `matlab/lib/thirdparty/FJ_codes/`, `jar/.../lib/fjcodes/`, `python/.../lib/thirdparty/fj/`, `cpp/include/line/api/fj/` |
| aoi-fluid toolbox | O. Dogan, N. Akar, E. U. Atay | BSD-2-Clause, source headers | `matlab/lib/aoi/`, `matlab/lib/thirdparty/aoi/`, `python/.../lib/thirdparty/aoi/` |
| rmf_tool | N. Gast, E. Rodriguez (Inria) | MIT, `python/.../lib/rmf_tool/LICENSE` | `python/.../lib/rmf_tool/`, `jar/.../lib/rmf/` |
| ILT-CME inverse Laplace coefficients | CME inverse-Laplace project | Not stated `ACTION`; **vendored 2026-07-25 with terms unstated, by explicit user decision** | `matlab/lib/thirdparty/iltcme/`, `python/.../lib/thirdparty/iltcme/`, `cpp/src/api/mam/iltcme_table.cpp` |
| libQBD | S. Astaf'ev, IAMR Karelian Research Centre RAS | BSD-3-Clause, `matlab/lib/thirdparty/libQBD/LICENSE` | `matlab/lib/thirdparty/libQBD/` |
| LSODA (ODE integrator port) | S. Frost | MIT, `matlab/lib/thirdparty/lsoda/LICENSE` | `matlab/lib/thirdparty/lsoda/`, `python/.../lib/lsoda.py`; see also `lsoda-java` in section 4 |
| freeLYAP and Lyapunov solvers | A. Townsend; Q. Fang, Nedialko, The MathWorks Inc., F. Glineur | MIT, `.../lyap/lyap.LICENSE.txt`; BSD-3-Clause, `.../lyap/license.txt` | `matlab/lib/thirdparty/lyap/` |
| QP (quadratic programming) | A. Barraud | BSD-2-Clause, `.../QP/license.txt` | `matlab/lib/m3a/lib/thirdparty/QP/` |
| MAPMsG | O. Gursoy | BSD-3-Clause, `matlab/lib/thirdparty/MAPMsG/license.txt` | `matlab/lib/thirdparty/MAPMsG/`, `python/.../lib/thirdparty/mapmsg/` |
| Hurst estimators | C. Chen | BSD-2-Clause, `.../hurst_estimators/license.txt` | `matlab/lib/thirdparty/hurst_estimators/`, `python/.../lib/thirdparty/hurst_estimators/` |
| graph_cc | M. Vedenyov | BSD-3-Clause, `.../graph_cc/license.txt` | `matlab/lib/thirdparty/graph_cc/` |
| xml_read / xml_write | J. Tuszynski | BSD-2-Clause, `.../xml_read/LICENSE.txt` | `matlab/lib/thirdparty/xml_read/` |

`ACTION` notes: for **BuTools, SMCSolver, MAMSolver** the
in-tree copies carry no license grant (BuTools' upstream declares none either,
checked 2026-08-16), so terms are genuinely UNKNOWN — obtain a written statement
upstream and add the license file; nothing in the tree may assert otherwise. For
**ILT-CME** (vendored 2026-07-25 on explicit user instruction with terms known to
be unstated): a dead gitlink was replaced with the three real files, and only the
127 reachable table entries were vendored into the C++ port; obtaining written
terms remains open. **Q-MAM** is RESOLVED (Apache-2.0); the stray BSD-2-Clause
`QMAM/license.txt` names an unrelated author (C. Greene) and is not Q-MAM's terms.

## 4. Java dependencies

Declared in `jar/pom.xml`; under `-P b` they are shaded into `common/jline.jar`
(redistributed in binary form). Licenses read from POM metadata / `META-INF`.

| Artifact | Version | License |
|---|---|---|
| org.apache.commons:commons-math3 | 3.6.1 | Apache-2.0 |
| org.apache.commons:commons-lang3 | 3.14.0 | Apache-2.0 |
| commons-io:commons-io | 2.11.0 | Apache-2.0 |
| commons-cli:commons-cli | 1.8.0 | Apache-2.0 |
| us.hebi.matlab.mat:mfl-core, mfl-ejml | 0.5.15 | Apache-2.0 |
| org.ejml:ejml-all | 0.41 | Apache-2.0 |
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
| lsoda-java:lsoda | 1.0 | MIT, see section 3; https://github.com/imperial-qore/lsoda-java |
| org.junit.jupiter:junit-jupiter-* | 5.10.1 | EPL-2.0, test scope, not redistributed |

## 5. Python dependencies

Declared in `python/pyproject.toml`, installed by the user's package manager (not
redistributed here). Every entry is permissive; keep it that way (see section 0).

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

## 6. Citing the tools LINE relies on

Where a paper's result was produced through a wrapper solver or embedded library,
cite the upstream work as well as LINE. Each wrapper prints its canonical
citation on first use, also available machine-readable: `line_citation('JMT')`
(MATLAB), `jline.io.InputOutput.line_citation("JMT")` (Java),
`line_solver.api.io.line_citation('JMT')` (Python). `'JMT'`, `'LQNS'`, `'QNS'`
are recognised (`'QNS'` returns the LQNS reference); keys match
`doc/latex/biblio.bib` and `BIBLIOGRAPHY.md`, where the section-3 library
references are collected.

## 7. Reporting an attribution error

If you author software listed here and the attribution, license, or
redistribution status is wrong, this is a defect and will be corrected. Open an
issue on the LINE repository or contact the maintainers via http://line-solver.sf.net.
