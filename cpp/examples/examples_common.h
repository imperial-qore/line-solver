/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_EXAMPLES_COMMON_H
#define LINE_EXAMPLES_COMMON_H

/**
 * The example gallery of `matlab/examples` and `python/examples`, in C++.
 *
 * ONE BINARY, NOT ONE PER EXAMPLE. A MATLAB or Python example is a script the
 * user runs by name; the C++ counterpart is a registered function the
 * `line-examples` binary runs by the SAME name, so `line-examples gallery_mm1`
 * is `gallery_mm1.py`. Each source file carries one directory of the reference
 * tree rather than one example, because a translation unit here instantiates
 * the solver templates it calls and 300 of them would cost more to compile than
 * the whole library.
 *
 * WHAT AN EXAMPLE PRINTS is the reference's table, column for column, and the
 * rendering is `src/cli/line_cli.cpp`'s so the two agree: an example and
 * `line-cli -a avg` on the same model must not differ by a space.
 *
 * WHERE C++ HAS NO SOLVER the example says so BY NAME through `na()` and moves
 * on. LDES, JMT, LQNS and QNS are the four the port does not carry, and an
 * example that answered a `SolverJMT` block with the MVA number would be a
 * silent solver substitution -- the defect the parity harness's `cpp_ported`
 * list exists to prevent (`line-test.git/parity/test_example.sh:126-140`).
 */

#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>
#include <string>
#include <vector>

#include "line/lang/qn/network_builder.h"
#include "line/num/number.h"
#include "line/util/error.h"

#include "example_util.h"
#include "parity_recorder.h"

/**
 * The repository root, so an example reads the SAME data file its MATLAB and
 * Python counterparts read instead of a copy that can drift from it. Defined by
 * the build; the fallback keeps the header usable in a hand-built translation
 * unit.
 */
#ifndef LINE_EXAMPLES_REPO_ROOT
#define LINE_EXAMPLES_REPO_ROOT "."
#endif

namespace line {
namespace examples {

using lang::Distrib;
using lang::SchedStrategy;
using lang::RoutingStrategy;

/**
 * The path of a data file under `cpp/examples/data/`, which is where every
 * repo-root read of a C++ example resolves.
 *
 * Its entries are SYMLINKS to the MATLAB and Python originals, so the tree
 * still holds exactly ONE copy of each trace and each `.lqnx` and no drift
 * between the codebases is possible, which is the rule the repo-root read was
 * introduced for. `zip` dereferences a symlink, so the C++-only release archive
 * gets the files as regular files under `cpp/examples/data/` and is
 * self-contained without shipping the `python/` and `matlab/` trees it has no
 * other use for.
 *
 * A new repo-root read in a C++ example means a new symlink here, and nothing
 * else: the release no longer carries a per-file include list that has to be
 * kept in step.
 */
inline std::string example_data_file(const std::string& name) {
    return std::string(LINE_EXAMPLES_REPO_ROOT) + "/cpp/examples/data/" + name;
}

/** The arithmetic every example runs in. `line-cli --arith` is the CLI's knob. */
using D = Distrib<double>;
using Net = qn::Network<double>;
using Routing = qn::RoutingMatrix<double>;
using Sn = qn::NetworkStruct<double>;

// ---------------------------------------------------------------------------
// Registry
// ---------------------------------------------------------------------------

/** One runnable example: the directory it came from, its name, its body. */
struct Example {
    std::string group;
    std::string name;
    std::function<void()> run;
};

inline std::vector<Example>& registry() {
    static std::vector<Example> reg;
    return reg;
}

struct Registrar {
    Registrar(const char* group, const char* name, std::function<void()> fn) {
        Example e;
        e.group = group;
        e.name = name;
        e.run = fn;
        registry().push_back(e);
    }
};

/** Register `fn` under its own name, in the reference directory `grp`. */
#define LINE_EXAMPLE(grp, fn) \
    static const ::line::examples::Registrar line_example_##fn(grp, #fn, fn)

// ---------------------------------------------------------------------------
// Printing
// ---------------------------------------------------------------------------

inline double to_d(double v) { return v; }

/** The banner a reference script prints before each solver block. */
/**
 * Declare which solver the tables below belong to.
 *
 * The printed banner is what the standalone scraper read; the recorder call
 * beside it is the same declaration made structurally, so a twin needs no edit
 * for its tables to become attributable. See parity_recorder.h.
 */
inline void section(const std::string& title) {
    std::printf("\nSOLVER: %s\n", title.c_str());
    line::examples::parity::set_solver(title);
}

/**
 * The same declaration where the BANNER and the GOLDEN'S KEY differ.
 *
 * A twin that solves the same model twice tells its reader which run is which
 * -- `CTMC (cutoff = 5)`, `NC (slotted)`, `MVA (method=sqd)` -- and that text
 * belongs in the output, because the reference prints it. It is not a solver
 * name: the goldens key that table under `CTMC`, and a record filed under the
 * decorated spelling aligns with nothing and reads as a solver that produced no
 * output. The Python recorder never has this problem because it labels from the
 * solver CLASS and not from anything printed.
 *
 * Use it ONLY where the decoration is presentation. Where the two runs are
 * genuinely different computations that the golden distinguishes -- `LN(NC)`
 * against `LN(moment3)` -- the qualified name IS the key and `section` is right.
 */
inline void section(const std::string& title, const std::string& golden_key) {
    std::printf("\nSOLVER: %s\n", title.c_str());
    line::examples::parity::set_solver(golden_key);
}

inline void note(const std::string& text) { std::printf("%s\n", text.c_str()); }

/**
 * Declare the solver the following tables belong to WITHOUT printing a banner.
 *
 * `section()` is the declaration a twin makes FOR ITS READER, and most of the
 * corpus has one. Some references print no banner at all -- each of the 57
 * `test_gallery_*` scripts is four lines that solve with MVA and print the
 * table -- and there the C++ twin must not invent one, because an example's
 * output is the reference's output. The recorder still needs the attribution:
 * without it `begin_table` has no solver to file the rows under and drops them,
 * so the row reads downstream as "solver MVA missing from the recorded results"
 * -- an unmeasured example rather than a measured one. This is the C++ twin of
 * what the Python recorder gets for free by wrapping the getter it was called
 * on.
 */
inline void attribute(const std::string& solver) {
    line::examples::parity::set_solver(solver);
}

/**
 * The solver this port does not carry, refused by name.
 *
 * A refusal is information; a substituted solver is not. `why` names what the
 * reference example asked for so a reader can tell a missing engine from a
 * model the C++ side cannot express.
 */
inline void na(const std::string& solver, const std::string& why) {
    std::printf("\nSOLVER: %s\nN/A: %s\n", solver.c_str(), why.c_str());
    // A REFUSAL IS RECORDED, NOT ONLY PRINTED. "this port does not carry MAM"
    // is a fact about the port and is a named skip downstream; a table that
    // simply never arrived is a failure. The two must not look alike.
    line::examples::parity::note_refusal(solver, why);
}

// ---------------------------------------------------------------------------
// Layered (LQN) result tables
// ---------------------------------------------------------------------------

/** One row of a layered result table: the element and its five measures. */
struct LnRow {
    std::string name;
    double q, u, r, w, t;  ///< NaN where the measure does not apply to it
};

/**
 * The golden's key for a layered solve: `LN(<layer engine>)`.
 *
 * AN ENSEMBLE'S MEMBER IS PART OF ITS IDENTITY, not decoration: MVA layers and
 * NC layers are different fixed points (`lcq_threehosts` converges to a cache
 * hit of 0.5 against 0.48331), and `lqn_twotasks` is goldened under the bare
 * `NC` for exactly that reason. The comparator reconciles `LN(NC)` against a
 * golden keyed `NC`, `LN(NC)` or plain `LN`; recording the member is what gives
 * it the evidence to do so safely.
 *
 * THE ENSEMBLE'S OWN METHOD WINS WHEN IT HAS ONE. `lqn_moment3` solves the same
 * model with the same NC layers twice, default and `moment3`, and its golden
 * holds the FIRST -- the two differ by 360x on T1, so labelling both `LN(NC)`
 * would let the second overwrite the first and report the default solve as a
 * 99.7% error.
 */
inline std::string ln_solver_key(const std::string& layer_engine, const std::string& method) {
    if (!method.empty() && method != "default") return "LN(" + method + ")";
    return "LN(" + layer_engine + ")";
}

/**
 * Record a layered result table under the golden's own (Station, JobClass) key.
 *
 * A LAYERED TABLE DOES NOT GO THROUGH `avg_rows`, which walks the stations of a
 * `NetworkStruct` that a layered solve does not have: it is keyed by ELEMENT
 * NAME, with the single class spelled `Jobs`, which is what `getAvgTable` prints
 * for a layered model and what the `LN` goldens carry.
 *
 * AN UNDEFINED CELL IS RECORDED AS NaN rather than dropped. NaN is what the
 * golden holds for a measure that does not apply to an element -- a Processor
 * has no response time -- and a dropped cell would instead read as a row the
 * twin failed to produce.
 */
inline void record_ln(const std::string& solver, const std::vector<LnRow>& rows) {
    namespace parity = line::examples::parity;
    if (!parity::enabled()) return;
    parity::set_solver(solver);
    parity::begin_table("avg", "Station", "JobClass");
    for (std::size_t i = 0; i < rows.size(); ++i) {
        std::vector<parity::Cell> cells;
        cells.push_back(parity::Cell{"QLen", rows[i].q});
        cells.push_back(parity::Cell{"Util", rows[i].u});
        cells.push_back(parity::Cell{"RespT", rows[i].r});
        cells.push_back(parity::Cell{"ResidT", rows[i].w});
        cells.push_back(parity::Cell{"Tput", rows[i].t});
        parity::add_row(rows[i].name, "Jobs", cells);
    }
}

/**
 * Record a quantity the EXAMPLE derived, under the golden's own key.
 *
 * The short spelling of `parity::add_derived`; see parity_recorder.h for why a
 * derived golden needs a key the example chooses rather than the (station,
 * class) pair a result table has. Costs nothing when `--record` was not asked
 * for.
 */
inline void derived(const std::string& solver, const std::string& row,
                    const std::string& col, double value,
                    const std::string& metric = "QLen") {
    line::examples::parity::add_derived(solver, row, col, value, metric);
}

/**
 * `derived` filed under the solver `section()` last declared.
 *
 * A quantity an example computes FROM a solver's answer belongs to that solver
 * -- a reward expectation is the CTMC's -- and naming it by hand in the printer
 * would be a second place for the two to disagree.
 */
inline void derived_here(const std::string& row, const std::string& col, double value,
                         const std::string& metric = "QLen") {
    const std::string& solver = line::examples::parity::current_solver();
    if (!solver.empty())
        line::examples::parity::add_derived(solver, row, col, value, metric);
}

/**
 * Record a CDF-derived statistic matrix under the golden's own keys.
 *
 * THE ROW NUMBER IS NOT THE STATION INDEX. These goldens key their rows
 * `Station1_AvgRespT`, `Station2_SCV` and so on, numbering only the stations
 * that MEASURED something: the reference loops over every station and stores a
 * zero where the law is empty (a Source has no response time of its own), and
 * the reader that generated the golden dropped those rows and renumbered what
 * was left. `values` is therefore the matrix AS PRINTED -- every station, zeros
 * included -- and the drop happens here, so a twin cannot shift `Station2` onto
 * `Station1` by filtering earlier.
 *
 * A cell that is zero or NaN is skipped for the same reason: it is a law the
 * run did not measure, not a measurement of zero.
 */
inline void record_cdf_matrix(const std::string& solver, bool scv,
                              const std::vector<std::vector<double> >& values) {
    if (!line::examples::parity::enabled()) return;
    const std::string metric = scv ? "_SCV" : "_AvgRespT";
    std::size_t row = 0;
    for (std::size_t i = 0; i < values.size(); ++i) {
        bool measured = false;
        for (std::size_t c = 0; c < values[i].size(); ++c)
            if (std::fabs(values[i][c]) > 1e-10) measured = true;
        if (!measured) continue;
        ++row;
        for (std::size_t c = 0; c < values[i].size(); ++c) {
            const double v = values[i][c];
            if (std::fabs(v) < 1e-10 || std::isnan(v)) continue;
            line::examples::parity::add_derived(
                solver, "Station" + std::to_string(row) + metric,
                "Class" + std::to_string(c + 1), v);
        }
    }
}

inline void avg_header() {
    std::printf("%-16s %-14s %12s %12s %12s %12s %12s %12s\n", "Station", "JobClass", "QLen",
                "Util", "RespT", "ResidT", "ArvR", "Tput");
}

/**
 * The AvgTable rows, dropping the (station, class) pairs the class never
 * visits exactly as `getAvgTable` does and as `emit_avg_table` does in the CLI.
 */
template <class T, class Get>
void avg_rows(const qn::NetworkStruct<T>& sn, Get get) {
    avg_header();
    // THE ONE PLACE EVERY AVERAGE TABLE IS MATERIALISED, which is why the
    // recorder hooks here and not at each of the hundred-odd call sites: the
    // printed row and the recorded row are then the same row by construction,
    // and a twin that gains a table gains a recorded table with it.
    namespace parity = line::examples::parity;
    if (parity::enabled()) parity::begin_table("avg", "Station", "JobClass");
    for (std::size_t i = 0; i < sn.nstations; ++i)
        for (std::size_t c = 0; c < sn.nclasses; ++c) {
            const double q = get(i, c, 0), u = get(i, c, 1), r = get(i, c, 2);
            const double w = get(i, c, 3), a = get(i, c, 4), t = get(i, c, 5);
            if (q == 0.0 && u == 0.0 && r == 0.0 && w == 0.0 && a == 0.0 && t == 0.0) continue;
            // FIVE SIGNIFICANT DIGITS, WHICH IS WHAT EVERY OTHER CODEBASE PRINTS.
            // MATLAB's table, pandas' repr and the JAR's printer all render a cell
            // at five, so that is the precision the shared goldens carry (0.76923,
            // 0.45505, 9.6499) and the precision a harness that PARSES this table
            // compares against. Printing six was not more agreement but less:
            // every cell then differed from its golden in the sixth digit, i.e.
            // by up to 6.5e-6 relative, which is above the 1e-6 a deterministic
            // solver is held to -- so a correct row failed on rounding alone
            // (ag_gnetwork: 11 AG cells "mismatched", every one agreeing to 10
            // digits). The CLI row already avoided this by re-rendering its
            // full-precision `-o json` at five (runner.CPP_TABLE_SIGFIGS); the
            // twin has no JSON to re-render, so it prints the shared precision.
            std::printf("%-16s %-14s %12.5g %12.5g %12.5g %12.5g %12.5g %12.5g\n",
                        sn.stations[i].name.c_str(), sn.classes[c].name.c_str(), q, u, r, w, a, t);
            if (!parity::enabled()) continue;
            // Recorded at FULL precision, unlike the printed row: quantizing to
            // the golden's precision is the comparator's job, and doing it here
            // would throw away the digits a full-precision golden would need.
            std::vector<parity::Cell> cells;
            cells.push_back(parity::Cell{"QLen", q});
            cells.push_back(parity::Cell{"Util", u});
            cells.push_back(parity::Cell{"RespT", r});
            cells.push_back(parity::Cell{"ResidT", w});
            cells.push_back(parity::Cell{"ArvR", a});
            cells.push_back(parity::Cell{"Tput", t});
            parity::add_row(sn.stations[i].name, sn.classes[c].name, cells);
        }
}

/**
 * `getAvgTable` for any solver returning `mva::AvgResult` -- MVA, NC, MAM, BA,
 * CTMC. Templated on the solution type so this header needs no solver include.
 */
template <class T, class Res>
void print_avg(const qn::NetworkStruct<T>& sn, const Res& r) {
    if (!r.warning.empty()) std::fprintf(stderr, "Warning: %s\n", r.warning.c_str());
    avg_rows(sn, [&](std::size_t i, std::size_t c, int k) {
        switch (k) {
            case 0: return num_traits<T>::to_double(r.QN(i, c));
            case 1: return num_traits<T>::to_double(r.UN(i, c));
            case 2: return num_traits<T>::to_double(r.RN(i, c));
            case 3: return num_traits<T>::to_double(r.WN(i, c));
            case 4: return num_traits<T>::to_double(r.AN(i, c));
            default: return num_traits<T>::to_double(r.TN(i, c));
        }
    });
}

/**
 * `getAvgTable` for the fluid and simulation solutions, which carry no separate
 * residence-time or arrival-rate matrix.
 *
 * NEITHER COLUMN IS A COPY OF ITS NEIGHBOUR, which is what this printer used to
 * make them. ResidT is the per-JOB residence time and RespT the per-VISIT
 * response time; they agree only where every station is visited once per cycle,
 * and reporting one for the other was a factor of 3 out on `sdroute_closed` and
 * 17 on Queue1 of `init_state_ps`. ArvR is a FLOW and parts from the throughput
 * at any station a job leaves by another route. `residt_from_respt` and
 * `arvr_from_tput` are the reference's own conversions and are what the CLI's
 * `-s fluid` and `-s ssa` arms apply, so an example and `line-cli -a avg` on the
 * same model still agree column for column. A Source keeps ArvR 0: it has no
 * arrivals to itself.
 *
 * `wsn` and `asn` override the struct each conversion reads its VISIT RATIOS
 * from, and default to `sn`. They exist for the models whose routing is a
 * RESULT: a cache's hit/miss split is what the analyzer converged to, not what
 * `link()` offered, and taking the offered one on `tut06_cache_lru_zipf` gave
 * both classes exactly half their response time -- a number that is not the
 * residence time of any model. This is the same distinction the CLI draws
 * between `sn`, its `refreshed` struct and `sn_with_ssa_cache_split`.
 */
template <class T, class Res>
void print_avg_sim(const qn::NetworkStruct<T>& sn, const Res& r,
                   const qn::NetworkStruct<T>* wsn = nullptr,
                   const qn::NetworkStruct<T>* asn = nullptr) {
    const std::size_t M = sn.nstations, K = sn.nclasses;
    Matrix<double> RN(M, K), TN(M, K);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t c = 0; c < K; ++c) {
            RN(i, c) = num_traits<T>::to_double(r.RN(i, c));
            TN(i, c) = num_traits<T>::to_double(r.TN(i, c));
        }
    const Matrix<double> WN = residt_from_respt(wsn ? *wsn : sn, RN);
    const Matrix<double> AN = arvr_from_tput(asn ? *asn : sn, TN);
    avg_rows(sn, [&](std::size_t i, std::size_t c, int k) {
        switch (k) {
            case 0: return num_traits<T>::to_double(r.QN(i, c));
            case 1: return num_traits<T>::to_double(r.UN(i, c));
            case 2: return RN(i, c);
            case 3: return WN(i, c);
            case 4: return sn.stations[i].sched == SchedStrategy::EXT ? 0.0 : AN(i, c);
            default: return TN(i, c);
        }
    });
}

/**
 * The samples of a trace file, one per line, as `Replayer(filename)` reads
 * them. A missing file is an error and not an empty trace: an empty Replayer
 * would silently become a zero-mean arrival process.
 */
inline std::vector<double> read_trace(const std::string& path) {
    std::ifstream in(path.c_str());
    if (!in) throw InputError("Replayer: cannot open the trace file '" + path + "'");
    std::vector<double> out;
    double v = 0.0;
    while (in >> v) out.push_back(v);
    if (out.empty()) throw InputError("Replayer: the trace file '" + path + "' has no samples");
    return out;
}

/** One labelled scalar, the shape a reference script prints a system metric in. */
inline void kv(const std::string& key, double value) {
    std::printf("%-28s %.6g\n", (key + ":").c_str(), value);
}

inline void kv(const std::string& key, const std::string& value) {
    std::printf("%-28s %s\n", (key + ":").c_str(), value.c_str());
}

/**
 * The scalar ON A LINE OF ITS OWN, beside the labelled one.
 *
 * A reference script that has no table to print puts its result bare --
 * `print(pr_ctmc)` in the Python twin, `Pmarg_ctmc =` and then the value in the
 * MATLAB one -- and the shared parity parser scrapes exactly that: a labelled
 * `key: value` line is neither a table nor a bare scalar, so it is read by
 * nothing. Without this the C++ twin computes the right number and the row
 * still compares no cell. Call it once per example, for the quantity the
 * reference prints bare, since the parser keeps the FIRST bare scalar it sees.
 *
 * Printed at 17 significant digits, the round-trip width of a double: the
 * golden holds the reference's full precision (0.0003484356916108198) and the
 * 6-digit `kv` form would compare a rounded number against it.
 */
inline void bare(double value) {
    std::printf("%.17g\n", value);
    // THE BARE SCALAR IS WHAT AN `OPT` GOLDEN HOLDS. Six `statepr_*` goldens
    // carry one row, ('OptResult', 'Value'), and it is this number: the reader
    // that generated them took the first bare value an example printed and
    // stopped, which is why the golden holds ONE probability where the example
    // computes two or three. Recording it here makes it attributable instead of
    // scraped, and the once-per-example rule above is what keeps it the same
    // value the golden was written from.
    line::examples::parity::add_derived("OPT", "OptResult", "Value", value);
}

// ---------------------------------------------------------------------------
// Model-building helpers
// ---------------------------------------------------------------------------

/**
 * `Network.serialRouting(n1, n2, ...)` for one class: each node routes to the
 * next with probability one. The class-indexed form is the `(r, r)` block, which
 * is what `P[jobclass] = Network.serial_routing(...)` writes in Python.
 */
inline void serial(Routing& P, const std::vector<std::size_t>& nodes) {
    for (std::size_t i = 0; i + 1 < nodes.size(); ++i) P.set(nodes[i], nodes[i + 1], 1.0);
}

inline void serial(Routing& P, std::size_t cls, const std::vector<std::size_t>& nodes) {
    for (std::size_t i = 0; i + 1 < nodes.size(); ++i) P.set(cls, cls, nodes[i], nodes[i + 1], 1.0);
}

/** The same, closed into a cycle: the last node routes back to the first. */
inline void cyclic(Routing& P, const std::vector<std::size_t>& nodes) {
    serial(P, nodes);
    if (nodes.size() > 1) P.set(nodes.back(), nodes.front(), 1.0);
}

inline void cyclic(Routing& P, std::size_t cls, const std::vector<std::size_t>& nodes) {
    serial(P, cls, nodes);
    if (nodes.size() > 1) P.set(cls, cls, nodes.back(), nodes.front(), 1.0);
}

}  // namespace examples
}  // namespace line

#endif  // LINE_EXAMPLES_COMMON_H
