/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `python/examples/basic/layeredModel/`: the layered-network gallery.
 *
 * Nineteen scripts, every one of them a LayeredNetwork solved by SolverLN. Two
 * of them read a `.lqnx`/`.xml` from the repository (`lqn_init`, `lqn_ofbiz`);
 * the rest are built with `LqnBuilder`, whose declaration order is the
 * reference's so that the element indices agree and the AvgTable rows line up.
 *
 * THE TABLE IS THE CLI'S. `print_ln` below prints exactly what
 * `src/cli/line_cli.cpp`'s `-i lqnx` arm prints -- the same banner naming the
 * layer engine and the method, the same seven columns, the same `ln_sanitize`
 * clean-up -- so an example and `line-cli model.lqnx` cannot disagree by a
 * space. The one field the banner drops is the wall clock, which would make the
 * example output irreproducible; where the reference itself times the solve
 * (`lqn_ofbiz`, `lqn_init`) the elapsed seconds are printed as the reference
 * prints them.
 *
 * WHAT IS REFUSED, and why each refusal is not a substitution. Every refusal
 * carries the reference's own call as a `TODO(cpp)` comment beside it, so the
 * source records what the reference runs and the gap is greppable, AND prints
 * `na()` at run time, so a reader of the output cannot mistake silence for
 * coverage (`examples/README.md`, "The four solvers this port does not carry"):
 *
 *  - Every `LQNS(...)` block. `lqns` and `lqsim` are external binaries this port
 *    does not carry. SolverLN is a different algorithm, not a stand-in.
 *  - `lqn_rrobin` and `lqn_jsq`. Their calls are ROUTED call groups
 *    (`synch_call_rrobin`, `synch_call_jsq`) solved under the squashed 'flat'
 *    layering; `LqnBuilder` has neither construct and `LnOptions` has no
 *    layering knob, so the model cannot be expressed at all. Building the
 *    probabilistic twin instead would answer a round-robin question with a
 *    random-branching number.
 *  - The monkey-patched `solver.pre` / `solver.post` / `solver.analyze` hooks of
 *    the three `*_debug` scripts. `SolverLN` exposes its converged state and its
 *    per-iteration layer results, and those are printed; it does not let a
 *    caller run code BETWEEN two layer solves, so the mid-iteration snapshots
 *    are refused rather than approximated by the converged values.
 *
 * `SolverLN(SolverNC)` is NOT refused: `LnOptions::layer_solver` takes `nc`
 * (also `fluid` and `ssa`), so the scripts that ask for NC layers get NC layers.
 *
 * TWO KNOWN GAPS, left to fail loudly rather than papered over:
 *
 *  - `lqn_setup` throws. A setup layer is routed to `solve_layer_mam`, and
 *    `solver_mam_basic` refuses a class-switching chain at an FCFS station, so
 *    the model only solves when the setup task's HOST is PS (which is what
 *    tests/test_ln_setup.cpp builds). The reference's P2 is FCFS. Python solves
 *    it over MVA layers and reports P1 Util 0.43577, P2 Util 0.036314, E1 RespT
 *    2.2948, E2 RespT 1.2948, Tput 0.43577.
 *  - `lqn_bpmn` reaches `iter_max` without converging. The
 *    Python reference does not get an answer on this model at all (its LN(MVA)
 *    run raises from the layer solver's product-form check), so there is no
 *    reference row to compare the fixed point against yet.
 */

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <exception>
#include <limits>
#include <string>
#include <vector>

#include "examples_common.h"
#include "line/api/infer/infer_lqn_ekf.h"
#include "line/api/infer/infer_lqn_getobs.h"
#include "line/lang/dist_fitters.h"
#include "line/lang/lqn/lqn_builder.h"
#include "line/lang/lqn/lqn_reader.h"
#include "line/solvers/fluid/fluid_runner.h"
#include "line/solvers/ln/solver_ln.h"
#include "line/util/method_type.h"

namespace line {
namespace examples {

namespace {

using Lqn = lqn::LqnBuilder<double>;
using LqnS = lqn::LqnStruct<double>;

/** `Exp(rate)` and `Exp.fitMean(mean)`, the two spellings the scripts use. */
D Er(double rate) { return Exp(rate); }
D Em(double mean) { return D::exp_mean(mean); }

/** The path of a data file of this reference directory. */
std::string layered_file(const std::string& name) {
    return example_data_file("layeredModel/" + name);
}

// ---------------------------------------------------------------------------
// The LQN AvgTable, in the CLI's layout
// ---------------------------------------------------------------------------

/** getAvgTable's numeric clean-up, verbatim from the CLI's `ln_sanitize`. */
double ln_sanitize(double x) {
    const double r = std::round(x * 10.0);
    if (std::fabs(x * 10.0 - r) < lang::GlobalConstants::CoarseTol * x * 10.0) x = r / 10.0;
    if (x <= lang::GlobalConstants::FineTol) x = 0.0;
    return x;
}

const char* ln_element_kind(const LqnS& l, std::size_t i) {
    switch (l.type[i]) {
        case lang::LqnElement::HOST: return "Processor";
        case lang::LqnElement::TASK: return l.isref[i] ? "RefTask" : "Task";
        case lang::LqnElement::ENTRY: return "Entry";
        default: return "Activity";
    }
}

const char* ln_layer_engine_name(const std::string& layer_solver) {
    if (layer_solver == "fluid") return "Fluid";
    if (layer_solver == "nc") return "NC";
    if (layer_solver == "ssa") return "SSA";
    return "MVA";
}

/** The golden's key for the layer engine `opt.layer_solver` names. */
const char* ln_layer_golden_key(const std::string& layer_solver) {
    if (layer_solver == "fluid") return "FLD";
    if (layer_solver == "nc") return "NC";
    if (layer_solver == "ssa") return "SSA";
    return "MVA";
}

/** The banner and the seven columns of `line-cli -i lqnx`, minus the wall clock. */
void print_ln(const LqnS& l, const ln::LnSolution<double>& sol, const ln::LnOptions& opt,
              std::size_t nlayers) {
    std::printf("SolverLN(Solver%s) arith=double type=%s layers=%zu iterations=%d converged=%d\n",
                ln_layer_engine_name(opt.layer_solver),
                ::line::util::method_type("LN", opt.method).c_str(), nlayers, sol.iterations,
                int(sol.converged));
    std::printf("%-62s %-10s %12s %12s %12s %12s %12s\n", "Node", "NodeType", "QLen", "Util",
                "RespT", "ResidT", "Tput");
    // THE PRINTED ROW AND THE RECORDED ROW ARE THE SAME ROW, sanitized value for
    // sanitized value: the goldens were generated from this table, so recording
    // the raw solution instead would compare a number the reference never
    // showed. An undefined measure is NaN in both, which is what the golden
    // holds for it.
    std::vector<LnRow> recorded;
    for (std::size_t i = 1; i <= l.nidx; ++i) {
        auto cell = [&](const std::vector<double>& v, const std::vector<bool>& d) {
            return d[i] ? ln_sanitize(v[i]) : std::numeric_limits<double>::quiet_NaN();
        };
        LnRow row;
        row.name = l.names[i];
        row.q = cell(sol.QN, sol.defined_Q);
        row.u = cell(sol.UN, sol.defined_U);
        row.r = cell(sol.RN, sol.defined_R);
        row.w = cell(sol.WN, sol.defined_W);
        row.t = cell(sol.TN, sol.defined_T);
        auto fmt = [](double v, char* buf) {
            if (std::isnan(v)) std::snprintf(buf, 24, "%12s", "NaN");
            else std::snprintf(buf, 24, "%12.6g", v);
        };
        char q[24], u[24], rr[24], w[24], t[24];
        fmt(row.q, q);
        fmt(row.u, u);
        fmt(row.r, rr);
        fmt(row.w, w);
        fmt(row.t, t);
        std::printf("%-62s %-10s %s %s %s %s %s\n", l.names[i].c_str(), ln_element_kind(l, i), q, u,
                    rr, w, t);
        recorded.push_back(row);
    }
    record_ln(ln_solver_key(ln_layer_golden_key(opt.layer_solver), opt.method), recorded);
}

/** Solve and print in one step, the shape most of these scripts have. */
void solve_ln(const LqnS& l, const ln::LnOptions& opt = ln::LnOptions()) {
    ln::SolverLN<double> s(l, opt);
    const ln::LnSolution<double> sol = s.get_ensemble_avg();
    print_ln(l, sol, opt, s.nlayers());
}

/** `LN(model, @(x) NC(x))`: the same fixed point over NC layers. */
ln::LnOptions nc_layers() {
    ln::LnOptions opt;
    opt.layer_solver = "nc";
    return opt;
}

/** The layer-structure dump of `line-cli -o layers`, so the two agree. */
void dump_layers(const ln::SolverLN<double>& solver) {
    const std::vector<qn::Layer<double> >& ens = solver.layers();
    for (std::size_t k = 0; k < ens.size(); ++k) {
        const qn::Layer<double>& L = ens[k];
        std::printf("LAYER %zu %s nstations=%zu nclasses=%zu nchains=%zu\n", k + 1, L.name.c_str(),
                    L.stations.size(), L.classes.size(), L.nchains);
        for (std::size_t i = 0; i < L.stations.size(); ++i)
            std::printf("  STATION %zu %s sched=%s nservers=%g\n", i + 1, L.stations[i].name.c_str(),
                        lang::sched_to_text(L.stations[i].sched), L.stations[i].nservers);
        for (std::size_t r = 0; r < L.classes.size(); ++r)
            std::printf("  CLASS %zu %s pop=%.17g refstat=%zu completes=%d\n", r + 1,
                        L.classes[r].name.c_str(), L.classes[r].population, L.classes[r].refstat,
                        int(L.classes[r].completes));
        for (std::size_t i = 0; i < L.stations.size(); ++i)
            for (std::size_t r = 0; r < L.classes.size(); ++r) {
                if (L.disabled.empty() || L.disabled[i][r]) continue;
                std::printf("  RATE %s %s %.17g scv=%.17g\n", L.stations[i].name.c_str(),
                            L.classes[r].name.c_str(), L.rates(i, r), L.scv(i, r));
            }
        for (const auto& kv : L.P) {
            const Matrix<double>& B = kv.second;
            for (std::size_t i = 0; i < B.rows(); ++i)
                for (std::size_t j = 0; j < B.cols(); ++j) {
                    if (B(i, j) == 0.0) continue;
                    std::printf("  ROUTE %s->%s %s->%s %.17g\n",
                                L.classes[kv.first.first - 1].name.c_str(),
                                L.classes[kv.first.second - 1].name.c_str(), L.nodes[i].name.c_str(),
                                L.nodes[j].name.c_str(), B(i, j));
                }
        }
    }
}

/** 1-based element index of a hashname (`T:T1`, `A:AS1`, ...); 0 when absent. */
std::size_t idx_of(const LqnS& l, const std::string& hashname) {
    for (std::size_t i = 1; i <= l.nidx; ++i)
        if (l.hashnames[i] == hashname) return i;
    return 0;
}

/** The largest entry of a per-iteration layer matrix, NaN-free as numpy's nanmax. */
double mat_max(const Matrix<double>& M) {
    double best = 0.0;
    bool seen = false;
    for (std::size_t i = 0; i < M.rows(); ++i)
        for (std::size_t j = 0; j < M.cols(); ++j) {
            const double v = M(i, j);
            if (std::isnan(v)) continue;
            if (!seen || v > best) {
                best = v;
                seen = true;
            }
        }
    return seen ? best : 0.0;
}

// ---------------------------------------------------------------------------
// Model factories shared by an example and its debug twin
// ---------------------------------------------------------------------------

/** `test_LQN_4`: two hosts, three tasks, a two-deep synchronous call chain. */
LqnS lqn_basic_model() {
    Lqn b;
    b.processor("P1", 2, SchedStrategy::PS);
    b.processor("P2", 3, SchedStrategy::PS);
    b.task("T1", 50, SchedStrategy::REF, "P1");
    b.think_time("T1", Er(1.0 / 2.0));
    b.task("T2", 50, SchedStrategy::FCFS, "P1");
    b.think_time("T2", Er(1.0 / 3.0));
    b.task("T3", 25, SchedStrategy::FCFS, "P2");
    b.think_time("T3", Er(1.0 / 4.0));
    b.entry("E1", "T1");
    b.entry("E2", "T2");
    b.entry("E3", "T3");
    b.activity("AS1", Er(10), "T1");
    b.bound_to("AS1", "E1");
    b.sync_call("AS1", "E2", 1.0);
    b.activity("AS2", Er(20), "T2");
    b.bound_to("AS2", "E2");
    b.sync_call("AS2", "E3", 5.0);
    b.replies_to("AS2", "E2");
    b.activity("AS3", Er(50), "T3");
    b.bound_to("AS3", "E3");
    b.replies_to("AS3", "E3");
    return b.build();
}

/** The BPMN-derived model: 7 hosts, 7 tasks, 16 entries, 29 activities. */
LqnS lqn_bpmn_model() {
    Lqn b;
    b.processor("R1_Processor", 100, SchedStrategy::FCFS);
    b.processor("R2_Processor", lang::GlobalConstants::MaxInt, SchedStrategy::INF);
    b.processor("R3_Processor", 2, SchedStrategy::FCFS);
    b.processor("R1A_Processor", 7, SchedStrategy::FCFS);
    b.processor("R1B_Processor", 3, SchedStrategy::FCFS);
    b.processor("R2A_Processor", 4, SchedStrategy::FCFS);
    b.processor("R2B_Processor", 5, SchedStrategy::FCFS);

    b.task("R1_Task", 100, SchedStrategy::REF, "R1_Processor");
    b.think_time("R1_Task", Em(20));
    b.task("R2_Task", lang::GlobalConstants::MaxInt, SchedStrategy::INF, "R2_Processor");
    b.think_time("R2_Task", Immediate());
    b.task("R3_Task", 2, SchedStrategy::FCFS, "R3_Processor");
    b.think_time("R3_Task", Immediate());
    b.task("R1A_Task", 7, SchedStrategy::FCFS, "R1A_Processor");
    b.think_time("R1A_Task", Immediate());
    b.task("R1B_Task", 3, SchedStrategy::FCFS, "R1B_Processor");
    b.think_time("R1B_Task", Immediate());
    b.task("R2A_Task", 4, SchedStrategy::FCFS, "R2A_Processor");
    b.think_time("R2A_Task", Immediate());
    b.task("R2B_Task", 5, SchedStrategy::FCFS, "R2B_Processor");
    b.think_time("R2B_Task", Immediate());

    b.entry("R1_Ref_Entry", "R1_Task");
    b.entry("R2_Synch_A2_Entry", "R2_Task");
    b.entry("R2_Synch_A5_Entry", "R2_Task");
    b.entry("R3_Synch_A9_Entry", "R3_Task");
    b.entry("R1A_Synch_A1_Entry", "R1A_Task");
    b.entry("R1A_Synch_A2_Entry", "R1A_Task");
    b.entry("R1A_Synch_A3_Entry", "R1A_Task");
    b.entry("R1B_Synch_A4_Entry", "R1B_Task");
    b.entry("R1B_Synch_A5_Entry", "R1B_Task");
    b.entry("R1B_Synch_A6_Entry", "R1B_Task");
    b.entry("R2A_Synch_A7_Entry", "R2A_Task");
    b.entry("R2A_Synch_A8_Entry", "R2A_Task");
    b.entry("R2A_Synch_A11_Entry", "R2A_Task");
    b.entry("R2B_Synch_A9_Entry", "R2B_Task");
    b.entry("R2B_Synch_A10_Entry", "R2B_Task");
    b.entry("R2B_Synch_A12_Entry", "R2B_Task");

    // T1
    b.activity("A1_Empty", Immediate(), "R1_Task");
    b.bound_to("A1_Empty", "R1_Ref_Entry");
    b.sync_call("A1_Empty", "R1A_Synch_A1_Entry", 1.0);
    b.activity("A2_Empty", Immediate(), "R1_Task");
    b.sync_call("A2_Empty", "R1A_Synch_A2_Entry", 1.0);
    b.activity("A5_Empty", Immediate(), "R1_Task");
    b.sync_call("A5_Empty", "R1B_Synch_A5_Entry", 1.0);
    b.activity("A6_Empty", Immediate(), "R1_Task");
    b.sync_call("A6_Empty", "R1B_Synch_A6_Entry", 1.0);
    b.activity("A3_Empty", Immediate(), "R1_Task");
    b.sync_call("A3_Empty", "R1A_Synch_A3_Entry", 1.0);
    b.activity("A4_Empty", Immediate(), "R1_Task");
    b.sync_call("A4_Empty", "R1B_Synch_A4_Entry", 1.0);

    // T2
    b.activity("E4_Empty", Immediate(), "R2_Task");
    b.bound_to("E4_Empty", "R2_Synch_A2_Entry");
    b.activity("A7_Empty", Immediate(), "R2_Task");
    b.sync_call("A7_Empty", "R2A_Synch_A7_Entry", 1.0);
    b.activity("A8_Empty", Immediate(), "R2_Task");
    b.sync_call("A8_Empty", "R2A_Synch_A8_Entry", 1.0);
    b.activity("A9_Empty", Immediate(), "R2_Task");
    b.sync_call("A9_Empty", "R2B_Synch_A9_Entry", 1.0);
    b.activity("A11_Empty", Immediate(), "R2_Task");
    b.sync_call("A11_Empty", "R2A_Synch_A11_Entry", 1.0);
    b.replies_to("A11_Empty", "R2_Synch_A2_Entry");
    b.activity("A12_Empty", Immediate(), "R2_Task");
    b.bound_to("A12_Empty", "R2_Synch_A5_Entry");
    b.sync_call("A12_Empty", "R2B_Synch_A12_Entry", 1.0);
    b.replies_to("A12_Empty", "R2_Synch_A5_Entry");
    b.activity("A10_Empty", Immediate(), "R2_Task");
    b.sync_call("A10_Empty", "R2B_Synch_A10_Entry", 1.0);

    // T3
    b.activity("A13", Em(10), "R3_Task");
    b.bound_to("A13", "R3_Synch_A9_Entry");
    b.replies_to("A13", "R3_Synch_A9_Entry");

    // T4
    b.activity("A1", Em(7), "R1A_Task");
    b.bound_to("A1", "R1A_Synch_A1_Entry");
    b.replies_to("A1", "R1A_Synch_A1_Entry");
    b.activity("A2", Em(4), "R1A_Task");
    b.bound_to("A2", "R1A_Synch_A2_Entry");
    b.activity("A3", Em(5), "R1A_Task");
    b.bound_to("A3", "R1A_Synch_A3_Entry");
    b.replies_to("A3", "R1A_Synch_A3_Entry");
    b.activity("A2_Res_Empty", Immediate(), "R1A_Task");
    b.sync_call("A2_Res_Empty", "R2_Synch_A2_Entry", 1.0);
    b.replies_to("A2_Res_Empty", "R1A_Synch_A2_Entry");

    // T5
    b.activity("A4", Em(8), "R1B_Task");
    b.bound_to("A4", "R1B_Synch_A4_Entry");
    b.replies_to("A4", "R1B_Synch_A4_Entry");
    b.activity("A5", Em(4), "R1B_Task");
    b.bound_to("A5", "R1B_Synch_A5_Entry");
    b.activity("A6", Em(6), "R1B_Task");
    b.bound_to("A6", "R1B_Synch_A6_Entry");
    b.replies_to("A6", "R1B_Synch_A6_Entry");
    b.activity("A5_Res_Empty", Immediate(), "R1B_Task");
    b.sync_call("A5_Res_Empty", "R2_Synch_A5_Entry", 1.0);
    b.replies_to("A5_Res_Empty", "R1B_Synch_A5_Entry");

    // T6
    b.activity("A7", Em(6), "R2A_Task");
    b.bound_to("A7", "R2A_Synch_A7_Entry");
    b.replies_to("A7", "R2A_Synch_A7_Entry");
    b.activity("A8", Em(8), "R2A_Task");
    b.bound_to("A8", "R2A_Synch_A8_Entry");
    b.replies_to("A8", "R2A_Synch_A8_Entry");
    b.activity("A11", Em(4), "R2A_Task");
    b.bound_to("A11", "R2A_Synch_A11_Entry");
    b.replies_to("A11", "R2A_Synch_A11_Entry");

    // T7
    b.activity("A9", Em(4), "R2B_Task");
    b.bound_to("A9", "R2B_Synch_A9_Entry");
    b.activity("A10", Em(6), "R2B_Task");
    b.bound_to("A10", "R2B_Synch_A10_Entry");
    b.replies_to("A10", "R2B_Synch_A10_Entry");
    b.activity("A12", Em(8), "R2B_Task");
    b.bound_to("A12", "R2B_Synch_A12_Entry");
    b.replies_to("A12", "R2B_Synch_A12_Entry");
    b.activity("A9_Res_Empty", Immediate(), "R2B_Task");
    b.sync_call("A9_Res_Empty", "R3_Synch_A9_Entry", 1.0);
    b.replies_to("A9_Res_Empty", "R2B_Synch_A9_Entry");

    b.serial("A1_Empty", "A2_Empty");
    b.serial("A5_Empty", "A6_Empty");
    b.serial("E4_Empty", "A7_Empty");
    b.serial("A9_Empty", "A10_Empty");
    b.serial("A2", "A2_Res_Empty");
    b.serial("A5", "A5_Res_Empty");
    b.serial("A9", "A9_Res_Empty");
    b.or_fork("A2_Empty", {"A3_Empty", "A4_Empty"}, {0.6, 0.4});
    b.and_fork("A7_Empty", {"A8_Empty", "A9_Empty"});
    b.or_join({"A3_Empty", "A4_Empty"}, "A5_Empty");
    b.and_join({"A8_Empty", "A10_Empty"}, "A11_Empty");
    return b.build();
}

/** `myLayeredModel` of lqn_twotasks / lqn_moment3: one caller, two callee entries. */
LqnS lqn_twotasks_model() {
    Lqn b;
    b.processor("P1", 1, SchedStrategy::PS);
    b.task("T1", 100, SchedStrategy::REF, "P1");
    b.entry("E1", "T1");
    b.processor("P2", 1, SchedStrategy::PS);
    b.task("T2", 1, SchedStrategy::INF, "P2");
    b.entry("E2", "T2");
    b.entry("E3", "T2");
    b.think_time("T1", lang::erlang_fit_mean_order<double>(10.0, 1));
    b.activity("A1", Er(1), "T1");
    b.bound_to("A1", "E1");
    b.sync_call("A1", "E2", 1.0);
    b.sync_call("A1", "E3", 1.0);
    b.activity("A20", Er(1), "T2");
    b.bound_to("A20", "E2");
    b.activity("A21", Er(1), "T2");
    b.activity("A22", Er(1), "T2");
    b.replies_to("A22", "E2");
    b.serial("A20", "A21");
    b.serial("A21", "A22");
    b.activity("A3", Er(1), "T2");
    b.bound_to("A3", "E3");
    b.replies_to("A3", "E3");
    return b.build();
}

/**
 * A processor whose servers are not interchangeable.
 *
 * P1 declares three servers, but they are not a homogeneous pool: S1 is
 * dedicated to task T2, S3 to task T3, and only S2 can take either. Neither
 * task can therefore reach more than two of the three servers, and the model is
 * a different system from a plain multiplicity-3 processor even though it holds
 * the same number of servers.
 *
 *   P1 (3 servers)     S1 --- T2
 *                      S2 --< T2, T3
 *                      S3 --- T3
 *
 * SolverLN lowers the declaration to the activated-server rate of
 * api::sn_compat_rate, carried onto the layer station as a joint dependence, so
 * a compatibility declaration is an APPROXIMATION inside a layer and is
 * admitted only under the class-switching layerings ("srvn.cs", "flat.cs").
 */
LqnS lqn_server_pools_model(bool compatibility) {
    Lqn b;
    b.processor("P0", 1, SchedStrategy::PS);
    b.processor("P1", 3, SchedStrategy::PS);
    b.task("T1", 3, SchedStrategy::REF, "P0");
    b.think_time("T1", Er(1.0));
    b.task("T2", 3, SchedStrategy::FCFS, "P1");
    b.task("T3", 3, SchedStrategy::FCFS, "P1");
    b.entry("E1", "T1");
    b.entry("E2", "T2");
    b.entry("E3", "T3");
    b.activity("A1", Er(2.0), "T1");
    b.bound_to("A1", "E1");
    b.sync_call("A1", "E2", 1.0);
    b.sync_call("A1", "E3", 1.0);
    b.activity("A2", Er(3.0), "T2");
    b.bound_to("A2", "E2");
    b.replies_to("A2", "E2");
    b.activity("A3", Er(2.0), "T3");
    b.bound_to("A3", "E3");
    b.replies_to("A3", "E3");
    if (compatibility) {
        b.add_server_type("P1", "S1", 1, {"T2"}, 1.0);        // dedicated to T2
        b.add_server_type("P1", "S2", 1, {"T2", "T3"}, 1.0);  // shared
        b.add_server_type("P1", "S3", 1, {"T3"}, 1.0);        // dedicated to T3
    } else {
        // one pool of three, every task eligible on every server: the neutral
        // declaration, which reproduces the plain multiserver
        b.add_server_type("P1", "All", 3, {"T2", "T3"}, 1.0);
    }
    return b.build();
}

}  // namespace

// ---------------------------------------------------------------------------
// Programmatic models
// ---------------------------------------------------------------------------

/** The two-processor three-task chain, solved by LN over MVA layers. */
void lqn_basic() { solve_ln(lqn_basic_model()); }

/** The compatibility-pool processor, neutral pool first, then the graph. */
void lqn_server_pools() {
    ln::LnOptions opt;
    // a compatibility declaration is a station rate law, which only the
    // class-switching layerings can carry
    opt.method = "srvn.cs";
    std::cout << "--- homogeneous pool on P1 ---" << std::endl;
    solve_ln(lqn_server_pools_model(false), opt);
    std::cout << "--- compatibility pool on P1 ---" << std::endl;
    solve_ln(lqn_server_pools_model(true), opt);
}

/**
 * The `lqn_basic` model with the LN state exposed.
 *
 * The reference replaces `solver.pre` and `solver.post` with tracing closures,
 * which prints the internal vectors BETWEEN two layer solves. `SolverLN` has no
 * such hook, so what is printed here is the converged state of the same vectors
 * plus the per-iteration layer results the solver does retain.
 */
void lqn_basic_debug() {
    const LqnS l = lqn_basic_model();
    ln::LnOptions opt;
    ln::SolverLN<double> s(l, opt);
    const ln::LnSolution<double> sol = s.get_ensemble_avg();

    // TODO(cpp): solver.pre = debug_pre; solver.post = debug_post; solver.get_avg_table()
    na("SolverLN.pre/post hooks",
       "the reference traces the update maps (servt_classes_updmap, thinkt_classes_updmap, "
       "call_classes_updmap, idxhash) from inside the iteration; SolverLN keeps those private and "
       "offers no per-iteration callback, so the mid-iteration snapshots have no counterpart");

    std::printf("\n=== LAYERS ===\n");
    std::printf("nlayers: %zu\n", s.nlayers());
    for (std::size_t k = 0; k < s.nlayers(); ++k)
        std::printf("  Layer %zu: %s\n", k, s.layers()[k].name.c_str());

    const char* tasks[] = {"T:T1", "T:T2", "T:T3"};
    const char* entries[] = {"E:E1", "E:E2", "E:E3"};
    const char* acts[] = {"A:AS1", "A:AS2", "A:AS3"};

    std::printf("\n=== CONVERGED STATE ===\n");
    for (int k = 0; k < 3; ++k) {
        const std::size_t t = idx_of(l, tasks[k]);
        std::printf("tput[%s]=%.6e util[%s]=%.6e thinkt[%s]=%.6e thinktproc mean=%.6e\n", tasks[k],
                    s.state_tput()[t], tasks[k], s.state_util()[t], tasks[k], s.state_thinkt()[t],
                    s.state_thinktproc()[t].mean);
    }
    for (int k = 0; k < 3; ++k) {
        const std::size_t a = idx_of(l, acts[k]);
        std::printf("servt[%s]=%.6e servtproc mean=%.6e\n", acts[k], s.state_servt()[a],
                    s.state_servtproc()[a].mean);
    }
    for (int k = 0; k < 3; ++k) {
        const std::size_t e = idx_of(l, entries[k]);
        std::printf("servt[%s]=%.6e\n", entries[k], s.state_servt()[e]);
    }
    for (std::size_t c = 1; c <= l.ncalls; ++c)
        std::printf("callservt[%s]=%.6e callservtproc mean=%.6e callproc_mean=%.6g\n",
                    l.callhashnames[c].c_str(), s.state_callservt()[c],
                    s.state_callservtproc()[c].mean, l.callproc_mean[c]);

    const std::vector<std::vector<ln::LayerResult<double> > >& hist = s.iteration_results();
    std::printf("\n=== PER-ITERATION LAYER RESULTS (%zu iterations) ===\n", hist.size());
    for (std::size_t it = 0; it < hist.size(); ++it)
        for (std::size_t e = 0; e < hist[it].size(); ++e)
            std::printf("  Iter %zu, Layer %zu: max(TN)=%.6f max(QN)=%.6f max(UN)=%.6f\n", it + 1, e,
                        mat_max(hist[it][e].TN), mat_max(hist[it][e].QN), mat_max(hist[it][e].UN));

    std::printf("\n=== FINAL RESULTS ===\n");
    print_ln(l, sol, opt, s.nlayers());
}

/**
 * The `lqn_basic` model with its layer structure dumped.
 *
 * The reference walks `solver.ensemble` printing nodes, classes, service times,
 * routing and visits, then solves one layer at a time and checks that no
 * station's total utilization exceeds its server count. The layer dump here is
 * the CLI's `-o layers` form, which carries the same information in the layout
 * the rest of this port already uses.
 */
void lqn_debug2() {
    const LqnS l = lqn_basic_model();
    ln::LnOptions opt;
    ln::SolverLN<double> s(l, opt);

    std::printf("\n=== LAYER STRUCTURE ===\n");
    dump_layers(s);

    const ln::LnSolution<double> sol = s.get_ensemble_avg();
    const std::vector<std::vector<ln::LayerResult<double> > >& hist = s.iteration_results();
    if (!hist.empty()) {
        std::printf("\n=== FIRST ITERATION ===\n");
        for (std::size_t e = 0; e < hist[0].size(); ++e) {
            const qn::Layer<double>& L = s.layers()[e];
            const Matrix<double>& UN = hist[0][e].UN;
            std::printf("Layer %zu (%s):\n", e, L.name.c_str());
            for (std::size_t i = 0; i < L.stations.size(); ++i) {
                double total = 0.0;
                for (std::size_t r = 0; r < UN.cols(); ++r) total += UN(i, r);
                std::printf("  Total utilization at %s: %.4f (servers: %g)\n",
                            L.stations[i].name.c_str(), total, L.stations[i].nservers);
                if (total > L.stations[i].nservers)
                    std::printf("    WARNING: Unstable! Utilization exceeds capacity!\n");
            }
        }
    }

    std::printf("\n=== FINAL RESULTS ===\n");
    print_ln(l, sol, opt, s.nlayers());
}

/** Two processors, serial activities inside each task, one synchronous call. */
void lqn_serial() {
    Lqn b;
    b.processor("P1", 1, SchedStrategy::PS);
    b.processor("P2", 1, SchedStrategy::PS);
    b.task("T1", 10, SchedStrategy::REF, "P1");
    b.think_time("T1", Em(100));
    b.task("T2", 1, SchedStrategy::FCFS, "P2");
    b.think_time("T2", Immediate());
    b.entry("E1", "T1");
    b.entry("E2", "T2");
    b.activity("AS1", Em(1.6), "T1");
    b.bound_to("AS1", "E1");
    b.activity("AS2", Immediate(), "T1");
    b.sync_call("AS2", "E2", 1.0);
    b.activity("AS3", Em(5), "T2");
    b.bound_to("AS3", "E2");
    b.activity("AS4", Em(1), "T2");
    b.replies_to("AS4", "E2");
    b.serial("AS1", "AS2");
    b.serial("AS3", "AS4");
    const LqnS l = b.build();

    // TODO(cpp): solver = LQNS(model, keep=True); print(solver.get_avg_table())
    na("LQNS", "the reference solves this model with the external lqns binary, which this port "
               "does not carry; SolverLN below is a different algorithm, not a stand-in");
    solve_ln(l);
}

/** One reference task calling two entries of one infinite-server task. */
void lqn_twotasks() {
    const LqnS l = lqn_twotasks_model();
    // TODO(cpp): solver_lqns = LQNS(model, keep=True, verbose=False)
    // TODO(cpp): print(solver_lqns.get_avg_table())
    na("LQNS", "the reference solves this model with the external lqns binary, which this port "
               "does not carry");
    solve_ln(l, nc_layers());
}

/** The same model under the default update and under the APH moment update. */
void lqn_moment3() {
    const LqnS l = lqn_twotasks_model();
    solve_ln(l, nc_layers());

    ln::LnOptions opt = nc_layers();
    opt.method = "moment3";
    try {
        solve_ln(l, opt);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "\nLN(moment3) failed: %s\n", e.what());
        solve_ln(l, nc_layers());
    }
}

/** One call of multiplicity three into an APH-served entry, over MVA and NC layers. */
void lqn_multi_solvers() {
    Lqn b;
    b.processor("P1", std::numeric_limits<double>::infinity(), SchedStrategy::INF);
    b.task("T1", 1, SchedStrategy::REF, "P1");
    b.entry("E1", "T1");
    b.processor("P2", std::numeric_limits<double>::infinity(), SchedStrategy::INF);
    b.task("T2", std::numeric_limits<double>::infinity(), SchedStrategy::INF, "P2");
    b.entry("E2", "T2");
    b.think_time("T1", lang::erlang_fit_mean_order<double>(0.0001, 2));
    b.activity("A1", Er(1.0), "T1");
    b.bound_to("A1", "E1");
    b.sync_call("A1", "E2", 3.0);
    b.activity("A2", lang::aph_fit_mean_scv<double>(1.0, 10.0), "T2");
    b.bound_to("A2", "E2");
    b.replies_to("A2", "E2");
    const LqnS l = b.build();

    // TODO(cpp): options = LQNS.default_options(); options.keep = True; options.verbose = 1
    // TODO(cpp): print(LQNS(model, options).avg_table())
    na("LQNS", "the reference's first solver is the external lqns binary, which this port does not "
               "carry");
    // The reference also sets lnoptions.seed = 2300; it selects nothing here,
    // since neither MVA nor NC layers draw random numbers.
    solve_ln(l);
    solve_ln(l, nc_layers());
}

/** A loop, an AND fork/join and an OR fork/join, one per task. */
void lqn_workflows() {
    std::printf("This example illustrates a layered network with a loop.\n");
    Lqn b;
    b.processor("P1", lang::GlobalConstants::MaxInt, SchedStrategy::INF);
    b.task("T1", 1, SchedStrategy::REF, "P1");
    b.think_time("T1", Immediate());
    b.entry("Entry", "T1");
    b.processor("P2", lang::GlobalConstants::MaxInt, SchedStrategy::INF);
    b.task("T2", 1, SchedStrategy::INF, "P2");
    b.think_time("T2", Immediate());
    b.entry("E2", "T2");
    b.processor("P3", 5, SchedStrategy::PS);
    b.task("T3", 20, SchedStrategy::INF, "P3");
    b.think_time("T3", Em(10));
    b.entry("E1", "T3");

    b.activity("A1", Em(1), "T1");
    b.bound_to("A1", "Entry");
    b.activity("A2", Em(2), "T1");
    b.activity("A3", Em(3), "T1");
    b.sync_call("A3", "E2", 1.0);

    b.activity("B1", Em(0.1), "T2");
    b.bound_to("B1", "E2");
    b.activity("B2", Em(0.2), "T2");
    b.activity("B3", Em(0.3), "T2");
    b.activity("B4", Em(0.4), "T2");
    b.activity("B5", Em(0.5), "T2");
    b.activity("B6", Em(0.6), "T2");
    b.sync_call("B6", "E1", 1.0);
    b.replies_to("B6", "E2");

    b.activity("C1", Em(0.1), "T3");
    b.bound_to("C1", "E1");
    b.activity("C2", Em(0.2), "T3");
    b.activity("C3", Em(0.3), "T3");
    b.activity("C4", Em(0.4), "T3");
    b.activity("C5", Em(0.5), "T3");
    b.replies_to("C5", "E1");

    // Loop(A1, {A2, A3}, 3): the last element of the MATLAB post list is the
    // loop EXIT, the ones before it the body.
    b.loop("A1", {"A2"}, "A3", 3.0);
    b.serial("B4", "B5");
    b.and_fork("B1", {"B2", "B3", "B4"});
    b.and_join({"B2", "B3", "B5"}, "B6");
    b.or_fork("C1", {"C2", "C3", "C4"}, {0.3, 0.3, 0.4});
    b.or_join({"C2", "C3", "C4"}, "C5");
    const LqnS l = b.build();

    // TODO(cpp): if LQNS.isAvailable(): print(LQNS(model).avg_table())
    na("LQNS", "the reference solves this model with the external lqns binary when it is "
               "installed; this port does not carry it");
    solve_ln(l);
}

/** The BPMN-derived network, over MVA layers. */
void lqn_bpmn() {
    const LqnS l = lqn_bpmn_model();
    std::printf("This example illustrates the solution of a complex layered queueing network "
                "extracted from a BPMN model.\n");
    // TODO(cpp): options = LQNS.default_options(); options.keep = True
    // TODO(cpp): if LQNS.isAvailable(): print(LQNS(model).avg_table())
    na("LQNS", "the reference solves this model with the external lqns binary when it is "
               "installed; this port does not carry it");
    solve_ln(l);
}

/** A task that pays a setup time after an idle period and a delay-off timer. */
void lqn_setup() {
    std::printf("Example of layered model with a task paying setup and delay-off times\n");
    Lqn b;
    b.processor("P1", lang::GlobalConstants::MaxInt, SchedStrategy::INF);
    b.task("T1", 1, SchedStrategy::REF, "P1");
    b.entry("E1", "T1");
    b.processor("P2", 4, SchedStrategy::FCFS);
    b.task("F2", 6, SchedStrategy::FCFS, "P2");
    b.think_time("F2", Em(8.0));
    b.setup_time("F2", Er(1.0), Er(2.0));
    b.entry("E2", "F2");
    b.activity("A1", Er(1.0), "T1");
    b.bound_to("A1", "E1");
    b.sync_call("A1", "E2", 1.0);
    b.activity("A2", Er(3.0), "F2");
    b.bound_to("A2", "E2");
    b.replies_to("A2", "E2");
    // lnoptions.seed = 23000 in the reference selects nothing over MVA layers.
    solve_ln(b.build());
}

/** An entry fed by an exogenous Poisson stream instead of a rendezvous. */
void lqn_open_arrival() {
    Lqn b;
    b.processor("P1", 1, SchedStrategy::PS);
    b.task("T1", 1, SchedStrategy::FCFS, "P1");
    b.think_time("T1", Immediate());
    b.entry("E1", "T1");
    b.open_arrival("E1", Er(0.2));
    b.activity("A1", Em(1.6), "T1");
    b.bound_to("A1", "E1");
    b.replies_to("A1", "E1");
    solve_ln(b.build());
}

/**
 * An AND fork/join on a task that ALSO takes an entry-level open arrival.
 *
 * SolverLN is expected to refuse the combination: the fork-join transform mints
 * its own Source and collides with the Source the open stream is routed
 * through. The refusal is the reference's own output, so it is caught and
 * printed rather than allowed to fail the example.
 */
void lqn_fork_open_arrival() {
    Lqn b;
    b.processor("P1", std::numeric_limits<double>::infinity(), SchedStrategy::INF);
    b.task("Client", 1, SchedStrategy::REF, "P1");
    b.think_time("Client", Em(1.0));
    b.entry("CE", "Client");
    b.processor("P2", std::numeric_limits<double>::infinity(), SchedStrategy::INF);
    b.task("Server", 1, SchedStrategy::FCFS, "P2");
    b.think_time("Server", Immediate());
    b.entry("SE", "Server");
    b.entry("OE", "Server");
    b.open_arrival("OE", Er(0.1));

    b.activity("CA", Em(0.5), "Client");
    b.bound_to("CA", "CE");
    b.sync_call("CA", "SE", 1.0);

    b.activity("RA1", Em(0.2), "Server");
    b.bound_to("RA1", "SE");
    b.activity("RA2", Em(0.3), "Server");
    b.activity("RA3", Em(0.4), "Server");
    b.activity("RA4", Em(0.1), "Server");
    b.replies_to("RA4", "SE");
    b.and_fork("RA1", {"RA2", "RA3"});
    b.and_join({"RA2", "RA3"}, "RA4");

    b.activity("OA1", Em(0.2), "Server");
    b.bound_to("OA1", "OE");
    b.activity("OA2", Em(0.3), "Server");
    b.activity("OA3", Em(0.4), "Server");
    b.activity("OA4", Em(0.1), "Server");
    b.and_fork("OA1", {"OA2", "OA3"});
    b.and_join({"OA2", "OA3"}, "OA4");
    const LqnS l = b.build();

    try {
        solve_ln(l);
    } catch (const std::exception& e) {
        std::printf("LN refuses this model: %s\n", e.what());
    }
}

/**
 * The Sock Shop microservice benchmark: 7 hosts, 7 tasks, 12 entries.
 *
 * `LqnBuilder` has no fan-out/fan-in setter, so the declarations are written
 * onto the raw model before `lqn_finalize` -- dropping them silently would be a
 * different model, and P2_1 carries a replication of 2 that fan-out interacts
 * with (SolverLN's layer-pooling test reads `fanout_at`). The processor quantum
 * the reference sets is an lqns scheduling attribute that no analyzer in any
 * codebase reads, so it has no counterpart here.
 */
void lqn_sockshop() {
    Lqn b;
    b.processor("P1", 1, SchedStrategy::INF);
    b.processor("P2_1", 1, SchedStrategy::PS, 2.0);
    b.processor("P2_2", 1, SchedStrategy::PS);
    b.processor("P2_3", 1, SchedStrategy::PS);
    b.processor("P3_1", 1, SchedStrategy::PS);
    b.processor("P3_2", 1, SchedStrategy::PS);
    b.processor("P3_3", 1, SchedStrategy::PS);

    b.task("T0", 1000, SchedStrategy::REF, "P1");
    b.think_time("T0", Em(7.0));
    b.task("T1", 24, SchedStrategy::FCFS, "P2_1");
    b.think_time("T1", Immediate());
    b.task("T2", 21, SchedStrategy::FCFS, "P2_2");
    b.think_time("T2", Immediate());
    b.task("T6", 100, SchedStrategy::FCFS, "P2_3");
    b.think_time("T6", Immediate());
    b.task("T3", 139, SchedStrategy::FCFS, "P3_1");
    b.think_time("T3", Immediate());
    b.task("T4", 16, SchedStrategy::FCFS, "P3_2");
    b.think_time("T4", Immediate());
    b.task("T5", 151, SchedStrategy::FCFS, "P3_3");
    b.think_time("T5", Immediate());

    b.entry("E0", "T0");
    b.entry("E1", "T1");
    b.entry("E2", "T2");
    b.entry("E3", "T2");
    b.entry("E4", "T2");
    b.entry("E11", "T6");
    b.entry("E5", "T3");
    b.entry("E6", "T3");
    b.entry("E7", "T3");
    b.entry("E8", "T4");
    b.entry("E9", "T4");
    b.entry("E10", "T5");

    b.activity("AS", Em(0.00000005), "T0");
    b.bound_to("AS", "E0");
    b.activity("AS0", Em(0.00000005), "T0");
    b.sync_call("AS0", "E1", 1.0);

    b.activity("AS1", Em(0.00000005), "T1");
    b.bound_to("AS1", "E1");
    b.activity("AS2", Em(0.0022216), "T1");
    b.sync_call("AS2", "E2", 0.33);
    b.sync_call("AS2", "E3", 0.17);
    b.sync_call("AS2", "E4", 0.50);
    b.replies_to("AS2", "E1");

    b.activity("AH1", Em(0.00000005), "T2");
    b.bound_to("AH1", "E2");
    b.activity("AH2", Em(0.0021319), "T2");
    b.replies_to("AH2", "E2");
    b.activity("AH3", Em(0.00000005), "T2");
    b.bound_to("AH3", "E3");
    b.activity("AH4", Em(0.0037561), "T2");
    b.sync_call("AH4", "E8", 0.5);
    b.sync_call("AH4", "E9", 0.5);
    b.replies_to("AH4", "E3");
    b.activity("AH5", Em(0.00000005), "T2");
    b.bound_to("AH5", "E4");
    b.activity("AH6", Em(0.0051774), "T2");
    b.sync_call("AH6", "E5", 0.33);
    b.sync_call("AH6", "E6", 0.33);
    b.sync_call("AH6", "E7", 0.33);
    b.replies_to("AH6", "E4");

    b.activity("AH15", Em(0.0000000005), "T6");
    b.bound_to("AH15", "E11");
    b.activity("AH16", Em(0.0040355), "T6");
    b.replies_to("AH16", "E11");

    b.activity("AH7", Em(0.0000000005), "T3");
    b.bound_to("AH7", "E5");
    b.activity("AH8", Em(0.0029469), "T3");
    b.sync_call("AH8", "E11", 1.0);
    b.replies_to("AH8", "E5");
    b.activity("AH9", Em(0.0000000005), "T3");
    b.bound_to("AH9", "E6");
    b.activity("AH10", Em(0.012323), "T3");
    b.sync_call("AH10", "E11", 1.0);
    b.replies_to("AH10", "E6");
    b.activity("AH11", Em(0.0000000005), "T3");
    b.bound_to("AH11", "E7");
    b.activity("AH12", Em(0.0033488), "T3");
    b.sync_call("AH12", "E11", 1.0);
    b.replies_to("AH12", "E7");

    b.activity("AS3", Em(0.0000000005), "T4");
    b.bound_to("AS3", "E8");
    b.activity("AS4", Em(0.0034925), "T4");
    b.sync_call("AS4", "E10", 1.0);
    b.replies_to("AS4", "E8");
    b.activity("AS5", Em(0.0000000005), "T4");
    b.bound_to("AS5", "E9");
    b.activity("AS6", Em(0.0030162), "T4");
    b.sync_call("AS6", "E10", 1.0);
    b.replies_to("AS6", "E9");

    b.activity("AH13", Em(0.0000000005), "T5");
    b.bound_to("AH13", "E10");
    b.activity("AH14", Em(0.0032434), "T5");
    b.replies_to("AH14", "E10");

    b.serial("AS", "AS0");
    b.serial("AS1", "AS2");
    b.serial("AH1", "AH2");
    b.serial("AH3", "AH4");
    b.serial("AH5", "AH6");
    b.serial("AH15", "AH16");
    b.serial("AH7", "AH8");
    b.serial("AH9", "AH10");
    b.serial("AH11", "AH12");
    b.serial("AS3", "AS4");
    b.serial("AS5", "AS6");
    b.serial("AH13", "AH14");

    lqn::LqnModel<double> m = b.model();
    auto task_slot = [&](const std::string& name) {
        for (std::size_t i = 0; i < m.tasks.size(); ++i)
            if (m.tasks[i].name == name) return i;
        throw InputError("lqn_sockshop: unknown task '" + name + "'");
    };
    m.tasks[task_slot("T1")].fanout.push_back(std::make_pair(std::string("T2"), 1.0));
    m.tasks[task_slot("T2")].fanout.push_back(std::make_pair(std::string("T3"), 1.0));
    m.tasks[task_slot("T2")].fanout.push_back(std::make_pair(std::string("T4"), 1.0));
    m.tasks[task_slot("T2")].fanin.push_back(std::make_pair(std::string("T1"), 1.0));
    m.tasks[task_slot("T6")].fanin.push_back(std::make_pair(std::string("T3"), 1.0));
    m.tasks[task_slot("T3")].fanout.push_back(std::make_pair(std::string("T6"), 1.0));
    m.tasks[task_slot("T3")].fanin.push_back(std::make_pair(std::string("T2"), 1.0));
    m.tasks[task_slot("T4")].fanout.push_back(std::make_pair(std::string("T5"), 1.0));
    m.tasks[task_slot("T4")].fanin.push_back(std::make_pair(std::string("T2"), 1.0));
    m.tasks[task_slot("T5")].fanin.push_back(std::make_pair(std::string("T4"), 1.0));

    solve_ln(lqn::lqn_finalize(m));
}

// ---------------------------------------------------------------------------
// Models read from a file
// ---------------------------------------------------------------------------

/**
 * `lqn_serial.xml`, solved without a warm start, then the response-time CDF.
 *
 * The reference reads the CDF off `FLD(ensemble[2])`, the third layer of the
 * converged ensemble; `qn::Layer` derives from `NetworkStruct`, so that layer
 * goes straight into `solver_fluid_cdf_respt`.
 */
void lqn_init() {
    std::printf("This example illustrates the initialization of LN using the output of LQNS.\n");
    const LqnS l = lqn::read_lqnx<double>(layered_file("lqn_serial.xml"));

    // TODO(cpp): options = LQNS.default_options(); options.keep = True
    // TODO(cpp): if LQNS.isAvailable(): print(LQNS(model, options).avg_table())
    na("LQNS", "the reference initializes LN from the external lqns binary's answer; this port "
               "carries neither lqns nor its output, so there is nothing to initialize from");

    std::printf("\nSolve with LN without initialization:\n");
    ln::LnOptions opt;
    const auto t0 = std::chrono::steady_clock::now();
    ln::SolverLN<double> s(l, opt);
    const ln::LnSolution<double> sol = s.get_ensemble_avg();
    const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0)
                               .count();
    print_ln(l, sol, opt, s.nlayers());
    std::printf("Time elapsed: %.3fs\n", elapsed);

    std::printf("\nWe now obtain the CDF of response times:\n");
    if (s.nlayers() < 3) {
        std::printf("Model ensemble not available or insufficient layers\n");
        return;
    }
    const qn::Layer<double>& L = s.layers()[2];
    fluid::FluidOptions fopt;
    const std::vector<std::vector<fluid::FluidPassage> > RD =
        fluid::solver_fluid_cdf_respt<double>(L, fopt);
    for (std::size_t i = 0; i < RD.size(); ++i)
        for (std::size_t c = 0; c < RD[i].size(); ++c) {
            if (RD[i][c].t.empty()) continue;
            const fluid::FluidPassage& p = RD[i][c];
            std::printf("RD[%s][%s]: %zu points, t in [%.6g, %.6g], F in [%.6g, %.6g]\n",
                        L.stations[i].name.c_str(), L.classes[c].name.c_str(), p.t.size(),
                        p.t.front(), p.t.back(), p.cdf.front(), p.cdf.back());
        }
}

/** A moderately large LQN read from `lqn_ofbiz.xml`, over NC layers. */
void lqn_ofbiz() {
    std::printf("This example illustrates the solution of a moderately large LQN.\n");
    const LqnS l = lqn::read_lqnx<double>(layered_file("lqn_ofbiz.xml"));
    const ln::LnOptions opt = nc_layers();
    const auto t0 = std::chrono::steady_clock::now();
    ln::SolverLN<double> s(l, opt);
    const ln::LnSolution<double> sol = s.get_ensemble_avg();
    const double tnoinit = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0)
                               .count();
    print_ln(l, sol, opt, s.nlayers());
    std::printf("Tnoinit = %.6f\n", tnoinit);
}

// ---------------------------------------------------------------------------
// Routed call groups
// ---------------------------------------------------------------------------

namespace {

/**
 * The client-and-three-servers model both dispatch examples share.
 *
 * ONE call per invocation, its destination chosen by `strategy` over the three
 * interchangeable servers. A mean of 1 over 3 targets is where the policy
 * actually bites: the probabilistic twin makes 0..3 calls per invocation with
 * the same mean, a routed group makes exactly one.
 */
LqnS routed_group_model(bool jsq) {
    Lqn b;
    b.processor("PC", 1, SchedStrategy::INF);
    b.processor("PS", 1, SchedStrategy::PS);
    b.task("TC", 10, SchedStrategy::REF, "PC");
    b.think_time("TC", Exp(1.0 / 5.0));
    b.task("TS1", 5, SchedStrategy::FCFS, "PS");
    b.task("TS2", 5, SchedStrategy::FCFS, "PS");
    b.task("TS3", 5, SchedStrategy::FCFS, "PS");
    b.entry("EC", "TC");
    b.entry("ES1", "TS1");
    b.entry("ES2", "TS2");
    b.entry("ES3", "TS3");

    b.activity("AC", Exp(2.0), "TC");
    b.bound_to("AC", "EC");
    std::vector<std::string> targets;
    targets.push_back("ES1");
    targets.push_back("ES2");
    targets.push_back("ES3");
    if (jsq) {
        b.sync_call_jsq("AC", targets, 1.0);
    } else {
        b.sync_call_round_robin("AC", targets, 1.0);
    }

    b.activity("AS1", Exp(1.0), "TS1");
    b.bound_to("AS1", "ES1");
    b.replies_to("AS1", "ES1");
    b.activity("AS2", Exp(1.0), "TS2");
    b.bound_to("AS2", "ES2");
    b.replies_to("AS2", "ES2");
    b.activity("AS3", Exp(1.0), "TS3");
    b.bound_to("AS3", "ES3");
    b.replies_to("AS3", "ES3");
    return b.build();
}

/**
 * The only configuration a routed call group admits.
 *
 * `SolverLN::assert_call_groups` refuses anything else, and is right to: under
 * `srvn` each server task lives in its own submodel and is replaced, in the
 * client's submodel, by a surrogate delay, so no node ever has arcs to more than
 * one of them; and MVA, NC and FLD read the routing matrix, into which
 * `refresh_routing` has already expanded the strategy as a uniform split, so
 * they would report a deterministic policy as a coin.
 */
ln::LnOptions routed_group_options() {
    ln::LnOptions opt;
    opt.method = "flat";
    opt.layer_solver = "ssa";
    return opt;
}

}  // namespace

/** Round-robin call dispatch over three interchangeable server tasks. */
void lqn_rrobin() {
    solve_ln(routed_group_model(false), routed_group_options());
}

/** Join-the-shortest-queue call dispatch over three interchangeable server tasks. */
void lqn_jsq() {
    solve_ln(routed_group_model(true), routed_group_options());
}

// ---------------------------------------------------------------------------
// Parameter identification
// ---------------------------------------------------------------------------

namespace {

/** The 16 x 3 standard normal draws of `np.random.seed(12345)` in the reference. */
const double kNoise[16][3] = {
    {-0.20470765948471295, 0.47894333805754824, -0.51943871505673811},
    {-0.55573030434749005, 1.9657805725027142, 1.3934058329729904},
    {0.092907876743717671, 0.28174615283020249, 0.76902256761183874},
    {1.2464347363862822, 1.0071893575830049, -1.2962211091122635},
    {0.27499163343212402, 0.22891287893531592, 1.3529168351654497},
    {0.88642934059158884, -2.0016373096603974, -0.37184253714025439},
    {1.6690253095248706, -0.43856973583557191, -0.53974144552166281},
    {0.47698501041229951, 3.2489439194307548, -1.0212275243555968},
    {-0.57708730304076716, 0.12412127567340774, 0.30261356191251138},
    {0.52377206815041655, 0.00094027777533288513, 1.3438097936141322},
    {-0.71354398509638317, -0.83115353885391396, -2.3702316539567447},
    {-1.8607607885507347, -0.86075739843174859, 0.56014529302803417},
    {-1.2659344916936925, 0.11982712466055132, -1.0635124480535407},
    {0.3328827156076713, -2.3594188073836815, -0.19954295533566749},
    {-1.5419955278741118, -0.97073591225918254, -1.3070302509708969},
    {0.28634974701415511, 0.37798411093737727, -0.75388653478989709}};

/** The light-population LQN whose think time and AS3 host demand are hidden. */
LqnS paramident_model() {
    Lqn b;
    b.processor("P1", 1, SchedStrategy::PS);
    b.processor("P2", 1, SchedStrategy::PS);
    b.task("T1", 5, SchedStrategy::REF, "P1");
    b.think_time("T1", Em(0.5));
    b.task("T2", 5, SchedStrategy::FCFS, "P1");
    b.think_time("T2", Em(1.0 / 3.0));
    b.task("T3", 3, SchedStrategy::FCFS, "P2");
    b.think_time("T3", Em(0.25));
    b.entry("E1", "T1");
    b.entry("E2", "T2");
    b.entry("E3", "T3");
    b.activity("AS1", Em(0.1), "T1");
    b.bound_to("AS1", "E1");
    b.sync_call("AS1", "E2", 1.0);
    b.activity("AS2", Em(0.05), "T2");
    b.bound_to("AS2", "E2");
    b.sync_call("AS2", "E3", 5.0);
    b.replies_to("AS2", "E2");
    b.activity("AS3", Em(0.02), "T3");
    b.bound_to("AS3", "E3");
    b.replies_to("AS3", "E3");
    return b.build();
}

}  // namespace

/**
 * Track two hidden LQN parameters with an Extended Kalman Filter.
 *
 * Zheng, Yang, Woodside, Litoiu, Iszlai, "Tracking Time-Varying Parameters in
 * Software Systems with Extended Kalman Filters", CASCON 2005. The C++ port
 * carries the filter (`infer_lqn_ekf`) and the observation selector
 * (`infer_lqn_getobs`) but not the `infer_lqn` wrapper, whose job is to build
 * Q, R and P0 from the scale factors and to close over the model; that wrapper
 * is written out here, in the reference's own equations (9a) and (9b).
 *
 * The measurement noise is Python's `np.random.randn` under seed 12345, which
 * no C++ generator reproduces, so the realized draws are the table above.
 */
void lqn_paramident() {
    const std::size_t nsteps = 16, np = 2, no = 3;
    LqnS l = paramident_model();

    const std::size_t t1 = infer::infer_lqn_findbyname(l.names, "T1");
    const std::size_t as3 = infer::infer_lqn_findbyname(l.names, "AS3");
    if (t1 == infer::INFER_LQN_NOT_FOUND || as3 == infer::INFER_LQN_NOT_FOUND)
        throw InputError("lqn_paramident: the hidden parameters name no element of the LQN");

    std::vector<infer::LqnObsSpec> obs(3);
    obs[0].metric = infer::LqnMetric::RespT;
    obs[0].name = "E1";
    obs[1].metric = infer::LqnMetric::Util;
    obs[1].name = "P1";
    obs[2].metric = infer::LqnMetric::Util;
    obs[2].name = "P2";

    // infer_lqn_setparams: a mean is injected as an exponential, i.e. SCV = 1
    auto hfun = [&](const std::vector<double>& a) {
        l.think[t1] = Em(a[0]);
        l.hostdem[as3] = Em(a[1]);
        ln::LnOptions opt;
        ln::SolverLN<double> s(l, opt);
        const ln::LnSolution<double> sol = s.get_ensemble_avg();
        infer::LqnMetrics<double> mt;
        mt.QLen = sol.QN;
        mt.Util = sol.UN;
        mt.RespT = sol.RN;
        mt.Tput = sol.TN;
        return infer::infer_lqn_getobs(l.names, mt, obs);
    };

    Matrix<double> a_true(np, nsteps, 0.0);
    for (std::size_t k = 0; k < nsteps; ++k) {
        a_true(0, k) = k < nsteps / 2 ? 0.5 : 1.0;
        a_true(1, k) = (k >= 4 && k < 12) ? 1.0 / 25.0 : 1.0 / 50.0;
    }

    Matrix<double> Z(no, nsteps, 0.0);
    for (std::size_t k = 0; k < nsteps; ++k) {
        std::vector<double> a(np);
        for (std::size_t i = 0; i < np; ++i) a[i] = a_true(i, k);
        const std::vector<double> z = hfun(a);
        for (std::size_t i = 0; i < no; ++i) Z(i, k) = z[i] * (1.0 + 0.02 * kNoise[k][i]);
    }

    // eq (9a) and (9b) of infer_lqn: the covariances are scale factors on |a0|
    // and on the mean measurement, not free parameters
    const double eps = std::numeric_limits<double>::epsilon();
    const std::vector<double> a0 = {0.7, 1.0 / 40.0};
    Matrix<double> Q(np, np, 0.0), P0(np, np, 0.0), R(no, no, 0.0);
    for (std::size_t i = 0; i < np; ++i) {
        Q(i, i) = std::max(std::pow(0.1 * std::fabs(a0[i]), 2.0), eps);
        P0(i, i) = std::max(std::pow(0.5 * std::fabs(a0[i]), 2.0), eps);
    }
    for (std::size_t i = 0; i < no; ++i) {
        double zbar = 0.0;
        for (std::size_t k = 0; k < nsteps; ++k) zbar += Z(i, k);
        zbar /= double(nsteps);
        R(i, i) = std::max(std::pow(0.2 * std::fabs(zbar) / 1.96, 2.0), eps);
    }

    infer::EkfOptions<double> eopt;
    eopt.a_true.push_back(a_true(0, nsteps - 1));
    eopt.a_true.push_back(a_true(1, nsteps - 1));
    const infer::EkfResult<double> info = infer::infer_lqn_ekf<double>(hfun, a0, P0, Z, Q, R, eopt);

    std::printf("\nStep |  Z_true  Z_hat |  Sd_true  Sd_hat | ||e||\n");
    for (std::size_t k = 0; k < nsteps; ++k) {
        double en = 0.0;
        for (std::size_t i = 0; i < no; ++i) en += info.e(i, k) * info.e(i, k);
        std::printf("%4zu | %6.3f  %6.3f | %7.4f  %7.4f | %.3g\n", k + 1, a_true(0, k),
                    info.ahat(0, k), a_true(1, k), info.ahat(1, k), std::sqrt(en));
    }
    std::printf("\nFinal: Z = %.4f (true %.4f), Sd = %.5f (true %.5f)\n",
                info.ahat(0, nsteps - 1), a_true(0, nsteps - 1), info.ahat(1, nsteps - 1),
                a_true(1, nsteps - 1));
    std::printf("RMS tracking Ea = %.4g, prediction Er = %.4g\n", info.has_Ea ? info.Ea : 0.0,
                info.Er);
}

// ---------------------------------------------------------------------------
// lqn_flatph / lqn_srvnph: the layering encodings compared
// ---------------------------------------------------------------------------

namespace {

/** A three-tier serial call chain, one entry per tier, all demands Exp(5). */
LqnS lqn_flatph_model() {
    Lqn b;
    b.processor("P1", 1, SchedStrategy::PS);
    b.processor("P2", 1, SchedStrategy::PS);
    b.processor("P3", 1, SchedStrategy::PS);
    b.task("T1", 4, SchedStrategy::REF, "P1");
    b.think_time("T1", Er(1.0));
    b.task("T2", 2, SchedStrategy::FCFS, "P2");
    b.task("T3", 1, SchedStrategy::FCFS, "P3");
    b.entry("E1", "T1");
    b.entry("E2", "T2");
    b.entry("E3", "T3");
    b.activity("A1", Er(5.0), "T1");
    b.bound_to("A1", "E1");
    b.sync_call("A1", "E2", 1.0);
    b.activity("A2", Er(5.0), "T2");
    b.bound_to("A2", "E2");
    b.sync_call("A2", "E3", 1.0);
    b.replies_to("A2", "E2");
    b.activity("A3", Er(5.0), "T3");
    b.bound_to("A3", "E3");
    b.replies_to("A3", "E3");
    return b.build();
}

/** The nested fork/join/loop model whose shape is what the encodings differ on. */
LqnS lqn_srvnph_model() {
    Lqn b;
    b.processor("P1", 1, SchedStrategy::INF);
    b.task("T1", 20, SchedStrategy::REF, "P1");
    b.think_time("T1", Em(1.0));
    b.entry("E1", "T1");
    b.processor("P2", 1, SchedStrategy::PS);
    b.task("T2", 5, SchedStrategy::FCFS, "P2");
    b.entry("E2", "T2");
    b.entry("E3", "T2");
    b.processor("P3", 1, SchedStrategy::PS);
    b.task("T3", 3, SchedStrategy::FCFS, "P3");
    b.entry("E4", "T3");

    b.activity("A1", Em(0.1), "T1");
    b.bound_to("A1", "E1");
    b.sync_call("A1", "E2", 1.0);
    b.sync_call("A1", "E3", 1.0);

    b.activity("A20", Em(0.2), "T2");
    b.bound_to("A20", "E2");
    b.activity("A21", Em(0.3), "T2");
    b.activity("A22", Em(0.2), "T2");
    b.sync_call("A22", "E4", 1.0);
    b.activity("A23", Em(0.1), "T2");
    b.replies_to("A23", "E2");

    b.activity("A30", Em(0.1), "T2");
    b.bound_to("A30", "E3");
    b.activity("A31", Em(0.2), "T2");
    b.activity("A32", Em(0.1), "T2");
    b.activity("A33", Em(0.1), "T2");
    b.replies_to("A33", "E3");

    b.activity("A4", Em(0.15), "T3");
    b.bound_to("A4", "E4");
    b.replies_to("A4", "E4");

    b.and_fork("A20", {"A21", "A22"});
    b.and_join({"A21", "A22"}, "A23");
    // The reference's Loop(pre, {A31, A32, A33}, 3) packs the body and the loop
    // EXIT into one list; this builder takes them apart, and the LAST element of
    // the reference's post list is the exit.
    b.loop("A30", {"A31", "A32"}, "A33", 3.0);
    return b.build();
}

/** The four encodings, solved and reported the way the reference reports them. */
void solve_ln_labelled(const LqnS& l, const std::string& method) {
    ln::LnOptions opt;
    opt.method = method;
    ln::SolverLN<double> s(l, opt);
    const ln::LnSolution<double> sol = s.get_ensemble_avg();
    std::printf("\n--- method=%s (built as %s, %zu submodel(s))\n", method.c_str(),
                s.state_lnmethod().c_str(), s.nlayers());
    print_ln(l, sol, opt, s.nlayers());
}

}  // namespace

/**
 * The four layering encodings on one serial call chain.
 *
 * `srvn.*` gives every task its own submodel and replaces a callee by a surrogate
 * delay; `flat.*` squashes the whole graph into ONE submodel. The `.cs` and `.ph`
 * suffixes then say how a task's internal sequence is carried: class switching,
 * or a phase-type service. The four answers differ because the approximations
 * differ, not because the model does.
 */
void lqn_flatph() {
    const LqnS l = lqn_flatph_model();
    for (const std::string& m : {std::string("srvn.cs"), std::string("srvn.ph"),
                                 std::string("flat.cs"), std::string("flat.ph")})
        solve_ln_labelled(l, m);
}

/**
 * The routing encoding against the phase-type one, timed.
 *
 * The reference used to ask for this with the bare token `srvnph`, which
 * `ln_requested_method` does not recognise -- it accepts `srvn.ph` and `ph` --
 * and every unrecognised token falls through to `srvn.cs`. So the comparison was
 * silently routing-against-routing, while the header claimed it was
 * routing-against-phase-type and the class counts below were read as evidence
 * for it. The token is now spelled the way it is resolved, in all four
 * codebases; the fallback itself is deliberate and stays.
 */
void lqn_srvnph() {
    const LqnS l = lqn_srvnph_model();
    double secs[2] = {0.0, 0.0};
    std::size_t nclasses[2] = {0, 0};
    const char* label[2] = {"default", "srvn.ph"};
    const char* method[2] = {"default", "srvn.ph"};

    for (int k = 0; k < 2; ++k) {
        ln::LnOptions opt;
        opt.method = method[k];
        const std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
        ln::SolverLN<double> s(l, opt);
        const ln::LnSolution<double> sol = s.get_ensemble_avg();
        secs[k] = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        for (std::size_t e = 0; e < s.layers().size(); ++e)
            nclasses[k] += s.layers()[e].nclasses;
        std::printf("\nLN(%s) Results [%.3f s]:\n", label[k], secs[k]);
        print_ln(l, sol, opt, s.nlayers());
    }
    std::printf("\nLayer classes: default %zu, srvn.ph %zu. Runtime: %.3fs vs %.3fs (%.2fx).\n",
                nclasses[0], nclasses[1], secs[0], secs[1],
                secs[1] > 0.0 ? secs[0] / secs[1] : 0.0);
}

// ---------------------------------------------------------------------------
// lqn_bpmn_trace
// ---------------------------------------------------------------------------

/**
 * The BPMN model with the fixed point traced iteration by iteration.
 *
 * The reference drives the loop by hand through `solver.pre/analyze/post`, which
 * this port has no counterpart for: `SolverLN` exposes no mid-iteration hooks.
 * It does RETAIN every iteration's per-layer result, though, and that is all the
 * reference prints -- so the trace below is the same numbers, read back after the
 * solve rather than during it. The reference caps its manual loop at five
 * iterations, and iterations one to five of a ten-iteration solve are the same
 * five iterations, so nothing is approximated here.
 */
void lqn_bpmn_trace() {
    note("=== Running LN solver with tracing ===");
    const LqnS l = lqn_bpmn_model();
    ln::LnOptions opt;
    opt.iter_max = 10;
    ln::SolverLN<double> s(l, opt);
    const ln::LnSolution<double> sol = s.get_ensemble_avg();

    const std::vector<std::vector<ln::LayerResult<double> > >& hist = s.iteration_results();
    const std::size_t shown = std::min<std::size_t>(5, hist.size());
    for (std::size_t it = 0; it < shown; ++it) {
        std::printf("=== Iteration %zu ===\n", it + 1);
        for (std::size_t e = 0; e < hist[it].size(); ++e) {
            const double tn = mat_max(hist[it][e].TN), qn = mat_max(hist[it][e].QN);
            if (tn > 1e-10 || qn > 1e-10)
                std::printf("  Layer %zu: max(TN)=%.4e, max(QN)=%.4e\n", e + 1, tn, qn);
        }
    }

    note("=== Final check: Getting AvgTable ===");
    const char* tasks[7] = {"R:R1_Task", "T:R2_Task",  "T:R3_Task", "T:R1A_Task",
                            "T:R1B_Task", "T:R2A_Task", "T:R2B_Task"};
    for (int k = 0; k < 7; ++k) {
        const std::size_t i = idx_of(l, tasks[k]);
        std::printf("%s: Tput=%.6f\n", l.names[i].c_str(), ln_sanitize(sol.TN[i]));
    }
}

// ---------------------------------------------------------------------------
// lqn_transient
// ---------------------------------------------------------------------------

/**
 * The per-layer transient of a two-task chain under fluid layers.
 *
 * NO HORIZON IS NAMED, as the reference names none: each layer's fluid solver
 * picks its own by the analyzer's rule, and the coupled path defers to the
 * decoupled one when there is no shared grid to relax over. So this reports the
 * DECOUPLED transient, with the inter-layer demands frozen at the fixed point --
 * which is exactly what the reference reports for the same input.
 *
 * The reference then PLOTS the traces. What is printed here instead is the table
 * `line-cli -a tran` prints, which carries the same series, so the example and
 * the CLI cannot disagree about the same model.
 */
void lqn_transient() {
    Lqn b;
    b.processor("P1", 1, SchedStrategy::PS);
    b.processor("P2", 1, SchedStrategy::PS);
    b.task("T1", 5, SchedStrategy::REF, "P1");
    b.think_time("T1", Er(1.0));
    b.task("T2", 5, SchedStrategy::FCFS, "P2");
    b.think_time("T2", Er(1.0));
    b.entry("E1", "T1");
    b.entry("E2", "T2");
    b.activity("A1", Er(2.0), "T1");
    b.bound_to("A1", "E1");
    b.sync_call("A1", "E2", 1.0);
    b.activity("A2", Er(3.0), "T2");
    b.bound_to("A2", "E2");
    b.replies_to("A2", "E2");
    const LqnS l = b.build();

    ln::LnOptions opt;
    opt.layer_solver = "fluid";
    ln::SolverLN<double> s(l, opt);
    const ln::LnTranSolution tr = s.get_tran_avg();
    const std::vector<qn::Layer<double> >& layers = s.layers();

    std::printf("LN.getTranAvg returned transient traces for %zu layers.\n", tr.layers.size());
    std::printf("SolverLN(SolverFluid) getTranAvg arith=double mode=%s layers=%zu "
                "iterations=%ld gap=%.3e\n",
                tr.mode.c_str(), tr.layers.size(), tr.iterations, tr.gap);
    for (std::size_t e = 0; e < tr.layers.size(); ++e) {
        const std::vector<double>& t = tr.layers[e].t;
        if (t.empty()) continue;
        std::printf("\nLayer %s  (%zu points on [%.6g, %.6g])\n", layers[e].name.c_str(), t.size(),
                    t.front(), t.back());
        std::printf("%-30s %-24s %12s %12s %12s %12s\n", "Station", "JobClass", "QLen(0)",
                    "QLen(end)", "Util(end)", "Tput(end)");
        for (std::size_t i = 0; i < tr.layers[e].QN.size(); ++i)
            for (std::size_t r = 0; r < tr.layers[e].QN[i].size(); ++r) {
                const std::vector<double>& q = tr.layers[e].QN[i][r];
                if (q.empty()) continue;
                std::printf("%-30s %-24s %12.6g %12.6g %12.6g %12.6g\n",
                            layers[e].stations[i].name.c_str(), layers[e].classes[r].name.c_str(),
                            ln_sanitize(q.front()), ln_sanitize(q.back()),
                            ln_sanitize(tr.layers[e].UN[i][r].back()),
                            ln_sanitize(tr.layers[e].TN[i][r].back()));
            }
    }
}

LINE_EXAMPLE("basic/layeredModel", lqn_basic);
LINE_EXAMPLE("basic/layeredModel", lqn_basic_debug);
LINE_EXAMPLE("basic/layeredModel", lqn_bpmn);
LINE_EXAMPLE("basic/layeredModel", lqn_bpmn_trace);
LINE_EXAMPLE("basic/layeredModel", lqn_flatph);
LINE_EXAMPLE("basic/layeredModel", lqn_srvnph);
LINE_EXAMPLE("basic/layeredModel", lqn_transient);
LINE_EXAMPLE("basic/layeredModel", lqn_debug2);
LINE_EXAMPLE("basic/layeredModel", lqn_fork_open_arrival);
LINE_EXAMPLE("basic/layeredModel", lqn_init);
LINE_EXAMPLE("basic/layeredModel", lqn_jsq);
LINE_EXAMPLE("basic/layeredModel", lqn_moment3);
LINE_EXAMPLE("basic/layeredModel", lqn_multi_solvers);
LINE_EXAMPLE("basic/layeredModel", lqn_ofbiz);
LINE_EXAMPLE("basic/layeredModel", lqn_open_arrival);
LINE_EXAMPLE("basic/layeredModel", lqn_paramident);
LINE_EXAMPLE("basic/layeredModel", lqn_rrobin);
LINE_EXAMPLE("basic/layeredModel", lqn_serial);
LINE_EXAMPLE("basic/layeredModel", lqn_server_pools);
LINE_EXAMPLE("basic/layeredModel", lqn_setup);
LINE_EXAMPLE("basic/layeredModel", lqn_sockshop);
LINE_EXAMPLE("basic/layeredModel", lqn_twotasks);
LINE_EXAMPLE("basic/layeredModel", lqn_workflows);

}  // namespace examples
}  // namespace line
