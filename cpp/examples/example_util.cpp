/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * The one translation unit that instantiates the solver stack for the examples.
 *
 * Every entry point of example_util.h forwards to the runner the CLI calls, so
 * an example and `line-cli -s <solver>` on the same model produce the same
 * numbers by construction. Nothing is computed here beyond the reshaping the
 * table needs.
 */

#include "example_util.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>

#include "line/api/mam/aph_fit.h"
#include "line/api/mam/aph_fit_moments.h"
#include "line/solvers/auto/solver_auto.h"
#include "line/solvers/ba/solver_ba_runner.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/ctmc/solver_ctmc_cdf.h"
#include "line/solvers/ctmc/solver_ctmc_getters.h"
#include "line/solvers/ctmc/solver_ctmc_prob.h"
#include "line/solvers/ctmc/solver_ctmc_waitq.h"
#include "line/solvers/fluid/fluid_runner.h"
#include "line/solvers/fluid/solver_fluid.h"
#include "line/solvers/mam/solver_mam_runner.h"
#include "line/solvers/mva/solver_mva_runner.h"
#include "line/solvers/nc/solver_nc_runner.h"
#include "line/solvers/ssa/ssa_dispatch.h"
#include "line/solvers/wrappers/jmt/jmt_logs.h"
#include "line/solvers/wrappers/jmt/solver_jmt.h"
#include "line/solvers/wrappers/ldes/solver_ldes.h"
#include "parity_recorder.h"

namespace line {
namespace examples {

namespace {

/** MATLAB prints a table with these columns; the widths are the CLI's. */
void print_header() {
    std::printf("%-16s %-14s %12s %12s %12s %12s %12s %12s\n", "Station", "JobClass", "QLen",
                "Util", "RespT", "ResidT", "ArvR", "Tput");
}

/** Fill the label columns and the per-class system metrics of a table. */
void fill_table(AvgTable& t, const qn::NetworkStruct<double>& sn, const mva::AvgResult<double>& r) {
    for (std::size_t i = 0; i < sn.nstations; ++i)
        for (std::size_t c = 0; c < sn.nclasses; ++c) {
            const double q = r.QN(i, c), u = r.UN(i, c), rt = r.RN(i, c);
            const double w = r.WN(i, c), a = r.AN(i, c), x = r.TN(i, c);
            if (q == 0.0 && u == 0.0 && rt == 0.0 && w == 0.0 && a == 0.0 && x == 0.0) continue;
            t.Station.push_back(sn.stations[i].name);
            t.JobClass.push_back(sn.classes[c].name);
            t.QLen.push_back(q);
            t.Util.push_back(u);
            t.RespT.push_back(rt);
            t.ResidT.push_back(w);
            t.ArvR.push_back(a);
            t.Tput.push_back(x);
        }
    for (std::size_t c = 0; c < sn.nclasses && c < r.CN.size(); ++c) {
        t.SysClass.push_back(sn.classes[c].name);
        t.SysRespT.push_back(r.CN[c]);
        t.SysTput.push_back(c < r.XN.size() ? r.XN[c] : 0.0);
    }
    t.iter = r.iter;
    t.warning = r.warning;
    t.ListCost = r.listcost;
    if (r.lognormconst.has_value()) {
        t.has_lognormconst = true;
        t.lognormconst = r.lognormconst.value();
    }
}

/**
 * The same fill for a solver whose solution type is not `mva::AvgResult`.
 *
 * THE SSA AND FLUID ARMS OWE BOTH DERIVED COLUMNS. Neither analyzer reports a
 * residence time or an arrival rate of its own, and neither is a copy of its
 * neighbour: ResidT is the per-JOB time and RespT the per-VISIT one, and they
 * agree only where every station is visited once per cycle -- reporting one for
 * the other was a factor of 3 out on `sdroute_closed` and 17 on Queue1 of
 * `init_state_ps`. ArvR is a FLOW and parts from the throughput at any station a
 * job leaves by another route. `sn_get_residt_from_respt` and
 * `sn_get_arvr_from_tput` are the reference's own conversions and are what the
 * CLI's `-s fluid` and `-s ssa` arms apply. A Source's ArvR stays zero: nothing
 * arrives TO it.
 */
template <class Sol>
void fill_table_sim(AvgTable& t, const qn::NetworkStruct<double>& sn, const Sol& r) {
    const std::size_t M = sn.nstations, K = sn.nclasses;
    Matrix<double> RN(M, K), TN(M, K);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t c = 0; c < K; ++c) {
            RN(i, c) = r.RN(i, c);
            TN(i, c) = r.TN(i, c);
        }
    const Matrix<double> WN = mva::sn_get_residt_from_respt<double>(sn, RN);
    const Matrix<double> AN = mva::sn_get_arvr_from_tput<double>(sn, TN);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t c = 0; c < K; ++c) {
            const double q = r.QN(i, c), u = r.UN(i, c), rt = r.RN(i, c), x = r.TN(i, c);
            const bool is_source = sn.stations[i].sched == lang::SchedStrategy::EXT;
            const double w = WN(i, c), a = is_source ? 0.0 : AN(i, c);
            // The all-zero row test over ALL SIX columns, as `avg_rows` and the
            // reference make it: a station only ARRIVALS reach still has a row.
            if (q == 0.0 && u == 0.0 && rt == 0.0 && w == 0.0 && a == 0.0 && x == 0.0) continue;
            t.Station.push_back(sn.stations[i].name);
            t.JobClass.push_back(sn.classes[c].name);
            t.QLen.push_back(q);
            t.Util.push_back(u);
            t.RespT.push_back(rt);
            t.ResidT.push_back(w);
            t.ArvR.push_back(a);
            t.Tput.push_back(x);
        }
    for (std::size_t c = 0; c < sn.nclasses && c < r.CN.size(); ++c) {
        t.SysClass.push_back(sn.classes[c].name);
        t.SysRespT.push_back(r.CN[c]);
        t.SysTput.push_back(c < r.XN.size() ? r.XN[c] : 0.0);
    }
}

/**
 * The LDES fill. Its solution type carries the six metrics separately and its
 * per-class CN/XN are (1 x nclasses) MATRICES rather than vectors, so neither
 * `fill_table` nor `fill_table_sim` can read it.
 *
 * ResidT IS RECOMPUTED RATHER THAN TAKEN, which is what the CLI's `-a avg` LDES
 * arm does and for the same reason: the engine counts ONE VISIT PER STATION, so
 * the residence time it reports is its response time wherever a visit ratio is
 * not 1. On `sdroute_jsq`, whose three queues each take a third of the flow,
 * that was a factor of three on every one of them.
 */
void fill_table_ldes(AvgTable& t, const qn::NetworkStruct<double>& sn,
                     const ldes::LdesResult& r) {
    const Matrix<double> WNfix = mva::sn_get_residt_from_respt<double>(sn, r.RN);
    for (std::size_t i = 0; i < sn.nstations; ++i)
        for (std::size_t c = 0; c < sn.nclasses; ++c) {
            const double q = r.QN(i, c), u = r.UN(i, c), rt = r.RN(i, c);
            const double w = WNfix(i, c), a = r.AN(i, c), x = r.TN(i, c);
            if (q == 0.0 && u == 0.0 && rt == 0.0 && w == 0.0 && a == 0.0 && x == 0.0) continue;
            t.Station.push_back(sn.stations[i].name);
            t.JobClass.push_back(sn.classes[c].name);
            t.QLen.push_back(q);
            t.Util.push_back(u);
            t.RespT.push_back(rt);
            t.ResidT.push_back(w);
            t.ArvR.push_back(a);
            t.Tput.push_back(x);
        }
    for (std::size_t c = 0; c < sn.nclasses && c < r.CN.cols(); ++c) {
        t.SysClass.push_back(sn.classes[c].name);
        t.SysRespT.push_back(r.CN(0, c));
        t.SysTput.push_back(c < r.XN.cols() ? r.XN(0, c) : 0.0);
    }
}

ctmc::CtmcOptions ctmc_options(const SolverOpts& o) {
    ctmc::CtmcOptions c;
    if (!o.method.empty()) c.method = o.method;
    if (o.cutoff > 0.0) c.cutoff = o.cutoff;
    if (!o.cutoff_vec.empty()) c.cutoff_vec = o.cutoff_vec;
    if (o.state_max) c.state_max = o.state_max;
    return c;
}

fluid::FluidOptions fluid_options(const SolverOpts& o) {
    fluid::FluidOptions f;
    if (!o.method.empty()) f.method = o.method;
    if (o.tol >= 0.0) f.tol = o.tol;
    if (o.iter_tol >= 0.0) f.iter_tol = o.iter_tol;
    if (o.iter_max >= 0) f.iter_max = static_cast<std::size_t>(o.iter_max);
    if (o.timespan_end > 0.0) f.timespan_end = o.timespan_end;
    if (!o.init_sol.empty()) f.init_sol = o.init_sol;
    if (o.seed) f.seed = o.seed;
    if (!o.highvar.empty()) f.highvar = o.highvar;
    f.stiff = o.stiff;
    return f;
}

}  // namespace

// ---------------------------------------------------------------------------
// Model-building helpers
// ---------------------------------------------------------------------------

void serial_routing(Routing& P, std::size_t r, std::size_t s,
                    const std::vector<std::size_t>& nodes) {
    for (std::size_t k = 0; k + 1 < nodes.size(); ++k) P.set(r, s, nodes[k], nodes[k + 1], 1.0);
}

void serial_routing(Routing& P, const std::vector<std::size_t>& nodes) {
    serial_routing(P, 1, 1, nodes);
}

void link_serial(Net& m, const std::vector<std::size_t>& nodes) {
    Routing P;
    const std::size_t K = m.raw_struct().classes.size();
    for (std::size_t r = 1; r <= K; ++r) {
        for (std::size_t k = 0; k + 1 < nodes.size(); ++k) P.set(r, r, nodes[k], nodes[k + 1], 1.0);
        // MATLAB's serialRouting CLOSES the ring unless the chain ends at a Sink,
        // which is how a closed model built from a node list circulates.
        const qn::NetworkStruct<double>& sn = m.raw_struct();
        if (!nodes.empty() && sn.nodes[nodes.back() - 1].nodetype != lang::NodeType::Sink)
            P.set(r, r, nodes.back(), nodes.front(), 1.0);
    }
    m.link(P);
}

D hyperexp_fit(double mean, double scv) {
    // HyperExp.fitMeanAndSCV: the balanced-mean two-phase fit of MATLAB's
    // Coxian/HyperExp family, p mu1^-1 = (1-p) mu2^-1 is NOT imposed; the
    // reference solves for (p, mu1, mu2) from the first two moments with the
    // third degree of freedom fixed by p = 0.5 * (1 - sqrt((scv-1)/(scv+1))).
    if (!(scv > 1.0))
        throw InputError("HyperExp.fitMeanAndSCV: the SCV of a hyperexponential exceeds one");
    const double p = 0.5 * (1.0 - std::sqrt((scv - 1.0) / (scv + 1.0)));
    const double mu1 = 2.0 * p / mean;
    const double mu2 = 2.0 * (1.0 - p) / mean;
    return D::hyperexp(p, mu1, mu2);
}

D erlang_fit(double mean, double scv) { return D::erlang_fit(mean, scv); }

D erlang_fit_order(double mean, std::size_t k) {
    if (k == 0) throw InputError("Erlang.fitMeanAndOrder: the order must be at least one");
    return D::erlang(static_cast<double>(k) / mean, k);
}

D aph_fit(double mean, double scv) {
    const mam::Map<double> M = mam::aph_fit_mean_scv(mean, scv);
    return D::phase_type(mam::map_pie(M), M.D0, true);
}

D replayer_from_file(const std::string& path) {
    std::ifstream in(path.c_str());
    if (!in) throw InputError("Replayer: cannot open the trace file '" + path + "'");
    std::vector<double> samples;
    double v = 0.0;
    while (in >> v) samples.push_back(v);
    if (samples.empty()) throw InputError("Replayer: the trace file '" + path + "' is empty");
    // The PATH travels with the samples: JMT is handed a `ReplayerPar` naming a
    // file and cannot take a trace inline.
    return D::replayer_from(samples, path);
}

std::vector<double> zipf(double alpha, std::size_t n) {
    std::vector<double> p(n, 0.0);
    double z = 0.0;
    for (std::size_t i = 1; i <= n; ++i) z += std::pow(static_cast<double>(i), -alpha);
    for (std::size_t i = 1; i <= n; ++i) p[i - 1] = std::pow(static_cast<double>(i), -alpha) / z;
    return p;
}

std::string line_root_folder() {
#ifdef LINE_MP_REPO_ROOT
    return std::string(LINE_MP_REPO_ROOT);
#else
    return std::string(".");
#endif
}

// ---------------------------------------------------------------------------
// The average table
// ---------------------------------------------------------------------------

double AvgTable::get(const std::string& col, const std::string& station,
                     const std::string& jobclass) const {
    for (std::size_t i = 0; i < Station.size(); ++i) {
        if (Station[i] != station || JobClass[i] != jobclass) continue;
        if (col == "QLen") return QLen[i];
        if (col == "Util") return Util[i];
        if (col == "RespT") return RespT[i];
        if (col == "ResidT") return ResidT[i];
        if (col == "ArvR") return ArvR[i];
        if (col == "Tput") return Tput[i];
        throw InputError("AvgTable::get: unknown column '" + col + "'");
    }
    return std::numeric_limits<double>::quiet_NaN();
}

std::vector<double> AvgTable::column(const std::string& col) const {
    if (col == "QLen") return QLen;
    if (col == "Util") return Util;
    if (col == "RespT") return RespT;
    if (col == "ResidT") return ResidT;
    if (col == "ArvR") return ArvR;
    if (col == "Tput") return Tput;
    throw InputError("AvgTable::column: unknown column '" + col + "'");
}

AvgTable solve_avg(const std::string& solver, Net& m, const SolverOpts& o) {
    AvgTable t;
    t.solver = solver;
    const qn::NetworkStruct<double>& sn = m.get_struct();

    if (solver == "LQNS" || solver == "QNS")
        throw UnsupportedError("solve_avg: '" + solver +
                               "' analyses a LayeredNetwork rather than a Network; the example "
                               "carries the call commented out");

    std::string name = solver;
    if (name == "AUTO") {
        // `SolverAuto` picks by feature set, exactly as the CLI's own auto arm.
        name = autosolver::auto_solver_name(autosolver::auto_choose_avg_solver(sn));
        std::transform(name.begin(), name.end(), name.begin(), ::toupper);
        if (name == "FLUID") name = "FLD";
    }

    if (name == "MVA") {
        mva::MvaOptions opt;
        if (!o.method.empty()) opt.method = o.method;
        if (o.tol >= 0.0) opt.tol = o.tol;
        if (o.iter_tol >= 0.0) opt.iter_tol = o.iter_tol;
        if (o.iter_max >= 0) opt.iter_max = o.iter_max;
        if (!o.multiserver.empty()) opt.multiserver = o.multiserver;
        if (!o.highvar.empty()) opt.highvar = o.highvar;
        if (!o.np_priority.empty()) opt.np_priority = o.np_priority;
        if (!o.fork_join.empty()) opt.fork_join = o.fork_join;
        Matrix<double> init;
        const mva::AvgResult<double> r = mva::solver_mva_run_analyzer(sn, opt, init);
        t.method = r.actualmethod;
        fill_table(t, sn, r);
    } else if (name == "NC") {
        nc::NcSolverOptions opt;
        if (!o.method.empty()) opt.method = o.method;
        if (o.tol >= 0.0) opt.tol = o.tol;
        if (o.iter_tol >= 0.0) opt.iter_tol = o.iter_tol;
        if (o.iter_max >= 0) opt.iter_max = o.iter_max;
        if (!o.highvar.empty()) opt.highvar = o.highvar;
        if (!o.multiserver.empty()) opt.multiserver = o.multiserver;
        if (!o.fork_join.empty()) opt.fork_join = o.fork_join;
        if (o.samples) opt.samples = o.samples;
        if (o.seed) opt.seed = o.seed;
        const mva::AvgResult<double> r = nc::solver_nc_run_analyzer(sn, opt);
        t.method = r.actualmethod;
        fill_table(t, sn, r);
    } else if (name == "CTMC") {
        const ctmc::CtmcOptions opt = ctmc_options(o);
        const mva::AvgResult<double> r = ctmc::solver_ctmc_run_analyzer_any(sn, opt);
        t.method = r.actualmethod;
        fill_table(t, sn, r);
    } else if (name == "MAM") {
        mam::MamOptions opt;
        if (!o.method.empty()) opt.method = o.method;
        if (o.tol >= 0.0) opt.tol = o.tol;
        if (o.iter_max >= 0) opt.iter_max = o.iter_max;
        const mva::AvgResult<double> r = mam::solver_mam_run_analyzer(sn, opt);
        t.method = r.actualmethod;
        fill_table(t, sn, r);
    } else if (name == "BA") {
        ba::BaOptions opt;
        if (!o.method.empty()) opt.method = o.method;
        opt.level = o.level;
        const mva::AvgResult<double> r = ba::solver_ba_run_analyzer(sn, opt);
        t.method = r.actualmethod;
        fill_table(t, sn, r);
    } else if (name == "SSA") {
        ssa::SsaOptions opt;
        if (!o.method.empty()) opt.method = o.method;
        if (o.samples) opt.samples = o.samples;
        if (o.seed) opt.seed = o.seed;
        if (!o.state_space_gen.empty()) opt.state_space_gen = o.state_space_gen;
        opt.verbose = o.verbose;
        const ssa::SsaSolution r = ssa::solver_ssa(sn, opt);
        t.method = r.method;
        fill_table_sim(t, sn, r);
    } else if (name == "FLD" || name == "FLUID") {
        const fluid::FluidOptions opt = fluid_options(o);
        const fluid::FluidSolution r = fluid::solver_fluid_run_analyzer(sn, opt);
        t.method = r.method;
        t.iter = static_cast<int>(r.iters);
        fill_table_sim(t, sn, r);
    } else if (name == "JMT") {
        jmt::JmtOptions opt;
        if (!o.method.empty()) opt.method = o.method;
        if (o.samples) opt.samples = static_cast<double>(o.samples);
        if (o.seed) opt.seed = static_cast<long>(o.seed);
        opt.keep = o.keep;
        opt.verbose = o.verbose;
        const jmt::JmtResult<double> r = jmt::solver_jmt_run_analyzer(sn, opt);
        t.method = r.avg.actualmethod;
        fill_table(t, sn, r.avg);
    } else if (name == "LDES") {
        ldes::LdesOptions opt;
        if (!o.method.empty()) opt.method = o.method;
        if (o.samples) opt.samples = o.samples;
        if (o.seed) opt.seed = static_cast<long>(o.seed);
        opt.verbose = o.verbose;
        const ldes::LdesResult r = ldes::solver_ldes(sn, opt);
        t.method = opt.method;
        fill_table_ldes(t, sn, r);
    } else {
        throw UnsupportedError("solve_avg: unknown solver '" + solver + "'");
    }
    return t;
}

SolverOpts sim_opts(unsigned long seed, std::size_t samples) {
    SolverOpts o;
    o.seed = seed;
    o.samples = samples;
    return o;
}

WrapperAvg jmt_avg(Net& m, const SolverOpts& o) {
    const qn::NetworkStruct<double>& sn = m.get_struct();
    jmt::JmtOptions opt;
    if (!o.method.empty()) opt.method = o.method;
    if (o.samples) opt.samples = static_cast<double>(o.samples);
    if (o.seed) opt.seed = static_cast<long>(o.seed);
    opt.keep = o.keep;
    opt.verbose = o.verbose;
    const jmt::JmtResult<double> r = jmt::solver_jmt_run_analyzer(sn, opt);
    WrapperAvg w;
    w.QN = r.avg.QN;
    w.UN = r.avg.UN;
    w.RN = r.avg.RN;
    w.WN = r.avg.WN;
    w.AN = r.avg.AN;
    w.TN = r.avg.TN;
    w.method = r.avg.actualmethod;
    return w;
}

WrapperAvg ldes_avg(Net& m, const SolverOpts& o) {
    const qn::NetworkStruct<double>& sn = m.get_struct();
    ldes::LdesOptions opt;
    if (!o.method.empty()) opt.method = o.method;
    if (o.samples) opt.samples = o.samples;
    if (o.seed) opt.seed = static_cast<long>(o.seed);
    opt.verbose = o.verbose;
    const ldes::LdesResult r = ldes::solver_ldes(sn, opt);
    WrapperAvg w;
    w.QN = r.QN;
    w.UN = r.UN;
    w.RN = r.RN;
    // Recomputed, not taken: the engine counts one visit per station, so its own
    // residence time is its response time wherever a visit ratio is not 1. The
    // CLI's `-a avg` LDES arm makes the same substitution.
    w.WN = mva::sn_get_residt_from_respt<double>(sn, r.RN);
    w.AN = r.AN;
    w.TN = r.TN;
    w.method = opt.method;
    return w;
}

void print_avg(const AvgTable& t, const std::string& caption) {
    namespace parity = line::examples::parity;
    if (!caption.empty()) std::printf("%s\n", caption.c_str());
    std::printf("Solver%s method=%s\n", t.solver.c_str(), t.method.c_str());
    if (!t.warning.empty()) std::printf("Warning: %s\n", t.warning.c_str());
    print_header();
    // This printer does NOT go through avg_rows, so it records for itself. The
    // METHOD is taken from the table, which knows which one the solver actually
    // resolved -- `MAM(dec.source)` and `MAM(inap)` differ by 236% on the same
    // model, so it is not decoration.
    if (parity::enabled()) {
        // THE TABLE NAMES ITS OWN SOLVER, which is attribution a reference that
        // prints no banner cannot otherwise give. It yields to `section()`; see
        // parity_recorder.h.
        parity::set_solver_from_table(t.solver);
        parity::set_method(t.method);
        parity::begin_table("avg", "Station", "JobClass");
    }
    for (std::size_t i = 0; i < t.Station.size(); ++i) {
        std::printf("%-16s %-14s %12.6g %12.6g %12.6g %12.6g %12.6g %12.6g\n", t.Station[i].c_str(),
                    t.JobClass[i].c_str(), t.QLen[i], t.Util[i], t.RespT[i], t.ResidT[i], t.ArvR[i],
                    t.Tput[i]);
        if (!parity::enabled()) continue;
        std::vector<parity::Cell> cells;
        cells.push_back(parity::Cell{"QLen", t.QLen[i]});
        cells.push_back(parity::Cell{"Util", t.Util[i]});
        cells.push_back(parity::Cell{"RespT", t.RespT[i]});
        cells.push_back(parity::Cell{"ResidT", t.ResidT[i]});
        cells.push_back(parity::Cell{"ArvR", t.ArvR[i]});
        cells.push_back(parity::Cell{"Tput", t.Tput[i]});
        parity::add_row(t.Station[i], t.JobClass[i], cells);
    }
    if (!t.ListCost.empty()) {
        std::printf("%-16s", "ListCost");
        for (std::size_t j = 0; j < t.ListCost.size(); ++j) std::printf(" %12.6g", t.ListCost[j]);
        std::printf("\n");
    }
}

void print_sys(const AvgTable& t, const std::string& caption) {
    namespace parity = line::examples::parity;
    if (!caption.empty()) std::printf("%s\n", caption.c_str());
    std::printf("%-16s %14s %14s\n", "JobClass", "SysRespT", "SysTput");
    if (parity::enabled()) parity::begin_table("sys", "Chain", "JobClass");
    for (std::size_t c = 0; c < t.SysClass.size(); ++c) {
        std::printf("%-16s %14.6g %14.6g\n", t.SysClass[c].c_str(), t.SysRespT[c], t.SysTput[c]);
        if (!parity::enabled()) continue;
        std::vector<parity::Cell> cells;
        cells.push_back(parity::Cell{"RespT", t.SysRespT[c]});
        cells.push_back(parity::Cell{"Tput", t.SysTput[c]});
        parity::add_row(t.SysClass[c], std::string(), cells);
    }
}

void banner(const std::string& text) { std::printf("== %s ==\n", text.c_str()); }

void print_scalar(const std::string& label, double v) {
    std::printf("%s = %.10g\n", label.c_str(), v);
    // The derived goldens -- a state probability, a phase-type moment, a cache
    // hit rate -- reach the record through here, because a twin computes them
    // itself and no table carries them.
    line::examples::parity::add_scalar(label, v);
}

void print_vector(const std::string& label, const std::vector<double>& v) {
    std::printf("%s =", label.c_str());
    for (std::size_t i = 0; i < v.size(); ++i) std::printf(" %.6g", v[i]);
    std::printf("\n");
}

void print_matrix(const std::string& label, const Matrix<double>& M) {
    std::printf("%s = [%zux%zu]\n", label.c_str(), M.rows(), M.cols());
    for (std::size_t i = 0; i < M.rows(); ++i) {
        for (std::size_t j = 0; j < M.cols(); ++j) std::printf(" %10.6g", M(i, j));
        std::printf("\n");
    }
}

Matrix<double> residt_from_respt(const qn::NetworkStruct<double>& sn,
                                 const Matrix<double>& RN) {
    return mva::sn_get_residt_from_respt<double>(sn, RN);
}

Matrix<double> arvr_from_tput(const qn::NetworkStruct<double>& sn,
                              const Matrix<double>& TN) {
    return mva::sn_get_arvr_from_tput<double>(sn, TN);
}

// ---------------------------------------------------------------------------
// The other analyses
// ---------------------------------------------------------------------------

double ctmc_prob_aggr(Net& m, std::size_t node, const std::vector<double>& state,
                      const SolverOpts& o) {
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const ctmc::CtmcOptions opt = ctmc_options(o);
    const ctmc::CtmcSolution<double> d = ctmc::solver_ctmc_analyzer(sn, opt);
    const Matrix<double> A = ctmc::ctmc_get_state_space_aggr(sn, d);
    const std::size_t ist = sn.nodes[node - 1].station;
    if (ist == 0) throw InputError("getProbAggr: the node is not a station");
    const std::size_t K = sn.nclasses;
    if (state.size() != K)
        throw InputError("getProbAggr: the state must have one entry per class");
    double p = 0.0;
    for (std::size_t s = 0; s < A.rows(); ++s) {
        bool hit = true;
        for (std::size_t k = 0; k < K && hit; ++k) hit = A(s, (ist - 1) * K + k) == state[k];
        if (hit) p += d.pi[s];
    }
    return p;
}

std::vector<double> ctmc_marg_aggr(Net& m, std::size_t node, const SolverOpts& o) {
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const ctmc::CtmcOptions opt = ctmc_options(o);
    const ctmc::CtmcSolution<double> d = ctmc::solver_ctmc_analyzer(sn, opt);
    const Matrix<double> A = ctmc::ctmc_get_state_space_aggr(sn, d);
    const std::size_t ist = sn.nodes[node - 1].station;
    if (ist == 0) throw InputError("getProbStateAggr: the node is not a station");
    const std::size_t K = sn.nclasses;
    std::size_t nmax = 0;
    for (std::size_t s = 0; s < A.rows(); ++s) {
        double tot = 0.0;
        for (std::size_t k = 0; k < K; ++k) tot += A(s, (ist - 1) * K + k);
        nmax = std::max(nmax, static_cast<std::size_t>(tot));
    }
    std::vector<double> pmf(nmax + 1, 0.0);
    for (std::size_t s = 0; s < A.rows(); ++s) {
        double tot = 0.0;
        for (std::size_t k = 0; k < K; ++k) tot += A(s, (ist - 1) * K + k);
        pmf[static_cast<std::size_t>(tot)] += d.pi[s];
    }
    return pmf;
}

Matrix<double> ctmc_state_space(Net& m, const SolverOpts& o) {
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const ctmc::CtmcOptions opt = ctmc_options(o);
    const ctmc::CtmcSolution<double> d = ctmc::solver_ctmc_analyzer(sn, opt);
    return ctmc::ctmc_get_state_space(sn, d).flat;
}

Matrix<double> ctmc_generator(Net& m, const SolverOpts& o) {
    const qn::NetworkStruct<double>& sn = m.get_struct();
    ctmc::CtmcOptions opt = ctmc_options(o);
    opt.keep_filtration = true;
    const ctmc::CtmcSolution<double> d = ctmc::solver_ctmc_analyzer(sn, opt);
    return ctmc::ctmc_get_infgen(sn, d).Q;
}

void print_inf_gen(const Matrix<double>& Q, const Matrix<double>& space) {
    std::printf("%8s  %-28s %s\n", "State", "Marking", "Rates (to: rate)");
    for (std::size_t i = 0; i < Q.rows(); ++i) {
        std::string mark;
        if (i < space.rows()) {
            std::ostringstream os;
            os << "[";
            for (std::size_t c = 0; c < space.cols(); ++c) os << (c ? " " : "") << space(i, c);
            os << "]";
            mark = os.str();
        }
        std::printf("%8zu  %-28s", i + 1, mark.c_str());
        for (std::size_t j = 0; j < Q.cols(); ++j)
            if (i != j && Q(i, j) != 0.0) std::printf("  %zu: %.6g", j + 1, Q(i, j));
        std::printf("\n");
    }
}

std::vector<std::vector<CdfCurve> > cdf_respt(const std::string& solver, Net& m,
                                              const SolverOpts& o) {
    const qn::NetworkStruct<double>& sn = m.get_struct();
    std::vector<std::vector<CdfCurve> > out(sn.nstations, std::vector<CdfCurve>(sn.nclasses));
    if (solver == "CTMC") {
        const std::vector<std::vector<ctmc::CdfCurve<double> > > R =
            ctmc::solver_ctmc_cdf_respt(sn, ctmc_options(o));
        for (std::size_t i = 0; i < R.size() && i < out.size(); ++i)
            for (std::size_t c = 0; c < R[i].size() && c < out[i].size(); ++c) {
                out[i][c].t = R[i][c].t;
                out[i][c].F = R[i][c].F;
            }
    } else if (solver == "FLD" || solver == "FLUID") {
        const std::vector<std::vector<fluid::FluidPassage> > R =
            fluid::solver_fluid_cdf_respt(sn, fluid_options(o));
        for (std::size_t i = 0; i < R.size() && i < out.size(); ++i)
            for (std::size_t c = 0; c < R[i].size() && c < out[i].size(); ++c) {
                out[i][c].t = R[i][c].t;
                out[i][c].F = R[i][c].cdf;
            }
    } else if (solver == "JMT") {
        // The EMPIRICAL law, read back out of the JMT arrival/departure logs --
        // the same measurement `line-cli -s jmt -a cdf` reports. The reference
        // examples cross-check the fluid response-time law against this one, so
        // without it the twin refuses the solver its golden was produced by.
        jmt::JmtOptions jopt;
        if (!o.method.empty()) jopt.method = o.method;
        if (o.samples) jopt.samples = static_cast<double>(o.samples);
        if (o.seed) jopt.seed = static_cast<long>(o.seed);
        jopt.keep = o.keep;
        jopt.verbose = o.verbose;
        const std::map<std::pair<std::size_t, std::size_t>,
                       std::vector<std::pair<double, double> > >
            R = jmt::jmt_get_cdf_resp_t(sn, jopt);
        for (std::map<std::pair<std::size_t, std::size_t>,
                      std::vector<std::pair<double, double> > >::const_iterator it = R.begin();
             it != R.end(); ++it) {
            const std::size_t i = it->first.first, c = it->first.second;
            if (i == 0 || c == 0 || i > out.size() || c > out[i - 1].size()) continue;
            // The map holds `ecdf`'s (F, x) pairs, in that order.
            for (std::size_t k = 0; k < it->second.size(); ++k) {
                out[i - 1][c - 1].F.push_back(it->second[k].first);
                out[i - 1][c - 1].t.push_back(it->second[k].second);
            }
        }
    } else {
        throw UnsupportedError("getCdfRespT: no C++ port for solver '" + solver + "'");
    }
    return out;
}

void print_cdf_summary(const std::string& label, const CdfCurve& c) {
    if (c.t.empty()) {
        std::printf("%s: (no law)\n", label.c_str());
        return;
    }
    const double pct[5] = {0.5, 0.75, 0.9, 0.95, 0.99};
    std::printf("%s: points=%zu", label.c_str(), c.t.size());
    for (std::size_t p = 0; p < 5; ++p) {
        double q = c.t.back();
        for (std::size_t i = 0; i < c.F.size(); ++i)
            if (c.F[i] >= pct[p]) {
                q = c.t[i];
                break;
            }
        std::printf("  p%02d=%.6g", static_cast<int>(pct[p] * 100), q);
    }
    std::printf("\n");
}

TranAvg fluid_tran_avg(Net& m, const SolverOpts& o) {
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::vector<fluid::FluidTranPoint> pts =
        fluid::solver_fluid_tran_avg(sn, fluid_options(o), 100);
    TranAvg out;
    for (std::size_t i = 0; i < sn.nstations; ++i)
        for (std::size_t c = 0; c < sn.nclasses; ++c)
            out.label.push_back(sn.stations[i].name + "/" + sn.classes[c].name);
    for (std::size_t s = 0; s < pts.size(); ++s) {
        out.t.push_back(pts[s].t);
        std::vector<double> row;
        for (std::size_t i = 0; i < sn.nstations; ++i)
            for (std::size_t c = 0; c < sn.nclasses; ++c) row.push_back(pts[s].QN(i, c));
        out.QNt.push_back(row);
    }
    return out;
}

std::vector<double> fluid_initsol_from_marginal(Net& m,
                                                const std::vector<std::vector<double> >& n) {
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const fluid::FluidLayout L = fluid::fluid_layout(sn);
    std::vector<double> y0(L.nstates, 0.0);
    for (std::size_t i = 0; i < sn.nstations && i < n.size(); ++i)
        for (std::size_t r = 0; r < sn.nclasses && r < n[i].size(); ++r)
            // `L.enabled` is load-bearing: a disabled pair's `qidx` is the NEXT
            // pair's start, so writing it would corrupt a neighbouring block.
            if (L.enabled[i][r] && n[i][r] > 0.0) y0[L.qidx[i][r]] = n[i][r];
    return y0;
}

std::vector<double> fluid_initsol_from_state_prior(Net& m,
                                                   const std::vector<std::vector<double> >& n,
                                                   std::size_t prior_station) {
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const fluid::FluidLayout L = fluid::fluid_layout(sn);
    std::vector<double> y0 = fluid_initsol_from_marginal(m, n);
    const std::size_t K = sn.nclasses;
    if (prior_station == 0 || prior_station > sn.nstations || prior_station - 1 >= n.size())
        return y0;
    const std::size_t i = prior_station - 1;

    std::vector<std::size_t> ni(K, 0), ph(K, 1), shift(K, 0);
    std::size_t w = 0;
    for (std::size_t r = 0; r < K; ++r) {
        if (r < n[i].size() && n[i][r] > 0.0) ni[r] = static_cast<std::size_t>(n[i][r] + 0.5);
        ph[r] = sn.phases_of(prior_station, r + 1);
        shift[r] = w;
        w += ph[r];
    }
    const std::vector<std::vector<double> > rows = qn::from_marginal(sn, prior_station, ni, ph);
    if (rows.empty()) return y0;

    // The UNIFORM prior's first moment: every state carrying the marginal is
    // equally likely, so the ODE starts at the AVERAGE per-phase occupancy over
    // them rather than at the first row's. That is the whole difference between
    // this and initFromMarginal, and it is small (1.5e-4 on init_state_fcfs_nonexp)
    // precisely because only the phase assignment moves.
    std::vector<std::vector<double> > kbar(K);
    for (std::size_t r = 0; r < K; ++r) kbar[r].assign(ph[r], 0.0);
    std::size_t used = 0;
    for (std::size_t s = 0; s < rows.size(); ++s) {
        qn::Marginal<double> mg;
        try {
            mg = qn::to_marginal(sn, prior_station, rows[s], ph, shift,
                                 sn.nvars_of(sn.node_of_station(prior_station)));
        } catch (const Error&) {
            continue;  // a row the encoding cannot decode is not a state
        }
        if (mg.kir.size() != K) continue;
        for (std::size_t r = 0; r < K; ++r)
            for (std::size_t k = 0; k < kbar[r].size() && k < mg.kir[r].size(); ++k)
                kbar[r][k] += mg.kir[r][k];
        ++used;
    }
    if (used == 0) return y0;

    for (std::size_t r = 0; r < K; ++r) {
        if (!L.enabled[i][r]) continue;
        const std::size_t np = L.kic[i][r];
        double served = 0.0;
        for (std::size_t k = 1; k < kbar[r].size(); ++k) served += kbar[r][k] / used;
        for (std::size_t k = 0; k < np; ++k) {
            const double v = (k == 0) ? (r < n[i].size() ? n[i][r] : 0.0) - served
                                      : (k < kbar[r].size() ? kbar[r][k] / used : 0.0);
            y0[L.qidx[i][r] + k] = v;
        }
    }
    return y0;
}

BoundsTable ba_bounds(Net& m, const SolverOpts& o) {
    const qn::NetworkStruct<double>& sn = m.get_struct();
    ba::BaOptions opt;
    if (!o.method.empty()) opt.method = o.method;
    opt.level = o.level;
    const ba::BaBounds<double> b = ba::ba_bounds(sn, opt);
    BoundsTable t;
    t.method = opt.method;
    for (std::size_t i = 0; i < sn.nstations; ++i)
        for (std::size_t c = 0; c < sn.nclasses; ++c) {
            if (i < b.keep.size() && c < b.keep[i].size() && !b.keep[i][c]) continue;
            t.Station.push_back(sn.stations[i].name);
            t.JobClass.push_back(sn.classes[c].name);
            t.Qlower.push_back(b.Qlower(i, c));
            t.Qupper.push_back(b.Qupper(i, c));
            t.Tlower.push_back(b.Tlower(i, c));
            t.Tupper.push_back(b.Tupper(i, c));
        }
    return t;
}

void print_bounds(const BoundsTable& b, const std::string& caption) {
    if (!caption.empty()) std::printf("%s\n", caption.c_str());
    std::printf("%-16s %-14s %12s %12s %12s %12s\n", "Station", "JobClass", "Qlower", "Qupper",
                "Tlower", "Tupper");
    for (std::size_t i = 0; i < b.Station.size(); ++i)
        std::printf("%-16s %-14s %12.6g %12.6g %12.6g %12.6g\n", b.Station[i].c_str(),
                    b.JobClass[i].c_str(), b.Qlower[i], b.Qupper[i], b.Tlower[i], b.Tupper[i]);
}

std::vector<std::string> list_valid_methods(const std::string& solver, Net& m) {
    const qn::NetworkStruct<double>& sn = m.get_struct();
    if (solver == "MVA") return mva::list_valid_methods(sn);
    if (solver == "NC") return nc::list_valid_methods();
    if (solver == "CTMC") return ctmc::list_valid_methods();
    if (solver == "MAM") return mam::list_valid_methods();
    if (solver == "BA") return ba::list_valid_methods();
    throw UnsupportedError("listValidMethods: unknown solver '" + solver + "'");
}

}  // namespace examples
}  // namespace line
