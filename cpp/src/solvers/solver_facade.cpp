/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * @file
 * The one translation unit that instantiates the solver stack for the facade.
 *
 * `line/solvers/solver.h` declares a `double`-only, non-template solver API so
 * that this file can carry every body: a caller including that header pays for
 * none of the template instantiation below, and `line_mp_api` pays for it once
 * for the CLI, the tests and the examples together.
 *
 * Every entry point forwards to the runner the CLI calls, so a program written
 * against the facade and `line-cli -s <solver>` on the same model produce the
 * same numbers by construction. Nothing is computed here beyond the reshaping
 * the tables need.
 */

#include "line/solvers/solver.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <utility>

#include "line/solvers/auto/auto_methods.h"
#include "line/solvers/auto/solver_auto.h"
#include "line/lang/qn/nodes.h"
#include "line/solvers/ba/solver_ba_runner.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/ctmc/solver_ctmc_cdf.h"
#include "line/solvers/ctmc/solver_ctmc_getters.h"
#include "line/solvers/ctmc/solver_ctmc_prob.h"
#include "line/solvers/ctmc/solver_ctmc_symbolic.h"
#include "line/solvers/ctmc/solver_ctmc_waitq.h"
#include "line/solvers/fluid/fluid_runner.h"
#include "line/solvers/fluid/solver_fluid.h"
#include "line/solvers/map_env.h"
#include "line/solvers/map_env_stages.h"
#include "line/solvers/mam/solver_mam_runner.h"
#include "line/solvers/mva/solver_mva_runner.h"
#include "line/solvers/nc/solver_nc_runner.h"
#include "line/solvers/ssa/ssa_dispatch.h"
#include "line/solvers/wrappers/jmt/jmt_logs.h"
#include "line/solvers/wrappers/jmt/solver_jmt.h"
#include "line/solvers/wrappers/ldes/solver_ldes.h"

namespace line {

// ---------------------------------------------------------------------------
// AvgTable
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

namespace {

template <class Keep>
AvgTable filter_avg_table(const AvgTable& source, Keep keep) {
    AvgTable out;
    out.solver = source.solver;
    out.method = source.method;
    out.iter = source.iter;
    out.has_lognormconst = source.has_lognormconst;
    out.lognormconst = source.lognormconst;
    out.ListCost = source.ListCost;
    out.warning = source.warning;
    out.SysClass = source.SysClass;
    out.SysRespT = source.SysRespT;
    out.SysTput = source.SysTput;
    for (std::size_t i = 0; i < source.Station.size(); ++i) {
        if (!keep(i)) continue;
        out.Station.push_back(source.Station[i]);
        out.JobClass.push_back(source.JobClass[i]);
        out.QLen.push_back(source.QLen[i]);
        out.Util.push_back(source.Util[i]);
        out.RespT.push_back(source.RespT[i]);
        out.ResidT.push_back(source.ResidT[i]);
        out.ArvR.push_back(source.ArvR[i]);
        out.Tput.push_back(source.Tput[i]);
    }
    return out;
}

}  // namespace

AvgTable AvgTable::filter_by(const std::string& name) const {
    return filter_avg_table(*this, [&](std::size_t i) {
        return Station[i] == name || JobClass[i] == name;
    });
}

AvgTable AvgTable::filter_by(const std::string& station, const std::string& jobclass) const {
    return filter_avg_table(*this, [&](std::size_t i) {
        return Station[i] == station && JobClass[i] == jobclass;
    });
}

AvgTable AvgTable::filter_by(const Node& node) const {
    return filter_avg_table(*this, [&](std::size_t i) {
        return Station[i] == node.get_name();
    });
}

AvgTable AvgTable::filter_by(const ::line::JobClass& jobclass) const {
    return filter_avg_table(*this, [&](std::size_t i) {
        return JobClass[i] == jobclass.get_name();
    });
}

AvgTable AvgTable::filter_by(const Node& node, const ::line::JobClass& jobclass) const {
    return filter_by(node.get_name(), jobclass.get_name());
}

AvgTable AvgTable::filter_by(const ::line::JobClass& jobclass, const Node& node) const {
    return filter_by(node, jobclass);
}

AvgTable AvgTable::get(const Node& node) const { return filter_by(node); }
AvgTable AvgTable::get(const ::line::JobClass& jobclass) const { return filter_by(jobclass); }
AvgTable AvgTable::get(const Node& node, const ::line::JobClass& jobclass) const {
    return filter_by(node, jobclass);
}
AvgTable AvgTable::get(const ::line::JobClass& jobclass, const Node& node) const {
    return filter_by(node, jobclass);
}
AvgTable AvgTable::tget(const Node& node) const { return filter_by(node); }
AvgTable AvgTable::tget(const ::line::JobClass& jobclass) const { return filter_by(jobclass); }
AvgTable AvgTable::tget(const Node& node, const ::line::JobClass& jobclass) const {
    return filter_by(node, jobclass);
}
AvgTable AvgTable::tget(const ::line::JobClass& jobclass, const Node& node) const {
    return filter_by(node, jobclass);
}
AvgTable AvgTable::operator()(const Node& node) const { return filter_by(node); }
AvgTable AvgTable::operator()(const ::line::JobClass& jobclass) const {
    return filter_by(jobclass);
}
AvgTable AvgTable::operator()(const Node& node, const ::line::JobClass& jobclass) const {
    return filter_by(node, jobclass);
}
AvgTable AvgTable::operator()(const ::line::JobClass& jobclass, const Node& node) const {
    return filter_by(node, jobclass);
}

void AvgTable::print(std::ostream& out) const {
    const std::ios::fmtflags flags = out.flags();
    const std::streamsize precision = out.precision();
    out << std::left << std::setw(16) << "Station" << std::setw(14) << "JobClass"
        << std::right << std::setw(12) << "QLen" << std::setw(12) << "Util"
        << std::setw(12) << "RespT" << std::setw(12) << "ResidT"
        << std::setw(12) << "ArvR" << std::setw(12) << "Tput" << '\n';
    out << std::setprecision(6);
    for (std::size_t i = 0; i < Station.size(); ++i)
        out << std::left << std::setw(16) << Station[i] << std::setw(14) << JobClass[i]
            << std::right << std::setw(12) << QLen[i] << std::setw(12) << Util[i]
            << std::setw(12) << RespT[i] << std::setw(12) << ResidT[i]
            << std::setw(12) << ArvR[i] << std::setw(12) << Tput[i] << '\n';
    out.flags(flags);
    out.precision(precision);
}

void AvgTable::print() const { print(std::cout); }

std::ostream& operator<<(std::ostream& out, const AvgTable& table) {
    table.print(out);
    return out;
}

namespace {

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
 * agree only where every station is visited once per cycle. ArvR is a FLOW and
 * parts from the throughput at any station a job leaves by another route.
 * `sn_get_residt_from_respt` and `sn_get_arvr_from_tput` are the reference's own
 * conversions and are what the CLI's `-s fluid` and `-s ssa` arms apply. A
 * Source's ArvR stays zero: nothing arrives TO it.
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
 * not 1.
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

/**
 * The caller-facing `map_env` knobs, as the gate wants them.
 *
 * They live on `SolverOptions` and NOT on `MvaOptions`/`NcSolverOptions`/
 * `FluidOptions`, because no `runAnalyzer` in the reference reads them: the
 * fallback is decided one level ABOVE the runner, in `getAvg`, which is what
 * this facade's `avg_table` plays the part of here.
 */
solvers::MapEnvConfig map_env_config(const SolverOptions& o) {
    solvers::MapEnvConfig c;
    if (!o.map_env.empty()) c.mode = o.map_env;
    if (!o.map_env_method.empty()) c.method = o.map_env_method;
    if (o.map_env_maxstages) c.max_stages = o.map_env_maxstages;
    return c;
}

ctmc::CtmcOptions ctmc_options(const SolverOptions& o) {
    ctmc::CtmcOptions c;
    if (!o.method.empty()) c.method = o.method;
    if (o.cutoff > 0.0) c.cutoff = o.cutoff;
    if (!o.cutoff_vec.empty()) c.cutoff_vec = o.cutoff_vec;
    if (o.state_max) c.state_max = o.state_max;
    return c;
}

/**
 * The symbolic knobs, which live on their own struct rather than on
 * CtmcOptions: the backend and its timeout configure the computer-algebra
 * request, not the generator build. An empty or non-positive field keeps the
 * engine default, as everywhere else here.
 */
ctmc::CtmcSymbolicOptions ctmc_symbolic_options(const SolverOptions& o) {
    ctmc::CtmcSymbolicOptions s;
    if (!o.symbolic.empty()) s.backend = o.symbolic;
    if (o.symbolic_timeout > 0) s.timeout_s = o.symbolic_timeout;
    return s;
}

fluid::FluidOptions fluid_options(const SolverOptions& o) {
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

/**
 * The JMT client's options.
 *
 * `timespan_end` becomes `max_simulated_time`, which is the finite horizon the
 * replicated arms REQUIRE, and `replications` is `options.config.replications`, the
 * size of their ensemble. `iter_max` means nothing to JMT and is not read.
 */
jmt::JmtOptions jmt_options(const SolverOptions& o) {
    jmt::JmtOptions j;
    if (!o.method.empty()) j.method = o.method;
    if (o.samples) j.samples = static_cast<double>(o.samples);
    if (o.seed) j.seed = static_cast<long>(o.seed);
    if (o.replications > 0) j.replications = o.replications;
    if (o.timespan_end > 0.0) j.max_simulated_time = o.timespan_end;
    j.keep = o.keep;
    j.verbose = o.verbose;
    return j;
}

/** `AUTO` resolves to the solver its feature set picks, as the CLI's auto arm does. */
std::string resolve_auto(const std::string& name, const qn::NetworkStruct<double>& sn) {
    if (name != "AUTO") return name;
    std::string picked = autosolver::auto_solver_name(autosolver::auto_choose_avg_solver(sn));
    std::transform(picked.begin(), picked.end(), picked.begin(), ::toupper);
    if (picked == "FLUID") picked = "FLD";
    return picked;
}

}  // namespace

// ---------------------------------------------------------------------------
// NetworkSolver
// ---------------------------------------------------------------------------

const AvgTable& NetworkSolver::avg_table() {
    if (solved_) return table_;
    AvgTable t;
    t.solver = name_;
    const qn::NetworkStruct<double>& sn = model_->get_struct();
    const SolverOptions& o = opts_;

    if (name_ == "LQNS")
        throw UnsupportedError("avg_table: '" + name_ +
                               "' is not wrapped by this facade; its Network path is "
                               "lqns::solve_network_run_analyzer");

    const std::string name = resolve_auto(name_, sn);

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
        const mva::AvgResult<double> r = solvers::run_avg<double>(
            sn, "SolverMVA", qn::mva_feature_set(mva::resolve_method(sn, opt.method), sn),
            map_env_config(o),
            [&opt, &init](const qn::NetworkStruct<double>& m) {
                return mva::solver_mva_run_analyzer(m, opt, init);
            },
            solvers::mva_stage_fn(opt), opt.method);
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
        const mva::AvgResult<double> r = solvers::run_avg<double>(
            sn, "SolverNC", qn::nc_feature_set(opt.method), map_env_config(o),
            [&opt](const qn::NetworkStruct<double>& m) { return nc::solver_nc_run_analyzer(m, opt); },
            solvers::nc_stage_fn(opt), opt.method);
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
        // DELIBERATELY NOT WRAPPED IN `run_avg`. A bound request must be answered
        // with a bound, and the environment image is an APPROXIMATION of the
        // model, so bounds computed on it do not bracket the original. The
        // reference refuses SolverBA inside `needsMapEnv` itself; here the
        // exclusion is the absence of the wrapper, so do not "fix" the
        // inconsistency with the arms above.
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
        // The fluid arm cannot go through `run_avg`: its runner reports a
        // `FluidSolution` and the driver an `AvgResult`, so the two halves fill
        // the table through different helpers. The DECISION is the same one
        // `run_avg` makes, asked here directly rather than duplicated.
        const solvers::MapEnvConfig mecfg = map_env_config(o);
        const solvers::MapEnvDecision d =
            solvers::needs_map_env(qn::fluid_feature_set(opt.method), sn, mecfg);
        if (d.needed) {
            const mva::AvgResult<double> r = solvers::map_env_approx<double>(
                sn, "SolverFLD", mecfg, solvers::fluid_stage_fn(opt), opt.method);
            t.method = r.actualmethod;
            fill_table(t, sn, r);
        } else {
            const fluid::FluidSolution r = fluid::solver_fluid_run_analyzer(sn, opt);
            t.method = r.method;
            t.iter = static_cast<int>(r.iters);
            fill_table_sim(t, sn, r);
        }
    } else if (name == "JMT") {
        const jmt::JmtOptions opt = jmt_options(o);
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
        throw UnsupportedError("avg_table: unknown solver '" + name_ + "'");
    }
    table_ = t;
    solved_ = true;
    return table_;
}

std::string NetworkSolver::method_used() { return avg_table().method; }

std::vector<std::string> NetworkSolver::list_valid_methods() const {
    const qn::NetworkStruct<double>& sn = model_->get_struct();
    if (name_ == "MVA") return mva::list_valid_methods(sn);
    if (name_ == "NC") return nc::list_valid_methods();
    if (name_ == "CTMC") return ctmc::list_valid_methods();
    if (name_ == "MAM") return mam::list_valid_methods();
    // The model-aware overload, not the bare one: the reduction bounds are
    // derived for a single-class closed network of single servers and the three
    // open-network bounds for its mirror image, so the list a caller may act on
    // depends on the model. `SolverBA` gates on exactly this list.
    if (name_ == "BA") return ba::list_valid_methods(sn);
    if (name_ == "FLD") return fluid::fluid_list_valid_methods();
    if (name_ == "SSA") return ssa::list_valid_methods();
    if (name_ == "JMT") return jmt::jmt_list_valid_methods();
    if (name_ == "LDES") return ldes::list_valid_methods();
    // AUTO answers about every family it can delegate to, gated by each one's
    // feature set on this model; see `auto_methods.h`.
    if (name_ == "AUTO") return autosolver::auto_list_valid_methods(sn);
    throw UnsupportedError("list_valid_methods: unknown solver '" + name_ + "'");
}

// ---------------------------------------------------------------------------
// SolverCTMC
// ---------------------------------------------------------------------------

Matrix<double> SolverCTMC::state_space() {
    const qn::NetworkStruct<double>& sn = model_->get_struct();
    const ctmc::CtmcSolution<double> d = ctmc::solver_ctmc_analyzer(sn, ctmc_options(opts_));
    return ctmc::ctmc_get_state_space(sn, d).flat;
}

Matrix<double> SolverCTMC::generator() {
    const qn::NetworkStruct<double>& sn = model_->get_struct();
    ctmc::CtmcOptions opt = ctmc_options(opts_);
    opt.keep_filtration = true;
    const ctmc::CtmcSolution<double> d = ctmc::solver_ctmc_analyzer(sn, opt);
    return ctmc::ctmc_get_infgen(sn, d).Q;
}

double SolverCTMC::prob_aggr(std::size_t node, const std::vector<double>& state) {
    const qn::NetworkStruct<double>& sn = model_->get_struct();
    const ctmc::CtmcSolution<double> d = ctmc::solver_ctmc_analyzer(sn, ctmc_options(opts_));
    const Matrix<double> A = ctmc::ctmc_get_state_space_aggr(sn, d);
    const std::size_t ist = sn.nodes[node - 1].station;
    if (ist == 0) throw InputError("prob_aggr: the node is not a station");
    const std::size_t K = sn.nclasses;
    if (state.size() != K) throw InputError("prob_aggr: the state must have one entry per class");
    double p = 0.0;
    for (std::size_t s = 0; s < A.rows(); ++s) {
        bool hit = true;
        for (std::size_t k = 0; k < K && hit; ++k) hit = A(s, (ist - 1) * K + k) == state[k];
        if (hit) p += d.pi[s];
    }
    return p;
}

std::vector<double> SolverCTMC::marg_aggr(std::size_t node) {
    const qn::NetworkStruct<double>& sn = model_->get_struct();
    const ctmc::CtmcSolution<double> d = ctmc::solver_ctmc_analyzer(sn, ctmc_options(opts_));
    const Matrix<double> A = ctmc::ctmc_get_state_space_aggr(sn, d);
    const std::size_t ist = sn.nodes[node - 1].station;
    if (ist == 0) throw InputError("marg_aggr: the node is not a station");
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

std::vector<std::vector<CdfCurve> > SolverCTMC::cdf_respt() {
    const qn::NetworkStruct<double>& sn = model_->get_struct();
    std::vector<std::vector<CdfCurve> > out(sn.nstations, std::vector<CdfCurve>(sn.nclasses));
    const std::vector<std::vector<ctmc::CdfCurve<double> > > R =
        ctmc::solver_ctmc_cdf_respt(sn, ctmc_options(opts_));
    for (std::size_t i = 0; i < R.size() && i < out.size(); ++i)
        for (std::size_t c = 0; c < R[i].size() && c < out[i].size(); ++c) {
            out[i][c].t = R[i][c].t;
            out[i][c].F = R[i][c].F;
        }
    return out;
}

SymbolicSolution SolverCTMC::symbolic_solution() {
    const qn::NetworkStruct<double>& sn = model_->get_struct();
    const ctmc::CtmcSymbolicSolution<double> r =
        ctmc::ctmc_symbolic_solution(sn, ctmc_options(opts_), ctmc_symbolic_options(opts_));
    SymbolicSolution out;
    out.pi = r.pi;
    out.num = r.num;
    out.den = r.den;
    out.engine = r.engine;
    out.symbols = r.symbols;
    out.rate0 = r.rate0;
    return out;
}

void SolverCTMC::print_inf_gen(const Matrix<double>& Q, const Matrix<double>& space) {
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

// ---------------------------------------------------------------------------
// SolverFLD
// ---------------------------------------------------------------------------

TranAvg SolverFLD::tran_avg() {
    const qn::NetworkStruct<double>& sn = model_->get_struct();
    const std::vector<fluid::FluidTranPoint> pts =
        fluid::solver_fluid_tran_avg(sn, fluid_options(opts_), 100);
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

std::vector<std::vector<CdfCurve> > SolverFLD::cdf_respt() {
    const qn::NetworkStruct<double>& sn = model_->get_struct();
    std::vector<std::vector<CdfCurve> > out(sn.nstations, std::vector<CdfCurve>(sn.nclasses));
    const std::vector<std::vector<fluid::FluidPassage> > R =
        fluid::solver_fluid_cdf_respt(sn, fluid_options(opts_));
    for (std::size_t i = 0; i < R.size() && i < out.size(); ++i)
        for (std::size_t c = 0; c < R[i].size() && c < out[i].size(); ++c) {
            out[i][c].t = R[i][c].t;
            out[i][c].F = R[i][c].cdf;
        }
    return out;
}

// ---------------------------------------------------------------------------
// SolverJMT
// ---------------------------------------------------------------------------

std::vector<std::vector<CdfCurve> > SolverJMT::cdf_respt() {
    const qn::NetworkStruct<double>& sn = model_->get_struct();
    std::vector<std::vector<CdfCurve> > out(sn.nstations, std::vector<CdfCurve>(sn.nclasses));
    // `true` is the getCdfRespT contract: the model is first solved for its
    // steady-state queue lengths and the logged run starts preloaded at their
    // rounded values, which is what makes the seeded curve comparable to the
    // reference's for the same seed.
    const std::map<std::pair<std::size_t, std::size_t>, std::vector<std::pair<double, double> > >
        rd = jmt::jmt_get_cdf_resp_t(sn, jmt_options(opts_), true);
    for (std::map<std::pair<std::size_t, std::size_t>,
                  std::vector<std::pair<double, double> > >::const_iterator it = rd.begin();
         it != rd.end(); ++it) {
        const std::size_t i = it->first.first, c = it->first.second;
        if (i == 0 || i > out.size() || c == 0 || c > sn.nclasses) continue;
        // THE PAIRS ARE (F, X), as `ecdf` returns them, and the curve stores
        // the two apart: reading them the other way round would report a
        // probability as a time.
        for (std::size_t j = 0; j < it->second.size(); ++j) {
            out[i - 1][c - 1].F.push_back(it->second[j].first);
            out[i - 1][c - 1].t.push_back(it->second[j].second);
        }
    }
    return out;
}

TranAvg SolverJMT::tran_avg() {
    const qn::NetworkStruct<double>& sn = model_->get_struct();
    const jmt::JmtReplication<double> r = jmt::jmt_transient_replications(sn, jmt_options(opts_));
    TranAvg out;
    for (std::size_t i = 0; i < sn.nstations; ++i)
        for (std::size_t c = 0; c < sn.nclasses; ++c)
            out.label.push_back(sn.stations[i].name + "/" + sn.classes[c].name);
    out.t = r.t;
    for (std::size_t g = 0; g < r.t.size(); ++g) {
        std::vector<double> row;
        for (std::size_t i = 0; i < sn.nstations; ++i)
            for (std::size_t c = 0; c < sn.nclasses; ++c)
                row.push_back(i < r.QNt.size() && c < r.QNt[i].size() && g < r.QNt[i][c].size()
                                  ? r.QNt[i][c][g]
                                  : 0.0);
        out.QNt.push_back(row);
    }
    return out;
}

double SolverJMT::prob_aggr(std::size_t node, const std::vector<double>& state) {
    const qn::NetworkStruct<double>& sn = model_->get_struct();
    if (node == 0 || node > sn.nodes.size()) throw InputError("prob_aggr: node index out of range");
    const std::size_t ist = sn.nodes[node - 1].station;
    if (ist == 0) throw InputError("prob_aggr: the node is not a station");
    const jmt::JmtProbAggr r = jmt::jmt_prob_aggr(sn, jmt_options(opts_), ist, state);
    return r.station[ist - 1];
}

SamplePath SolverJMT::sample(std::size_t events, std::size_t node) {
    const qn::NetworkStruct<double>& sn = model_->get_struct();
    const jmt::JmtOptions opt = jmt_options(opts_);
    SamplePath out;
    if (node) {
        if (node > sn.nodes.size() || sn.nodes[node - 1].station == 0)
            throw InputError("sample: the node is not a station, so nothing is logged for it");
        const jmt::JmtNodeTrace<double> tr = jmt::jmt_sample_aggr(sn, node, events, opt);
        for (std::size_t c = 0; c < sn.nclasses; ++c) out.label.push_back(sn.classes[c].name);
        out.t = tr.t;
        out.state = tr.qlen;
        return out;
    }
    const jmt::JmtSysTrace<double> tr = jmt::jmt_sample_sys_aggr(sn, events, opt);
    for (std::size_t i = 0; i < sn.nstations; ++i)
        for (std::size_t c = 0; c < sn.nclasses; ++c)
            out.label.push_back(sn.stations[i].name + "/" + sn.classes[c].name);
    out.t = tr.t;
    for (std::size_t g = 0; g < tr.t.size(); ++g) {
        std::vector<double> row;
        for (std::size_t i = 0; i < sn.nstations; ++i)
            for (std::size_t c = 0; c < sn.nclasses; ++c)
                row.push_back(i < tr.state.size() && g < tr.state[i].size() &&
                                      c < tr.state[i][g].size()
                                  ? tr.state[i][g][c]
                                  : 0.0);
        out.state.push_back(row);
    }
    return out;
}

// ---------------------------------------------------------------------------
// SolverBA
// ---------------------------------------------------------------------------

BoundsTable SolverBA::bounds_table() {
    const qn::NetworkStruct<double>& sn = model_->get_struct();
    ba::BaOptions opt;
    if (!opts_.method.empty()) opt.method = opts_.method;
    opt.level = opts_.level;
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

}  // namespace line
