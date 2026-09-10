/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_EXAMPLES_EXAMPLE_UTIL_H
#define LINE_EXAMPLES_EXAMPLE_UTIL_H

/**
 * What a ported example script is allowed to call.
 *
 * A MATLAB example is two things: a model built with the `Network` API, and a
 * few solver calls whose tables it prints. The first is already here --
 * `qn::Network<double>` IS that API -- and the second is what this header adds:
 * one `solve_avg` over every native solver, plus the non-average analyses the
 * examples reach for (`getProb`, `getCdfRespT`, `getTranAvg`, the CTMC state
 * space and generator, the bounds table).
 *
 * IT IS A FACADE, NOT A SOLVER. Every entry point forwards to the same runner
 * the CLI calls, so an example and `line-cli` on the same model answer with the
 * same numbers, and nothing is computed here. The solver headers are included
 * by example_util.cpp ALONE: they instantiate the whole template stack, and
 * paying that in each of the ~30 example translation units would dominate the
 * build.
 *
 * WRAPPER SOLVERS. JMT and LDES go through the same `solve_avg` facade as the
 * native solvers, forwarding to `jmt::solver_jmt_run_analyzer` and `ldes::solver_ldes`,
 * so a ported example runs the simulation its reference ran rather than
 * declaring it absent. Both need an external artefact -- `common/JMT.jar` plus
 * a JVM, and `common/ldes` or `common/ldes.jar` -- and throw by name when it is
 * missing, which is a different statement from "not ported". LQNS and QNS stay
 * out: they analyse a LayeredNetwork, not the `Net` this facade takes.
 */

#include <cstddef>
#include <string>
#include <vector>

#include "line/lang/distribution.h"
#include "line/lang/distributions.h"
#include "line/lang/qn/network_builder.h"
#include "line/lang/qn/nodes.h"
#include "line/solvers/solver.h"

namespace line {
namespace examples {

// The spellings a ported script uses, so a line of C++ reads like the line of
// MATLAB it came from.
typedef lang::Distrib<double> D;
typedef qn::Network<double> Net;
typedef qn::RoutingMatrix<double> Routing;
using lang::DropStrategy;
using lang::NodeType;
using lang::PollingType;
using lang::ReplacementStrategy;
using lang::RoutingStrategy;
using lang::SchedStrategy;
typedef qn::CacheParam<double> CacheParam;
typedef qn::TransitionParam<double> TransitionParam;

// ---------------------------------------------------------------------------
// Model-building helpers that MATLAB has as static methods
// ---------------------------------------------------------------------------

/** `Network.serialRouting(nodes)` for one class pair: 1 -> 2 -> ... -> n. */
void serial_routing(Routing& P, std::size_t r, std::size_t s,
                    const std::vector<std::size_t>& nodes);

/** `Network.serialRouting(nodes)` on the (only) class of a single-class model. */
void serial_routing(Routing& P, const std::vector<std::size_t>& nodes);

/**
 * `model.link(Network.serialRouting(nodes))` for a single-class CLOSED chain or
 * an open one: the last node is wired back to the first unless it is a Sink.
 */
void link_serial(Net& m, const std::vector<std::size_t>& nodes);

/** `HyperExp.fitMeanAndSCV(mean, scv)`, the balanced-mean two-phase fit. */
D hyperexp_fit(double mean, double scv);

/** `Erlang.fitMeanAndSCV(mean, scv)`. */
D erlang_fit(double mean, double scv);

/** `Erlang.fitMeanAndOrder(mean, k)`. */
D erlang_fit_order(double mean, std::size_t k);

/** `APH.fitMeanAndSCV(mean, scv)` -- the acyclic phase-type moment fit. */
D aph_fit(double mean, double scv);

/** `Replayer(filename)`: the trace file, one sample per line. */
D replayer_from_file(const std::string& path);

/** `Zipf(alpha, n)` as the read popularity of a Cache node: p_i ~ i^-alpha. */
std::vector<double> zipf(double alpha, std::size_t n);

/** The repository root, so an example can name a trace file the suite ships. */
std::string line_root_folder();

// ---------------------------------------------------------------------------
// The average table
// ---------------------------------------------------------------------------

/** The knobs the examples set; a negative or empty field keeps the default. */
struct SolverOpts {
    std::string method = "default";
    double tol = -1.0;
    double iter_tol = -1.0;
    int iter_max = -1;
    std::size_t samples = 0;
    unsigned long seed = 0;
    /** CTMC `options.cutoff`, scalar or per class. */
    double cutoff = -1.0;
    std::vector<std::size_t> cutoff_vec;
    /** CTMC refusal threshold on the state-space size. */
    std::size_t state_max = 0;
    /** Fluid `options.timespan(2)` and warm start. */
    double timespan_end = -1.0;
    std::vector<double> init_sol;
    bool stiff = false;
    /** `options.config.*` of the MVA / NC / fluid families. */
    std::string multiserver;
    std::string highvar;
    std::string np_priority;
    /** MVA / NC `options.config.fork_join`: 'default'/'mmt'/'fjt' or 'ht'. */
    std::string fork_join;
    /** SSA `options.config.state_space_gen`. */
    std::string state_space_gen;
    /** SolverBA `options.level`. */
    int level = 2;
    /** JMT `options.keep`: leave the scratch directory in place after the solve. */
    bool keep = false;
    bool verbose = false;
};

/** `getAvgTable`, one row per (station, class) that carries a metric. */
struct AvgTable {
    std::vector<std::string> Station, JobClass;
    std::vector<double> QLen, Util, RespT, ResidT, ArvR, Tput;
    /** `getAvgSysTable`: system response time and throughput, per class. */
    std::vector<std::string> SysClass;
    std::vector<double> SysRespT, SysTput;
    std::string solver, method;
    int iter = 0;
    bool has_lognormconst = false;
    double lognormconst = 0.0;
    /** `getAvgCacheTable`'s ListCost column; empty on a model without item sizes. */
    std::vector<double> ListCost;
    /** The reference's own warning text, empty when it did not warn. */
    std::string warning;

    /** One cell of the table, by station and class NAME; NaN when absent. */
    double get(const std::string& column, const std::string& station,
               const std::string& jobclass) const;
    /** The column of a station over every class it has a row for. */
    std::vector<double> column(const std::string& column) const;
};

/**
 * Solve and return the average table.
 *
 * `solver` is one of MVA, NC, CTMC, FLD, SSA, MAM, BA, AUTO, JMT, LDES -- the
 * names MATLAB's `SolverX` classes carry, so a ported line names the same
 * solver its source did. LQNS and QNS are refused BY NAME rather than mapped
 * onto a native solver.
 */
AvgTable solve_avg(const std::string& solver, Net& m, const SolverOpts& opt = SolverOpts());

/** Print an AvgTable the way the CLI prints it, under a caption. */
void print_avg(const AvgTable& t, const std::string& caption = std::string());

/**
 * The metric matrices of a wrapper solve, in the layout `print_avg(sn, r)` of
 * examples_common.h reads.
 *
 * It exists so that an example written against that header can call JMT or LDES
 * without including the wrapper -- and, more to the point, without instantiating
 * the whole subprocess and JSON stack in its own translation unit. The two
 * `*_avg` entry points below are the only place either wrapper is instantiated.
 */
struct WrapperAvg {
    Matrix<double> QN, UN, RN, WN, AN, TN;
    /** Read by `print_avg`; a wrapper reports none, so it stays empty. */
    std::string warning;
    std::string method;
};

/** `JMT(model, 'seed', opt.seed, 'samples', opt.samples).getAvgTable()`. */
WrapperAvg jmt_avg(Net& m, const SolverOpts& opt = SolverOpts());

/** `LDES(model, 'seed', opt.seed, 'samples', opt.samples).getAvgTable()`. */
WrapperAvg ldes_avg(Net& m, const SolverOpts& opt = SolverOpts());

/** The usual `SolverOpts` of a reference JMT/LDES call: a seed and a budget. */
SolverOpts sim_opts(unsigned long seed, std::size_t samples = 0);

/** Print the system table (`getAvgSysTable`). */
void print_sys(const AvgTable& t, const std::string& caption = std::string());

/** The banner an example prints before its blocks, `MODEL: <name>`. */
void banner(const std::string& text);

/** Print a labelled scalar, vector or matrix, so a ported `disp` has a home. */
void print_scalar(const std::string& label, double v);
void print_vector(const std::string& label, const std::vector<double>& v);
void print_matrix(const std::string& label, const Matrix<double>& M);

// ---------------------------------------------------------------------------
// The two columns a fluid or simulation analyzer does not report itself
// ---------------------------------------------------------------------------

/**
 * `sn_get_residt_from_respt`: the per-JOB residence time RN implies.
 *
 * ResidT IS NOT RespT UNLESS EVERY STATION IS VISITED ONCE PER CYCLE. The fluid
 * and simulation analyzers report a per-VISIT response time and no residence
 * time of their own, so a caller owes the reference's own conversion; reporting
 * RespT in its place was a factor of 3 out on `sdroute_closed` and 17 on Queue1
 * of `init_state_ps`, both multi-visit closed models. This is the conversion the
 * CLI's `-s fluid` and `-s ssa` arms apply, and it is a pure function of `sn`
 * and RN.
 */
Matrix<double> residt_from_respt(const qn::NetworkStruct<double>& sn,
                                 const Matrix<double>& RN);

/**
 * `sn_get_arvr_from_tput`: the arrival rate the routing and TN imply.
 *
 * THE ARRIVAL RATE IS A FLOW, NOT A COPY OF THE THROUGHPUT: the two agree only
 * where every job a station serves it also completes, and they part on a station
 * a job LEAVES by another route -- a Join, or the Delay stations of
 * `cache_replc_routing`, which take arrivals the fluid solution gives zero
 * throughput.
 */
Matrix<double> arvr_from_tput(const qn::NetworkStruct<double>& sn,
                              const Matrix<double>& TN);


// ---------------------------------------------------------------------------
// The analyses that are not the average table
// ---------------------------------------------------------------------------

/** `getProbAggr(node, state)` and `getProb(node, state)` under SolverCTMC. */
double ctmc_prob_aggr(Net& m, std::size_t node, const std::vector<double>& state,
                      const SolverOpts& opt = SolverOpts());

/** `getProbStateAggr` for the whole system: the marginal over each station. */
std::vector<double> ctmc_marg_aggr(Net& m, std::size_t node,
                                   const SolverOpts& opt = SolverOpts());

/** `getStateSpace`: the aggregate state space, one row per state. */
Matrix<double> ctmc_state_space(Net& m, const SolverOpts& opt = SolverOpts());

/** `getGenerator`: the infinitesimal generator over that space. */
Matrix<double> ctmc_generator(Net& m, const SolverOpts& opt = SolverOpts());

/** `CTMC.printInfGen(Q, SS)`: the generator beside the state it belongs to. */
void print_inf_gen(const Matrix<double>& Q, const Matrix<double>& space);

/** One response-time CDF curve: column 0 the CDF value, column 1 the time. */
struct CdfCurve {
    std::vector<double> F, t;
};

/** `getCdfRespT` under SolverCTMC or SolverFLD, per (station, class). */
std::vector<std::vector<CdfCurve> > cdf_respt(const std::string& solver, Net& m,
                                              const SolverOpts& opt = SolverOpts());

/** Print the percentiles of a CDF curve, which is what the examples show. */
void print_cdf_summary(const std::string& label, const CdfCurve& c);

/** `getTranAvg` under SolverFLD: the transient mean queue length per station. */
struct TranAvg {
    std::vector<double> t;                    ///< the time axis
    std::vector<std::vector<double> > QNt;    ///< [step][station*class]
    std::vector<std::string> label;           ///< the (station, class) of each column
};
TranAvg fluid_tran_avg(Net& m, const SolverOpts& opt = SolverOpts());

/**
 * `Network.initFromMarginal(n)` FOR A FLUID TRANSIENT: the ODE initial
 * condition that places `n[station][class]` jobs, all in phase one.
 *
 * The port carries no `sn.state`, so a marginal cannot be turned into a
 * declared network state -- but the fluid ODE never needs one: it needs the
 * FIRST MOMENT, and `State.fromMarginal` writes `init(1) = si(r)`, i.e. every
 * placed job in phase one. That is exactly what this writes into
 * `FluidOptions.init_sol`, so a transient started here is the reference's
 * transient started from the same point. Priors that spread jobs ACROSS phases
 * (`setStatePrior`) are still out of reach and stay refused by name.
 *
 * @param n (nstations x nclasses) job counts; a short row is zero-filled
 */
std::vector<double> fluid_initsol_from_marginal(Net& m,
                                                const std::vector<std::vector<double> >& n);

/**
 * `StatefulNode.setStatePrior(uniform)`: the same placement, with the PHASE
 * assignment at `prior_station` (1-based) averaged uniformly over every state
 * carrying that station's marginal.
 *
 * `from_marginal` enumerates those states and `to_marginal` decodes each back
 * to `kir`, so the first moment the ODE needs is their mean -- which is all a
 * fluid transient can see of a prior. Only the phase split moves, so this and
 * `fluid_initsol_from_marginal` differ by very little (1.5e-4 of queue length on
 * `init_state_fcfs_nonexp`), and that difference is exactly what the reference's
 * Prior 3 block reports.
 */
std::vector<double> fluid_initsol_from_state_prior(Net& m,
                                                   const std::vector<std::vector<double> >& n,
                                                   std::size_t prior_station);

/** `SolverBA(model, method).getBoundsTable()`. */
struct BoundsTable {
    std::vector<std::string> Station, JobClass;
    std::vector<double> Qlower, Qupper, Tlower, Tupper;
    std::string method;
};
BoundsTable ba_bounds(Net& m, const SolverOpts& opt = SolverOpts());
void print_bounds(const BoundsTable& b, const std::string& caption = std::string());

/** The method names a solver advertises on this model (`listValidMethods`). */
std::vector<std::string> list_valid_methods(const std::string& solver, Net& m);

}  // namespace examples
}  // namespace line

#endif  // LINE_EXAMPLES_EXAMPLE_UTIL_H
