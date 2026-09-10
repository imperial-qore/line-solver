/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_SOLVER_H
#define LINE_SOLVERS_SOLVER_H

/**
 * The solver API a user writes, spelled as its Python twin.
 *
 *     SolverMVA  s(model, "lin");
 *     AvgTable   t = s.avg_table();
 *
 * against Python's
 *
 *     s = MVA(model, method='lin')
 *     t = s.avg_table()
 *
 * NAMING. A multi-word getter drops the `get_` prefix (`avg_table`,
 * `tran_avg`, `cdf_respt`), because `line_solver/_aliasing.py` resolves exactly
 * that spelling onto Python's own `getAvgTable`/`getTranAvg`/`getCdfRespT`; a
 * SINGLE-word getter keeps it (`get_name`), because that same module
 * deliberately refuses to expand a bare lowercase word and Python has no
 * `name()`. So every name here is a name that also resolves in Python.
 *
 * IT IS A FACADE, NOT A SOLVER. Every method forwards to the runner the CLI
 * calls, so a program written against this header and `line-cli` on the same
 * model answer with the same numbers.
 *
 * THESE CLASSES ARE `double`-ONLY AND NOT TEMPLATES, deliberately. The solver
 * templates are a heavy instantiation and the multiprecision arithmetic has no
 * Python counterpart, so the bodies are compiled ONCE into `line_mp_api` and a
 * translation unit including this header pays for none of it. Reach for
 * `line::mva::solver_mva_run_analyzer` and its siblings directly when you want
 * another arithmetic.
 */

#include <cstddef>
#include <string>
#include <vector>

#include "line/lang/qn/network_builder.h"
#include "line/solvers/avg_table.h"
#include "line/solvers/solver_options.h"
#include "line/util/matrix.h"

namespace line {

/** `Network` as a user names it: the `double` model, Python's `Network`. */
typedef qn::Network<double> Network;
/** `model.init_routing_matrix()`'s type, Python's `RoutingMatrix`. */
typedef qn::RoutingMatrix<double> RoutingMatrix;

/**
 * The shared surface of every solver, Python's `NetworkSolver`.
 *
 * The solve is performed on first demand and cached, which is what
 * `hasResults`/`is_solved` reports; `reset()` discards it.
 */
class NetworkSolver {
  public:
    virtual ~NetworkSolver() {}

    /** `getName()`: the solver's own name, e.g. "MVA". */
    const std::string& get_name() const { return name_; }
    /** The options this solver was built with. */
    const SolverOptions& options() const { return opts_; }
    /** The model this solver was built on. */
    Network& model() const { return *model_; }

    /** `getAvgTable()`: the average table, solving on first demand. */
    const AvgTable& avg_table();
    /** `getAvgSysTable()`: the per-class system columns of the same solve. */
    const AvgTable& avg_sys_table() { return avg_table(); }
    /** `runAnalyzer()`: force the solve, returning the table it produced. */
    const AvgTable& run_analyzer() { return avg_table(); }
    /** `listValidMethods()`: the methods this solver advertises on this model. */
    std::vector<std::string> list_valid_methods() const;
    /** `isSolved()` / `hasResults()`: whether a solve has been performed. */
    bool is_solved() const { return solved_; }
    /** `reset()`: discard the cached solve. */
    void reset() { solved_ = false; table_ = AvgTable(); }
    /** `getMethodUsed()`: the method the solve actually resolved to. */
    std::string method_used();

  protected:
    NetworkSolver(Network& m, const std::string& name, const SolverOptions& o)
        : model_(&m), name_(name), opts_(o) {}

    Network* model_;
    std::string name_;
    SolverOptions opts_;
    AvgTable table_;
    bool solved_ = false;
};

/**
 * Every solver takes the model and either a method name or a full option set,
 * mirroring Python's `SolverX(model, method_or_options=None, **kwargs)`.
 */
#define LINE_DECLARE_SOLVER(Cls, tag)                                                  \
    class Cls : public NetworkSolver {                                                 \
      public:                                                                          \
        explicit Cls(Network& m, const SolverOptions& o = SolverOptions())              \
            : NetworkSolver(m, tag, o) {}                                              \
        Cls(Network& m, const std::string& method)                                     \
            : NetworkSolver(m, tag, SolverOptions().set_method(method)) {}             \
    }

LINE_DECLARE_SOLVER(SolverMVA, "MVA");
LINE_DECLARE_SOLVER(SolverNC, "NC");
LINE_DECLARE_SOLVER(SolverMAM, "MAM");
LINE_DECLARE_SOLVER(SolverSSA, "SSA");
LINE_DECLARE_SOLVER(SolverJMT, "JMT");
LINE_DECLARE_SOLVER(SolverLDES, "LDES");
LINE_DECLARE_SOLVER(SolverAUTO, "AUTO");

#undef LINE_DECLARE_SOLVER

/** `SolverCTMC`: the average table plus the chain it was computed from. */
class SolverCTMC : public NetworkSolver {
  public:
    explicit SolverCTMC(Network& m, const SolverOptions& o = SolverOptions())
        : NetworkSolver(m, "CTMC", o) {}
    SolverCTMC(Network& m, const std::string& method)
        : NetworkSolver(m, "CTMC", SolverOptions().set_method(method)) {}

    /** `getStateSpace()`: the aggregate state space, one row per state. */
    Matrix<double> state_space();
    /** `getGenerator()`: the infinitesimal generator over that space. */
    Matrix<double> generator();
    /** `getProbAggr(node, state)`: the aggregate marginal of one state. */
    double prob_aggr(std::size_t node, const std::vector<double>& state);
    /** `getProbStateAggr(node)`: the marginal over every state of one station. */
    std::vector<double> marg_aggr(std::size_t node);
    /** `getCdfRespT()`: the response-time CDF per (station, class). */
    std::vector<std::vector<CdfCurve> > cdf_respt();

    /** `CTMC.printInfGen(Q, SS)`: the generator beside the state it belongs to. */
    static void print_inf_gen(const Matrix<double>& Q, const Matrix<double>& space);
};

/** `SolverFLD`: the fluid solver, which also answers a transient. */
class SolverFLD : public NetworkSolver {
  public:
    explicit SolverFLD(Network& m, const SolverOptions& o = SolverOptions())
        : NetworkSolver(m, "FLD", o) {}
    SolverFLD(Network& m, const std::string& method)
        : NetworkSolver(m, "FLD", SolverOptions().set_method(method)) {}

    /** `getTranAvg()`: the transient mean queue length per station. */
    TranAvg tran_avg();
    /** `getCdfRespT()`: the response-time CDF per (station, class). */
    std::vector<std::vector<CdfCurve> > cdf_respt();
};

/** `SolverBA`: the bounding solver, whose result is a bounds table. */
class SolverBA : public NetworkSolver {
  public:
    explicit SolverBA(Network& m, const SolverOptions& o = SolverOptions())
        : NetworkSolver(m, "BA", o) {}
    SolverBA(Network& m, const std::string& method)
        : NetworkSolver(m, "BA", SolverOptions().set_method(method)) {}

    /** `getBoundsTable()`: the per-class queue-length and throughput bounds. */
    BoundsTable bounds_table();
};

// Python's short aliases: `MVA(model)` is `SolverMVA(model)`.
typedef SolverMVA MVA;
typedef SolverNC NC;
typedef SolverCTMC CTMC;
typedef SolverSSA SSA;
typedef SolverFLD FLD;
typedef SolverFLD SolverFluid;
typedef SolverMAM MAM;
typedef SolverJMT JMT;
typedef SolverLDES LDES;
typedef SolverAUTO AUTO;
typedef SolverBA BA;

}  // namespace line

#endif  // LINE_SOLVERS_SOLVER_H
