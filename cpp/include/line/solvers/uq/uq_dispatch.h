/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_UQ_UQ_DISPATCH_H
#define LINE_SOLVERS_UQ_UQ_DISPATCH_H

/**
 * The stage solver of SolverUQ, named rather than passed.
 *
 * `UQ(model, @@SolverMVA)` gives the reference a factory; a caller who has a
 * solver NAME -- the CLI, a host bridge -- needs that name turned into one, and
 * this is where the turning happens. It is separate from `solver_uq.h` for the
 * reason `env_dispatch.h` is separate from `solver_env.h`: the UQ machinery
 * itself depends on nothing but the `AvgResult` contract, while this file pulls
 * in every Network solver in the port, and a caller that already has a functor
 * should not pay for that.
 *
 * WHAT EACH TOKEN COSTS, since UQ multiplies it by the design size: `mva`, `nc`
 * and `ba` are closed forms or short iterations; `mam`, `ctmc`, `ssa` and
 * `fluid` are not, and a 121-point design over two continuous Priors under
 * `-s ctmc` enumerates the state space 121 times. That is the intended
 * behaviour -- each design point IS a different model -- and it is why the
 * design cap exists.
 *
 * THE ARITHMETIC RESTRICTIONS ARE THE INNER SOLVER'S, and they are enforced
 * here at COMPILE time by `if constexpr` plus a run-time refusal: `mam` fits
 * phase-type representations, `ssa` draws exponential clocks and `fluid`
 * integrates with LSODA, so all three are double-only and an exact or
 * high-precision instantiation of them would fail to compile rather than refuse.
 *
 * ARVR AND RESIDT ON THE SIMULATION AND FLUID PATHS are filled exactly as
 * `line_cli.cpp` fills them for `-s ssa` and `-s fluid` -- residence time equals
 * response time, arrival rate equals throughput except at a Source, which has no
 * arrivals to itself. Those two solution types carry no separate AN/WN matrix,
 * and inventing a different convention here would make the UQ expectation of a
 * column disagree with the same column printed by the solver alone.
 */

#include <cstddef>
#include <string>
#include <type_traits>
#include <vector>

#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/solvers/ba/solver_ba_runner.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/ctmc/solver_ctmc_waitq.h"
#include "line/solvers/fluid/fluid_runner.h"
#include "line/solvers/fluid/solver_fluid.h"
#include "line/solvers/mam/solver_mam_runner.h"
#include "line/solvers/mva/solver_mva_runner.h"
#include "line/solvers/nc/solver_nc_runner.h"
#include "line/solvers/ssa/ssa_dispatch.h"
#include "line/solvers/uq/solver_uq.h"
#include "line/util/error.h"

namespace line {
namespace uq {

/** The inner solver's knobs, carried through untranslated. */
struct UqStageOptions {
    /** `mva` | `nc` | `mam` | `ba` | `ctmc` | `fluid` | `ssa`. */
    std::string solver;
    /** The method WITHIN that solver; empty or `default` leaves its own default. */
    std::string method;
    double tol = -1.0;       ///< < 0 = not given
    double iter_tol = -1.0;
    int iter_max = -1;
    std::size_t samples = 0;  ///< ssa run length; 0 = not given
    unsigned long seed = 0;   ///< ssa stream; 0 = not given
    double cutoff = -1.0;     ///< ctmc open-population cutoff; < 0 = not given
};

/** The solver method names `uq_stage_solver` accepts, for a caller that lists them. */
inline std::vector<std::string> uq_list_stage_solvers() {
    return std::vector<std::string>{"mva", "nc", "mam", "ba", "ctmc", "fluid", "ssa"};
}

namespace detail {

/** The double-matrix solutions of the fluid and simulation paths, as an AvgResult. */
template <class T, class Sol>
mva::AvgResult<T> from_double_solution(const qn::NetworkStruct<T>& sn, const Sol& r) {
    mva::AvgResult<T> a;
    const std::size_t M = r.QN.rows(), K = r.QN.cols();
    const T zero = num_traits<T>::from_int(0);
    a.QN = Matrix<T>(M, K, zero);
    a.UN = Matrix<T>(M, K, zero);
    a.RN = Matrix<T>(M, K, zero);
    a.TN = Matrix<T>(M, K, zero);
    a.AN = Matrix<T>(M, K, zero);
    a.WN = Matrix<T>(M, K, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t c = 0; c < K; ++c) {
            a.QN(i, c) = num_traits<T>::from_double(r.QN(i, c));
            a.UN(i, c) = num_traits<T>::from_double(r.UN(i, c));
            a.RN(i, c) = num_traits<T>::from_double(r.RN(i, c));
            a.TN(i, c) = num_traits<T>::from_double(r.TN(i, c));
            // Residence time IS the response time and the arrival rate IS the
            // throughput on these paths, except at a Source, which does not
            // arrive at itself; see the header note.
            a.WN(i, c) = a.RN(i, c);
            a.AN(i, c) = sn.stations[i].sched == lang::SchedStrategy::EXT ? zero : a.TN(i, c);
        }
    for (double v : r.CN) a.CN.push_back(num_traits<T>::from_double(v));
    for (double v : r.XN) a.XN.push_back(num_traits<T>::from_double(v));
    a.method = r.method;
    a.actualmethod = r.method;
    return a;
}

}  // namespace detail

/**
 * The stage solver named by `o.solver`.
 *
 * An unknown or unported name is refused BY NAME rather than answered with a
 * default engine: which solver ran is a property of every number UQ reports.
 */
template <class T>
UqStageSolver<T> uq_stage_solver(const UqStageOptions& o) {
    const std::string s = o.solver;
    if (s == "mva") {
        return [o](const qn::NetworkStruct<T>& sn) {
            mva::MvaOptions opt;
            if (!o.method.empty() && o.method != "default") opt.method = o.method;
            if (o.tol >= 0.0) opt.tol = o.tol;
            if (o.iter_tol >= 0.0) opt.iter_tol = o.iter_tol;
            if (o.iter_max >= 0) opt.iter_max = o.iter_max;
            Matrix<T> init;
            return mva::solver_mva_run_analyzer(sn, opt, init);
        };
    }
    if (s == "nc") {
        return [o](const qn::NetworkStruct<T>& sn) {
            nc::NcSolverOptions opt;
            if (!o.method.empty() && o.method != "default") opt.method = o.method;
            if (o.tol >= 0.0) opt.tol = o.tol;
            if (o.iter_tol >= 0.0) opt.iter_tol = o.iter_tol;
            if (o.iter_max >= 0) opt.iter_max = o.iter_max;
            return nc::solver_nc_run_analyzer(sn, opt);
        };
    }
    if (s == "ba") {
        return [o](const qn::NetworkStruct<T>& sn) {
            ba::BaOptions opt;
            if (!o.method.empty() && o.method != "default") opt.method = o.method;
            return ba::solver_ba_run_analyzer(sn, opt);
        };
    }
    if (s == "ctmc") {
        return [o](const qn::NetworkStruct<T>& sn) {
            ctmc::CtmcOptions opt;
            if (!o.method.empty() && o.method != "default") opt.method = o.method;
            if (o.cutoff >= 0.0) opt.cutoff = o.cutoff;
            return ctmc::solver_ctmc_run_analyzer_any(sn, opt);
        };
    }
    if (s == "mam") {
        if constexpr (std::is_same<T, double>::value) {
            return [o](const qn::NetworkStruct<T>& sn) {
                mam::MamOptions opt;
                if (!o.method.empty() && o.method != "default") opt.method = o.method;
                if (o.tol >= 0.0) opt.tol = o.tol;
                if (o.iter_max >= 0) opt.iter_max = o.iter_max;
                return mam::solver_mam_run_analyzer(sn, opt);
            };
        } else {
            throw UnsupportedError(
                "SolverUQ: the MAM stage solver fits phase-type representations, whose fitter "
                "requires transcendental arithmetic; run UQ under the double backend");
        }
    }
    if (s == "fluid" || s == "fld") {
        if constexpr (std::is_same<T, double>::value) {
            return [o](const qn::NetworkStruct<T>& sn) {
                fluid::FluidOptions opt;
                if (!o.method.empty()) opt.method = o.method;
                if (o.tol >= 0.0) opt.tol = o.tol;
                if (o.iter_tol >= 0.0) opt.iter_tol = o.iter_tol;
                if (o.iter_max >= 0) opt.iter_max = static_cast<std::size_t>(o.iter_max);
                return detail::from_double_solution<T>(sn, fluid::solver_fluid_run_analyzer(sn, opt));
            };
        } else {
            throw UnsupportedError(
                "SolverUQ: the fluid stage solver integrates its drift with LSODA, which is double "
                "precision by construction; run UQ under the double backend");
        }
    }
    if (s == "ssa") {
        if constexpr (std::is_same<T, double>::value) {
            return [o](const qn::NetworkStruct<T>& sn) {
                ssa::SsaOptions opt;
                if (!o.method.empty() && o.method != "default") opt.method = o.method;
                if (o.samples) opt.samples = o.samples;
                if (o.seed) opt.seed = o.seed;
                // EVERY DESIGN POINT GETS THE SAME STREAM, deliberately: the
                // points differ by the model, so a per-point seed would mix
                // Monte Carlo error into the epistemic spread the design is
                // there to measure. Common random numbers is the standard
                // variance-reduction pairing for exactly this comparison.
                return detail::from_double_solution<T>(sn, ssa::solver_ssa(sn, opt));
            };
        } else {
            throw UnsupportedError(
                "SolverUQ: an SSA sample path is generated from exponential clocks, which are "
                "transcendental; run UQ under the double backend");
        }
    }
    if (s.empty())
        throw InputError(
            "SolverUQ: no stage solver was named. UQ solves nothing itself; name the solver that "
            "runs at each design point (mva, nc, mam, ba, ctmc, fluid or ssa)");
    throw UnsupportedError("SolverUQ: '" + s +
                           "' is not a stage solver this port carries; the ported names are mva, "
                           "nc, mam, ba, ctmc, fluid and ssa");
}

}  // namespace uq
}  // namespace line

#endif  // LINE_SOLVERS_UQ_UQ_DISPATCH_H
