/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_MVA_SOLVER_MVA_SJN_H
#define LINE_SOLVERS_MVA_SOLVER_MVA_SJN_H

/**
 * Closed networks with shortest-job-next (SJF) stations, ladder branch 0.
 *
 * Templated port of `matlab/src/solvers/MVA/solver_mva_sjn_analyzer.m`,
 * cross-checked against
 * `jar/src/main/java/jline/solvers/mva/analyzers/Solver_mva_sjn_analyzer.java`.
 *
 * A non-preemptive SJF station serves the shortest queued job first, the job
 * size being known on arrival. It is not product form and no ordinary MVA
 * equation covers it, so the station is modelled by the conditional waiting
 * time equation of K. Kant, "MVA approximations for SJN scheduling",
 * Performance Evaluation 15(1):41-61, 1992, which `pfqn_sjn.h` carries.
 *
 * WHY THE LADDER PUTS THIS FIRST. `mvaDispatch.m` tests `hasSJN` ahead of every
 * other branch, and the branches are NOT disjoint: a closed model with an SJF
 * station and a delay would otherwise fall through to the generic AMVA path,
 * which reads only the mean service time and would report the SJF station as if
 * it scheduled size-blind. The answer would be a plausible number for a
 * different model, which is worse than a refusal.
 *
 * CLOSED MODELS ONLY, BY NAME. The conditional waiting time equation is a
 * POPULATION recursion, so the open case has nothing to recur over. MATLAB's
 * getFeatureSet declares SchedStrategy_SJF unconditionally and the dispatch
 * then refuses the open case by name; the boolean registry cannot express
 * "this discipline, but only in a closed model", so the C++ gate declares the
 * feature and this file carries the same imperative refusal.
 *
 * LATTICE OR FIXED POINT. The exact recursion steps over prod(N+1) states, so
 * `default` switches to the Schweitzer closure once the lattice exceeds
 * `kLatticeMax`. Ask for `exact` or `mva` to force the lattice, `amva` or `bs`
 * to force the fixed point; `sjn.mva` and `sjn.amva` name them directly. The
 * reported `actualmethod` is the one actually run.
 *
 * OTHER STATIONS. The remaining stations are solved with the single-server MVA
 * equation, so the reference restricts them to INF (folded into the think time)
 * and single-server PS, LCFS-PR, FCFS or SIRO. Anything else, and any
 * multi-server queueing station, is refused by name.
 *
 * WARNINGS, AND HOW FAR THEY GET. The reference calls line_warning in two
 * places that still return a usable answer: the utilization cap binding, and
 * the fixed point exhausting its iterations. `pfqn_sjn.h` turns both into
 * result flags because the api layer has no warning channel; this file turns
 * the flags back into the reference's own text on
 * `SjnAnalyzerResult::warning`, which `mva_dispatch.h` copies onto
 * `DispatchResult::warning`. That is as far as it goes: this port has NO
 * general warning facility, and `solver_mva_runner.h` does not read the field,
 * so a user going through SolverMVA still does not see it. Recorded here as a
 * known divergence rather than left silent -- the capped answer is stable and
 * satisfies the population law exactly, but its accuracy is not warranted, and
 * a number that looks authoritative while carrying an accuracy claim the
 * reference declines to make is worse than a visibly wrong one.
 *
 * Arithmetic: TRANSCENDENTAL. The recursion evaluates regularized incomplete
 * gammas on a Simpson grid; `pfqn_sjn.h` is a double-precision API for that
 * reason, and the exact (Rational) path is refused by name below.
 */

#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "line/api/pfqn/pfqn_sjn.h"
#include "line/lang/lang_types.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/mva/mva_types.h"
#include "line/solvers/mva/sn_chain.h"

namespace line {
namespace mva {

/**
 * Lattice size above which `default` prefers the fixed point.
 *
 * The reference reads it, and the grid parameters ns / Lfactor / umax, from
 * `options.config.sjn_*`. MvaOptions carries no config map, so the four sit at
 * the reference's own defaults here; the method name still selects the route.
 */
inline constexpr double kSjnLatticeMax = 1e5;

/** What the analyzer returns, the reference's metrics plus its `actualmethod`. */
template <class T>
struct SjnAnalyzerResult {
    MvaSolution<T> sol;
    std::string actualmethod;
    /**
     * Non-empty when the reference would have called line_warning AND RETURNED,
     * carrying its text. Today: the utilization cap binding, and the fixed
     * point running out of iterations.
     *
     * This is the shape solver_nc_cdf.h:83 established for a reference warning
     * that accompanies a usable answer, as opposed to solver_mam_retrial.h:437,
     * which throws because its warning accompanies an answer the port declines
     * to stand behind. The SJN cap is the first kind: the numbers are stable
     * and the population law still holds exactly, but their ACCURACY is not
     * warranted, and MATLAB says so out loud. Dropping it would leave the
     * caller an authoritative-looking number with a silent accuracy claim the
     * reference explicitly declines to make.
     */
    std::string warning;
};

/** True when the layer has an SJF station, the reference's `any(sn.sched == SchedStrategy.SJF)`. */
template <class T>
bool sn_has_sjn(const qn::NetworkStruct<T>& L) {
    for (const auto& st : L.stations)
        if (st.sched == lang::SchedStrategy::SJF) return true;
    return false;
}

/** Port of `solver_mva_sjn_analyzer.m`. */
template <class T>
SjnAnalyzerResult<T> solver_mva_sjn_analyzer(const qn::NetworkStruct<T>& L, const MvaOptions& opt) {
    using lang::SchedStrategy;
    if constexpr (!num_traits<T>::has_transcendental) {
        (void)L;
        (void)opt;
        throw UnsupportedError(
            "solver_mva_sjn_analyzer: the SJN conditional waiting time recursion needs "
            "transcendental arithmetic (regularized incomplete gammas on a quadrature grid); "
            "rerun this model with --arith double or --arith real");
    } else {
    const std::size_t M = L.nstations, C = L.nchains;
    const ChainDemands<T> d = sn_get_demands_chain(L);

    for (std::size_t c = 0; c < C; ++c)
        if (std::isinf(d.Nchain[c]))
            throw UnsupportedError(
                "SolverMVA supports shortest-job-next (SJF) scheduling only in closed models, the "
                "conditional waiting time equation being a population recursion. Use SolverLDES, "
                "or SolverMVA with SRPT or PSJF for the preemptive size-based open queue");

    // rows: the queueing stations handed to the SJN recursion; infrows: the
    // delay stations, folded into the think time; sjnrows: positions of the SJN
    // stations WITHIN rows, which is the indexing pfqn_sjn expects
    std::vector<std::size_t> rows, infrows, sjnrows;
    for (std::size_t i = 0; i < M; ++i) {
        const SchedStrategy s = L.stations[i].sched;
        const double k = L.stations[i].nservers;
        if (s == SchedStrategy::EXT) continue;  // no arrival stream in a closed model
        if (s == SchedStrategy::INF) {
            infrows.push_back(i);
            continue;
        }
        if (s == SchedStrategy::SJF) {
            if (k != 1.0)
                throw UnsupportedError(
                    "solver_mva_sjn_analyzer: SJN scheduling at station " +
                    std::to_string(i + 1) +
                    " requires a single server, the response time equation is a single-server one");
            rows.push_back(i);
            sjnrows.push_back(rows.size() - 1);
            continue;
        }
        if (s == SchedStrategy::PS || s == SchedStrategy::LCFSPR || s == SchedStrategy::FCFS ||
            s == SchedStrategy::SIRO) {
            if (k != 1.0)
                throw UnsupportedError(
                    "solver_mva_sjn_analyzer: station " + std::to_string(i + 1) + " has " +
                    (std::isfinite(k) ? std::to_string(static_cast<long long>(k))
                                      : std::string("infinitely many")) +
                    " servers, the SJN analyzer solves the remaining stations with the "
                    "single-server MVA equation");
            rows.push_back(i);
            continue;
        }
        throw UnsupportedError("solver_mva_sjn_analyzer: the SJN analyzer does not support " +
                               std::string(lang::sched_to_text(s)) +
                               " scheduling at the other stations");
    }

    const std::size_t Mq = rows.size();
    Matrix<double> Ld(Mq, C, 0.0), Vd(Mq, C, 0.0), scvd(Mq, C, 1.0);
    for (std::size_t j = 0; j < Mq; ++j) {
        const std::size_t i = rows[j];
        for (std::size_t c = 0; c < C; ++c) {
            Ld(j, c) = num_traits<T>::to_double(T(d.STchain(i, c) * d.Vchain(i, c)));
            Vd(j, c) = num_traits<T>::to_double(d.Vchain(i, c));
        }
    }
    // only the SJN stations read an SCV: the size distribution is reconstructed
    // from it, and the reference leaves the others at one
    for (std::size_t j : sjnrows) {
        const std::size_t i = rows[j];
        for (std::size_t c = 0; c < C; ++c) {
            const double v = num_traits<T>::to_double(d.SCVchain(i, c));
            if (std::isfinite(v) && v > 0.0) scvd(j, c) = v;
        }
    }
    std::vector<double> Zd(C, 0.0), Nd(C, 0.0);
    for (std::size_t c = 0; c < C; ++c) {
        for (std::size_t i : infrows)
            Zd[c] += num_traits<T>::to_double(T(d.STchain(i, c) * d.Vchain(i, c)));
        Nd[c] = d.Nchain[c];
    }

    pfqn::SjnOptions sjnopt;
    sjnopt.tol = opt.iter_tol;
    sjnopt.iter_max = static_cast<std::size_t>(opt.iter_max);
    // SJN applies within a class and the classes are then non-preemptively
    // prioritised; without distinct priorities the jobs of every class are
    // compared by size directly (the pooled reading, options.prio empty).
    // The reference indexes classprio by CLASS and hands it to a per-CHAIN
    // argument, which is well defined only because the two counts agree here.
    if (L.nchains == L.nclasses) {
        bool distinct = true;
        for (std::size_t i = 0; i < L.nclasses && distinct; ++i)
            for (std::size_t j = i + 1; j < L.nclasses; ++j)
                if (L.classes[i].prio == L.classes[j].prio) {
                    distinct = false;
                    break;
                }
        if (distinct) {
            sjnopt.prio.assign(C, 0);
            for (std::size_t c = 0; c < C; ++c) sjnopt.prio[c] = L.classes[c].prio;
        }
    }

    double lattice = 1.0;
    for (std::size_t c = 0; c < C; ++c) lattice *= Nd[c] + 1.0;
    bool uselattice;
    if (opt.method == "amva" || opt.method == "bs" || opt.method == "sjn.amva")
        uselattice = false;
    else if (opt.method == "exact" || opt.method == "mva" || opt.method == "sjn.mva")
        uselattice = true;
    else
        uselattice = lattice <= kSjnLatticeMax;

    pfqn::SjnResult sjn;
    if (uselattice) {
        try {
            sjn = pfqn::pfqn_mvasjn(Ld, Nd, Zd, scvd, sjnrows, Vd, sjnopt);
        } catch (const pfqn::SjnStarvationError&) {
            if (opt.method != "default") throw;
            sjn = pfqn::pfqn_amvasjn(Ld, Nd, Zd, scvd, sjnrows, Vd, sjnopt);
            uselattice = false;
        }
    } else {
        sjn = pfqn::pfqn_amvasjn(Ld, Nd, Zd, scvd, sjnrows, Vd, sjnopt);
    }

    Matrix<T> Qchain(M, C, num_traits<T>::from_int(0)), Uchain(M, C, num_traits<T>::from_int(0));
    Matrix<T> Rchain(M, C, num_traits<T>::from_int(0)), Tchain(M, C, num_traits<T>::from_int(0));
    std::vector<T> Xchain(C, num_traits<T>::from_int(0));
    for (std::size_t c = 0; c < C; ++c) Xchain[c] = num_traits<T>::from_double(sjn.XN[c]);
    for (std::size_t j = 0; j < Mq; ++j) {
        const std::size_t i = rows[j];
        for (std::size_t c = 0; c < C; ++c) {
            Qchain(i, c) = num_traits<T>::from_double(sjn.QN(j, c));
            Uchain(i, c) = num_traits<T>::from_double(sjn.UN(j, c));
        }
    }
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t c = 0; c < C; ++c) Tchain(i, c) = T(Xchain[c] * d.Vchain(i, c));
    // a delay station holds T S jobs, all of them in service
    for (std::size_t i : infrows)
        for (std::size_t c = 0; c < C; ++c) {
            Qchain(i, c) = T(Tchain(i, c) * d.STchain(i, c));
            Uchain(i, c) = Qchain(i, c);
        }
    const T zero = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t c = 0; c < C; ++c)
            if (Tchain(i, c) != zero) Rchain(i, c) = T(Qchain(i, c) / Tchain(i, c));
    // The reference also runs an isfinite sweep here. It has nothing to do in
    // this port: every division above is guarded, so no non-finite value can be
    // written in the first place. What does remain is the empty-chain columns.
    for (std::size_t c = 0; c < C; ++c) {
        if (d.Nchain[c] != 0.0) continue;
        Xchain[c] = zero;
        for (std::size_t i = 0; i < M; ++i) {
            Qchain(i, c) = zero;
            Uchain(i, c) = zero;
            Rchain(i, c) = zero;
            Tchain(i, c) = zero;
        }
    }

    // Qchain and Uchain are DELIBERATELY not handed to the deaggregation: the
    // reference passes [] in both slots (solver_mva_sjn_analyzer.m, last line)
    // and lets it rebuild them from Rchain, Tchain and the per-class alpha. The
    // Java port passes them and so disagrees; MATLAB is ground truth.
    const ClassResults<T> cr = sn_deaggregate_chain_results(L, d, Matrix<T>(), Matrix<T>(), Rchain,
                                                            Tchain, Xchain);
    SjnAnalyzerResult<T> out;
    out.sol.Q = cr.Q;
    out.sol.U = cr.U;
    out.sol.R = cr.R;
    out.sol.Tp = cr.Tp;
    out.sol.C = cr.C;
    out.sol.X = cr.X;
    out.sol.method = opt.method;
    out.sol.iter = static_cast<int>(sjn.iter);
    // the reference returns lG = NaN: an SJN solve carries no normalizing constant
    out.sol.lG = std::numeric_limits<double>::quiet_NaN();
    out.actualmethod = uselattice ? "sjn.mva" : "sjn.amva";
    // Carry the reference's two line_warning texts rather than drop them.
    if (sjn.capped)
        out.warning =
            "the utilization cap of 0.999 was binding at an SJN station: the station is in the "
            "starvation regime, where long jobs are held back and the arrival theorem is badly "
            "violated. The results are stable but their accuracy is not warranted, use SolverCTMC "
            "or SolverLDES there";
    if (!sjn.converged) {
        if (!out.warning.empty()) out.warning += "; ";
        out.warning += "the SJN fixed point did not converge in " + std::to_string(sjn.iter) +
                       " iterations";
    }
    return out;
    }  // if constexpr has_transcendental
}

}  // namespace mva
}  // namespace line

#endif  // LINE_SOLVERS_MVA_SOLVER_MVA_SJN_H
