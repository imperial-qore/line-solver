/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_CTMC_SOLVER_CTMC_CFTP_H
#define LINE_SOLVERS_CTMC_SOLVER_CTMC_CFTP_H

/**
 * The `cftp` and `cftp.approx` methods of SolverCTMC: stationary analysis of a
 * closed single-class product-form network by PERFECT SAMPLING rather than
 * state-space enumeration.
 *
 * Port of matlab/src/solvers/CTMC/solver_ctmc_cftp.m. Reference: S. Kijima and
 * T. Matsui, "Approximate/Perfect Samplers for Closed Jackson Networks", Winter
 * Simulation Conference 2005; the sampler itself is `pfqn::pfqn_cftp`.
 *
 * WHAT KIND OF NUMBER THIS IS. States are drawn iid from the EXACT stationary
 * distribution, so there is no truncation and no cutoff, but the means are
 * sample averages and carry Monte Carlo error O(samples^(-1/2)). Diffing a
 * `cftp` row against an exact solver at solver tolerance therefore reads as a
 * defect and is not one -- it is the same contract as the SSA row. `cftp.approx`
 * additionally drops exactness of the DRAW, running the paper's rapidly-mixing
 * sampler M_A for a deterministic number of updates instead of coupling from the
 * past; its running time is bounded where perfect sampling's is not.
 *
 * WHY THE MODEL CLASS IS GATED SO NARROWLY. The sampler's balance function
 * encodes the closed single-class product form and nothing else, so a model
 * outside that class is REFUSED rather than approximated: the sampler would
 * return states of a different network and the estimator would converge, with
 * shrinking error bars, to the wrong answer.
 *
 * THE TWO ESTIMATORS ARE NOT INTERCHANGEABLE, and the split is deliberate.
 * Throughput is taken at the REFERENCE station and propagated through the visit
 * ratios, which matches the CTMC convention (XN is the arrival rate at the
 * reference station) and keeps flow balance, Little's law and C = N/X exact in
 * the reported table. Utilization instead keeps its own estimator
 * E[min(n_i,c_i)]/c_i, which is unbiased and confined to [0,1] by construction,
 * whereas deriving it from the reference-station throughput lets Monte Carlo
 * error push a saturated station above one.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include "line/api/pfqn/pfqn_cftp.h"
#include "line/api/pfqn/pfqn_mc_common.h"
#include "line/lang/lang_types.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/solvers/ctmc/solver_ctmc.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/mva/sn_chain.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace ctmc {

// NetworkStruct, NodeType and SchedStrategy come in from solver_ctmc.h.

/** The knobs of one perfect-sampling run. */
struct CtmcCftpOptions {
    /** Number of iid stationary draws; the reference has no default here. */
    std::size_t samples = 0;
    /** Stream seed, so a row is reproducible within this port. */
    unsigned long seed = 23000;
};

/** What one `cftp` solve produces beside the means. */
template <class T>
struct CtmcCftpSolution {
    CtmcAvg<T> avg;
    /** (samples x stations) the sampled states, one draw per row. */
    Matrix<int> states;
    /** (samples) per-draw coalescence horizon, or the mixing steps of M_A. */
    std::vector<long> horizon;
    /** The distinct sampled states, aligned with `paggr`. */
    std::vector<std::vector<int>> distinct_states;
    /** Empirical probability of each distinct sampled state. */
    std::vector<T> paggr;
    std::string actualmethod = "cftp";
};

namespace cftp_detail {

/** Port of the reference's method-string to sampler map. */
inline pfqn::CftpMethod cftp_sampler(const std::string& method) {
    if (method == "cftp" || method == "cftp.exact" || method.empty())
        return pfqn::CftpMethod::Cftp;
    if (method == "cftp.approx") return pfqn::CftpMethod::Approx;
    throw InputError("SolverCTMC(cftp): unknown cftp variant '" + method +
                     "'. Use 'cftp' or 'cftp.approx'");
}

/**
 * Can the cftp perfect sampler be asked for this model?
 *
 * The model-class gate asked as a predicate rather than thrown. `cftp_assert`
 * below refuses with it, and the feature-set / AUTO report reaches the same
 * call, so that a caller sees the verdict before paying for a run. A second
 * copy of the rules is how the report and the run drift into two answers.
 *
 * The sampler is exact only on the closed single-class product form its balance
 * function encodes; anything else must be refused, not approximated. What the
 * feature registry CAN name is also declared in `qn::ctmc_feature_set("cftp")`;
 * this predicate carries the structural rules the registry has no name for --
 * the class count, the station count and the phase count.
 *
 * @param sn the refreshed struct of the model
 * @return an empty string when the sampler may run, else the refusal
 */
template <class T>
std::string cftp_supports_reason(const NetworkStruct<T>& sn) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t M = sn.nstations;
    if (sn.nclasses != 1)
        return (
            "SolverCTMC(cftp): the cftp method supports single-class models only, this model has " +
            std::to_string(sn.nclasses) + " classes");
    const std::vector<double> njobs = sn.njobs();
    if (!std::isfinite(njobs[0]) || njobs[0] < 1.0)
        return ("SolverCTMC(cftp): the cftp method supports closed models only, "
                               "with a finite positive population");
    if (M < 2)
        return ("SolverCTMC(cftp): the cftp method requires at least two stations");
    for (std::size_t ind = 0; ind < sn.nodes.size(); ++ind) {
        const NodeType t = sn.nodes[ind].nodetype;
        if (t != NodeType::Queue && t != NodeType::Delay && t != NodeType::Router)
            return (
                "SolverCTMC(cftp): the cftp method supports Queue, Delay and Router nodes only, "
                "node " +
                std::to_string(ind + 1) + " is of a different type");
    }
    for (std::size_t i = 0; i < M; ++i) {
        const SchedStrategy s = sn.stations[i].sched;
        if (!(s == SchedStrategy::INF || s == SchedStrategy::PS || s == SchedStrategy::FCFS ||
              s == SchedStrategy::SIRO || s == SchedStrategy::LCFSPR))
            return (
                "SolverCTMC(cftp): the cftp method requires a product-form scheduling strategy "
                "(INF, PS, FCFS, SIRO, LCFSPR) at station " +
                std::to_string(i + 1));
        if (sn.phases_of(i + 1, 1) > 1)
            return (
                "SolverCTMC(cftp): the cftp method requires exponential service times, station " +
                std::to_string(i + 1) + " has " + std::to_string(sn.phases_of(i + 1, 1)) +
                " phases");
        if (std::isfinite(sn.cap[i]) && sn.cap[i] < njobs[0])
            return (
                "SolverCTMC(cftp): the cftp method requires infinite buffers, station " +
                std::to_string(i + 1) + " has capacity " + std::to_string(sn.cap[i]));
        if (!(sn.rates(i, 0) > zero))
            return ("SolverCTMC(cftp): the cftp method requires a finite positive "
                                   "service rate at station " +
                                   std::to_string(i + 1));
        if (!sn.stations[i].lldscaling.empty() || sn.stations[i].cdscaling ||
            sn.stations[i].jdscaling)
            return ("SolverCTMC(cftp): the cftp method does not support "
                                   "load-dependent, class-dependent or joint-dependent service "
                                   "rates");
    }
    if (!sn.regions.empty())
        return (
            "SolverCTMC(cftp): the cftp method does not support finite capacity regions");
    for (std::size_t ind = 0; ind < sn.nodes.size(); ++ind)
        for (std::size_t r = 0; r < sn.nodes[ind].routing.size(); ++r) {
            const lang::RoutingStrategy rs = sn.nodes[ind].routing[r];
            if (rs != lang::RoutingStrategy::PROB && rs != lang::RoutingStrategy::RAND &&
                rs != lang::RoutingStrategy::DISABLED)
                return (
                    "SolverCTMC(cftp): the cftp method requires Markovian routing (PROB, RAND), "
                    "node " +
                    std::to_string(ind + 1) + " uses a state-dependent strategy");
        }
    return std::string();
}

/**
 * Port of the reference's model-class gate.
 *
 * The sampler is exact only on the closed single-class product form its balance
 * function encodes; anything else must be refused, not approximated.
 */
template <class T>
void cftp_assert(const NetworkStruct<T>& sn) {
    const std::string reason = cftp_supports_reason(sn);
    if (!reason.empty()) throw UnsupportedError(reason);
}

}  // namespace cftp_detail

/**
 * The `cftp` model-class gate as a public predicate.
 *
 * The one call a REPORT can make: `auto_family_refusal` reaches the structural
 * rules the feature set has no name for (the class count, the station count and
 * the phase count) through this, and `solver_ctmc_cftp` refuses through the
 * same body, so a pair the report offers is a pair the sampler runs.
 *
 * @param sn the refreshed struct of the model
 * @return an empty string when the sampler may run, else the refusal
 */
template <class T>
std::string solver_ctmc_cftp_supports(const NetworkStruct<T>& sn) {
    return cftp_detail::cftp_supports_reason(sn);
}

/**
 * Solve with the `cftp` / `cftp.approx` method.
 *
 * @param sn the refreshed struct of a CLOSED single-class product-form network
 * @param opt the SolverCTMC knobs; `method` picks the sampler
 * @param cftpopt the run length and the stream
 */
template <class T>
CtmcCftpSolution<T> solver_ctmc_cftp(const NetworkStruct<T>& sn, const CtmcOptions& opt,
                                     const CtmcCftpOptions& cftpopt) {
    static_assert(num_traits<T>::has_transcendental,
                  "solver_ctmc_cftp draws random states and forms the station balance functions "
                  "in the log domain, neither of which exists in exact rational arithmetic");
    const T zero = num_traits<T>::from_int(0);
    const std::size_t M = sn.nstations, R = sn.nclasses;

    const pfqn::CftpMethod sampler = cftp_detail::cftp_sampler(opt.method);
    cftp_detail::cftp_assert(sn);
    if (cftpopt.samples < 1)
        throw InputError("SolverCTMC(cftp): the cftp method requires a finite positive number of "
                         "samples; pass --samples");

    const mva::ChainDemands<T> ch = mva::sn_get_demands_chain(sn);
    std::vector<T> L(M, zero);
    for (std::size_t i = 0; i < M; ++i) L[i] = ch.Lchain(i, 0);
    const int N = static_cast<int>(ch.Nchain[0]);
    std::vector<int> S(M, 1);
    for (std::size_t i = 0; i < M; ++i)
        S[i] = sn.stations[i].sched == SchedStrategy::INF
                   ? pfqn::cftp_inf_servers
                   : static_cast<int>(sn.stations[i].nservers);

    pfqn::McRng rng(static_cast<std::uint64_t>(cftpopt.seed));
    const pfqn::CftpResult<T> s = pfqn::pfqn_cftp(L, N, S, cftpopt.samples, sampler, rng);

    CtmcCftpSolution<T> o;
    o.avg.QN = Matrix<T>(M, R, zero);
    o.avg.UN = Matrix<T>(M, R, zero);
    o.avg.RN = Matrix<T>(M, R, zero);
    o.avg.TN = Matrix<T>(M, R, zero);
    o.avg.XN.assign(R, zero);
    o.avg.CN.assign(R, zero);

    // E[min(n_i, c_i)]: the busy-server count, which is the utilization law's
    // numerator and, at the reference station, the throughput law's.
    std::vector<T> busy(M, zero);
    for (std::size_t i = 0; i < M; ++i) {
        T acc = zero;
        for (std::size_t r = 0; r < cftpopt.samples; ++r) {
            const int n = s.X(r, i);
            const int cap = S[i] == pfqn::cftp_inf_servers ? n : (n < S[i] ? n : S[i]);
            acc += num_traits<T>::from_int(cap);
        }
        busy[i] = T(acc / num_traits<T>::from_int(static_cast<int>(cftpopt.samples)));
    }

    const std::size_t iref = ch.refstatchain[0] - 1;
    if (ch.STchain(iref, 0) > zero && ch.Vchain(iref, 0) > zero)
        o.avg.XN[0] = T(busy[iref] / ch.STchain(iref, 0) / ch.Vchain(iref, 0));
    if (o.avg.XN[0] > zero)
        o.avg.CN[0] = T(num_traits<T>::from_int(N) / o.avg.XN[0]);

    for (std::size_t i = 0; i < M; ++i) {
        o.avg.QN(i, 0) = s.Q[i];
        o.avg.TN(i, 0) = T(ch.Vchain(i, 0) * o.avg.XN[0]);
        o.avg.UN(i, 0) = sn.stations[i].sched == SchedStrategy::INF
                             ? o.avg.QN(i, 0)
                             : T(busy[i] / num_traits<T>::from_double(sn.stations[i].nservers));
        if (o.avg.TN(i, 0) > zero) o.avg.RN(i, 0) = T(o.avg.QN(i, 0) / o.avg.TN(i, 0));
    }

    // The empirical law over the DISTINCT draws, the reference's [SSq, pAggr].
    // Ordered by the state tuple so two runs of the same tape list them alike.
    std::map<std::vector<int>, std::size_t> counts;
    for (std::size_t r = 0; r < cftpopt.samples; ++r) {
        std::vector<int> row(M, 0);
        for (std::size_t i = 0; i < M; ++i) row[i] = s.X(r, i);
        ++counts[row];
    }
    for (std::map<std::vector<int>, std::size_t>::const_iterator it = counts.begin();
         it != counts.end(); ++it) {
        o.distinct_states.push_back(it->first);
        o.paggr.push_back(T(num_traits<T>::from_int(static_cast<int>(it->second)) /
                            num_traits<T>::from_int(static_cast<int>(cftpopt.samples))));
    }

    o.states = s.X;
    o.horizon = s.horizon;
    o.actualmethod = opt.method.empty() ? std::string("cftp") : opt.method;
    return o;
}

/** Solve with `cftp` and format, for a caller with no use for the sampled states. */
template <class T>
mva::AvgResult<T> solver_ctmc_cftp_run_analyzer(const NetworkStruct<T>& sn, const CtmcOptions& opt,
                                       const CtmcCftpOptions& cftpopt) {
    const CtmcCftpSolution<T> s = solver_ctmc_cftp(sn, opt, cftpopt);
    CtmcSolution<T> d;
    d.avg = s.avg;
    d.actualmethod = s.actualmethod;
    return solver_ctmc_avg_table(sn, d, opt.method);
}

}  // namespace ctmc
}  // namespace line

#endif  // LINE_SOLVERS_CTMC_SOLVER_CTMC_CFTP_H
