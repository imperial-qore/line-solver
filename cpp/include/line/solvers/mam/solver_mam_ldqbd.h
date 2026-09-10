/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_MAM_SOLVER_MAM_LDQBD_H
#define LINE_SOLVERS_MAM_SOLVER_MAM_LDQBD_H

/**
 * Port of `solver_mam_ldqbd.m`: the level-dependent QBD analyzer for a
 * single-class network of one infinite server and one FCFS queue.
 *
 * WHY THIS MATTERS MORE THAN ITS SIZE SUGGESTS. It is the one branch of the MAM
 * ladder where the reference prefers an EXACT method over the `dec.source`
 * decomposition, and `solver_mam_analyzer.m` routes `default` here for a
 * single-class closed Delay+Queue. The level-dependent arrival rate
 * `(N - n) lambda` captures the population constraint that `dec.source` can only
 * approximate through its throughput fixed point, so on this shape the two
 * answers are not close: on the test model below the LD-QBD queue length is
 * 1.36 against dec.source's 1.42, and the LD-QBD one is right.
 *
 * EXACTNESS: exact for exponential service at any number of servers, and for PH
 * service at any number of servers. The multiserver PH chain comes from
 * `ldqbd_mphc`, whose level coordinate is the MULTISET of the phases the
 * min(n,c) busy servers sit in. The collapsed single-phase approximation the
 * reference carried until 2026-08-18 -- one PH process run at min(n,c) times its
 * speed, ~1e-2 relative against SolverCTMC -- is gone from every codebase.
 *
 * TWO REGIMES, one generator. Closed: level n is the number at the queue,
 * 0 <= n <= N, and the arrival rate out of level n is `(N-n) lambda_eff`, which
 * vanishes at n = N and closes the chain by itself. Open: Poisson arrivals at a
 * constant `lambda_eff`, truncated at a level chosen so the tail probability is
 * below 1e-10 (or at `options.cutoff`). Only the per-level arrival rate and the
 * top level differ.
 *
 * ARITHMETIC. `ldqbd_R` runs a backward recursion of matrix inverses and its
 * `pinv` fallback needs singular vectors, so the whole path is gated on
 * transcendental arithmetic and refuses by name under exact/Rational.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "line/api/mam/ldqbd.h"
#include "line/api/mam/ldqbd_mphc.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/qbd_setupdelayoff.h"
#include "line/lang/distribution.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/mam/mam_types.h"
#include "line/solvers/mva/mva_types.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/**
 * The LD-QBD blocks and parameters, the reference's optional eighth output.
 *
 * Built unconditionally here rather than behind a `nargout` test: C++ has no
 * such thing, the cost is a few pointers, and the SolverENV state-vector
 * analyzer is the consumer the reference built it for.
 */
template <class T>
struct LdqbdBlocks {
    std::vector<Matrix<T>> Q0, Q1, Q2;
    std::size_t Nlev = 0;
    std::size_t nPhases = 1;
    bool isPH = false;
    bool isOpen = false;
    std::size_t queueIdx = 0, refIdx = 0, M = 0;
    double nServers = 1.0;
    T mean_service = num_traits<T>::from_int(0);
    bool hasLLD = false;
    /** Per-level service factor sf(n), with sf[0] unused so it lines up by level. */
    std::vector<T> sf;
    /**
     * The capacity that normalizes the utilization: max(c, max(alpha)), the
     * LARGEST factor the load-dependence table declares rather than the
     * saturated one, since a non-monotone alpha peaks in the middle. Same rule
     * as CTMC's ceff, which is what makes the two report the same number.
     */
    double utilPeak = 1.0;
    T lambda_eff = num_traits<T>::from_int(0);
    T delayRate = num_traits<T>::from_int(0);
    double N = 0.0;
    /**
     * The station alternates OFF -> setup -> busy -> delay-off around the
     * service, so the chain carries phases the block builder has no place for.
     * When this is set the blocks above describe a server that is ALWAYS warm
     * and must not be used: the closed regime hands the whole chain to
     * `qbd_setupdelayoff_closed` instead.
     */
    bool hasSetup = false;
    T alpharate = num_traits<T>::from_int(0), alphascv = num_traits<T>::from_int(1);
    T betarate = num_traits<T>::from_int(0), betascv = num_traits<T>::from_int(1);
};

/** What the analyzer returns: the metrics plus the blocks it built them from. */
template <class T>
struct LdqbdSolution {
    mva::MvaSolution<T> sol;
    LdqbdBlocks<T> ld;
};

/**
 * Port of `solver_mam_ldqbd.m`.
 *
 * @param L   the refreshed struct; must be single-class, two stations, and
 *            either Delay+Queue (closed) or Source+Queue (open)
 * @param opt the MAM options; `cutoff` bounds the open truncation
 */
template <class T>
LdqbdSolution<T> solver_mam_ldqbd(const qn::NetworkStruct<T>& L, const MamOptions& opt) {
    if constexpr (!num_traits<T>::has_transcendental) {
        throw UnsupportedError(
            "solver_mam_ldqbd: the level-dependent QBD recursion inverts a matrix per level and "
            "falls back to a pseudo-inverse (singular vectors) when a level is singular, neither "
            "of which is exact arithmetic; rerun with --arith double or --arith real");
    } else {
    using lang::SchedStrategy;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t M = L.nstations, K = L.nclasses;

    if (K != 1)
        throw UnsupportedError("solver_mam_ldqbd: the LDQBD method requires a single-class model");

    std::size_t nDelay = 0, nQueue = 0, nSource = 0;
    std::size_t delayIdx = 0, queueIdx = 0, srcIdx = 0;
    for (std::size_t i = 1; i <= M; ++i) {
        const SchedStrategy s = L.stations[i - 1].sched;
        if (s == SchedStrategy::INF) { ++nDelay; delayIdx = i; }
        else if (s == SchedStrategy::FCFS) { ++nQueue; queueIdx = i; }
        else if (s == SchedStrategy::EXT) { ++nSource; srcIdx = i; }
    }
    const bool isOpen = std::isinf(L.classes[0].population);
    if (isOpen) {
        if (nSource != 1 || nQueue != 1 || M != 2)
            throw UnsupportedError(
                "solver_mam_ldqbd: the open LDQBD method requires exactly one Source and one "
                "Queue station");
    } else {
        if (nDelay != 1 || nQueue != 1 || M != 2)
            throw UnsupportedError(
                "solver_mam_ldqbd: the closed LDQBD method requires exactly one Delay and one "
                "Queue station");
    }

    // ---- the service process at the queue --------------------------------
    const Map<T> PHq = lang::dist_to_map(L.service[queueIdx - 1][0]);
    const double nServers = L.stations[queueIdx - 1].nservers;
    const std::size_t nPhases = PHq.D0.rows();
    const bool isPH = nPhases > 1;
    T mu = zero, mean_service = zero;
    std::vector<T> alpha;
    if (!isPH) {
        mu = T(-PHq.D0(0, 0));
        if (!(num_traits<T>::to_double(mu) > 0.0))
            throw InputError("solver_mam_ldqbd: the queue has a non-positive service rate");
        mean_service = T(one / mu);
    } else {
        alpha = map_pie(PHq);
        mean_service = map_mean(PHq);
    }

    // ---- setup and delay-off at the queue --------------------------------
    // Refused BY NAME outside the closed, single-server, exponential,
    // load-independent case rather than answered as if the server were always
    // warm, which is what this solver did until 2026-09 and is BUG-78.
    bool hasSetup = false;
    T alpharate = zero, alphascv = one, betarate = zero, betascv = one;
    {
        const typename std::map<std::size_t, qn::SetupDelayOffParam<T>>::const_iterator sit =
            L.setupparam.find(queueIdx);
        if (sit != L.setupparam.end()) {
            lang::Distrib<T> su, doff;
            if (sit->second.last(su, doff) && !doff.disabled) {
                if (isOpen)
                    throw InputError(
                        "solver_mam_ldqbd: open LDQBD does not model a setup/delay-off server; "
                        "use method 'dec.source', whose qbd_setupdelayoff covers the open case");
                if (isPH || nServers > 1)
                    throw InputError(
                        "solver_mam_ldqbd: closed LDQBD models a setup/delay-off server with "
                        "exponential service at a single server only; this station has "
                        "phase-type service or several servers");
                hasSetup = true;
                alpharate = T(one / su.mean);
                alphascv = su.scv;
                betarate = T(one / doff.mean);
                betascv = doff.scv;
            }
        }
    }

    // ---- the per-level service factor ------------------------------------
    // sn.lldscaling when present, else min(n, c).
    const std::vector<T>& lld = L.stations[queueIdx - 1].lldscaling;
    bool hasLLD = false;
    for (const T& v : lld)
        if (v != one) hasLLD = true;
    const double sfMax = hasLLD ? num_traits<T>::to_double(lld.back()) : nServers;
    // The capacity that normalizes the utilization is the LARGEST factor the
    // table declares, not the saturated one: a non-monotone alpha peaks in the
    // middle. Same rule as CTMC's ceff = max(nservers, max(lldscaling(ist,:))),
    // which is what makes the two report the same number.
    double utilPeak = nServers;
    if (hasLLD)
        for (const T& v : lld) utilPeak = std::max(utilPeak, num_traits<T>::to_double(v));

    // ---- the per-level arrival rate and the number of levels -------------
    const std::size_t Kc = K;
    auto rt_at = [&](std::size_t from, std::size_t to) -> T {
        const std::size_t a = (L.stateful_of_station(from) - 1) * Kc;
        const std::size_t b = (L.stateful_of_station(to) - 1) * Kc;
        if (a >= L.rt.rows() || b >= L.rt.cols()) return zero;
        return L.rt(a, b);
    };

    std::size_t Nlev = 0;
    T lambda_eff = zero, delayRate = zero;
    std::vector<T> arrRate;
    if (isOpen) {
        const Map<T> arv = lang::dist_to_map(L.service[srcIdx - 1][0]);
        if (arv.D0.rows() > 1)
            throw UnsupportedError(
                "solver_mam_ldqbd: the open LDQBD method currently supports Poisson (exponential) "
                "arrivals only; the Source uses a MAP/MMPP process");
        const T lambda = L.rates(srcIdx - 1, 0);
        lambda_eff = T(lambda * rt_at(srcIdx, queueIdx));
        const double rho =
            num_traits<T>::to_double(T(lambda_eff * mean_service)) / (sfMax > 0.0 ? sfMax : 1.0);
        if (rho >= 1.0)
            throw NumericError(
                "solver_mam_ldqbd: the open LDQBD method requires a stable queue (rho = " +
                std::to_string(rho) +
                " >= 1). Increase service capacity or reduce the arrival rate");
        const std::size_t c = static_cast<std::size_t>(
            std::isfinite(nServers) ? std::llround(nServers) : 1);
        if (opt.cutoff > 0) {
            Nlev = std::max(c + 1, opt.cutoff);
        } else {
            const double tailTol = 1e-10;
            const long lv =
                static_cast<long>(c) + static_cast<long>(std::ceil(std::log(tailTol) / std::log(rho)));
            Nlev = static_cast<std::size_t>(
                std::min<long>(std::max<long>(lv, static_cast<long>(c) + 10), 100000));
        }
        arrRate.assign(Nlev + 1, lambda_eff);
        arrRate[Nlev] = zero;  // truncation: no arrivals above the top level
    } else {
        delayRate = L.rates(delayIdx - 1, 0);
        lambda_eff = T(delayRate * rt_at(delayIdx, queueIdx));
        const double Nd = L.classes[0].population;
        Nlev = static_cast<std::size_t>(std::llround(Nd));
        arrRate.assign(Nlev + 1, zero);
        // Finite-source rate (N - n) lambda_eff, which is zero at n = N and
        // closes the chain without a truncation.
        for (std::size_t n = 0; n <= Nlev; ++n)
            arrRate[n] = T(num_traits<T>::from_double(Nd - static_cast<double>(n)) * lambda_eff);
    }
    if (Nlev < 1)
        throw UnsupportedError(
            "solver_mam_ldqbd: the model has no levels to solve (a zero population)");

    std::vector<T> sf(Nlev + 1, zero);  // sf[n] for n = 1..Nlev
    for (std::size_t n = 1; n <= Nlev; ++n) {
        if (hasLLD)
            sf[n] = lld[std::min(n, lld.size()) - 1];
        else
            sf[n] = num_traits<T>::from_double(
                std::min(static_cast<double>(n), std::isfinite(nServers) ? nServers : 1.0));
    }

    // ---- the block-tridiagonal generator ---------------------------------
    // Q0[n] level n -> n+1 (arrival), Q1[n] local, Q2[n] level n -> n-1.
    // Q2 carries an unused entry at index 0 so the three line up by level, as
    // the C++ ldqbd takes them.
    std::vector<Matrix<T>> Q0(Nlev), Q1(Nlev + 1), Q2(Nlev + 1);
    if (!isPH) {
        for (std::size_t n = 0; n + 1 <= Nlev; ++n) Q0[n] = Matrix<T>(1, 1, arrRate[n]);
        for (std::size_t n = 0; n <= Nlev; ++n) {
            const T dep = (n > 0) ? T(sf[n] * mu) : zero;
            Q1[n] = Matrix<T>(1, 1, T(-(arrRate[n] + dep)));
        }
        Q2[0] = Matrix<T>(1, 1, zero);
        for (std::size_t n = 1; n <= Nlev; ++n) Q2[n] = Matrix<T>(1, 1, T(sf[n] * mu));
    } else {
        // PH service: the level carries the MULTISET of the phases the min(n,c)
        // busy servers sit in, which is exact at any number of servers. At c = 1
        // the multiset is just the phase, so this reproduces the single-server
        // blocks (sf(n) D0 - arr I, sf(n) D1) entry for entry. `sf` is indexed
        // from 1 here and from 0 there, hence the shifted copy.
        std::vector<T> sf1(Nlev, zero);
        for (std::size_t n = 1; n <= Nlev; ++n) sf1[n - 1] = sf[n];
        LdqbdMphcBlocks<T> blk =
            ldqbd_mphc(PHq.D0, PHq.D1, alpha, nServers, arrRate, sf1);
        Q0 = blk.Q0;
        Q1 = blk.Q1;
        Q2 = blk.Q2;
    }

    std::vector<T> p;
    T mean_queue = zero, x_setup = zero;
    if (hasSetup) {
        // SETUP AND DELAY-OFF, the closed vacation queue. The level-dependent
        // chain this needs is the one above with two extra phase families -- the
        // setup above level 0 and the delay-off at level 0 -- and
        // `qbd_setupdelayoff_closed` builds and solves exactly that, so it is
        // called rather than duplicated. Without it the blocks above describe a
        // server that is ALWAYS warm and the answer is byte-identical across any
        // setup mean (BUG-78).
        const mam::SetupDelayoffClosed<T> cr = mam::qbd_setupdelayoff_closed(
            num_traits<T>::from_double(L.classes[0].population), T(one / lambda_eff), mu,
            alpharate, alphascv, betarate, betascv);
        mean_queue = cr.QN;
        x_setup = cr.XN;
    } else {
        const LdqbdResult<T> res = ldqbd(Q0, Q1, Q2);
        p = res.pi.pi;
        if (p.size() != Nlev + 1)
            throw NumericError("solver_mam_ldqbd: the LD-QBD solve returned the wrong level count");

        for (std::size_t n = 0; n <= Nlev; ++n)
            mean_queue += num_traits<T>::from_int(static_cast<int>(n)) * p[n];
    }

    // Utilization is the fraction of the station's PEAK capacity in use,
    // sum_n p(n)*sf(n)/utilPeak, the work-based convention CTMC, MVA, NC and
    // serial SSA all report. Without load dependence sf(n) = min(n,c) and
    // utilPeak = c, giving the mean fraction of the c servers in use; at c = 1
    // that is sf(n) = 1 for every n >= 1, so the sum collapses to 1 - p(0).
    //
    // It used to report 1 - p(0) under load dependence, i.e. P(busy), which
    // reads a station running alpha(n) times faster as no busier than one at
    // its nominal rate: 0.9587 against CTMC's 0.6612 on a 4-job closed model
    // with alpha = [1 1.5 2 2.5].
    T util = zero;
    if (hasSetup) {
        // With a setup the server is DELIVERING work only in the busy phase, so
        // the level occupancy over-counts it: a level is occupied during the
        // setup too. The utilization law gives the same work-based number
        // without the per-phase vector, X*E[S]/peak, which is what the level sum
        // reduces to without a setup.
        util = T(x_setup * mean_service / num_traits<T>::from_double(utilPeak));
    } else {
        const T peak = num_traits<T>::from_double(utilPeak);
        for (std::size_t n = 1; n <= Nlev; ++n) util += T(sf[n] / peak * p[n]);
    }

    LdqbdSolution<T> out;
    mva::MvaSolution<T>& s = out.sol;
    s.Q = Matrix<T>(M, K, zero);
    s.U = Matrix<T>(M, K, zero);
    s.R = Matrix<T>(M, K, zero);
    s.Tp = Matrix<T>(M, K, zero);
    s.C.assign(K, zero);
    s.X.assign(K, zero);
    s.iter = 1;  // LDQBD is a direct method

    if (isOpen) {
        // Served throughput = arrival rate less the truncation blocking.
        const T X = T(lambda_eff * T(one - p[Nlev]));
        const T Rq = (X > zero) ? T(mean_queue / X) : zero;
        s.Tp(srcIdx - 1, 0) = X;
        s.Q(queueIdx - 1, 0) = mean_queue;
        s.U(queueIdx - 1, 0) = util;
        s.R(queueIdx - 1, 0) = Rq;
        s.Tp(queueIdx - 1, 0) = X;
        s.X[0] = X;
        s.C[0] = Rq;
    } else {
        const T mean_delay = T(num_traits<T>::from_double(L.classes[0].population) - mean_queue);
        const T X = T(mean_delay * lambda_eff);
        const T Rq = (X > zero) ? T(mean_queue / X) : zero;
        const T Rd = T(one / delayRate);
        // The delay completes at mean_delay * lambda_d, of which only the
        // fraction rt(delay, queue) proceeds to the queue, so its throughput is
        // NOT the queue flow X whenever the delay routes elsewhere.
        s.Q(delayIdx - 1, 0) = mean_delay;
        s.U(delayIdx - 1, 0) = mean_delay;  // infinite server: U = Q
        s.R(delayIdx - 1, 0) = Rd;
        s.Tp(delayIdx - 1, 0) = T(mean_delay * delayRate);
        s.Q(queueIdx - 1, 0) = mean_queue;
        s.U(queueIdx - 1, 0) = util;
        s.R(queueIdx - 1, 0) = Rq;
        s.Tp(queueIdx - 1, 0) = X;
        s.X[0] = X;
        s.C[0] = T(Rd + Rq);
    }

    LdqbdBlocks<T>& b = out.ld;
    b.Q0 = Q0;
    b.Q1 = Q1;
    b.Q2 = Q2;
    b.Nlev = Nlev;
    b.nPhases = nPhases;
    b.isPH = isPH;
    b.isOpen = isOpen;
    b.queueIdx = queueIdx;
    b.refIdx = isOpen ? srcIdx : delayIdx;
    b.M = M;
    b.nServers = nServers;
    b.mean_service = mean_service;
    b.hasLLD = hasLLD;
    b.sf = sf;
    b.utilPeak = utilPeak;
    b.lambda_eff = lambda_eff;
    b.hasSetup = hasSetup;
    b.alpharate = alpharate;
    b.alphascv = alphascv;
    b.betarate = betarate;
    b.betascv = betascv;
    b.delayRate = isOpen ? zero : delayRate;
    b.N = isOpen ? std::numeric_limits<double>::infinity() : L.classes[0].population;
    return out;
    }  // if constexpr has_transcendental
}

}  // namespace mam
}  // namespace line

#endif  // LINE_SOLVERS_MAM_SOLVER_MAM_LDQBD_H
