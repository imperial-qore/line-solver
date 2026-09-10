/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_NC_DPS_H
#define LINE_SOLVERS_NC_DPS_H

/**
 * Heavy-usage asymptotic analysis of the closed two-station network with one think
 * (infinite-server) station and one discriminatory processor-sharing station, by the
 * generating-function expansion of J.A. Morrison, "Asymptotic analysis of a large closed queueing
 * network with discriminatory processor sharing", Queueing Systems 9 (1991) 191-214.
 *
 * Port of matlab/src/solvers/NC/nc_is_dps_model.m and solver_nc_dps_analyzer.m, matching
 * jar/src/main/java/jline/solvers/nc/analyzers/Solver_nc_dps_analyzer.java and
 * python/line_solver/solvers/solver_nc/solver_nc_dps_analyzer.py. The kernel is
 * line/api/npfqn/npfqn_dps_morrison.h.
 *
 * THERE IS NO NORMALIZING CONSTANT HERE. A DPS station is not product-form -- that is the premise
 * of the paper -- so lG is NaN, as on the maximum-entropy route. NC hosts this method because NC is
 * where LINE keeps the asymptotic expansions of generating functions and normalizing-constant
 * integrals (pana, mmint2, le, ble, gleint, rayint), which is the family Morrison's expansion
 * belongs to, not because a constant is being computed.
 *
 * Response times come from Little's law on the queue-length result rather than from the expanded
 * RESULT 2 (eq. 4.17), so that Q = R*T holds exactly in the returned table; the two agree to the
 * order of the approximation, since Morrison derives (4.17) as the ratio (4.11)/(4.15).
 */

#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "line/api/npfqn/npfqn_dps_morrison.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/solvers/nc/nc_types.h"
#include "line/solvers/nc/solver_nc_oi.h"  // detail::oi_visits, shared with the OI route
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace nc {

/** True when any station of the network is scheduled DPS. */
template <class T>
bool sn_has_dps(const qn::NetworkStruct<T>& sn) {
    for (std::size_t i = 0; i < sn.nstations; ++i)
        if (sn.stations[i].sched == qn::SchedStrategy::DPS) return true;
    return false;
}

/**
 * True when the model is the closed two-station network Morrison's expansion is derived for: one
 * infinite-server (think) station and one single-server DPS station, exponential service, every
 * class alternating between the two. The shape is checked exactly, not approximately: outside it
 * the expansion has no derivation behind it.
 */
template <class T>
bool nc_is_dps_model(const qn::NetworkStruct<T>& sn) {
    using qn::SchedStrategy;
    if (sn.nstations != 2) return false;
    bool anyPositive = false;
    for (const qn::JobClass& c : sn.classes) {
        if (std::isinf(c.population)) return false;
        if (c.population > 0) anyPositive = true;
    }
    if (!anyPositive) return false;

    std::size_t iInf = sn.nstations, iDps = sn.nstations;
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        const SchedStrategy s = sn.stations[i].sched;
        if (s == SchedStrategy::INF) {
            if (iInf != sn.nstations) return false;
            iInf = i;
        } else if (s == SchedStrategy::DPS) {
            if (iDps != sn.nstations) return false;
            iDps = i;
        } else {
            return false;
        }
    }
    if (iInf == sn.nstations || iDps == sn.nstations) return false;
    const double c = sn.stations[iDps].nservers;
    if (std::isfinite(c) && c != 1.0) return false;  // multi-server DPS is not Morrison's share
    if (sn.nchains != sn.nclasses) return false;     // class switching
    for (std::size_t ch = 0; ch < sn.nchains; ++ch)
        if (sn.inchain[ch].size() > 1) return false;

    for (std::size_t r = 0; r < sn.nclasses; ++r) {
        if (!(sn.classes[r].population > 0)) return false;
        const std::size_t sts[2] = {iInf, iDps};
        for (std::size_t k = 0; k < 2; ++k) {
            const double rate = num_traits<T>::to_double(sn.rates(sts[k], r));
            if (!std::isfinite(rate) || rate <= 0) return false;
            const double scv = num_traits<T>::to_double(sn.scv(sts[k], r));
            if (std::isfinite(scv) && std::fabs(scv - 1.0) > 1e-6) return false;
        }
        if (sn.stations[iDps].schedparam.size() != sn.nclasses) return false;
        const double wgt = num_traits<T>::to_double(sn.stations[iDps].schedparam[r]);
        if (!std::isfinite(wgt) || wgt <= 0) return false;
    }

    const Matrix<T> V = detail::oi_visits(sn);
    for (std::size_t r = 0; r < sn.nclasses; ++r) {
        const double vi = num_traits<T>::to_double(V(iInf, r));
        const double vd = num_traits<T>::to_double(V(iDps, r));
        if (std::fabs(vi - vd) > 1e-9 * std::max(1.0, vi)) return false;  // unequal visits
    }
    return true;
}

/**
 * Analyzes the closed think+DPS network.
 *
 * @param sn the network structure, of the shape nc_is_dps_model accepts
 * @param opt solver options
 * @return the mean performance measures, with lG = NaN
 */
template <class T>
NcSolution<T> solver_nc_dps_analyzer(const qn::NetworkStruct<T>& sn, const NcSolverOptions& opt) {
    (void)opt;
    NcSolution<T> out;
    if constexpr (!num_traits<T>::has_transcendental) {
        (void)sn;
        throw UnsupportedError(
            "solver_nc_dps_analyzer: Morrison's W_m are erfc integrals; this backend has no "
            "transcendental arithmetic");
    } else {
        using qn::SchedStrategy;
        // The shape is re-checked here, not assumed from the caller: this analyzer is
        // reachable from the runner, the dispatch chain and directly from user code, and
        // every quantity below -- think time, DPS service time, weights, visit ratios --
        // is meaningless off the shape the expansion was derived for.
        if (!nc_is_dps_model(sn))
            throw InputError(
                "solver_nc_dps_analyzer: applies only to a CLOSED network of exactly two stations, "
                "one infinite-server (think) station and one single-server DPS station with "
                "exponential service and one visit each per cycle (see nc_is_dps_model).");
        const T zero = num_traits<T>::from_int(0);
        const std::size_t M = sn.nstations, K = sn.nclasses;
        std::size_t iInf = M, iDps = M;
        for (std::size_t i = 0; i < M; ++i) {
            if (sn.stations[i].sched == SchedStrategy::INF) iInf = i;
            if (sn.stations[i].sched == SchedStrategy::DPS) iDps = i;
        }
        if (iInf == M || iDps == M)
            throw InputError("solver_nc_dps_analyzer: requires one INF and one DPS station");

        std::vector<T> Npop(K, zero), Z(K, zero), S(K, zero), w(K, zero);
        for (std::size_t r = 0; r < K; ++r) {
            Npop[r] = num_traits<T>::from_double(sn.classes[r].population);
            Z[r] = T(num_traits<T>::from_int(1) / sn.rates(iInf, r));
            S[r] = T(num_traits<T>::from_int(1) / sn.rates(iDps, r));
            w[r] = sn.stations[iDps].schedparam[r];
        }

        const npfqn::DpsMorrisonResult<T> mor = npfqn::npfqn_dps_morrison(Npop, Z, S, w);

        const Matrix<T> V = detail::oi_visits(sn);
        Matrix<T> Q(M, K, zero), U(M, K, zero), R(M, K, zero), Tp(M, K, zero);
        std::vector<T> X(K, zero), C(K, zero);

        const double cd = sn.stations[iDps].nservers;
        const T c = num_traits<T>::from_double((std::isfinite(cd) && cd > 0.0) ? cd : 1.0);

        for (std::size_t r = 0; r < K; ++r) {
            T q = mor.Q[r];
            if (!(q >= zero)) q = zero;             // outside the regime the expansion can leave
            if (q > Npop[r]) q = Npop[r];           // [0,N]; clamp rather than report an impossible
            const T x = T((Npop[r] - q) / Z[r]);
            Q(iDps, r) = q;
            Q(iInf, r) = T(Npop[r] - q);            // population conservation (exact, closed)
            X[r] = x;
            for (std::size_t i = 0; i < M; ++i) Tp(i, r) = T(x * V(i, r));
            U(iInf, r) = Q(iInf, r);                // INF utilization convention
            U(iDps, r) = T(x * V(iDps, r) * S[r] / c);
            if (x > zero) C[r] = T(Npop[r] / x);
        }
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < K; ++r)
                if (Tp(i, r) > zero) R(i, r) = T(Q(i, r) / Tp(i, r));

        out.sol.Q = Q;
        out.sol.U = U;
        out.sol.R = R;
        out.sol.Tp = Tp;
        out.sol.X = X;
        out.sol.C = C;
        out.sol.iter = 1;
        out.sol.lG = std::numeric_limits<double>::quiet_NaN();  // not product-form
        out.sol.method = "morrison";
        out.actualmethod = "morrison";
    }
    return out;
}

}  // namespace nc
}  // namespace line

#endif  // LINE_SOLVERS_NC_DPS_H
