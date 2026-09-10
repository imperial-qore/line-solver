/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_NC_SOLVER_NC_SDR_H
#define LINE_SOLVERS_NC_SOLVER_NC_SDR_H

/**
 * Exact product-form analysis under state-dependent routing.
 *
 * Port of matlab/src/solvers/NC/solver_nc_sdr_analyzer.m, from A. E.
 * Krzesinski, "Multiclass Queueing Networks with State-Dependent Routing",
 * Performance Evaluation 7(2):125-143, 1987. The joint distribution is
 * eq. (16); the coefficients xi are those of Section 3.2, taken from the
 * state-independent part of the routing matrix.
 *
 * This is the ONLY exact route for an SDR model in this port. The generator
 * carries no per-state routing function, so the CTMC and SSA handlers cannot
 * represent eq. (10) and do not declare the feature; `refresh_routing` gives
 * the entry row the same uniform placeholder MATLAB's `getRoutingMatrix` does,
 * which this analyzer never reads because `pfqn_sdrvisits` overwrites it with
 * the collapsed e -> d arc.
 */

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "line/api/pfqn/pfqn_sdr.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/solvers/nc/nc_types.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace nc {

/**
 * Solves a network whose entry centre routes by state-dependent routing.
 *
 * @param sn  the refreshed struct, carrying a non-empty `sdr`
 * @param opt solver controls, unused beyond the method label
 */
template <class T>
NcSolution<T> solver_nc_sdr(const qn::NetworkStruct<T>& sn, const NcSolverOptions& opt) {
    const std::size_t M = sn.nstations, K = sn.nclasses;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);

    const std::vector<double> njobs = sn.njobs();
    for (std::size_t r = 0; r < K; ++r)
        if (std::isinf(njobs[r]))
            throw UnsupportedError("solver_nc_sdr: state-dependent routing is defined for closed "
                                   "networks only; the model has an open class");
    if (sn.nchains != K)
        throw UnsupportedError("solver_nc_sdr: state-dependent routing does not support class "
                               "switching; the product form of Krzesinski (1987) is stated over "
                               "closed chains whose customers keep their class");
    if (sn.stateful_nodes.size() != M)
        throw UnsupportedError("solver_nc_sdr: state-dependent routing requires every stateful "
                               "node to be a station; the product form is over queue lengths, and "
                               "a stateless node holds none");

    Matrix<T> S(M, K, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < K; ++r) {
            const T mu = sn.rates(i, r);
            if (mu > zero) S(i, r) = one / mu;
        }

    std::size_t Ntot = 0;
    std::vector<std::size_t> N(K, 0);
    for (std::size_t r = 0; r < K; ++r) {
        N[r] = static_cast<std::size_t>(njobs[r]);
        Ntot += N[r];
    }

    Matrix<T> alpha(M, Ntot > 0 ? Ntot : 1, one);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 1; k <= alpha.cols(); ++k) {
            T a = one;
            if (sn.stations[i].sched == qn::SchedStrategy::INF) {
                a = num_traits<T>::from_int(static_cast<long>(k));
            } else {
                const double c = sn.stations[i].nservers;
                if (std::isfinite(c) && c > 1.0)
                    a = num_traits<T>::from_int(
                        static_cast<long>(std::min<double>(static_cast<double>(k), c)));
            }
            alpha(i, k - 1) = a;
        }

    // A BCMP centre served FCFS must hold one rate for every chain; eq. (16)
    // admits chain-dependent rates only at the symmetric disciplines.
    for (std::size_t i = 0; i < M; ++i)
        if (sn.stations[i].sched == qn::SchedStrategy::FCFS) {
            bool haveRef = false;
            T ref = zero;
            for (std::size_t r = 0; r < K; ++r) {
                if (njobs[r] <= 0.0 || !(S(i, r) > zero)) continue;
                if (!haveRef) {
                    ref = S(i, r);
                    haveRef = true;
                } else if (!(ref == S(i, r))) {
                    throw UnsupportedError("solver_nc_sdr: station '" + sn.stations[i].name +
                                           "' is FCFS with chain-dependent service times, which "
                                           "has no BCMP product form; use PS, LCFSPR or INF, or "
                                           "equalize the service times");
                }
            }
        }

    std::vector<Matrix<T>> P;
    P.reserve(K);
    for (std::size_t r = 0; r < K; ++r) {
        Matrix<T> Pr(M, M, zero);
        for (std::size_t i = 0; i < M; ++i) {
            const std::size_t isf = sn.stateful_of_station(i + 1) - 1;
            for (std::size_t j = 0; j < M; ++j) {
                const std::size_t jsf = sn.stateful_of_station(j + 1) - 1;
                Pr(i, j) = sn.rt(isf * K + r, jsf * K + r);
            }
        }
        P.push_back(Pr);
    }
    const Matrix<T> xi = pfqn::pfqn_sdrvisits(sn.sdr, P);
    // 'sdr' evaluates the product form (16) exactly by state enumeration, which
    // is general in the branch topology; 'sdr.mva' runs the paper's Section 4
    // MVA and convolution, which costs O(J T M (V_1...V_J)^2) instead of the
    // state-space size but requires single-centre branches. Both are exact.
    const bool useMva = (opt.method == "sdr.mva");
    const pfqn::SdrResult<T> pf = useMva ? pfqn::pfqn_sdrmva(S, xi, N, sn.sdr, alpha)
                                         : pfqn::pfqn_sdr(S, xi, N, sn.sdr, alpha);

    NcSolution<T> out;
    out.sol.Q = pf.QN;
    out.sol.U = pf.UN;
    out.sol.R = pf.RN;
    out.sol.Tp = pf.XN;
    out.sol.X.assign(K, zero);
    out.sol.C.assign(K, zero);
    for (std::size_t r = 0; r < K; ++r) {
        const std::size_t ref = sn.classes[r].refstat - 1;
        out.sol.X[r] = pf.XN(ref, r);
        if (out.sol.X[r] > zero)
            out.sol.C[r] = num_traits<T>::from_double(njobs[r]) / out.sol.X[r];
    }
    out.sol.method = useMva ? "sdr.mva" : "sdr";
    out.sol.iter = 1;
    out.actualmethod = out.sol.method;
    out.STeff = S;
    return out;
}

}  // namespace nc
}  // namespace line

#endif  // LINE_SOLVERS_NC_SOLVER_NC_SDR_H
