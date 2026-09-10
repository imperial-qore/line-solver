/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_NC_SN_PF_PARAMS_H
#define LINE_SOLVERS_NC_SN_PF_PARAMS_H

/**
 * Port of `matlab/src/api/sn/sn_get_product_form_params.m`: the CLASS-level
 * product-form parameters.
 *
 * This is the sibling of `sn_get_product_form_chain_params` in
 * `solvers/mva/sn_chain.h`, and the difference is the whole point of having
 * both. The chain version aggregates classes into chains and is what the
 * mean-value solvers want; this one keeps the classes apart, which is what the
 * SOJOURN-TIME distribution needs -- `pfqn_stdf` conditions on a tagged job of
 * a named class, and a chain-aggregated demand has no such job in it.
 *
 * The demand is normalized by the visits of the chain's REFERENCE CLASS at the
 * reference station when the chain has one, and left unnormalized when it does
 * not, exactly as the reference does. That divisor is what makes D a demand per
 * system completion rather than per visit.
 *
 * `mu` is wider than the population on purpose: the reference sizes it
 * `ceil(sum N) + max(S)` because `pfqn_mvaldmx` indexes past |N| by the server
 * count. Narrowing it to |N| is a silent out-of-range read there.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/util/matrix.h"

namespace line {
namespace nc {

/** The `[lambda,D,N,Z,mu,S,V]` of the reference. */
template <class T>
struct PfParams {
    std::vector<T> lambda;   ///< (R) arrival rate, zero on a closed class
    Matrix<T> D;             ///< (Mq x R) demand at the queueing stations
    std::vector<double> N;   ///< (R) population, infinite on an open class
    Matrix<T> Z;             ///< (max(1,Mz) x R) demand at the delay stations
    Matrix<T> mu;            ///< (Mq x ceil(sum N)+max S) rate lattice
    std::vector<double> S;   ///< (Mq) server counts
    Matrix<T> V;             ///< (M x R) visits summed over chains
    std::vector<std::size_t> queue_stations;  ///< (Mq) 1-based station indices
    std::vector<std::size_t> delay_stations;  ///< (Mz) 1-based station indices
};

/**
 * Port of `sn_get_product_form_params`.
 *
 * @param sn the refreshed struct
 */
template <class T>
PfParams<T> sn_get_product_form_params(const qn::NetworkStruct<T>& sn) {
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const std::size_t R = sn.nclasses;

    PfParams<T> p;
    p.N.assign(R, 0.0);
    for (std::size_t r = 0; r < R; ++r) p.N[r] = sn.classes[r].population;

    // The reference selects on NODE type, not on the scheduling discipline: a
    // Delay is an infinite server and a Queue is not, whatever its sched.
    std::size_t sourceStation = 0;
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        const qn::NodeType nt = sn.stations[i].nodetype;
        if (nt == qn::NodeType::Queue) p.queue_stations.push_back(i + 1);
        else if (nt == qn::NodeType::Delay) p.delay_stations.push_back(i + 1);
        else if (nt == qn::NodeType::Source) sourceStation = i + 1;
    }
    const std::size_t Mq = p.queue_stations.size(), Mz = p.delay_stations.size();

    p.lambda.assign(R, zero);
    for (std::size_t r = 0; r < R; ++r)
        if (std::isinf(p.N[r]) && sourceStation > 0 && !sn.disabled[sourceStation - 1][r])
            p.lambda[r] = sn.rates(sourceStation - 1, r);

    p.S.assign(Mq, 1.0);
    double maxS = 1.0;
    for (std::size_t i = 0; i < Mq; ++i) {
        p.S[i] = sn.stations[p.queue_stations[i] - 1].nservers;
        if (std::isfinite(p.S[i]) && p.S[i] > maxS) maxS = p.S[i];
    }

    double Nct = 0.0;
    for (double n : p.N)
        if (std::isfinite(n)) Nct += n;
    const std::size_t width =
        static_cast<std::size_t>(std::ceil(Nct)) + static_cast<std::size_t>(maxS);

    // The per-class demand: visits at this station over the service rate, over
    // the reference class's visits at the reference station.
    const auto demand = [&](std::size_t ist, std::size_t r) {
        std::size_t c = sn.nchains;
        for (std::size_t cc = 0; cc < sn.nchains; ++cc)
            if (sn.chains[cc][r]) { c = cc; break; }
        if (c == sn.nchains || sn.disabled[ist - 1][r] || sn.rates(ist - 1, r) == zero)
            return zero;  // MATLAB divides by NaN here and clears with D(isnan(D))=0
        const std::size_t sf = sn.stateful_of_station(ist) - 1;
        T num = T(sn.visits[c](sf, r) / sn.rates(ist - 1, r));
        if (sn.refclass[c] > 0) {
            const std::size_t rsf = sn.stateful_of_station(sn.classes[r].refstat) - 1;
            const T den = sn.visits[c](rsf, sn.refclass[c] - 1);
            if (den == zero) return zero;
            num = T(num / den);
        }
        return num;
    };

    p.D = Matrix<T>(Mq, R, zero);
    p.mu = Matrix<T>(Mq, std::max<std::size_t>(1, width), one);
    for (std::size_t i = 0; i < Mq; ++i) {
        for (std::size_t r = 0; r < R; ++r) p.D(i, r) = demand(p.queue_stations[i], r);
        for (std::size_t n = 1; n <= p.mu.cols(); ++n)
            p.mu(i, n - 1) = num_traits<T>::from_double(
                std::min<double>(static_cast<double>(n), p.S[i]));
    }

    p.Z = Matrix<T>(std::max<std::size_t>(1, Mz), R, zero);
    for (std::size_t i = 0; i < Mz; ++i)
        for (std::size_t r = 0; r < R; ++r) p.Z(i, r) = demand(p.delay_stations[i], r);

    p.V = Matrix<T>(sn.nstations, R, zero);
    for (std::size_t c = 0; c < sn.nchains; ++c)
        for (std::size_t i = 0; i < sn.nstations; ++i) {
            const std::size_t sf = sn.stateful_of_station(i + 1) - 1;
            for (std::size_t r = 0; r < R; ++r) p.V(i, r) = T(p.V(i, r) + sn.visits[c](sf, r));
        }
    return p;
}

}  // namespace nc
}  // namespace line

#endif  // LINE_SOLVERS_NC_SN_PF_PARAMS_H
