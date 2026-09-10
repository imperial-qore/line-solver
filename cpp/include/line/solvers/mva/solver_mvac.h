/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_MVA_SOLVER_MVAC_H
#define LINE_SOLVERS_MVA_SOLVER_MVAC_H

/**
 * MVAC, exact mean value analysis BY CHAIN (Conway, de Souza e Silva and
 * Lavenberg, IEEE Trans. Computers 38(3):432-442, 1989). Port of
 * `solver_mvac.m`.
 *
 * MVAC IS EXACT for its model class, not an approximation: it returns what the
 * classic MVA recursion returns, by a different route. Where `pfqn_mva` recurs
 * on the population vector at a cost of prod(N+1), MVAC recurs on the CHAINS,
 * replacing each removed chain by self-looping single-customer chains, so its
 * cost is governed by the number of DISTINCT demand columns. That is the whole
 * reason to offer it beside `exact`, and it is why any disagreement with exact
 * MVA on a model both accept is a defect rather than a tolerance.
 *
 * ITS MODEL CLASS IS NARROWER THAN THE `exact` PATH's, and every restriction is
 * refused by name rather than approximated:
 *
 *   closed only          an open chain has no population to recur on
 *   product form only    the recursion is the BCMP one
 *   single server only   the api implements the SSFR arrival theorem, and has
 *                        no multiserver correction; `pfqn_mvacld` is the
 *                        load-dependent sibling, and the reference does NOT
 *                        route here to it -- it refuses, so that a multiserver
 *                        model takes the `exact` path it is already solved by
 *
 * TWO CHAIN-LEVEL QUANTITIES ARE COMPUTED AND THEN DISCARDED by the reference,
 * which passes an empty Qchain and Uchain to the deaggregation. The class-level
 * queue length is therefore reconstructed from the residence time by Little's
 * law and the utilization from the utilization law, NOT split from the chain
 * matrices by the visit share. Handing the chain matrices over instead would
 * change the numbers on any chain holding more than one class, which is why the
 * reference's utilization renormalization (its `sum(Uchain) > 1` rescaling) is
 * not reproduced here: it cannot reach the answer, and porting dead arithmetic
 * would suggest it does. The cycle time is discarded the same way, and survives
 * only as the divisor that turns the population into the throughput.
 *
 * Arithmetic: EXACT-CAPABLE, and therefore ungated. `pfqn_mvac` forms no
 * normalizing constant and evaluates no logarithm -- additions, multiplications
 * and divisions in the field of the inputs only -- so this analyzer runs under
 * exact arithmetic and returns exact fractions there. Contrast `solver_sqd`,
 * whose calibration is transcendental and is gated for that reason. `lG` is
 * NaN: MVAC forms no normalizing constant, and the reference says so rather
 * than reporting a zero that reads as one.
 */

#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "line/api/pfqn/pfqn_mvac.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/mva/mva_types.h"
#include "line/solvers/mva/sn_chain.h"
#include "line/util/error.h"

namespace line {
namespace mva {

/** Port of `solver_mvac.m`. */
template <class T>
MvaSolution<T> solver_mvac_analyzer(const qn::NetworkStruct<T>& L, const MvaOptions& opt) {
    using qn::SchedStrategy;
    // `options.tol` is read by the reference only inside the utilization
    // renormalization the deaggregation discards; nothing else here is tuned.
    (void)opt;
    const T zero = num_traits<T>::from_int(0);
    const ChainDemands<T> d = sn_get_demands_chain(L);
    const std::size_t M = L.nstations, K = L.nchains;

    // One predicate for the gate and the run: list_valid_methods asks
    // mva_mvac_reason before it advertises 'mvac', so a listed row is a row that
    // runs and the refusal reads the same either way.
    {
        const std::string mvac_reason = mva_mvac_reason(L, "mvac");
        if (!mvac_reason.empty()) throw UnsupportedError(mvac_reason);
    }
    for (std::size_t c = 0; c < K; ++c) {
        if (!std::isfinite(d.Nchain[c]))
            throw UnsupportedError(
                "solver_mvac_analyzer: MVAC supports closed models only; use method 'exact' for "
                "open or mixed networks");
        if (d.Nchain[c] != std::floor(d.Nchain[c]))
            throw UnsupportedError(
                "solver_mvac_analyzer: the MVAC recursion removes customers one at a time and has "
                "no fractional-population form");
    }

    std::vector<std::size_t> infSET, qSET;  // 0-based station indices
    for (std::size_t i = 0; i < M; ++i) {
        switch (L.stations[i].sched) {
            case SchedStrategy::EXT: break;  // no external world in a closed model
            case SchedStrategy::INF: infSET.push_back(i); break;
            case SchedStrategy::PS:
            case SchedStrategy::LCFSPR:
            case SchedStrategy::FCFS:
            case SchedStrategy::SIRO: qSET.push_back(i); break;
            default:
                throw UnsupportedError(std::string("solver_mvac_analyzer: MVAC does not support ") +
                                       lang::sched_to_text(L.stations[i].sched) + " scheduling");
        }
    }

    // A chain with no jobs carries no work. The reference solves on the active
    // set and re-expands, because an empty column would make the api recur on a
    // chain that has no customer to remove.
    std::vector<std::size_t> rset;
    for (std::size_t c = 0; c < K; ++c)
        if (d.Nchain[c] != 0.0) rset.push_back(c);

    std::vector<T> Xchain(K, zero);
    Matrix<T> Qchain(M, K, zero), Tchain(M, K, zero), Wchain(M, K, zero), Rchain(M, K, zero);
    std::vector<T> Ccycle(K, zero);

    if (!rset.empty()) {
        const std::size_t Mq = qSET.size(), Rr = rset.size();
        Matrix<T> Lq(Mq, Rr, zero), Zq(1, Rr, zero);
        std::vector<int> N(Rr, 0);
        for (std::size_t j = 0; j < Rr; ++j) {
            const std::size_t c = rset[j];
            N[j] = static_cast<int>(std::llround(d.Nchain[c]));
            for (std::size_t a = 0; a < Mq; ++a)
                Lq(a, j) = T(d.STchain(qSET[a], c) * d.Vchain(qSET[a], c));
            // every delay folds into one think time per chain
            for (std::size_t i : infSET) Zq(0, j) += T(d.STchain(i, c) * d.Vchain(i, c));
        }

        const pfqn::MvacResult<T> mr = pfqn::pfqn_mvac(Lq, N, Zq);
        for (std::size_t j = 0; j < Rr; ++j) {
            Xchain[rset[j]] = mr.X[j];
            for (std::size_t a = 0; a < Mq; ++a) Qchain(qSET[a], rset[j]) = mr.Q(a, j);
        }
        for (std::size_t i : infSET)
            for (std::size_t c = 0; c < K; ++c)
                Qchain(i, c) = T(Xchain[c] * d.STchain(i, c) * d.Vchain(i, c));

        for (std::size_t c : rset) {
            for (std::size_t i : infSET) Wchain(i, c) = d.STchain(i, c);
            // The isinf(nservers) branch the reference guards this with is
            // unreachable: the gate above already refused any queueing station
            // whose server count is not exactly one.
            for (std::size_t i : qSET) {
                if (d.Vchain(i, c) == zero || Xchain[c] == zero) continue;
                Wchain(i, c) = T(Qchain(i, c) / (Xchain[c] * d.Vchain(i, c)));
            }
        }

        // The throughput is re-derived from the residence times rather than
        // taken from the api, which is what makes the reported Q, T and X
        // consistent with one another at the model level.
        for (std::size_t c : rset) {
            T wsum = zero;
            for (std::size_t i = 0; i < M; ++i) wsum += Wchain(i, c);
            if (wsum == zero) {
                Xchain[c] = zero;
            } else {
                for (std::size_t i = 0; i < M; ++i)
                    Ccycle[c] += T(d.Vchain(i, c) * Wchain(i, c));
                Xchain[c] = T(num_traits<T>::from_double(d.Nchain[c]) / Ccycle[c]);
            }
            for (std::size_t i = 0; i < M; ++i) {
                Qchain(i, c) = T(Xchain[c] * d.Vchain(i, c) * Wchain(i, c));
                Tchain(i, c) = T(Xchain[c] * d.Vchain(i, c));
            }
        }
    }

    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t c = 0; c < K; ++c)
            if (Tchain(i, c) > zero) Rchain(i, c) = T(Qchain(i, c) / Tchain(i, c));

    const ClassResults<T> cr =
        sn_deaggregate_chain_results(L, d, Matrix<T>(), Matrix<T>(), Rchain, Tchain, Xchain);
    MvaSolution<T> out;
    out.Q = cr.Q;
    out.U = cr.U;
    out.R = cr.R;
    out.Tp = cr.Tp;
    out.C = cr.C;
    out.X = cr.X;
    out.method = "mvac";
    out.iter = 0;
    out.lG = std::numeric_limits<double>::quiet_NaN();
    return out;
}

}  // namespace mva
}  // namespace line

#endif  // LINE_SOLVERS_MVA_SOLVER_MVAC_H
