/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_MAM_SOLVER_MAM_LDQBD_FLATTEN_H
#define LINE_SOLVERS_MAM_SOLVER_MAM_LDQBD_FLATTEN_H

/**
 * The two reductions that let an LD-QBD stand in for an enumerated CTMC.
 *
 * Ports of `solver_mam_ldqbd_flatten.m` and `solver_mam_ldqbd_avg.m`. They
 * exist for one caller, the SolverENV state-vector analyzer, and the shape of
 * that caller is what explains them: it propagates a DISTRIBUTION over a stage's
 * state space through the stage's sojourn, so it needs the stage as a flat
 * generator it can exponentiate, and it needs to reduce an ARBITRARY vector over
 * that space to means -- not the stationary one, which `solver_mam_ldqbd` alone
 * would give it.
 *
 * WHY FLATTENING IS NOT A LOSS HERE. `solver_mam_ldqbd` produces the chain as
 * block-tridiagonal (Q0 up, Q1 within, Q2 down) and solves it by a level-by-
 * level recursion that never forms the whole matrix -- which is the point of a
 * QBD. The state-vector analyzer cannot use that recursion: it does not want a
 * stationary vector, it wants `exp(Q t)` applied to a vector it brings with it.
 * The blocks are finite here (`Nlev` is the closed population, or the open
 * truncation from `cutoff`), so the flat matrix EXISTS; it is O(Nlev^2 nPhases^2)
 * to store where the blocks are O(Nlev nPhases^2), and that is the cost of the
 * question being asked.
 *
 * LEVEL 0 IS ONE STATE AND THE REST ARE nPhases WIDE, so the offsets are not a
 * multiple of the level index. `levelOf` carries the level of each flat state
 * rather than making the caller recompute it, because getting that mapping wrong
 * misreports the queue length without making the generator invalid.
 */

#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/solvers/mam/solver_mam_ldqbd.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/** The flat generator of an LD-QBD, with the level each flat state belongs to. */
template <class T>
struct LdqbdFlat {
    Matrix<T> Q;                       ///< the dense generator over level/phase states
    std::vector<std::size_t> levelOf;  ///< levelOf[s] = the queue level of flat state s
};

/** Port of `solver_mam_ldqbd_flatten.m`. */
template <class T>
LdqbdFlat<T> solver_mam_ldqbd_flatten(const LdqbdBlocks<T>& ld) {
    const std::size_t Nlev = ld.Nlev;
    if (ld.Q1.size() != Nlev + 1)
        throw InputError(
            "solver_mam_ldqbd_flatten: the block set carries " + std::to_string(ld.Q1.size()) +
            " within-level blocks for " + std::to_string(Nlev + 1) + " levels");

    std::vector<std::size_t> levelSize(Nlev + 1, 0), levelStart(Nlev + 2, 0);
    for (std::size_t n = 0; n <= Nlev; ++n) {
        levelSize[n] = ld.Q1[n].rows();
        levelStart[n + 1] = levelStart[n] + levelSize[n];
    }
    const std::size_t dim = levelStart[Nlev + 1];

    LdqbdFlat<T> out;
    out.Q = Matrix<T>(dim, dim, num_traits<T>::from_int(0));
    out.levelOf.assign(dim, 0);
    for (std::size_t n = 0; n <= Nlev; ++n) {
        const std::size_t r0 = levelStart[n], rn = levelSize[n];
        for (std::size_t a = 0; a < rn; ++a) out.levelOf[r0 + a] = n;
        for (std::size_t a = 0; a < rn; ++a)
            for (std::size_t b = 0; b < rn; ++b) out.Q(r0 + a, r0 + b) = ld.Q1[n](a, b);
        if (n < Nlev) {
            const std::size_t c0 = levelStart[n + 1], cn = levelSize[n + 1];
            for (std::size_t a = 0; a < rn; ++a)
                for (std::size_t b = 0; b < cn && b < ld.Q0[n].cols(); ++b)
                    out.Q(r0 + a, c0 + b) = ld.Q0[n](a, b);
        }
        if (n >= 1) {
            // Q2 IS INDEXED BY LEVEL HERE, NOT BY THE MATLAB CELL POSITION.
            // `solver_mam_ldqbd.m` writes `Q2{n}` for the level n -> n-1 block,
            // 1-based, so the flatten there reads `Q2{n}`; this port keeps an
            // unused placeholder at index 0 so that Q0, Q1 and Q2 all line up by
            // level, which makes the same block `Q2[n]`. Reading `Q2[n-1]`
            // compiles, keeps the generator valid and conserves the population
            // -- it just shifts every departure rate down one level, which on a
            // three-job repairman model moved the mean queue from 1.42105 to
            // 1.98824 with the total still exactly 3.
            const std::size_t c0 = levelStart[n - 1], cn = levelSize[n - 1];
            for (std::size_t a = 0; a < rn && a < ld.Q2[n].rows(); ++a)
                for (std::size_t b = 0; b < cn && b < ld.Q2[n].cols(); ++b)
                    out.Q(r0 + a, c0 + b) = ld.Q2[n](a, b);
        }
    }
    return out;
}

/** Per-(station,class) means read off an arbitrary distribution over the LD-QBD. */
template <class T>
struct LdqbdAvg {
    Matrix<T> QN, UN, RN, TN;  ///< (M x 1), the model being single-class by construction
};

/**
 * Port of `solver_mam_ldqbd_avg.m`: map a distribution over the flat state
 * space to means.
 *
 * NOT THE STATIONARY DISTRIBUTION. This mirrors the metric formulas of
 * `solver_mam_ldqbd` but applies them to whatever vector it is handed -- in
 * practice the TIME-AVERAGE over an environment stage's sojourn, which is not a
 * stationary law of anything. That is why the formulas are written out again
 * here rather than shared: the stationary versions in `solver_mam_ldqbd` reach
 * for quantities (the level recursion's own R matrices) that only exist at the
 * fixed point.
 *
 * The vector is clipped at zero and renormalized first, as the reference does:
 * a transient vector that a quadrature has pushed a hair negative is a numerical
 * artifact of the propagation, not a signed measure to be propagated further.
 */
template <class T>
LdqbdAvg<T> solver_mam_ldqbd_avg(const LdqbdBlocks<T>& ld, const std::vector<T>& piflat_in,
                                 const std::vector<std::size_t>& levelOf) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t Nlev = ld.Nlev, M = ld.M;
    if (ld.queueIdx == 0 || ld.refIdx == 0 || M == 0)
        throw InputError("solver_mam_ldqbd_avg: the block set names no queue or reference station");
    if (piflat_in.size() != levelOf.size())
        throw InputError(
            "solver_mam_ldqbd_avg: the distribution and the level map disagree on the size of the "
            "state space");

    std::vector<T> piflat = piflat_in;
    T total = zero;
    for (std::size_t s = 0; s < piflat.size(); ++s) {
        if (num_traits<T>::to_double(piflat[s]) < 0.0) piflat[s] = zero;
        total = T(total + piflat[s]);
    }
    if (num_traits<T>::to_double(total) > 0.0)
        for (std::size_t s = 0; s < piflat.size(); ++s) piflat[s] = T(piflat[s] / total);

    std::vector<T> pLevel(Nlev + 1, zero);
    for (std::size_t s = 0; s < piflat.size(); ++s)
        if (levelOf[s] <= Nlev) pLevel[levelOf[s]] = T(pLevel[levelOf[s]] + piflat[s]);

    T mean_queue = zero;
    for (std::size_t n = 0; n <= Nlev; ++n)
        mean_queue = T(mean_queue + T(num_traits<T>::from_int(static_cast<long long>(n)) * pLevel[n]));

    // Utilization is the fraction of the station's PEAK capacity in use,
    // sum_n p(n)*sf(n)/utilPeak, the work-based convention CTMC, MVA, NC and
    // serial SSA all report. Without load dependence sf(n) = min(n,c) and
    // utilPeak = c, giving the mean fraction of the c servers in use; at c = 1
    // that is sf(n) = 1 for every n >= 1, so the sum collapses to 1 - p(0).
    // The `/utilPeak` is inside the sum, so the branches below must NOT divide
    // again.
    T util = zero;
    {
        const T peak = num_traits<T>::from_double(ld.utilPeak);
        for (std::size_t n = 1; n <= Nlev && n < ld.sf.size(); ++n)
            util = T(util + T(ld.sf[n] / peak * pLevel[n]));
    }

    LdqbdAvg<T> out;
    out.QN = Matrix<T>(M, 1, zero);
    out.UN = Matrix<T>(M, 1, zero);
    out.RN = Matrix<T>(M, 1, zero);
    out.TN = Matrix<T>(M, 1, zero);
    const std::size_t qi = ld.queueIdx - 1, ri = ld.refIdx - 1;

    if (ld.isOpen) {
        // The throughput is the arrival rate less the share LOST at the
        // truncation level, which is what makes the open answer depend on
        // `cutoff` rather than silently ignoring the loss.
        const T X = T(ld.lambda_eff * T(one - pLevel[Nlev]));
        const T Rq = num_traits<T>::to_double(X) > 0.0 ? T(mean_queue / X) : zero;
        out.TN(ri, 0) = X;
        out.QN(qi, 0) = mean_queue;
        out.UN(qi, 0) = util;
        out.RN(qi, 0) = Rq;
        out.TN(qi, 0) = X;
    } else {
        const T mean_delay = T(num_traits<T>::from_double(ld.N) - mean_queue);
        const T X = T(mean_delay * ld.lambda_eff);
        const T Rq = num_traits<T>::to_double(X) > 0.0 ? T(mean_queue / X) : zero;
        out.QN(ri, 0) = mean_delay;
        // A Delay's utilization IS its mean population: it has one server per
        // job, so "fraction busy" has no other meaning there.
        out.UN(ri, 0) = mean_delay;
        out.RN(ri, 0) = num_traits<T>::to_double(ld.delayRate) > 0.0 ? T(one / ld.delayRate) : zero;
        out.TN(ri, 0) = X;
        out.QN(qi, 0) = mean_queue;
        out.UN(qi, 0) = util;
        out.RN(qi, 0) = Rq;
        out.TN(qi, 0) = X;
    }
    return out;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_SOLVERS_MAM_SOLVER_MAM_LDQBD_FLATTEN_H
