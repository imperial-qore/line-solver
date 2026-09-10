/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_MAPDC_H
#define LINE_API_QSYS_QSYS_MAPDC_H

/**
 * The MAP/D/c FCFS queue: c servers, deterministic service of length s, fed by
 * a Markovian arrival process. The multiserver generalization of qsys_mapd1.
 *
 * ALGORITHM, AND HOW IT DIFFERS FROM THE MATLAB REFERENCE.
 * matlab/src/api/qsys/qsys_mapdc.m calls Q-MAM's Q_CT_MAP_D_C, which is not
 * transcribed here. This port computes the same quantities by Crommelin's
 * embedded lattice chain, which is exact for deterministic service and is NOT
 * an Erlang-k or a heavy-traffic approximation.
 *
 * THE LATTICE CHAIN IS EXACT, AND WHY. Sample the system at the epochs
 * t_n = n s and let N_n be the number in system at t_n, J_n the arrival phase.
 * Every job in service at t_n started somewhere in (t_n - s, t_n], so it
 * departs inside (t_n, t_{n+1}]; every job that starts service inside that
 * interval departs after t_{n+1}. Hence exactly min(N_n, c) departures occur
 * per interval -- all of them when the system is below capacity, since then no
 * job waits -- and
 *
 *     N_{n+1} = max(N_n - c, 0) + A_n,
 *
 * with A_n the number of MAP arrivals in the interval. (N_n, J_n) is therefore
 * a Markov chain with transition blocks A_k = P_k(s), the MAP counting
 * probabilities over an interval of exact length s. No residual service time
 * has to be carried, which is what makes the deterministic case tractable where
 * a general one would not be.
 *
 * The stationary law of that chain IS the time-stationary law of N. N(t)
 * converges in distribution as t -> infinity (the system is non-lattice: the
 * MAP has continuous interarrival times), so N(n s) converges to the same
 * limit, and a positive recurrent Markov chain converges only to its own
 * stationary law. No PASTA argument and no time-averaging step are involved --
 * unlike qsys_mapd1, whose chain is embedded at DEPARTURE epochs and therefore
 * needs the explicit inter-departure averaging of its step 4.
 *
 * FOUR STEPS.
 * 1. A_k = P_k(s), k = 0..K, as the top block row of exp(M s) with M block
 *    bidiagonal carrying D0 on the diagonal and D1 above it. Shared with
 *    qsys_mapd1 (line::qsys::detail::map_counting_at). K grows from a Poisson
 *    tail estimate until the missing mass max_i (1 - sum_k (A_k e)_i) is under
 *    the tolerance, so the truncation is a measured quantity.
 * 2. Grouping. The chain is skip-free to the left by c, not by 1, so c
 *    consecutive levels are grouped into one super-level: level n = L c + u
 *    becomes super-level L, sub-level u. Skip-freeness by one super-level then
 *    holds and the blocks are c m x c m with m the MAP order:
 *      repeating   A^(i)[(u,v)] = A_{i c + v - u},   for i c + v - u >= 0
 *      boundary    B^(i)[(u,v)] = A_{i c + v},       independent of u,
 *    the boundary row being the levels 0..c-1, from which max(n-c,0) = 0.
 * 3. Ramaswami's recursion, exactly as in qsys_mapd1: G the minimal solution of
 *    G = sum_i A^(i) G^i, the partial sums Ahat_i = sum_{k>=i} A^(k) G^(k-i)
 *    and Bhat_i likewise, then
 *      x_0 (B^(0) + Bhat_1 (I - Ahat_1)^-1 A^(0)) = x_0,
 *      x_L = [x_0 Bhat_L + sum_{k=1}^{L-1} x_k Ahat_{L-k+1}] (I - Ahat_1)^-1.
 * 4. Unfolding. Super-level L, sub-level u, phase j is the level L c + u, so
 *    P(N = n) is read off directly and E[N] = sum_n n P(N = n).
 *
 * meanWaitingTime is E[N]/lambda - s by Little's law and meanSojournTime is
 * E[N]/lambda. This is exactly the definition the MATLAB reference adopted when
 * its one-point quadrature was replaced; see the reference-defect note below.
 *
 * INVARIANTS THE TESTS ASSERT (both are identities of the construction, not
 * fitted quantities):
 *  - sum_n P(N = n) = 1 up to the level-tail truncation;
 *  - sum_n min(n, c) P(N = n) = lambda s. Taking expectations in the recursion
 *    gives E[min(N,c)] = E[A] = lambda s, i.e. the mean number of departures
 *    per interval equals the mean number of arrivals. Equivalently the
 *    utilization is rho = lambda s / c, which is what the result reports.
 *
 * ORACLES.
 *  - c = 1 collapses to MAP/D/1 and is cross-checked against qsys_mapd1, which
 *    reaches the same numbers by the departure-epoch chain, a genuinely
 *    different construction.
 *  - Poisson arrivals and c = 1 collapse to M/D/1, where Pollaczek-Khinchin
 *    gives L = rho + rho^2/(2(1-rho)) and Wq = rho s/(2(1-rho)) exactly.
 *  - Poisson arrivals and c > 1 are checked against the MATLAB reference and
 *    against the M/D/c heavy-traffic ordering Wq(M/D/c) < Wq(M/M/c).
 *
 * ARITHMETIC. Gated on num_traits<T>::has_transcendental for two independent
 * reasons, both inherited from qsys_mapd1: step 1 calls expm, which is a
 * tolerance-controlled Pade approximation and cannot be exact in any
 * arithmetic, and step 3 obtains G by a fixed-point iteration that does not
 * terminate in a finite number of field operations. Steps 2 and 4 are finite
 * exact matrix algebra and add no error of their own.
 *
 * REFERENCE DEFECTS.
 *  - FIXED IN MATLAB, and this port matches the fixed behaviour.
 *    qsys_mapdc.m used to return meanWaitingTime as a one-point left-rectangle
 *    quadrature of the Q-MAM waiting-time survival function (numSteps
 *    defaulting to 1). On the M/D/1 instance D0 = [-2], D1 = [2], s = 1/3 it
 *    returned 0.4444444443 against the exact Pollaczek-Khinchin value 1/3, an
 *    error of +33%. It now uses Little's law, meanQueueLength/lambda - s. This
 *    port reproduces the FIXED behaviour and is asserted against 1/3, never
 *    against the superseded number.
 *  - STILL OPEN. qsys_mapdc.m computes the utilization as rho = lambda s / c
 *    and returns the Q-MAM queue-length vector unchanged, so meanQueueLength is
 *    the mean number in SYSTEM; the Little's law line is consistent with that
 *    reading. No defect is claimed here, but note that the field is documented
 *    as "Mean number of customers in system" while the inline comment on the
 *    Q-MAM output says "ql(i) = Prob[(i-1) customers in the queue]". The two
 *    readings differ by rho c, and only the system reading makes
 *    meanWaitingTime non-negative at high load, so the system reading is the
 *    one implemented on both sides.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/api/mam/qbd_r.h"
#include "line/api/qsys/qsys_mapd1.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/expm.h"
#include "line/util/linalg.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace qsys {

/** Return value of qsys_mapdc, mirroring the MATLAB result struct. */
template <class T>
struct MapDcResult {
    T meanQueueLength;               ///< E[N], number in system, time-stationary
    T meanWaitingTime;               ///< Wq = E[N]/lambda - s
    T meanSojournTime;               ///< W = Wq + s = E[N]/lambda
    T utilization;                   ///< rho = lambda s / c, per server
    std::vector<T> queueLengthDist;  ///< P(N = n), n = 0, 1, ...
};

/**
 * MAP/D/c by Crommelin's exact embedded lattice chain.
 *
 * @param arrival      arrival MAP (D0, D1)
 * @param s            deterministic service time, s > 0
 * @param c            number of servers, c >= 1
 * @param dist_size    how many entries of queueLengthDist to materialize
 * @param max_arrivals cap on K, the number of arrivals per interval tracked
 * @param max_levels   cap on the number of super-levels generated
 * @param tol          tolerance on the counting-mass truncation, on the G
 *                     iteration and on the level-tail truncation
 */
template <class T>
MapDcResult<T> qsys_mapdc(const mam::Map<T>& arrival, const T& s, unsigned c,
                          std::size_t dist_size, unsigned max_arrivals, std::size_t max_levels,
                          const T& tol) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_mapdc requires transcendental arithmetic");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    if (s <= zero) throw InputError("qsys_mapdc: service time s must be positive");
    if (c == 0) throw InputError("qsys_mapdc: at least one server is required");
    if (dist_size == 0) throw InputError("qsys_mapdc: dist_size must be positive");
    if (max_levels == 0) throw InputError("qsys_mapdc: max_levels must be positive");
    const std::size_t m = arrival.D0.rows();
    if (arrival.D0.cols() != m || arrival.D1.rows() != m || arrival.D1.cols() != m)
        throw InputError("qsys_mapdc: D0 and D1 must be square and of equal order");

    const Matrix<T>& D0 = arrival.D0;
    const Matrix<T>& D1 = arrival.D1;
    const T lambda = mam::map_lambda(arrival);
    if (lambda <= zero) throw InputError("qsys_mapdc: non-positive arrival rate");
    const T cT = num_traits<T>::from_int(static_cast<long>(c));
    const T rho = lambda * s / cT;
    if (rho >= one) throw InputError("qsys_mapdc: load rho must be strictly less than 1");

    // ---- step 1: MAP counting probabilities over one interval of length s ----
    T numax = zero;
    for (std::size_t i = 0; i < m; ++i) {
        const T d = -D0(i, i);
        if (d > numax) numax = d;
    }
    const T est = numax * s;
    unsigned K = 8u;
    {
        const double e = num_traits<T>::to_double(est);
        const unsigned guess = static_cast<unsigned>(2.0 * e + 10.0 * std::sqrt(e + 1.0) + 8.0);
        if (guess > K) K = guess;
    }
    if (K > max_arrivals) K = max_arrivals;
    std::vector<Matrix<T>> A = detail::map_counting_at(D0, D1, s, K);
    while (true) {
        T missing = zero;
        for (std::size_t i = 0; i < m; ++i) {
            T row = zero;
            for (std::size_t k = 0; k < A.size(); ++k)
                for (std::size_t j = 0; j < m; ++j) row += A[k](i, j);
            const T lack = one - row;
            if (lack > missing) missing = lack;
        }
        if (missing <= tol || K >= max_arrivals) break;
        K = (2u * K < max_arrivals) ? 2u * K : max_arrivals;
        A = detail::map_counting_at(D0, D1, s, K);
    }
    const std::size_t Kmax = A.size() - 1;

    // ---- step 2: group c consecutive levels into one super-level ----
    // Ag[i][(u,j),(v,j')] = A_{i c + v - u}, Bg[i][(u,j),(v,j')] = A_{i c + v}.
    const std::size_t cs = static_cast<std::size_t>(c);
    const std::size_t nb = cs * m;  // super-level block order
    std::size_t Kg = (Kmax + cs - 1) / cs + 1;                 // top super-level index used
    if (Kg < 2) Kg = 2;  // the Horner sweeps below need at least A^(0), A^(1), A^(2)
    std::vector<Matrix<T>> Ag(Kg + 1, Matrix<T>(nb, nb, zero));
    std::vector<Matrix<T>> Bg(Kg + 1, Matrix<T>(nb, nb, zero));
    for (std::size_t i = 0; i <= Kg; ++i) {
        for (std::size_t u = 0; u < cs; ++u) {
            for (std::size_t v = 0; v < cs; ++v) {
                const long ka = static_cast<long>(i * cs + v) - static_cast<long>(u);
                if (ka >= 0 && static_cast<std::size_t>(ka) <= Kmax) {
                    const Matrix<T>& blk = A[static_cast<std::size_t>(ka)];
                    for (std::size_t p = 0; p < m; ++p)
                        for (std::size_t q = 0; q < m; ++q)
                            Ag[i](u * m + p, v * m + q) = blk(p, q);
                }
                const std::size_t kb = i * cs + v;
                if (kb <= Kmax) {
                    const Matrix<T>& blk = A[kb];
                    for (std::size_t p = 0; p < m; ++p)
                        for (std::size_t q = 0; q < m; ++q)
                            Bg[i](u * m + p, v * m + q) = blk(p, q);
                }
            }
        }
    }

    // ---- step 3: Ramaswami's recursion on the grouped chain ----
    Matrix<T> G(nb, nb, zero);
    const unsigned gmax = 100000u;
    for (unsigned it = 0; it < gmax; ++it) {
        // S = sum_{i=2}^{Kg} A^(i) G^(i-2) by Horner, so sum_{i>=2} A^(i) G^(i-1) = S G.
        Matrix<T> S = Ag[Kg];
        for (std::size_t i = Kg; i-- > 2;) S = mam::qbd_detail::madd(Ag[i], matmul(S, G));
        const Matrix<T> U = mam::qbd_detail::madd(Ag[1], matmul(S, G));
        const Matrix<T> ImU = mam::qbd_detail::msub(eye<T>(nb), U);
        const Matrix<T> Gn = matmul(inverse(ImU), Ag[0]);
        T gap = zero;
        for (std::size_t i = 0; i < nb; ++i)
            for (std::size_t j = 0; j < nb; ++j) {
                const T d = num_abs(T(Gn(i, j) - G(i, j)));
                if (d > gap) gap = d;
            }
        G = Gn;
        if (gap <= tol) break;
    }

    std::vector<Matrix<T>> Ahat(Kg + 2), Bhat(Kg + 2);
    Ahat[Kg + 1] = Matrix<T>(nb, nb, zero);
    Bhat[Kg + 1] = Matrix<T>(nb, nb, zero);
    for (std::size_t i = Kg + 1; i-- > 0;) {
        Ahat[i] = mam::qbd_detail::madd(Ag[i], matmul(Ahat[i + 1], G));
        Bhat[i] = mam::qbd_detail::madd(Bg[i], matmul(Bhat[i + 1], G));
    }

    const Matrix<T> W = inverse(mam::qbd_detail::msub(eye<T>(nb), Ahat[1]));
    Matrix<T> M0 = mam::qbd_detail::madd(Bg[0], matmul(matmul(Bhat[1], W), Ag[0]));
    for (std::size_t i = 0; i < nb; ++i) M0(i, i) -= one;
    std::vector<std::vector<T>> x;
    x.push_back(mam::qbd_detail::statvec(M0));

    T mass = zero;
    for (const T& v : x[0]) mass += v;
    for (std::size_t lvl = 1; lvl <= max_levels; ++lvl) {
        std::vector<T> acc(nb, zero);
        if (lvl <= Kg) {
            const std::vector<T> t = vecmul(x[0], Bhat[lvl]);
            for (std::size_t j = 0; j < nb; ++j) acc[j] += t[j];
        }
        const std::size_t kmin = (lvl + 1 > Kg) ? (lvl + 1 - Kg) : 1;
        for (std::size_t k = kmin; k + 1 <= lvl; ++k) {
            const std::vector<T> t = vecmul(x[k], Ahat[lvl - k + 1]);
            for (std::size_t j = 0; j < nb; ++j) acc[j] += t[j];
        }
        const std::vector<T> xn = vecmul(acc, W);
        T inc = zero;
        for (const T& v : xn) inc += v;
        x.push_back(xn);
        mass += inc;
        if (inc <= tol * mass) break;
    }
    for (std::vector<T>& row : x)
        for (T& v : row) v /= mass;

    // ---- step 4: unfold the super-levels into the level distribution ----
    MapDcResult<T> r;
    r.queueLengthDist.assign(dist_size, zero);
    T L = zero;
    for (std::size_t lvl = 0; lvl < x.size(); ++lvl) {
        for (std::size_t u = 0; u < cs; ++u) {
            T pn = zero;
            for (std::size_t j = 0; j < m; ++j) pn += x[lvl][u * m + j];
            const std::size_t n = lvl * cs + u;
            L += num_traits<T>::from_int(static_cast<long>(n)) * pn;
            if (n < dist_size) r.queueLengthDist[n] = pn;
        }
    }

    r.meanQueueLength = L;
    r.meanSojournTime = L / lambda;
    // waiting-time clamp rationale: see _kb/03-api-layer.md (cpp port notes: qsys)
    const T wq = r.meanSojournTime - s;
    r.meanWaitingTime = (wq > zero) ? wq : zero;
    r.utilization = rho;
    return r;
}

/**
 * qsys_mapdc with 100 materialized levels, an arrival-count cap of 4096, a
 * super-level cap of 20000 and tolerance 1e-14, matching the qsys_mapd1
 * defaults.
 */
template <class T>
MapDcResult<T> qsys_mapdc(const mam::Map<T>& arrival, const T& s, unsigned c) {
    return qsys_mapdc(arrival, s, c, static_cast<std::size_t>(100), 4096u,
                      static_cast<std::size_t>(20000), T(num_traits<T>::from_double(1e-14)));
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_MAPDC_H
