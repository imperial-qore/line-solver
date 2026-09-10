/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_MAPD1_H
#define LINE_API_QSYS_QSYS_MAPD1_H

/**
 * The MAP/D/1 FCFS queue: deterministic service of length s fed by a Markovian
 * arrival process.
 *
 * ALGORITHM, AND HOW IT DIFFERS FROM THE MATLAB REFERENCE.
 * matlab/src/api/qsys/qsys_mapd1.m delegates to qsys_mapdc, which calls
 * Q-MAM's Q_CT_MAP_D_C, which is not transcribed here. This port computes the
 * SAME quantities exactly, by the embedded-chain route, and NOT by an
 * Erlang-k approximation:
 * a deterministic service time has no phase-type representation, so no QBD can
 * carry it, but an M/G/1-type Markov chain can.
 *
 * The construction is in four exact steps.
 *
 * 1. MAP counting probabilities over one service. P_k(u) is the matrix whose
 *    (i,j) entry is P(k arrivals in [0,u], phase j at u | phase i at 0). They
 *    satisfy dP_k/du = P_k D0 + P_{k-1} D1, so the row (P_0(s) ... P_K(s)) is
 *    the top block row of exp(M s) with M the (K+1)-block bidiagonal matrix
 *    carrying D0 on the diagonal and D1 on the superdiagonal. One expm gives
 *    all of them. K is grown until the missing mass max_i (1 - sum_k P_k(s)e)
 *    is below the tolerance, so the truncation is a controlled quantity and not
 *    a modelling choice.
 *
 * 2. Time-in-level integrals. I_k = integral_0^s P_k(u) du follows from the
 *    same differential equation without a second expm: integrating gives
 *    P_k(s) - P_k(0) = I_k D0 + I_{k-1} D1, hence I_0 = (P_0(s) - I) D0^-1 and
 *    I_k = (P_k(s) - I_{k-1} D1) D0^-1. D0 is a non-singular sub-generator.
 *
 * 3. The embedded chain at departure epochs. With the level equal to the number
 *    left behind by a departure and the phase equal to the arrival phase, the
 *    chain is skip-free to the left with blocks A_k = P_k(s) from a busy level
 *    and B_k = (-D0)^-1 D1 A_k from the empty level (wait for the next arrival,
 *    then serve). Its stationary vector comes from Ramaswami's recursion: with
 *    G the minimal solution of G = sum_k A_k G^k, and the partial sums
 *    Ahat_i = sum_{k>=i} A_k G^(k-i), Bhat_i = sum_{k>=i} B_k G^(k-i),
 *
 *        x_0 (B_0 + Bhat_1 (I - Ahat_1)^-1 A_0) = x_0,
 *        x_n = [ x_0 Bhat_n + sum_{k=1}^{n-1} x_k Ahat_{n-k+1} ] (I - Ahat_1)^-1.
 *
 * 4. The time-stationary distribution. The departure-epoch vector x is NOT the
 *    time-stationary one unless the arrivals are Poisson, so the two are not
 *    interchanged here. Averaging over an inter-departure cycle with the
 *    integrals of step 2,
 *
 *        p_0   = lambda x_0 (-D0)^-1 e,
 *        p_n   = lambda [ y_0 I_{n-1} e + sum_{m=1}^{n} x_m I_{n-m} e ],  n >= 1,
 *
 *    with y_0 = x_0 (-D0)^-1 D1 the phase at the arrival that ends an idle
 *    period. Since I_k e sums to s e over k, this construction satisfies
 *    sum_n p_n = lambda (E[idle] + s) = 1 and p_0 = 1 - lambda s identically,
 *    which the test file asserts. The mean follows in closed form from
 *    sum_k I_k and sum_k k I_k without materializing the levels.
 *
 * meanWaitingTime is (L - rho)/lambda by Little's law and meanSojournTime is
 * that plus s.
 *
 * MEASURED AGREEMENT (MATLAB R2025a, T = double). meanQueueLength
 * agrees with the reference; meanWaitingTime does not, and the reference is
 * the one that is wrong -- see the defect note below.
 *  - M/D/1 collapse, qsys_mapd1(D0 = [-2], D1 = [2], s = 1/3), rho = 2/3.
 *    Textbook L = rho + rho^2/(2(1-rho)) = 4/3 and Wq = rho s/(2(1-rho)) = 1/3.
 *    The port returns 1.33333333333307 and 0.333333333333203, i.e. the
 *    textbook values to 2.0e-13 and 4.0e-13. MATLAB's meanQueueLength is
 *    1.333333330345585, itself 2.2e-9 below the textbook value (its own
 *    maxNumComp truncation), so the port and MATLAB differ by 2.24e-9.
 *  - Correlated MMPP2 arrivals D0 = [-2.5 0.2; 0.1 -0.7], D1 = diag(2.3, 0.6)
 *    (lambda = 7/6), s = 0.4, rho = 7/15: MATLAB meanQueueLength
 *    1.118300079597477, port 1.11830008358428, relative difference 3.57e-9.
 *  - Erlang-2 arrivals D0 = [-4 4; 0 -4], D1 = [0 0; 4 0] (lambda = 2),
 *    s = 0.3, rho = 0.6: MATLAB meanQueueLength 0.7758216464528508, port
 *    0.775821646757718, relative difference 3.93e-10.
 *
 * Both non-Poisson values were confirmed independently of Q-MAM by the
 * Erlang-k limit through LINE's own MATLAB qbd_mapmap1: replacing the
 * deterministic service by an Erlang-k of the same mean and Richardson-
 * extrapolating the O(1/k) convergence from k = 40 and k = 80 gives
 * 1.11832240806 and 0.775808057317, which match the port to 2.0e-5 and 1.7e-5,
 * the residual of the extrapolation itself. The port additionally satisfies
 * sum_n p_n = 1 and p_0 = 1 - rho to 1e-15 on all three instances, both being
 * identities of the construction rather than fitted quantities.
 *
 * At T = Real50 the port reproduces its own double results to 2e-14, so the
 * double values above are not precision-limited.
 *
 * DEFECT IN THE MATLAB REFERENCE (reported, not fixed here). qsys_mapdc, and
 * therefore qsys_mapd1, computes meanWaitingTime as a left-rectangle sum of
 * the survival function of the Q-MAM waiting-time CDF with step s/numSteps and
 * numSteps defaulting to 1. Two consequences:
 *   (a) even for Poisson arrivals the default is a one-point quadrature. On
 *       the M/D/1 instance above it returns 0.4444444443040794 against the
 *       exact 1/3, a +33% error, converging as O(1/numSteps): 0.3350694 at
 *       numSteps = 64 and 0.3334418 at numSteps = 1024.
 *   (b) for non-Poisson arrivals it does not converge to the right value at
 *       all. On the correlated instance it converges to 0.3540778458263887
 *       while Little's law applied to its own (correct) meanQueueLength gives
 *       0.5585429253692660; on the Erlang-2 instance it converges to
 *       0.1428343846295046 against 0.0879108232264254. Errors of -37% and
 *       +62%, in opposite directions, so this is not a quadrature artefact.
 * The reference's meanQueueLength is correct in all three cases -- that is what
 * the Erlang-k cross-check above establishes -- so this port takes Little's law
 * as the definition of meanWaitingTime and does not reproduce the reference's
 * waiting-time numbers. Reproduction: run qsys_mapd1([-2],[2],1/3,'numSteps',N)
 * for N = 1, 64, 1024 and compare meanWaitingTime against 1/3.
 *
 * ARITHMETIC. Gated on num_traits<T>::has_transcendental for two independent
 * reasons: step 1 calls expm, a scaling-and-squaring Pade approximation that is
 * tolerance-controlled and cannot be exact in any arithmetic, and step 3
 * computes G by a fixed-point iteration that does not terminate in a finite
 * number of field operations. Steps 2 and 4 are finite exact matrix algebra
 * given P_k(s) and G, and add no error of their own.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/api/mam/qbd_r.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/expm.h"
#include "line/util/linalg.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace qsys {

/** Return value of qsys_mapd1, mirroring the MATLAB result struct. */
template <class T>
struct MapD1Result {
    T meanQueueLength;               ///< L, number in system, time-stationary
    T meanWaitingTime;               ///< Wq = (L - rho)/lambda
    T meanSojournTime;               ///< W = Wq + s
    T utilization;                   ///< rho = lambda s
    std::vector<T> queueLengthDist;  ///< P(N = n), n = 0, 1, ...
};

namespace detail {

/** Right division X = C D^-1, formed from an explicit inverse of D. */
template <class T>
Matrix<T> rdivide(const Matrix<T>& C, const Matrix<T>& Dinv) {
    return matmul(C, Dinv);
}

/**
 * Counting probabilities P_k(s), k = 0..K, of a MAP over an interval of exact
 * length s, as the top block row of exp(M s) with M block bidiagonal. One expm
 * of order (K+1)n.
 */
template <class T>
std::vector<Matrix<T>> map_counting_at(const Matrix<T>& D0, const Matrix<T>& D1, const T& s,
                                       unsigned K) {
    const std::size_t n = D0.rows();
    const std::size_t dim = (static_cast<std::size_t>(K) + 1) * n;
    Matrix<T> M(dim, dim, num_traits<T>::from_int(0));
    for (unsigned b = 0; b <= K; ++b) {
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) M(b * n + i, b * n + j) = D0(i, j);
        if (b < K)
            for (std::size_t i = 0; i < n; ++i)
                for (std::size_t j = 0; j < n; ++j) M(b * n + i, (b + 1) * n + j) = D1(i, j);
    }
    const Matrix<T> E = expm(M, s);
    std::vector<Matrix<T>> P(static_cast<std::size_t>(K) + 1);
    for (unsigned b = 0; b <= K; ++b) {
        P[b] = Matrix<T>(n, n);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) P[b](i, j) = E(i, b * n + j);
    }
    return P;
}

}  // namespace detail

/**
 * MAP/D/1 by the exact embedded M/G/1-type chain.
 *
 * @param arrival      arrival MAP (D0, D1)
 * @param s            deterministic service time, s > 0
 * @param dist_size    how many entries of queueLengthDist to materialize
 * @param max_arrivals cap on K, the number of arrivals per service that is
 *                     tracked; K is grown from a Poisson-tail estimate until
 *                     the missing counting mass falls below tol
 * @param max_levels   cap on the number of embedded levels generated by
 *                     Ramaswami's recursion
 * @param tol          tolerance on the counting-mass truncation, on the G
 *                     iteration and on the level-tail truncation
 */
template <class T>
MapD1Result<T> qsys_mapd1(const mam::Map<T>& arrival, const T& s, std::size_t dist_size,
                          unsigned max_arrivals, std::size_t max_levels, const T& tol) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_mapd1 requires transcendental arithmetic");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    if (s <= zero) throw InputError("qsys_mapd1: service time s must be positive");
    if (dist_size == 0) throw InputError("qsys_mapd1: dist_size must be positive");
    if (max_levels == 0) throw InputError("qsys_mapd1: max_levels must be positive");
    const std::size_t n = arrival.D0.rows();
    if (arrival.D0.cols() != n || arrival.D1.rows() != n || arrival.D1.cols() != n)
        throw InputError("qsys_mapd1: D0 and D1 must be square and of equal order");

    const Matrix<T>& D0 = arrival.D0;
    const Matrix<T>& D1 = arrival.D1;
    const T lambda = mam::map_lambda(arrival);
    if (lambda <= zero) throw InputError("qsys_mapd1: non-positive arrival rate");
    const T rho = lambda * s;
    if (rho >= one) throw InputError("qsys_mapd1: load rho must be strictly less than 1");

    // counting-probability start rationale: see _kb/03-api-layer.md (cpp port notes: qsys)
    T numax = zero;
    for (std::size_t i = 0; i < n; ++i) {
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
        for (std::size_t i = 0; i < n; ++i) {
            T row = zero;
            for (std::size_t k = 0; k < A.size(); ++k)
                for (std::size_t j = 0; j < n; ++j) row += A[k](i, j);
            const T lack = one - row;
            if (lack > missing) missing = lack;
        }
        if (missing <= tol || K >= max_arrivals) break;
        K = (2u * K < max_arrivals) ? 2u * K : max_arrivals;
        A = detail::map_counting_at(D0, D1, s, K);
    }

    // ---- step 2: time-in-level integrals over one service ----
    const Matrix<T> invD0 = inverse(D0);
    std::vector<Matrix<T>> I(A.size());
    {
        Matrix<T> C = A[0];
        for (std::size_t i = 0; i < n; ++i) C(i, i) -= one;
        I[0] = detail::rdivide(C, invD0);
        for (std::size_t k = 1; k < A.size(); ++k) {
            Matrix<T> Ck = A[k];
            const Matrix<T> prev = matmul(I[k - 1], D1);
            for (std::size_t i = 0; i < n; ++i)
                for (std::size_t j = 0; j < n; ++j) Ck(i, j) -= prev(i, j);
            I[k] = detail::rdivide(Ck, invD0);
        }
    }

    // ---- step 3: embedded chain at departures, Ramaswami's recursion ----
    Matrix<T> negD0inv(n, n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) negD0inv(i, j) = -invD0(i, j);
    const Matrix<T> Pa = matmul(negD0inv, D1);  // phase at the arrival ending an idle period
    std::vector<Matrix<T>> B(A.size());
    for (std::size_t k = 0; k < A.size(); ++k) B[k] = matmul(Pa, A[k]);

    // G by the U-based natural iteration G <- (I - A_1 - sum_{k>=2} A_k G^(k-1))^-1 A_0.
    Matrix<T> G(n, n, zero);
    const unsigned gmax = 100000u;
    const std::size_t ktop = A.size() - 1;  // = K, and K >= 8 by construction
    for (unsigned it = 0; it < gmax; ++it) {
        // S = sum_{k=2}^{K} A_k G^(k-2) by Horner, so sum_{k>=2} A_k G^(k-1) = S G.
        Matrix<T> S = A[ktop];
        for (std::size_t k = ktop; k-- > 2;) S = mam::qbd_detail::madd(A[k], matmul(S, G));
        const Matrix<T> U = mam::qbd_detail::madd(A[1], matmul(S, G));
        const Matrix<T> ImU = mam::qbd_detail::msub(eye<T>(n), U);
        const Matrix<T> Gn = matmul(inverse(ImU), A[0]);
        T gap = zero;
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) {
                const T d = num_abs(T(Gn(i, j) - G(i, j)));
                if (d > gap) gap = d;
            }
        G = Gn;
        if (gap <= tol) break;
    }

    // Ahat_i = sum_{k>=i} A_k G^(k-i), Bhat_i likewise, by Horner from the top.
    const std::size_t K1 = A.size() - 1;
    std::vector<Matrix<T>> Ahat(K1 + 2), Bhat(K1 + 2);
    Ahat[K1 + 1] = Matrix<T>(n, n, zero);
    Bhat[K1 + 1] = Matrix<T>(n, n, zero);
    for (std::size_t i = K1 + 1; i-- > 0;) {
        Ahat[i] = mam::qbd_detail::madd(A[i], matmul(Ahat[i + 1], G));
        Bhat[i] = mam::qbd_detail::madd(B[i], matmul(Bhat[i + 1], G));
    }

    const Matrix<T> W = inverse(mam::qbd_detail::msub(eye<T>(n), Ahat[1]));
    // x_0 is the stationary vector of the level-0 censored chain.
    Matrix<T> M0 = mam::qbd_detail::madd(B[0], matmul(matmul(Bhat[1], W), A[0]));
    for (std::size_t i = 0; i < n; ++i) M0(i, i) -= one;
    std::vector<std::vector<T>> x;
    x.push_back(mam::qbd_detail::statvec(M0));

    T mass = zero;
    for (const T& v : x[0]) mass += v;
    for (std::size_t lvl = 1; lvl <= max_levels; ++lvl) {
        std::vector<T> acc(n, zero);
        if (lvl <= K1) {
            const std::vector<T> t = vecmul(x[0], Bhat[lvl]);
            for (std::size_t j = 0; j < n; ++j) acc[j] += t[j];
        }
        const std::size_t kmin = (lvl + 1 > K1) ? (lvl + 1 - K1) : 1;
        for (std::size_t k = kmin; k + 1 <= lvl; ++k) {
            const std::vector<T> t = vecmul(x[k], Ahat[lvl - k + 1]);
            for (std::size_t j = 0; j < n; ++j) acc[j] += t[j];
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

    // ---- step 4: time-stationary distribution and its mean ----
    const std::vector<T> e = ones<T>(n);
    const std::vector<T> y0 = vecmul(x[0], Pa);
    std::vector<std::vector<T>> Ie(I.size());
    for (std::size_t k = 0; k < I.size(); ++k) Ie[k] = mulvec(I[k], e);
    std::vector<T> Isum_e(n, zero), Iw_e(n, zero);
    for (std::size_t k = 0; k < I.size(); ++k)
        for (std::size_t j = 0; j < n; ++j) {
            Isum_e[j] += Ie[k][j];
            Iw_e[j] += num_traits<T>::from_int(static_cast<long>(k)) * Ie[k][j];
        }

    T L = zero;
    for (std::size_t m = 1; m < x.size(); ++m) {
        T a = zero, b = zero;
        for (std::size_t j = 0; j < n; ++j) {
            a += x[m][j] * Isum_e[j];
            b += x[m][j] * Iw_e[j];
        }
        L += num_traits<T>::from_int(static_cast<long>(m)) * a + b;
    }
    for (std::size_t j = 0; j < n; ++j) L += y0[j] * (Isum_e[j] + Iw_e[j]);
    L *= lambda;

    MapD1Result<T> r;
    r.meanQueueLength = L;
    r.meanWaitingTime = (L - rho) / lambda;
    r.meanSojournTime = r.meanWaitingTime + s;
    r.utilization = rho;

    r.queueLengthDist.assign(dist_size, zero);
    {
        const std::vector<T> idle = mulvec(negD0inv, e);
        T p0 = zero;
        for (std::size_t j = 0; j < n; ++j) p0 += x[0][j] * idle[j];
        r.queueLengthDist[0] = lambda * p0;
        for (std::size_t lvl = 1; lvl < dist_size; ++lvl) {
            T acc = zero;
            if (lvl - 1 < Ie.size())
                for (std::size_t j = 0; j < n; ++j) acc += y0[j] * Ie[lvl - 1][j];
            for (std::size_t m = 1; m <= lvl && m < x.size(); ++m) {
                const std::size_t kk = lvl - m;
                if (kk >= Ie.size()) continue;
                for (std::size_t j = 0; j < n; ++j) acc += x[m][j] * Ie[kk][j];
            }
            r.queueLengthDist[lvl] = lambda * acc;
        }
    }
    return r;
}

/**
 * qsys_mapd1 with 100 materialized levels, an arrival-count cap of 4096, a
 * level cap of 20000 and tolerance 1e-14.
 */
template <class T>
MapD1Result<T> qsys_mapd1(const mam::Map<T>& arrival, const T& s) {
    return qsys_mapd1(arrival, s, static_cast<std::size_t>(100), 4096u,
                      static_cast<std::size_t>(20000), T(num_traits<T>::from_double(1e-14)));
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_MAPD1_H
