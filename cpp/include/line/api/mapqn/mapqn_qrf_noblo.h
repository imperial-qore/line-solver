/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAPQN_MAPQN_QRF_NOBLO_H
#define LINE_API_MAPQN_MAPQN_QRF_NOBLO_H

/**
 * The QRF no-blocking nonlinear bounds: `qrf_noblo_mmi`, `qrf_noblo_mem` and
 * the load-dependent `qrf_noblo_mmi_ld`.
 *
 * Port of python/line_solver/api/mapqn/qrf_noblo_{mmi,mem,mmi_ld}.py and the
 * `sub_qrfcon_noblo` constraint inventory they share. `api/mapqn` has no MATLAB
 * implementation, so Python and the JAR are the references.
 *
 * THE ARITY OF q SELECTS THE MODEL, exactly as the AMPL skeletons do:
 *
 * - a 4D q[i][j][k][h] is the population-free form of `qrboundsbas_skel.mod`,
 *   used by MMI and MEM. THM1 is stated on the aggregated e[i,k], which is
 *   exact here because no alpha makes the rate population-dependent.
 * - a 5D q[i][j][k][h][n] is the load-dependent form of
 *   `qrboundsrsrd_skel.mod`, whose fifth index is the population of the
 *   EMITTING station. THM1 is then stated per population against the
 *   station-i marginal.
 *
 * THM30 AND THM3 ARE WHAT MAKE THE POLYTOPE DEPEND ON THE SERVICE RATES AT
 * ALL. They are the marginal-balance families; without them an LP over the
 * remaining constraints returns the vacuous [0,1] for every instance (the
 * reference verified that against glpsol). Both are emitted in either form and
 * differ only in the q lookup.
 *
 * ARITHMETIC: transcendental, inherited from the objectives.
 */

#include <algorithm>
#include <cstddef>
#include <string>
#include <vector>

#include "line/api/mapqn/mapqn_qrf_common.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mapqn {

/**
 * The transition rates the constraint inventory reads.
 *
 * `load_dependent` selects the arity: with it false only `q4` is read, with it
 * true only `q5`. They are kept as two members rather than one flattened array
 * because confusing the two is the failure mode the reference documents -- a
 * population-free q cannot distinguish the rate at which a station empties at
 * population n from the rate at n', so the polytope stops pinning the
 * utilization at all.
 */
template <class T>
struct QrfRates {
    bool load_dependent = false;
    std::size_t M = 0, Kmax = 0, N = 0;
    std::vector<T> q4;  ///< [i][j][k][h]
    std::vector<T> q5;  ///< [i][j][k][h][n], n the EMITTING station's population

    T at(std::size_t i, std::size_t j, std::size_t k, std::size_t h, std::size_t n) const {
        if (load_dependent)
            return q5[((((i * M + j) * Kmax + k) * Kmax + h) * (N + 1)) + n];
        return q4[((i * M + j) * Kmax + k) * Kmax + h];
    }
    /** The population-free lookup, for the arms that do not carry n. */
    T at4(std::size_t i, std::size_t j, std::size_t k, std::size_t h) const {
        return q4[((i * M + j) * Kmax + k) * Kmax + h];
    }
};

/** build_q_from_mu_v_rt: the population-free rates. */
template <class T>
QrfRates<T> qrf_build_q(std::size_t M, const std::vector<int>& K, const Matrix<T>& mu,
                        const Matrix<T>& v, const Matrix<T>& rt) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t Kmax = static_cast<std::size_t>(*std::max_element(K.begin(), K.end()));
    QrfRates<T> q;
    q.load_dependent = false;
    q.M = M;
    q.Kmax = Kmax;
    q.q4.assign(M * M * Kmax * Kmax, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t j = 0; j < M; ++j)
            for (std::size_t k = 0; k < static_cast<std::size_t>(K[i]); ++k)
                for (std::size_t h = 0; h < static_cast<std::size_t>(K[i]); ++h) {
                    const T mv = mu(i, k * Kmax + h);
                    q.q4[((i * M + j) * Kmax + k) * Kmax + h] =
                        (j != i) ? T(rt(i, j) * mv) : T(v(i, k * Kmax + h) + rt(i, i) * mv);
                }
    return q;
}

/**
 * build_q_ld: the load-dependent rates.
 *
 * `qrboundsrsrd_skel.mod:11` declares q with FIVE indices, the fifth being the
 * population n of the EMITTING station, with q[...,0] = 0, and the scaling
 * alpha[i,n] multiplying BOTH the background term v and the completion term
 * rt[i,i] mu. Dropping alpha from the v term is as wrong as dropping the index.
 *
 * @param alpha (M x N), 0-based in the population, so alpha(i, n-1) is the
 *              AMPL alpha[i,n]; an empty matrix means all ones
 */
template <class T>
QrfRates<T> qrf_build_q_ld(std::size_t M, const std::vector<int>& K, const Matrix<T>& mu,
                           const Matrix<T>& v, const Matrix<T>& rt, std::size_t N,
                           const Matrix<T>& alpha) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t Kmax = static_cast<std::size_t>(*std::max_element(K.begin(), K.end()));
    const bool unit = (alpha.rows() == 0 || alpha.cols() == 0);
    if (!unit && (alpha.rows() != M || alpha.cols() != N))
        throw InputError("qrf_build_q_ld: alpha must be (M x N)");
    QrfRates<T> q;
    q.load_dependent = true;
    q.M = M;
    q.Kmax = Kmax;
    q.N = N;
    q.q5.assign(M * M * Kmax * Kmax * (N + 1), zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t j = 0; j < M; ++j)
            for (std::size_t k = 0; k < static_cast<std::size_t>(K[i]); ++k)
                for (std::size_t h = 0; h < static_cast<std::size_t>(K[i]); ++h)
                    for (std::size_t n = 1; n <= N; ++n) {
                        const T a = unit ? one : alpha(i, n - 1);
                        const T mv = mu(i, k * Kmax + h);
                        const std::size_t at =
                            ((((i * M + j) * Kmax + k) * Kmax + h) * (N + 1)) + n;
                        q.q5[at] = (j != i) ? T(rt(i, j) * mv * a)
                                            : T(v(i, k * Kmax + h) * a + rt(i, i) * mv * a);
                    }
    return q;
}

/** extract_mu_v_from_maps: the completion and background rates of each MAP. */
template <class T>
void qrf_extract_mu_v(const std::vector<std::pair<Matrix<T>, Matrix<T>>>& MAPs, std::size_t M,
                      const std::vector<int>& K, Matrix<T>* mu, Matrix<T>* v) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t Kmax = static_cast<std::size_t>(*std::max_element(K.begin(), K.end()));
    *mu = Matrix<T>(M, Kmax * Kmax, zero);
    *v = Matrix<T>(M, Kmax * Kmax, zero);
    for (std::size_t i = 0; i < M; ++i) {
        const Matrix<T>& D0 = MAPs[i].first;
        const Matrix<T>& D1 = MAPs[i].second;
        for (std::size_t h = 0; h < static_cast<std::size_t>(K[i]); ++h)
            for (std::size_t k = 0; k < static_cast<std::size_t>(K[i]); ++k) {
                (*mu)(i, h * Kmax + k) = D1(h, k);
                // The DIAGONAL of D0 is the total exit rate, not a background
                // transition, so it is dropped here rather than carried.
                // (from, to), the order mu and the q assembly already use.
                (*v)(i, h * Kmax + k) = (h == k) ? zero : D0(h, k);
            }
    }
}

/** The two constraint blocks: g(x) <= 0 and h(x) = 0. */
template <class T>
struct QrfConstraints {
    std::vector<T> ineq;
    std::vector<T> eq;
};

/**
 * The full no-blocking constraint inventory, `sub_qrfcon_noblo`.
 *
 * @param x  decision vector
 * @param q  the rates, whose arity selects the model
 * @param BB (MR x M) blocking-state matrix
 * @param F  (M) capacity per queue
 */
template <class T>
QrfConstraints<T> sub_qrfcon_noblo(const std::vector<T>& x, const QrfRates<T>& q, std::size_t M,
                                   std::size_t MR, const Matrix<int>& BB,
                                   const std::vector<int>& F, std::size_t N,
                                   const std::vector<int>& K) {
    const T zero = num_traits<T>::from_int(0);
    const QrfVars<T> V = sub_qrfvar(x, M, N, K, MR);
    QrfConstraints<T> out;
    std::vector<T>& ceq = out.eq;
    std::vector<T>& c = out.ineq;
    auto Ki = [&K](std::size_t i) { return static_cast<std::size_t>(K[i]); };
    auto Fi = [&F](std::size_t i) { return static_cast<std::size_t>(F[i]); };

    // ONE: each station's diagonal marginal is a probability distribution.
    for (std::size_t j = 0; j < M; ++j) {
        T val = zero;
        for (std::size_t nj = 0; nj <= N; ++nj)
            for (std::size_t k = 0; k < Ki(j); ++k)
                for (std::size_t m = 0; m < MR; ++m) val += V.p(j, nj, k, j, nj, k, m);
        ceq.push_back(val - num_traits<T>::from_int(1));
    }

    // ZERO1/ZERO2/ZERO3: a pair entry that describes an impossible joint state.
    for (std::size_t j = 0; j < M; ++j)
        for (std::size_t k = 0; k < Ki(j); ++k)
            for (std::size_t nj = 0; nj <= N; ++nj)
                for (std::size_t i = 0; i < M; ++i)
                    for (std::size_t h = 0; h < Ki(i); ++h)
                        for (std::size_t ni = 0; ni <= N; ++ni)
                            for (std::size_t m = 0; m < MR; ++m) {
                                if (i == j && nj == ni && h != k)
                                    ceq.push_back(V.p(j, nj, k, i, ni, h, m));  // ZERO1
                                if (i == j && nj != ni)
                                    ceq.push_back(V.p(j, nj, k, i, ni, h, m));  // ZERO2
                                if (i != j && nj + ni > N)
                                    ceq.push_back(V.p(j, nj, k, i, ni, h, m));  // ZERO3
                            }

    // ZERO5: a blocked station holds no empty-population mass.
    for (std::size_t j = 0; j < M; ++j)
        for (std::size_t k = 0; k < Ki(j); ++k)
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t h = 0; h < Ki(i); ++h)
                    for (std::size_t ni = 0; ni <= Fi(i); ++ni)
                        for (std::size_t m = 1; m < MR; ++m)
                            if (BB(m, j) == 1) ceq.push_back(V.p(j, 0, k, i, ni, h, m));

    // ZERO6: above the buffer there is no mass.
    for (std::size_t j = 0; j < M; ++j)
        for (std::size_t k = 0; k < Ki(j); ++k)
            for (std::size_t nj = Fi(j) + 1; nj <= N; ++nj)
                for (std::size_t i = 0; i < M; ++i)
                    for (std::size_t h = 0; h < Ki(i); ++h)
                        for (std::size_t ni = 0; ni <= N; ++ni)
                            for (std::size_t m = 0; m < MR; ++m)
                                ceq.push_back(V.p(j, nj, k, i, ni, h, m));

    // ZERO7 is the blocking family and is empty at MR = 1, which is the whole
    // of the no-blocking model; it is not emitted rather than emitted as 0 = 0.

    // SYMMETRY: the pair tensor describes an unordered pair.
    for (std::size_t j = 0; j < M; ++j)
        for (std::size_t nj = 0; nj <= N; ++nj)
            for (std::size_t k = 0; k < Ki(j); ++k)
                for (std::size_t i = 0; i < M; ++i)
                    for (std::size_t ni = 0; ni <= N; ++ni)
                        for (std::size_t h = 0; h < Ki(i); ++h)
                            for (std::size_t m = 0; m < MR; ++m)
                                ceq.push_back(V.p(i, ni, h, j, nj, k, m) -
                                              V.p(j, nj, k, i, ni, h, m));

    // MARGINALS: the diagonal entry is the marginal of the pair over the other
    // station. The inner sum runs over the POPULATION range 0..N, not over the
    // 1-based index, which is where the reference's twin twice went wrong.
    for (std::size_t j = 0; j < M; ++j)
        for (std::size_t k = 0; k < Ki(j); ++k)
            for (std::size_t nj = 0; nj <= N; ++nj)
                for (std::size_t i = 0; i < M; ++i)
                    for (std::size_t m = 0; m < MR; ++m) {
                        if (i == j) continue;
                        T val = V.p(j, nj, k, j, nj, k, m);
                        for (std::size_t ni = 0; ni <= N; ++ni)
                            for (std::size_t h = 0; h < Ki(i); ++h)
                                val -= V.p(j, nj, k, i, ni, h, m);
                        ceq.push_back(val);
                    }

    // UEFF: the effective rate variable is the busy marginal of the pair.
    for (std::size_t j = 0; j < M; ++j)
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t ki = 0; ki < Ki(i); ++ki) {
                T val = V.e[i * V.Kmax + ki];
                for (std::size_t nj = 0; nj <= N; ++nj)
                    for (std::size_t kj = 0; kj < Ki(j); ++kj)
                        for (std::size_t m = 0; m < MR; ++m)
                            for (std::size_t ni = 1; ni <= N; ++ni)
                                if (BB(m, i) == 0) val -= V.p(j, nj, kj, i, ni, ki, m);
                ceq.push_back(val);
            }

    // THM1: phase balance at station i. Stated on the aggregated e in the
    // population-free form, and per population against the station-i marginal
    // in the load-dependent one, because there the rate depends on n.
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < Ki(i); ++k) {
            T val = zero;
            if (q.load_dependent) {
                for (std::size_t ni = 1; ni <= Fi(i); ++ni)
                    for (std::size_t m = 0; m < MR; ++m)
                        for (std::size_t j = 0; j < M; ++j)
                            for (std::size_t h = 0; h < Ki(i); ++h) {
                                val += q.at(i, j, k, h, ni) * V.p(i, ni, k, i, ni, k, m);
                                val -= q.at(i, j, h, k, ni) * V.p(i, ni, h, i, ni, h, m);
                            }
            } else {
                for (std::size_t j = 0; j < M; ++j)
                    for (std::size_t h = 0; h < Ki(i); ++h) {
                        val += q.at4(i, j, k, h) * V.e[i * V.Kmax + k];
                        val -= q.at4(i, j, h, k) * V.e[i * V.Kmax + h];
                    }
            }
            ceq.push_back(val);
        }

    // THM2: the population carried by the conditional distribution is N.
    for (std::size_t j = 0; j < M; ++j)
        for (std::size_t k = 0; k < Ki(j); ++k)
            for (std::size_t nj = 0; nj <= Fi(j); ++nj)
                for (std::size_t m = 0; m < MR; ++m) {
                    T val = zero;
                    for (std::size_t i = 0; i < M; ++i)
                        for (std::size_t ni = 1; ni <= Fi(i); ++ni)
                            for (std::size_t ki = 0; ki < Ki(i); ++ki)
                                val += num_traits<T>::from_int(static_cast<long>(ni)) *
                                       V.p(j, nj, k, i, ni, ki, m);
                    val -= num_traits<T>::from_int(static_cast<long>(N)) *
                           V.p(j, nj, k, j, nj, k, m);
                    ceq.push_back(val);
                }

    // COR1: the second moment of the total population is N^2.
    {
        T val = zero;
        for (std::size_t m = 0; m < MR; ++m)
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t j = 0; j < M; ++j)
                    for (std::size_t nj = 1; nj <= Fi(j); ++nj)
                        for (std::size_t ni = 1; ni <= Fi(i); ++ni)
                            for (std::size_t ki = 0; ki < Ki(i); ++ki)
                                for (std::size_t kj = 0; kj < Ki(j); ++kj)
                                    val += num_traits<T>::from_int(
                                               static_cast<long>(ni * nj)) *
                                           V.p(j, nj, kj, i, ni, ki, m);
        val -= num_traits<T>::from_int(static_cast<long>(N * N));
        ceq.push_back(val);
    }

    // THM30 {i, u}: balance across the ni = 0 boundary of station i.
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t u = 0; u < Ki(i); ++u) {
            T val = zero;
            for (std::size_t j = 0; j < M; ++j) {
                if (j == i) continue;
                for (std::size_t nj = 1; nj <= Fi(j); ++nj)
                    for (std::size_t k = 0; k < Ki(j); ++k) {
                        T coef = zero;
                        for (std::size_t h = 0; h < Ki(j); ++h) coef += q.at(j, i, k, h, nj);
                        if (coef == zero) continue;
                        for (std::size_t m = 0; m < MR; ++m)
                            val += coef * V.p(j, nj, k, i, 0, u, m);
                    }
            }
            for (std::size_t j = 0; j < M; ++j) {
                if (j == i) continue;
                for (std::size_t nj = 0; nj <= Fi(j); ++nj)
                    for (std::size_t k = 0; k < Ki(i); ++k) {
                        // The emitting population is 1: the transition leaves
                        // station i holding exactly one job.
                        const T coef = q.at(i, j, k, u, 1);
                        if (coef == zero) continue;
                        for (std::size_t h = 0; h < Ki(j); ++h)
                            for (std::size_t m = 0; m < MR; ++m)
                                val -= coef * V.p(j, nj, h, i, 1, k, m);
                    }
            }
            ceq.push_back(val);
        }

    // THM3 {i, ni in 0..F[i]-1}: balance across the ni -> ni+1 boundary.
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t ni = 0; ni + 1 <= Fi(i); ++ni) {
            T val = zero;
            for (std::size_t j = 0; j < M; ++j) {
                if (j == i) continue;
                for (std::size_t nj = 1; nj <= Fi(j); ++nj)
                    for (std::size_t k = 0; k < Ki(j); ++k) {
                        T coef = zero;
                        for (std::size_t h = 0; h < Ki(j); ++h) coef += q.at(j, i, k, h, nj);
                        if (coef == zero) continue;
                        for (std::size_t u = 0; u < Ki(i); ++u)
                            for (std::size_t m = 0; m < MR; ++m)
                                val += coef * V.p(j, nj, k, i, ni, u, m);
                    }
            }
            for (std::size_t j = 0; j < M; ++j) {
                if (j == i) continue;
                for (std::size_t nj = 0; nj <= Fi(j); ++nj)
                    for (std::size_t k = 0; k < Ki(i); ++k) {
                        T coef = zero;
                        for (std::size_t h = 0; h < Ki(i); ++h)
                            coef += q.at(i, j, k, h, ni + 1);
                        if (coef == zero) continue;
                        for (std::size_t u = 0; u < Ki(j); ++u)
                            for (std::size_t m = 0; m < MR; ++m)
                                val -= coef * V.p(j, nj, u, i, ni + 1, k, m);
                    }
            }
            ceq.push_back(val);
        }

    // THM4: an inequality, stated as >= upstream and stored with the sign
    // swapped so that every row of this block reads g(x) <= 0.
    for (std::size_t j = 0; j < M; ++j)
        for (std::size_t k = 0; k < Ki(j); ++k)
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t m = 0; m < MR; ++m) {
                    T val = zero;
                    for (std::size_t t = 0; t < M; ++t)
                        for (std::size_t h = 0; h < Ki(t); ++h)
                            for (std::size_t njx = 0; njx <= N; ++njx)
                                for (std::size_t nt = 0; nt <= N; ++nt)
                                    val -= num_traits<T>::from_int(static_cast<long>(nt)) *
                                           V.p(j, njx, k, t, nt, h, m);
                    for (std::size_t h = 0; h < Ki(i); ++h)
                        for (std::size_t njx = 0; njx <= N; ++njx)
                            for (std::size_t ni = 1; ni <= N; ++ni)
                                val += num_traits<T>::from_int(static_cast<long>(N)) *
                                       V.p(j, njx, k, i, ni, h, m);
                    c.push_back(val);
                }

    return out;
}

namespace qrfdetail {

/** Assemble, reduce, start and solve: the body every entry point shares. */
template <class T, class Obj, class Grad>
QrfMetrics<T> run_noblo(const QrfRates<T>& q, std::size_t M, std::size_t MR,
                        const Matrix<int>& BB, const std::vector<int>& F, std::size_t N,
                        const std::vector<int>& K, Obj objective, Grad gradient,
                        const std::string& name, const Matrix<T>* alpha = nullptr) {
    const std::size_t n = qrf_num_vars(M, N, K, MR);

    // Both blocks come from ONE residual callback each, so a family added to
    // the inventory reaches the equalities and the inequalities together.
    const QrfAffine<T> eq = qrf_affine_matrices<T>(
        [&](const std::vector<T>& z) { return sub_qrfcon_noblo(z, q, M, MR, BB, F, N, K).eq; }, n);
    const QrfAffine<T> ub = qrf_affine_matrices<T>(
        [&](const std::vector<T>& z) { return sub_qrfcon_noblo(z, q, M, MR, BB, F, N, K).ineq; },
        n);
    const QrfReduced<T> red = qrf_reduce_equalities(eq.A, eq.b);

    const std::vector<T> x0 = qrf_feasible_start(red.A, red.b, ub.A, ub.b, n);
    const std::vector<T> xopt =
        solve_qrf_nlp(objective, gradient, x0, red.A, red.b, ub.A, ub.b, name);
    return qrf_extract_results(sub_qrfvar(xopt, M, N, K, MR), M, K, F, MR, alpha);
}

}  // namespace qrfdetail

/**
 * `qrf_noblo_mmi`: the no-blocking bound under mutual-information minimization.
 *
 * @param M  number of queues
 * @param K  phases per queue
 * @param N  total population
 * @param mu (M x Kmax*Kmax) completion rates
 * @param v  (M x Kmax*Kmax) background rates
 * @param rt (M x M) routing matrix
 */
template <class T>
QrfMetrics<T> qrf_noblo_mmi(std::size_t M, const std::vector<int>& K, std::size_t N,
                            const Matrix<T>& mu, const Matrix<T>& v, const Matrix<T>& rt) {
    const std::size_t MR = 1;  // no blocking: one configuration, by definition
    Matrix<int> BB(1, M, 0);
    const std::vector<int> F(M, static_cast<int>(N));
    const QrfRates<T> q = qrf_build_q(M, K, mu, v, rt);
    const std::vector<long> idx = qrf_index_map(M, N, K, MR);
    return qrfdetail::run_noblo<T>(
        q, M, MR, BB, F, N, K,
        [&](const std::vector<T>& x) { return mmi_objective(x, M, N, K, F, MR); },
        [&](const std::vector<T>& x) { return mmi_gradient(x, M, N, K, F, MR, idx); },
        "qrf_noblo_mmi");
}

/**
 * `qrf_noblo_bethe`: the same polytope under the tree-reweighted free entropy.
 *
 * The polytope, the phase-1 feasible start and the NLP call are exactly those
 * of `qrf_noblo_mmi`; the objective is the only difference. See
 * `bethe_objective` in `mapqn_qrf_common.h` for what it is and why the uniform
 * edge weight is lambda = 1/M.
 *
 * ONE SOLVE, NO RESTARTS. The objective is convex on this polytope, so there
 * is no second local minimum for a restart to find; the single solve from the
 * phase-1 point returns the global optimum.
 *
 * @param M  number of queues
 * @param K  phases per queue
 * @param N  total population
 * @param mu (M x Kmax*Kmax) completion rates
 * @param v  (M x Kmax*Kmax) background rates
 * @param rt (M x M) routing matrix
 */
template <class T>
QrfMetrics<T> qrf_noblo_bethe(std::size_t M, const std::vector<int>& K, std::size_t N,
                              const Matrix<T>& mu, const Matrix<T>& v, const Matrix<T>& rt) {
    const std::size_t MR = 1;  // no blocking: one configuration, by definition
    Matrix<int> BB(1, M, 0);
    const std::vector<int> F(M, static_cast<int>(N));
    const QrfRates<T> q = qrf_build_q(M, K, mu, v, rt);
    const std::vector<long> idx = qrf_index_map(M, N, K, MR);
    return qrfdetail::run_noblo<T>(
        q, M, MR, BB, F, N, K,
        [&](const std::vector<T>& x) { return bethe_objective(x, M, N, K, F, MR); },
        [&](const std::vector<T>& x) { return bethe_gradient(x, M, N, K, F, MR, idx); },
        "qrf_noblo_bethe");
}

/** `qrf_noblo_mem`: the same polytope under maximum entropy. */
template <class T>
QrfMetrics<T> qrf_noblo_mem(const std::vector<std::pair<Matrix<T>, Matrix<T>>>& MAPs,
                            std::size_t N, const Matrix<T>& rt) {
    const std::size_t M = MAPs.size();
    std::vector<int> K(M, 0);
    for (std::size_t i = 0; i < M; ++i) K[i] = static_cast<int>(MAPs[i].first.rows());
    Matrix<T> mu, v;
    qrf_extract_mu_v(MAPs, M, K, &mu, &v);

    const std::size_t MR = 1;
    Matrix<int> BB(1, M, 0);
    const std::vector<int> F(M, static_cast<int>(N));
    const QrfRates<T> q = qrf_build_q(M, K, mu, v, rt);
    const std::vector<long> idx = qrf_index_map(M, N, K, MR);
    return qrfdetail::run_noblo<T>(
        q, M, MR, BB, F, N, K,
        [&](const std::vector<T>& x) { return mem_objective(x, M, N, K, F, MR); },
        [&](const std::vector<T>& x) { return mem_gradient(x, M, N, K, F, MR, idx); },
        "qrf_noblo_mem");
}

/**
 * `qrf_noblo_mmi_ld`: MMI on the LOAD-DEPENDENT polytope.
 *
 * The only difference from `qrf_noblo_mmi` is the arity of q, and that is the
 * whole point: with a population-free q the balance families cannot tell the
 * rate at which a station empties at population n from the rate at n', so the
 * polytope stops pinning the utilization.
 *
 * @param alpha (M x N) load-dependent scaling, alpha(i, n-1) being the AMPL
 *              alpha[i,n]; an empty matrix means all ones
 */
template <class T>
QrfMetrics<T> qrf_noblo_mmi_ld(std::size_t M, const std::vector<int>& K, std::size_t N,
                               const Matrix<T>& mu, const Matrix<T>& v, const Matrix<T>& rt,
                               const Matrix<T>& alpha) {
    const std::size_t MR = 1;
    Matrix<int> BB(1, M, 0);
    const std::vector<int> F(M, static_cast<int>(N));
    const QrfRates<T> q = qrf_build_q_ld(M, K, mu, v, rt, N, alpha);
    const std::vector<long> idx = qrf_index_map(M, N, K, MR);
    return qrfdetail::run_noblo<T>(
        q, M, MR, BB, F, N, K,
        [&](const std::vector<T>& x) { return mmi_objective(x, M, N, K, F, MR); },
        [&](const std::vector<T>& x) { return mmi_gradient(x, M, N, K, F, MR, idx); },
        "qrf_noblo_mmi_ld", &alpha);
}

/**
 * `qrf_noblo_mmi_linear`: the load-dependent no-blocking bound, under MMI.
 *
 * The `linear` in the name is about HOW the reference builds its constraints,
 * not about which constraints they are and not about the objective: it emits
 * the same inventory directly as sparse matrices instead of recovering it from
 * a residual callback, because scipy's SLSQP under-allocates its Fortran
 * workspace when the equality block outnumbers the variables and corrupts the
 * heap rather than refusing. This port recovers the matrices affinely for every
 * entry point and reduces the equality block before it reaches the optimizer,
 * so the distinction does not arise and the two spellings are one function
 * here. It is verified, not assumed: `test_mapqn_qrf_noblo.cpp` checks this
 * against the reference.
 *
 * Until 2026-08-29 the MATLAB reference called its own mem() here, with mmi()
 * surviving only in a commented-out line, and this port mirrored that: the
 * entry point named for mutual-information minimisation returned an entropy
 * extremum. The objective is now MMI in all four ports. MMI is not convex, so
 * unlike the MEM it replaces this entry point has no unique optimum and agrees
 * with the reference only where the polytope pins the answer.
 */
template <class T>
QrfMetrics<T> qrf_noblo_mmi_linear(const std::vector<std::pair<Matrix<T>, Matrix<T>>>& MAPs,
                                   std::size_t N, const Matrix<T>& rt, const Matrix<T>& alpha) {
    const std::size_t M = MAPs.size();
    std::vector<int> K(M, 0);
    for (std::size_t i = 0; i < M; ++i) K[i] = static_cast<int>(MAPs[i].first.rows());
    Matrix<T> mu, v;
    qrf_extract_mu_v(MAPs, M, K, &mu, &v);

    const std::size_t MR = 1;
    Matrix<int> BB(1, M, 0);
    const std::vector<int> F(M, static_cast<int>(N));
    const QrfRates<T> q = qrf_build_q_ld(M, K, mu, v, rt, N, alpha);
    const std::vector<long> idx = qrf_index_map(M, N, K, MR);
    return qrfdetail::run_noblo<T>(
        q, M, MR, BB, F, N, K,
        [&](const std::vector<T>& x) { return mmi_objective(x, M, N, K, F, MR); },
        [&](const std::vector<T>& x) { return mmi_gradient(x, M, N, K, F, MR, idx); },
        "qrf_noblo_mmi_linear", &alpha);
}

}  // namespace mapqn
}  // namespace line

#endif  // LINE_API_MAPQN_MAPQN_QRF_NOBLO_H
