/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_RETRIEVAL_RETRIEVAL_METRICS_H
#define LINE_API_RETRIEVAL_RETRIEVAL_METRICS_H

/**
 * Exact miss, hit and delayed-hit metrics of a delayed-hit (list-based) cache.
 *
 * Templated port of matlab/src/api/retrieval/retrieval_metrics.m,
 * cross-checked against
 * jar/src/main/java/jline/api/retrieval/Retrieval_metrics.java.
 *
 * With E(m) = E(0,m) the normalizing constant of the whole system and E_i the
 * constant of the system without item i (both from retrieval_nc):
 *
 *   miss ratio           pi_{i,0} = E_i(m)/E(m)
 *   hit ratio at list j  pi_{i,j} = m_j gamma_{i,j} E_i(m-1_j)/E(m)
 *   delayed hit at IS    phi_{0,i} = lambda_i eta_{0,i} E_i(m)/E(m)
 *   delayed hit at PS s  phi_{s,i} = lambda_i eta_{s,i} E_i(1_s,m)/E(m)
 *
 * and every item satisfies the balance
 * pi_{i,0} + sum_s phi_{s,i} + sum_j pi_{i,j} = 1.
 *
 * ARITHMETIC: ratios of exact normalizing constants, so a finite field
 * computation. In the exact instantiation the balance above holds as an
 * identity between rationals, not to within a tolerance.
 */

#include <cstddef>
#include <vector>

#include "line/api/retrieval/retrieval_nc.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace retrieval {

/** Mirrors the [pmiss, phit, pdh] return list of the MATLAB function. */
template <class T>
struct RetrievalMetricsResult {
    std::vector<T> pmiss;  ///< (n) miss ratios pi_{i,0}
    Matrix<T> phit;        ///< (h x n) hit ratios pi_{i,j}
    Matrix<T> pdh;         ///< ((r+1) x n) delayed-hit probabilities phi_{s,i}, s = 0..r
};

namespace detail {

/** Row i deleted from an (n x c) matrix, i.e. MATLAB's A([1:i-1,i+1:n],:). */
template <class T>
Matrix<T> drop_row(const Matrix<T>& A, std::size_t i) {
    Matrix<T> B(A.rows() - 1, A.cols());
    std::size_t r = 0;
    for (std::size_t a = 0; a < A.rows(); ++a) {
        if (a == i) continue;
        for (std::size_t b = 0; b < A.cols(); ++b) B(r, b) = A(a, b);
        ++r;
    }
    return B;
}

}  // namespace detail

/**
 * @param m      (h) cache list capacities
 * @param lambda (n) per-item arrival rates
 * @param eta    (n x (r+1)) fetching demands, column 0 = IS station, columns 1..r = PS
 * @param gamma  (n x h) access factors
 */
template <class T>
RetrievalMetricsResult<T> retrieval_metrics(const std::vector<int>& m, const std::vector<T>& lambda,
                                            const Matrix<T>& eta, const Matrix<T>& gamma) {
    const std::size_t n = lambda.size();
    const std::size_t h = m.size();
    if (eta.rows() != n || gamma.rows() != n)
        throw InputError("retrieval_metrics: eta/gamma and lambda disagree on the item count");
    if (eta.cols() == 0) throw InputError("retrieval_metrics: eta has no columns");
    const std::size_t r = eta.cols() - 1;
    const std::vector<int> v0(r, 0);

    const T E = retrieval_nc(v0, m, lambda, eta, gamma);
    if (E == num_traits<T>::from_int(0))
        throw NumericError("retrieval_metrics: the normalizing constant is zero");

    RetrievalMetricsResult<T> out;
    out.pmiss.assign(n, num_traits<T>::from_int(0));
    out.phit = Matrix<T>(h, n, num_traits<T>::from_int(0));
    out.pdh = Matrix<T>(r + 1, n, num_traits<T>::from_int(0));

    for (std::size_t i = 0; i < n; ++i) {
        std::vector<T> lambda_i;
        lambda_i.reserve(n - 1);
        for (std::size_t k = 0; k < n; ++k)
            if (k != i) lambda_i.push_back(lambda[k]);
        const Matrix<T> eta_i = detail::drop_row(eta, i);
        const Matrix<T> gamma_i = detail::drop_row(gamma, i);

        const T Ei = retrieval_nc(v0, m, lambda_i, eta_i, gamma_i);
        out.pmiss[i] = Ei / E;
        out.pdh(0, i) = lambda[i] * eta(i, 0) * Ei / E;

        for (std::size_t s = 0; s < r; ++s) {
            std::vector<int> vs(r, 0);
            vs[s] = 1;
            const T Eis = retrieval_nc(vs, m, lambda_i, eta_i, gamma_i);
            out.pdh(s + 1, i) = lambda[i] * eta(i, s + 1) * Eis / E;
        }

        for (std::size_t j = 0; j < h; ++j) {
            if (m[j] > 0) {
                std::vector<int> mj = m;
                mj[j] -= 1;
                const T Eij = retrieval_nc(v0, mj, lambda_i, eta_i, gamma_i);
                out.phit(j, i) =
                    num_traits<T>::from_int(static_cast<long>(m[j])) * gamma(i, j) * Eij / E;
            }
        }
    }
    return out;
}

}  // namespace retrieval
}  // namespace line

#endif  // LINE_API_RETRIEVAL_RETRIEVAL_METRICS_H
