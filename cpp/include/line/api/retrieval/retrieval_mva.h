/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_RETRIEVAL_RETRIEVAL_MVA_H
#define LINE_API_RETRIEVAL_RETRIEVAL_MVA_H

/**
 * Exact MVA-style recursion for delayed-hit (list-based) cache metrics.
 *
 * Templated port of matlab/src/api/retrieval/retrieval_mva.m. There is no JAR
 * counterpart: jline/api/retrieval/ ships nc, metrics, fpi and fpi_latency
 * only, so MATLAB is the sole reference for this function.
 *
 * This is the exact recursion that retrieval_fpi approximates. Writing
 * phi^{(k)} for the delayed-hit probability in the system WITHOUT item k,
 *
 *   theta_ij(m) = gamma_ij / (1 + lambda_i eta_{0,i}
 *                   + sum_s lambda_i eta_{s,i}(1 + sum_{k!=i} phi^{(i)}_{s,k}(m-1_j)))
 *   xi_j(m)     = m_j / sum_i theta_ij(m)(1 - pihit_i(m-1_j))
 *   pi_ij(m)    = theta_ij(m) xi_j(m) (1 - pihit_i(m-1_j))
 *   pi_i0(m)    = (1 - pihit_i(m)) / (1 + lambda_i eta_{0,i}
 *                   + sum_s lambda_i eta_{s,i}(1 + sum_{k!=i} phi^{(i)}_{s,k}(m)))
 *   phi_{s,k}(m)= lambda_k pi_k0(m) eta_{s,k}(1 + sum_{i!=k} phi^{(k)}_{s,i}(m))
 *   phi_{0,k}(m)= lambda_k eta_{0,k} pi_k0(m)
 *
 * memoized over (item subset, capacity vector). The recursion bottoms out at
 * the empty item set and at any capacity able to hold every remaining item,
 * where the items are permanently cached (pihit = 1, no fetching at all).
 * Cost is O(2^n n^2 h r prod_j (1+m_j)) time and memory, so it is a small-case
 * oracle: use retrieval_fpi beyond that.
 *
 * ARITHMETIC: additions, multiplications and divisions of the inputs only, so
 * a finite field computation, exact in the exact instantiation. It then agrees
 * with retrieval_metrics as an identity between rationals, which is the
 * strongest available check on either implementation.
 */

#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace retrieval {

/** Mirrors the [pmiss, phit, pdh] return list of the MATLAB function. */
template <class T>
struct RetrievalMvaResult {
    std::vector<T> pmiss;  ///< (n) miss ratios pi_{i,0}
    Matrix<T> phit;        ///< (h x n) hit ratios pi_{i,j}
    Matrix<T> pdh;         ///< ((r+1) x n) delayed-hit probabilities phi_{s,i}, s = 0..r
};

namespace detail {

/**
 * Memo tables and the recursion itself. Held in an object rather than in
 * statics so the function is reentrant and leaves no state behind.
 */
template <class T>
class RetrievalMvaSolver {
public:
    RetrievalMvaSolver(const std::vector<int>& m, const std::vector<T>& lambda,
                       const Matrix<T>& eta, const Matrix<T>& gamma)
        : m_(m), lambda_(lambda), eta_(eta), gamma_(gamma), n_(lambda.size()), h_(m.size()),
          r_(eta.cols() - 1) {
        radix_.resize(h_);
        ncap_ = 1;
        for (std::size_t j = 0; j < h_; ++j) {
            radix_[j] = static_cast<std::size_t>(m_[j]) + 1;
            ncap_ *= radix_[j];
        }
        nmask_ = static_cast<std::size_t>(1) << n_;
        const std::size_t cells = nmask_ * ncap_;
        done_.assign(cells, false);
        const T zero = num_traits<T>::from_int(0);
        pi0_.assign(cells * n_, zero);
        pihit_.assign(cells * n_, zero);
        pij_.assign(cells * n_ * h_, zero);
        phi_.assign(cells * n_ * (r_ + 1), zero);
    }

    RetrievalMvaResult<T> run() {
        const std::size_t full = nmask_ - 1;
        solve(full, m_);
        const std::size_t ci = capidx(m_);
        RetrievalMvaResult<T> out;
        out.pmiss.assign(n_, num_traits<T>::from_int(0));
        out.phit = Matrix<T>(h_, n_, num_traits<T>::from_int(0));
        out.pdh = Matrix<T>(r_ + 1, n_, num_traits<T>::from_int(0));
        for (std::size_t i = 0; i < n_; ++i) {
            out.pmiss[i] = pi0_[(full * ncap_ + ci) * n_ + i];
            for (std::size_t j = 0; j < h_; ++j)
                out.phit(j, i) = pij_[((full * ncap_ + ci) * n_ + i) * h_ + j];
            for (std::size_t s = 0; s <= r_; ++s)
                out.pdh(s, i) = phi_[((full * ncap_ + ci) * n_ + i) * (r_ + 1) + s];
        }
        return out;
    }

private:
    std::size_t capidx(const std::vector<int>& c) const {
        std::size_t idx = 0, mul = 1;
        for (std::size_t j = 0; j < h_; ++j) {
            idx += static_cast<std::size_t>(c[j]) * mul;
            mul *= radix_[j];
        }
        return idx;
    }

    T& PI0(std::size_t mask, std::size_t ci, std::size_t i) {
        return pi0_[(mask * ncap_ + ci) * n_ + i];
    }
    T& PIHIT(std::size_t mask, std::size_t ci, std::size_t i) {
        return pihit_[(mask * ncap_ + ci) * n_ + i];
    }
    T& PIJ(std::size_t mask, std::size_t ci, std::size_t i, std::size_t j) {
        return pij_[((mask * ncap_ + ci) * n_ + i) * h_ + j];
    }
    T& PHI(std::size_t mask, std::size_t ci, std::size_t i, std::size_t s) {
        return phi_[((mask * ncap_ + ci) * n_ + i) * (r_ + 1) + s];
    }

    void solve(std::size_t mask, const std::vector<int>& c) {
        const std::size_t ci = capidx(c);
        if (done_[mask * ncap_ + ci]) return;
        if (mask == 0) {
            done_[mask * ncap_ + ci] = true;
            return;
        }
        std::vector<std::size_t> active;
        for (std::size_t i = 0; i < n_; ++i)
            if (mask & (static_cast<std::size_t>(1) << i)) active.push_back(i);

        long csum = 0;
        for (int x : c) csum += x;
        if (csum >= static_cast<long>(active.size())) {
            // The cache holds every active item, so all of them are cached
            // permanently: no miss, no fetch, hit probability one.
            for (std::size_t i : active) PIHIT(mask, ci, i) = num_traits<T>::from_int(1);
            done_[mask * ncap_ + ci] = true;
            return;
        }

        const T one = num_traits<T>::from_int(1);

        // dependencies at m - 1_j, with and without each active item
        for (std::size_t j = 0; j < h_; ++j) {
            if (c[j] > 0) {
                std::vector<int> cj = c;
                cj[j] -= 1;
                solve(mask, cj);
                for (std::size_t i : active) solve(mask & ~(static_cast<std::size_t>(1) << i), cj);
            }
        }
        for (std::size_t k : active) solve(mask & ~(static_cast<std::size_t>(1) << k), c);

        // theta, xi and pi_ij, all evaluated on the cache recursion m - 1_j
        for (std::size_t j = 0; j < h_; ++j) {
            if (c[j] == 0) continue;
            std::vector<int> cj = c;
            cj[j] -= 1;
            const std::size_t cjx = capidx(cj);
            std::vector<T> theta(n_, num_traits<T>::from_int(0));
            for (std::size_t i : active) {
                const std::size_t maski = mask & ~(static_cast<std::size_t>(1) << i);
                T acc = num_traits<T>::from_int(0);
                for (std::size_t s = 0; s < r_; ++s) {
                    T sphi = num_traits<T>::from_int(0);
                    for (std::size_t k : active)
                        if (k != i) sphi += PHI(maski, cjx, k, s + 1);
                    acc += lambda_[i] * eta_(i, s + 1) * (one + sphi);
                }
                theta[i] = gamma_(i, j) / (one + lambda_[i] * eta_(i, 0) + acc);
            }
            T sden = num_traits<T>::from_int(0);
            for (std::size_t i : active) sden += theta[i] * (one - PIHIT(mask, cjx, i));
            if (sden == num_traits<T>::from_int(0))
                throw NumericError("retrieval_mva: degenerate list occupancy (zero denominator)");
            const T xi_j = num_traits<T>::from_int(static_cast<long>(c[j])) / sden;
            for (std::size_t i : active)
                PIJ(mask, ci, i, j) = theta[i] * xi_j * (one - PIHIT(mask, cjx, i));
        }

        for (std::size_t i : active) {
            T ph = num_traits<T>::from_int(0);
            for (std::size_t j = 0; j < h_; ++j) ph += PIJ(mask, ci, i, j);
            PIHIT(mask, ci, i) = ph;
        }

        for (std::size_t i : active) {
            const std::size_t maski = mask & ~(static_cast<std::size_t>(1) << i);
            T acc = num_traits<T>::from_int(0);
            for (std::size_t s = 0; s < r_; ++s) {
                T sphi = num_traits<T>::from_int(0);
                for (std::size_t k : active)
                    if (k != i) sphi += PHI(maski, ci, k, s + 1);
                acc += lambda_[i] * eta_(i, s + 1) * (one + sphi);
            }
            PI0(mask, ci, i) =
                (one - PIHIT(mask, ci, i)) / (one + lambda_[i] * eta_(i, 0) + acc);
        }

        for (std::size_t k : active) {
            const std::size_t maskk = mask & ~(static_cast<std::size_t>(1) << k);
            const T pi0k = PI0(mask, ci, k);
            PHI(mask, ci, k, 0) = lambda_[k] * eta_(k, 0) * pi0k;
            for (std::size_t s = 0; s < r_; ++s) {
                T sphi = num_traits<T>::from_int(0);
                for (std::size_t i : active)
                    if (i != k) sphi += PHI(maskk, ci, i, s + 1);
                PHI(mask, ci, k, s + 1) = lambda_[k] * pi0k * eta_(k, s + 1) * (one + sphi);
            }
        }

        done_[mask * ncap_ + ci] = true;
    }

    const std::vector<int>& m_;
    const std::vector<T>& lambda_;
    const Matrix<T>& eta_;
    const Matrix<T>& gamma_;
    std::size_t n_, h_, r_;
    std::vector<std::size_t> radix_;
    std::size_t ncap_ = 1, nmask_ = 1;
    std::vector<bool> done_;
    std::vector<T> pi0_, pihit_, pij_, phi_;
};

}  // namespace detail

/**
 * @param m      (h) cache list capacities
 * @param lambda (n) per-item arrival rates
 * @param eta    (n x (r+1)) fetching demands, column 0 = IS station, columns 1..r = PS
 * @param gamma  (n x h) access factors
 */
template <class T>
RetrievalMvaResult<T> retrieval_mva(const std::vector<int>& m, const std::vector<T>& lambda,
                                    const Matrix<T>& eta, const Matrix<T>& gamma) {
    const std::size_t n = lambda.size();
    if (eta.rows() != n || gamma.rows() != n)
        throw InputError("retrieval_mva: eta/gamma and lambda disagree on the item count");
    if (gamma.cols() != m.size())
        throw InputError("retrieval_mva: gamma and m disagree on the number of lists");
    if (eta.cols() == 0) throw InputError("retrieval_mva: eta has no columns");
    for (int x : m)
        if (x < 0) throw InputError("retrieval_mva: negative list capacity");
    if (n > 8 * sizeof(std::size_t) - 1)
        throw UnsupportedError("retrieval_mva: too many items for a bitmask subset enumeration");
    detail::RetrievalMvaSolver<T> solver(m, lambda, eta, gamma);
    return solver.run();
}

}  // namespace retrieval
}  // namespace line

#endif  // LINE_API_RETRIEVAL_RETRIEVAL_MVA_H
