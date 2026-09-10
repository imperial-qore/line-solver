/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MAP_COUNT_MOMENT_H
#define LINE_API_MAM_MAP_COUNT_MOMENT_H

/**
 * Power moments of the counts of a MAP in a window of length t.
 *
 * Templated port of matlab/lib/kpctoolbox/map/map_count_moment.m. The moment
 * generating function of the number of arrivals N(t) is
 *
 *   M(z) = theta exp(D0 t + e^z D1 t) e,   E[N(t)^k] = d^k M / dz^k |_{z=0},
 *
 * with theta the stationary phase vector.
 *
 * DIVERGENCE, deliberately: MATLAB evaluates those derivatives by NUMERICAL
 * differentiation (derivest, Richardson extrapolation) for orders up to 4 and
 * by symbolic differentiation beyond, so its accuracy degrades quickly with
 * the order and it needs the Symbolic Toolbox for order 5 and above. This port
 * takes the derivatives analytically, at the cost of one larger exponential.
 * Write the exponent as a polynomial in z,
 *
 *   B(z) = Q t + sum_{j>=1} (D1 t / j!) z^j,   Q = D0 + D1,
 *
 * truncated at z^K. Polynomials in z modulo z^(K+1) are represented faithfully
 * by block upper-triangular Toeplitz matrices, G[i][i+j] = B_j, and that
 * representation is a ring homomorphism, so exp(G) is the representation of
 * exp(B(z)) mod z^(K+1): its (0,j) block is exactly the z^j Taylor
 * coefficient, i.e. the j-th derivative divided by j!. Hence
 *
 *   E[N(t)^k] = k! theta [exp(G)]_{0,k} e,
 *
 * accurate to the tolerance of the exponential at every order, with no step
 * size to choose and no symbolic algebra. The cost is one exponential of an
 * (K+1)n square matrix.
 *
 * ARITHMETIC: transcendental (one matrix exponential).
 *
 * The JAR has no counterpart of this function.
 */

#include <cstddef>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/expm.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/**
 * @param m      the MAP (D0, D1)
 * @param t      window length
 * @param orders orders of the moments to compute (0 returns 1)
 * @return E[N(t)^k] for each requested order, in the order of orders
 */
template <class T>
std::vector<T> map_count_moment(const Map<T>& m, const T& t, const std::vector<unsigned>& orders) {
    static_assert(num_traits<T>::has_transcendental,
                  "map_count_moment requires transcendental arithmetic");
    if (t < num_traits<T>::from_int(0)) throw InputError("map_count_moment: negative window length");
    unsigned K = 0;
    for (std::size_t i = 0; i < orders.size(); ++i)
        if (orders[i] > K) K = orders[i];

    const std::size_t n = m.order();
    const std::vector<T> theta = map_prob(m);
    const Matrix<T> Q = map_infgen(m);

    // Block Toeplitz representation of B(z) = Qt + sum_j (D1 t / j!) z^j.
    const std::size_t N = (K + 1) * n;
    Matrix<T> G(N, N, num_traits<T>::from_int(0));
    for (std::size_t bi = 0; bi <= K; ++bi) {
        for (std::size_t bj = bi; bj <= K; ++bj) {
            const std::size_t j = bj - bi;
            const T fac = j == 0 ? num_traits<T>::from_int(1)
                                 : num_traits<T>::from_int(1) / num_factorial<T>(static_cast<unsigned>(j));
            for (std::size_t r = 0; r < n; ++r)
                for (std::size_t c = 0; c < n; ++c) {
                    const T coeff = j == 0 ? T(Q(r, c)) : T(m.D1(r, c) * fac);
                    G(bi * n + r, bj * n + c) = coeff * t;
                }
        }
    }
    const Matrix<T> E = expm(G);

    std::vector<T> out;
    out.reserve(orders.size());
    for (std::size_t i = 0; i < orders.size(); ++i) {
        const unsigned k = orders[i];
        T s = num_traits<T>::from_int(0);
        for (std::size_t r = 0; r < n; ++r)
            for (std::size_t c = 0; c < n; ++c) s += theta[r] * E(r, k * n + c);
        out.push_back(num_factorial<T>(k) * s);
    }
    return out;
}


/**
 * Per-class counting moments of a marked MAP, `mmap_count_moment`.
 *
 * Port of jar/src/main/java/jline/api/mam/Mmap_count_moment.java. Class c's own
 * counting process is the MAP whose arrivals are c's alone and whose hidden
 * transitions absorb every other class,
 *
 *   D0' = D0 + sum_{j != c} D1_j,   D1' = D1_c,
 *
 * so the per-class moments are `map_count_moment` on that marginal. That
 * marginalization is exact -- an arrival of another class IS a hidden phase
 * transition as far as class c's counter is concerned -- and it is why the
 * per-class counts are NOT independent: they share the phase process.
 *
 * @param m      the marked MAP
 * @param t      window length
 * @param orders moment orders
 * @return (orders x classes) matrix of counting moments
 */
template <class T>
Matrix<T> mmap_count_moment(const Mmap<T>& m, const T& t, const std::vector<unsigned>& orders) {
    const std::size_t K = m.classes();
    if (K == 0) throw InputError("mmap_count_moment: the MMAP has no classes");
    Matrix<T> out(orders.size(), K, num_traits<T>::from_int(0));
    for (std::size_t c = 0; c < K; ++c) {
        Map<T> marg;
        marg.D0 = m.D0;
        for (std::size_t j = 0; j < K; ++j) {
            if (j == c) continue;
            for (std::size_t a = 0; a < marg.D0.rows(); ++a)
                for (std::size_t b = 0; b < marg.D0.cols(); ++b) marg.D0(a, b) += m.Dc[j](a, b);
        }
        marg.D1 = m.Dc[c];
        const std::vector<T> mm = map_count_moment(marg, t, orders);
        for (std::size_t i = 0; i < orders.size(); ++i) out(i, c) = mm[i];
    }
    return out;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MAP_COUNT_MOMENT_H
