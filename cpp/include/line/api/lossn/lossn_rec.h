#ifndef LINE_API_LOSSN_LOSSN_REC_H
#define LINE_API_LOSSN_LOSSN_REC_H

/**
 * @file lossn_rec.h
 * @brief Exact analysis of a loss network by MDD-rec.
 *
 * The normalising constant is the sum of a product form over the admissible set
 * {n >= 0 : A n <= C}, which is what a decision diagram holding that set
 * computes in one memoised walk.
 *
 * A Kelly loss network carries offered load nu_r on route r and admits a call
 * only while the resource constraint A n <= C still holds after it. The
 * stationary law is the truncation of independent Poisson counts to that set,
 *
 *   P(n) = (1/G) prod_r nu_r^{n_r} / n_r!,   G = sum_{A n <= C} prod_r ...,
 *
 * so g_r(k) = nu_r^k/k! and `mdd_rec` returns G. By PASTA the acceptance
 * probability of a class-r call is the ratio of two such constants,
 *
 *   1 - B_r = G(C - A e_r) / G(C),
 *
 * which is one further diagram per class.
 *
 * WHY THIS EXISTS ALONGSIDE `lossn_manjunath` AND `lossn_erlangfp`. The
 * Manjunath-Sikdar transform evaluates G exactly as a multidimensional residue,
 * and the residue argument counts WHOLE UNITS: it needs an integral A and C.
 * This port's `lossn_erlangfp` needs integrality too, for its own reason -- it
 * raises (1-E_i) to an unsigned integer power -- so before MDD-rec a FRACTIONAL
 * region had no route here at all except the Monte Carlo `lossn_mci`, whose
 * answer is a random variable. MDD-rec needs only that the admissible set be
 * finite and bounded coordinate by coordinate, which a fractional constraint
 * still is, so it is exact there too and is the default the fractional case now
 * takes.
 *
 * ARITHMETIC. Everything here is a sum, a product and one factorial, so the
 * whole method is rational and available under exact arithmetic; only the
 * reported `lG` needs a logarithm, and it is a `double` diagnostic rather than a
 * `T`, exactly as the other analyzers treat it.
 *
 * References:
 *   F. P. Kelly, "Loss networks", Annals of Applied Probability 1(3), 1991.
 *   S. Balsamo, A. Marin, I. Stojic, "Computation of the normalising constant
 *   for product-form models of distributed systems with synchronisation",
 *   Future Generation Computer Systems 111 (2020) 475-490.
 *
 * @see lossn_manjunath, lossn_erlangfp, lossn_mci, mdd_rec
 */

#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "line/api/mdd/mdd_rec.h"
#include "line/api/mdd/mdd_reachset.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/num/number.h"

namespace line {
namespace lossn {

/** Carried load, blocking, log normalising constant and walk count. */
template <class T>
struct LossnRecResult {
    /** Mean number of class-r calls in progress, the carried load. */
    std::vector<T> QLen;
    /** Blocking probability per class. */
    std::vector<T> Loss;
    /** Normalising constant G(C). */
    T G;
    /** log G(C), a double diagnostic. */
    double lG = 0.0;
    /** Number of diagram walks performed, K + 1. */
    int iterations = 0;
};

namespace detail {

/**
 * The admissible set {n >= 0 : A n <= C}, generated one call at a time from the
 * empty network. Adding a call is the only move, so the breadth-first closure
 * visits exactly the admissible vectors.
 */
template <class T>
mdd::MddStruct lossn_rec_diagram(const Matrix<T>& A, const std::vector<T>& C,
                                 const std::vector<int>& bound) {
    const std::size_t K = bound.size();
    const std::size_t J = C.size();
    std::vector<int> domain(K);
    for (std::size_t r = 0; r < K; ++r) domain[r] = bound[r] + 1;

    const mdd::MddNextState nextfun = [&A, &C, &bound, K, J](const std::vector<int>& s) {
        std::vector<std::vector<int>> out;
        for (std::size_t r = 0; r < K; ++r) {
            if (s[r] >= bound[r]) continue;
            std::vector<int> t = s;
            ++t[r];
            bool ok = true;
            for (std::size_t j = 0; j < J && ok; ++j) {
                T sum = num_traits<T>::from_int(0);
                for (std::size_t q = 0; q < K; ++q)
                    sum += T(A(j, q) * num_traits<T>::from_int(t[q]));
                ok = !(sum > C[j]);
            }
            if (ok) out.push_back(t);
        }
        return out;
    };
    mdd::MDD diagram = mdd::mdd_reachset(domain, std::vector<int>(K, 0), nextfun);
    return diagram.to_struct();
}

/**
 * G over the admissible set at capacity C, keeping the per-class domains of the
 * FULL problem so that one set of factors g serves every reduced capacity.
 */
template <class T>
T lossn_rec_G(const Matrix<T>& A, const std::vector<T>& C, const std::vector<int>& bound,
              const std::vector<std::vector<T>>& g) {
    const T zero = num_traits<T>::from_int(0);
    for (std::size_t j = 0; j < C.size(); ++j)
        if (C[j] < zero) return zero;
    return mdd::mdd_rec<T>(lossn_rec_diagram<T>(A, C, bound), g);
}

}  // namespace detail

/**
 * Exact loss-network analysis by MDD-rec.
 *
 * @param nu offered load per class, length K
 * @param A  J x K non-negative resource requirement matrix
 * @param C  capacity vector, length J
 * @return the carried load, the blocking probabilities and the normalising constant
 */
template <class T>
LossnRecResult<T> lossn_rec(const std::vector<T>& nu, const Matrix<T>& A,
                            const std::vector<T>& C) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t K = nu.size();
    const std::size_t J = C.size();
    if (A.cols() != K)
        throw InputError("lossn_rec: A has " + std::to_string(A.cols()) +
                         " columns but there are " + std::to_string(K) + " classes");
    if (A.rows() != J)
        throw InputError("lossn_rec: A has " + std::to_string(A.rows()) + " rows but C has " +
                         std::to_string(J) + " entries");
    for (std::size_t j = 0; j < J; ++j)
        for (std::size_t r = 0; r < K; ++r)
            if (A(j, r) < zero)
                throw InputError("lossn_rec: the resource matrix A must be non-negative");

    // ---- per-class bound: the most calls the tightest constraint alone admits
    std::vector<int> bound(K, 0);
    for (std::size_t r = 0; r < K; ++r) {
        double b = std::numeric_limits<double>::infinity();
        for (std::size_t j = 0; j < J; ++j)
            if (A(j, r) > zero)
                b = std::min(b, std::floor(num_traits<T>::to_double(C[j]) /
                                           num_traits<T>::to_double(A(j, r))));
        if (!std::isfinite(b))
            throw InputError("lossn_rec: class " + std::to_string(r + 1) +
                             " consumes no resource, so the admissible set is unbounded in that "
                             "coordinate and its normalising constant diverges");
        bound[r] = static_cast<int>(std::max(0.0, b));
    }

    std::vector<std::vector<T>> g(K);
    for (std::size_t r = 0; r < K; ++r) {
        g[r].assign(bound[r] + 1, one);
        T fact = one, pw = one;
        for (int k = 0; k <= bound[r]; ++k) {
            if (k > 0) {
                fact = T(fact * num_traits<T>::from_int(k));
                pw = T(pw * nu[r]);
            }
            g[r][k] = T(pw / fact);
        }
    }

    const T G = detail::lossn_rec_G<T>(A, C, bound, g);
    if (!(G > zero))
        throw InputError("lossn_rec: the admissible set is empty: no call of any class fits "
                         "within C");

    // ---- carried load per class, from the marginals of the same diagram
    const mdd::MddStruct mdds = detail::lossn_rec_diagram<T>(A, C, bound);
    LossnRecResult<T> out;
    out.G = G;
    out.QLen.assign(K, zero);
    out.Loss.assign(K, zero);
    for (std::size_t r = 0; r < K; ++r) {
        const std::vector<T> pk = mdd::mdd_rec_marginal<T>(mdds, g, r);
        T s = zero;
        for (std::size_t k = 0; k < pk.size(); ++k)
            s += T(num_traits<T>::from_int(static_cast<int>(k)) * pk[k] / G);
        out.QLen[r] = s;
    }

    // ---- blocking: 1 - B_r = G(C - A e_r)/G(C), Kelly's ratio, by PASTA
    for (std::size_t r = 0; r < K; ++r) {
        std::vector<T> Cr(J, zero);
        bool fits = true;
        for (std::size_t j = 0; j < J; ++j) {
            Cr[j] = T(C[j] - A(j, r));
            if (Cr[j] < zero) fits = false;
        }
        if (!fits) {
            out.Loss[r] = one;                        // the call never fits
            continue;
        }
        const T Gr = detail::lossn_rec_G<T>(A, Cr, bound, g);
        out.Loss[r] = T(one - Gr / G);
        if (out.Loss[r] < zero) out.Loss[r] = zero;
        if (out.Loss[r] > one) out.Loss[r] = one;
    }

    if constexpr (num_traits<T>::has_transcendental) {
        out.lG = std::log(num_traits<T>::to_double(G));
    } else {
        out.lG = std::log(num_traits<T>::to_double(G));
    }
    out.iterations = static_cast<int>(K) + 1;
    return out;
}

}  // namespace lossn
}  // namespace line

#endif  // LINE_API_LOSSN_LOSSN_REC_H
