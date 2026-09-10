/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_LOSSN_LOSSN_ERLANGFP_H
#define LINE_API_LOSSN_LOSSN_ERLANGFP_H

/**
 * Erlang fixed-point (reduced-load) approximation for a loss network.
 *
 * Templated port of matlab/src/api/lossn/lossn_erlangfp.m. Each link j carries
 * an offered load
 *   rho_j = (1/(1-E_j)) sum_r nu_r A(j,r) prod_i (1-E_i)^A(i,r)
 * and blocks it with Erlang's loss formula B(rho_j, C_j); the vector E is the
 * fixed point of that map, reached by the shared damped iteration in
 * line::da::da_fpi. Carried traffic and per-class loss follow from E.
 *
 * MATLAB evaluates Erlang B through logs and exponentials to keep the
 * factorials in range. The port keeps that form, and consequently the whole
 * function is transcendental-gated: the fixed point is only defined to within
 * the iteration tolerance anyway, so an exact instantiation would promise more
 * than the algorithm delivers.
 *
 * The recursive form B_k = rho B_{k-1} / (k + rho B_{k-1}) IS rational, and a
 * future exact variant of the blocking formula alone could use it; the fixed
 * point around it would still be inexact.
 */

#include <cmath>
#include <cstddef>
#include <functional>
#include <vector>

#include "line/api/da/da_fpi.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace lossn {

template <class T>
struct ErlangFpResult {
    std::vector<T> QLen;  ///< carried traffic per class
    std::vector<T> Loss;  ///< blocking probability per class
    std::vector<T> E;     ///< per-link blocking probabilities (the fixed point)
    std::size_t iterations = 0;
    bool converged = false;
};

/**
 * Erlang's loss formula B(nu, C), evaluated through logs as MATLAB does so the
 * factorials stay in range for large C.
 */
template <class T>
T erlang_b(const T& nu, int C) {
    static_assert(num_traits<T>::has_transcendental,
                  "erlang_b as ported evaluates through log/exp; use the rational recursion "
                  "B_k = nu B_{k-1} / (k + nu B_{k-1}) for an exact variant");
    if (C < 0) throw InputError("erlang_b: negative capacity");
    using std::exp;
    using std::log;
    const T lnu = log(nu);
    T den = num_traits<T>::from_int(0);
    for (int i = 0; i <= C; ++i)
        den += exp(num_traits<T>::from_int(i) * lnu - log(num_factorial<T>(static_cast<unsigned>(i))));
    const T lb = num_traits<T>::from_int(C) * lnu -
                 log(num_factorial<T>(static_cast<unsigned>(C))) - log(den);
    return exp(lb);
}

/**
 * @param nu offered load per class (R)
 * @param A  (J x R) route matrix: A(j,r) is the number of circuits class r
 *           takes on link j
 * @param C  (J) link capacities
 * @param options fixed-point options (tolerance, iteration cap, damping)
 */
template <class T>
ErlangFpResult<T> lossn_erlangfp(const std::vector<T>& nu, const Matrix<T>& A,
                                 const std::vector<int>& C,
                                 const da::FpiOptions& options = da::FpiOptions()) {
    static_assert(num_traits<T>::has_transcendental,
                  "lossn_erlangfp requires transcendental arithmetic (Erlang B through logs, and "
                  "a tolerance-driven fixed point)");
    const std::size_t R = nu.size();
    const std::size_t J = C.size();
    if (A.rows() != J || A.cols() != R)
        throw InputError("lossn_erlangfp: the route matrix does not match nu and C");

    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);

    // One sweep of the reduced-load map, in the (xnew, xref) form da_fpi wants.
    const std::function<std::pair<std::vector<T>, std::vector<T>>(const std::vector<T>&, std::size_t)>
        sweep = [&](const std::vector<T>& Eprev, std::size_t) {
            std::vector<T> Enew = Eprev;
            for (std::size_t j = 0; j < J; ++j) {
                T rhoj = zero;
                for (std::size_t r = 0; r < R; ++r) {
                    if (A(j, r) <= zero) continue;
                    T term = nu[r] * A(j, r);
                    for (std::size_t i = 0; i < J; ++i) {
                        if (A(i, r) <= zero) continue;
                        const long e = static_cast<long>(num_traits<T>::to_double(A(i, r)));
                        term *= num_pow_int(T(one - Eprev[i]), static_cast<unsigned>(e));
                    }
                    rhoj += term;
                }
                const T avail = one - Eprev[j];
                if (avail == zero) throw NumericError("lossn_erlangfp: link blocked with probability 1");
                rhoj /= avail;
                Enew[j] = erlang_b(rhoj, C[j]);
            }
            return std::make_pair(Enew, Eprev);
        };

    std::vector<T> E0(J, num_traits<T>::from_rational(1, 2));
    da::FpiResult<T> fp = da::da_fpi<T>(sweep, E0, options);

    ErlangFpResult<T> r;
    r.E = fp.x;
    r.iterations = fp.iterations;
    r.converged = fp.converged;
    r.QLen.assign(R, zero);
    r.Loss.assign(R, zero);
    for (std::size_t cls = 0; cls < R; ++cls) {
        T q = nu[cls];
        for (std::size_t j = 0; j < J; ++j) {
            if (A(j, cls) <= zero) continue;
            const long e = static_cast<long>(num_traits<T>::to_double(A(j, cls)));
            q *= num_pow_int(T(one - r.E[j]), static_cast<unsigned>(e));
        }
        r.QLen[cls] = q;
        r.Loss[cls] = (nu[cls] == zero) ? zero : one - q / nu[cls];
    }
    return r;
}

}  // namespace lossn
}  // namespace line

#endif  // LINE_API_LOSSN_LOSSN_ERLANGFP_H
