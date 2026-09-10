/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MC_CTMC_TESTPF_KOLMOGOROV_H
#define LINE_API_MC_CTMC_TESTPF_KOLMOGOROV_H

/**
 * Kolmogorov reversibility criterion, used as a product-form test.
 *
 * Templated port of jar/src/main/java/jline/api/mc/Ctmc_testpf_kolmogorov.java,
 * which has no MATLAB twin. Kolmogorov's criterion states that an irreducible
 * chain is reversible iff for every cycle c0 -> c1 -> ... -> c0 the product of
 * the rates around it equals the product around the reverse cycle,
 *   prod_i q(c_i, c_i+1) = prod_i q(c_i+1, c_i),
 * both products taken on the SAME generator. A cycle whose reverse edges are
 * not all present fails outright.
 *
 * REFERENCE DEFECT (fixed here and in the JAR, 2026-08-01). The reference took
 * the reverse product on the time-reversed generator Qr instead of on Q. Since
 * Qr(a,b) = Q(b,a) pi_b / pi_a, that product is
 *   prod_i Q(c_i, c_i+1) pi_c_i / pi_c_i+1,
 * whose pi factors telescope to 1 around any cycle, so it equals the forward
 * product identically and the test returned true for every chain, reversible or
 * not. The two formulations agree exactly when the chain IS reversible, which is
 * why the defect never showed up as a wrong "false".
 *
 * COST. Every simple cycle through every edge is enumerated, so the work is
 * exponential in the number of states. That is the reference algorithm and is
 * kept: the callers apply it to small chains only.
 *
 * ARITHMETIC: the tolerance 1e-6 on the relative gap is a double comparison, so
 * the verdict is a floating-point one even at Rational.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/mc/ctmc_solve.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mc {

namespace detail {

/**
 * Every simple path from `start` to `target` in the adjacency matrix `adj`,
 * avoiding the states already marked used. A path is returned as the list of
 * its states, `start` first and `target` last.
 */
template <class T>
void kolmogorov_paths(const Matrix<T>& adj, std::vector<bool>& used, std::vector<std::size_t>& path,
                      std::size_t start, std::size_t target,
                      std::vector<std::vector<std::size_t>>& out) {
    const T zero = num_traits<T>::from_int(0);
    const bool wasUsed = used[start];
    used[start] = true;
    path.push_back(start);
    if (start == target) {
        out.push_back(path);
    } else {
        for (std::size_t j = 0; j < adj.cols(); ++j)
            if (adj(start, j) > zero && !used[j]) kolmogorov_paths(adj, used, path, j, target, out);
    }
    path.pop_back();
    used[start] = wasUsed;
}

}  // namespace detail

/**
 * @param Qin generator of an irreducible CTMC
 * @return true when Kolmogorov's criterion holds on every simple cycle
 */
template <class T>
bool ctmc_testpf_kolmogorov(const Matrix<T>& Qin) {
    const std::size_t n = Qin.rows();
    if (Qin.cols() != n) throw InputError("ctmc_testpf_kolmogorov: generator is not square");
    const T zero = num_traits<T>::from_int(0);

    // Clip the negative off-diagonals, then reset the diagonal to close the rows.
    Matrix<T> Q = Qin;
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            if (i != j && Q(i, j) < zero) Q(i, j) = zero;
    Q = ctmc_makeinfgen(Q);

    Matrix<T> A(n, n, zero);
    const T one = num_traits<T>::from_int(1);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            if (i != j && Q(i, j) > zero) A(i, j) = one;

    for (std::size_t start = 0; start < n; ++start) {
        for (std::size_t target = 0; target < n; ++target) {
            if (target == start || !(A(target, start) > zero)) continue;
            std::vector<bool> used(n, false);
            std::vector<std::size_t> path;
            std::vector<std::vector<std::size_t>> cycles;
            detail::kolmogorov_paths(A, used, path, start, target, cycles);

            for (std::size_t c = 0; c < cycles.size(); ++c) {
                std::vector<std::size_t> cyc = cycles[c];
                cyc.push_back(start);  // close the cycle back onto its origin
                T q = one, qr = one;
                for (std::size_t i = 0; i + 1 < cyc.size(); ++i) {
                    q = T(q * Q(cyc[i], cyc[i + 1]));
                    qr = T(qr * Q(cyc[i + 1], cyc[i]));
                }
                const double qd = num_traits<T>::to_double(q);
                const double qrd = num_traits<T>::to_double(qr);
                if (qd == 0.0) continue;
                if (std::fabs(qd - qrd) / std::fabs(qd) > 1e-6) return false;
            }
        }
    }
    return true;
}

}  // namespace mc
}  // namespace line

#endif  // LINE_API_MC_CTMC_TESTPF_KOLMOGOROV_H
