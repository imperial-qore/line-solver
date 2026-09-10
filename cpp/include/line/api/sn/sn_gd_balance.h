/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SN_SN_GD_BALANCE_H
#define LINE_API_SN_SN_GD_BALANCE_H

/**
 * Whittle balance check for a globally state-dependent rate scaling.
 *
 * For every state n of the lattice 0..cutoffs and every pair of stations (s,t)
 * populated in n, the balance property requires
 *
 *   phi_s(n) phi_t(n - e_s) = phi_t(n) phi_s(n - e_t).
 *
 * When it holds, the chain is reversible with pi(n) ~ Phi(n) prod rho^n for the
 * balance function Phi implied by phi, and the stationary law is insensitive to
 * the service-time distribution beyond its mean. When it fails the model is
 * still solvable by SolverCTMC, but it has no product form and IS sensitive --
 * which is exactly the distinction this routine exists to make checkable, since
 * nothing in a `set_global_dependence` declaration announces it.
 *
 * Twin of matlab/src/api/sn/sn_gd_balance.m, python
 * line_solver.api.sn.sn_gd_balance and jline.api.sn.SnGdBalance.
 *
 * Reference: P. Whittle, "Partial balance and insensitivity", J. Appl. Prob.
 * 22(1), 1985; T. Bonald, A. Proutiere, "Insensitivity in processor-sharing
 * networks", Perf. Eval. 49, 2002.
 */

#include <cmath>
#include <cstddef>
#include <functional>
#include <vector>

#include "line/util/error.h"
#include "line/num/number.h"

namespace line {
namespace sn {

/**
 * Worst relative violation of the balance property over the given lattice.
 *
 * @param phi     scaling evaluated on an (nstations) population vector,
 *                returning one scalar (broadcast) or one entry per station
 * @param cutoffs per-station lattice bound, one entry per station
 * @return the worst relative violation; 0 to rounding when phi is balanced
 */
template <class T>
T sn_gd_balance(const std::function<std::vector<T>(const std::vector<T>&)>& phi,
                const std::vector<std::size_t>& cutoffs) {
    if (!phi) throw InputError("sn_gd_balance: phi must be callable");
    const std::size_t S = cutoffs.size();
    if (S < 2)
        throw InputError(
            "sn_gd_balance: cutoffs must have one entry per station (at least two stations are "
            "needed for a balance pair)");

    std::size_t total = 1;
    for (std::size_t s = 0; s < S; ++s) total *= cutoffs[s] + 1;

    const T zero = num_traits<T>::from_int(0);
    T viol = zero;
    std::vector<T> n(S, zero);
    const std::function<std::vector<T>(const std::vector<T>&)> ev =
        [&phi, S](const std::vector<T>& x) {
            std::vector<T> v = phi(x);
            if (v.size() == 1) return std::vector<T>(S, v[0]);
            if (v.size() != S)
                throw InputError("sn_gd_balance: phi must return a scalar or one entry per station");
            return v;
        };

    for (std::size_t idx = 0; idx < total; ++idx) {
        std::size_t rem = idx;
        for (std::size_t s = 0; s < S; ++s) {
            n[s] = num_traits<T>::from_int(static_cast<long>(rem % (cutoffs[s] + 1)));
            rem /= cutoffs[s] + 1;
        }
        for (std::size_t s = 0; s < S; ++s) {
            if (num_traits<T>::to_double(n[s]) == 0) continue;
            for (std::size_t t = s + 1; t < S; ++t) {
                if (num_traits<T>::to_double(n[t]) == 0) continue;
                const std::vector<T> xn = ev(n);
                std::vector<T> m = n;
                m[s] = T(m[s] - num_traits<T>::from_int(1));
                const std::vector<T> xs = ev(m);
                m = n;
                m[t] = T(m[t] - num_traits<T>::from_int(1));
                const std::vector<T> xt = ev(m);
                const double lhs = num_traits<T>::to_double(T(xn[s] * xs[t]));
                const double rhs = num_traits<T>::to_double(T(xn[t] * xt[s]));
                const double scale = std::max(std::fabs(lhs), std::fabs(rhs));
                if (scale > 0) {
                    const double v = std::fabs(lhs - rhs) / scale;
                    if (v > num_traits<T>::to_double(viol)) viol = num_traits<T>::from_double(v);
                }
            }
        }
    }
    return viol;
}

}  // namespace sn
}  // namespace line

#endif  // LINE_API_SN_SN_GD_BALANCE_H
