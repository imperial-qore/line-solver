/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_OI_FNC_H
#define LINE_API_PFQN_PFQN_OI_FNC_H

/**
 * Order-independent (OI) functional server: the balance function Psi and the
 * rate mu_f of an auxiliary station whose insertion turns the mean of a
 * queue-dependent function into a ratio of normalizing constants.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_oi_fnc.m, the OI generalization
 * of pfqn_fnc (Casale, QEST 2006, Theorem 3 and Corollary 1). Psi is defined
 * by the lattice convolution identity
 *
 *   (Psi * Phi)(n) = (1 + f(n)) Phi(n),   (Psi * Phi)(n) = sum_{0<=k<=n} Psi(k) Phi(n-k)
 *
 * deconvolved triangularly in column-major order, and then inverted to the
 * rate by balanced fairness,
 *
 *   mu_f(n) = ( sum_{r: n_r > 0} Psi(n - e_r) ) / Psi(n).
 *
 * With that station in the model, E[f(n)] = G+/G - 1, with no probabilities
 * and no Little's law. For f(n) = sum(n) this is the exact total mean queue
 * length of the station.
 *
 * SIGNED BALANCE. Psi and mu_f may be negative or otherwise non-physical, as
 * the reference notes; this is immaterial because only the normalizing-constant
 * ratio is used. A state with Psi(n) = 0 has no defined rate and returns
 * infinity, exactly as in MATLAB.
 *
 * ARITHMETIC. Deconvolution and division only, no transcendental function
 * anywhere, so the routine is EXACT in rational arithmetic and is deliberately
 * left ungated. That matters here more than usual: the deconvolution is a
 * triangular solve with alternating signs, which is exactly the pattern where
 * floating point loses digits, and the exact instantiation is a real check on
 * a double one.
 */

#include <cmath>
#include <limits>
#include <cstddef>
#include <functional>
#include <vector>

#include "line/api/pfqn/pfqn_asympt_common.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace pfqn {

/** Return value of pfqn_oi_fnc, mirroring [muf, Psi, mu] flattened. */
template <class T>
struct OiFncResult {
    std::vector<T> Psi;              ///< column-major over the lattice
    std::vector<T> mu;               ///< column-major over the lattice
    std::vector<std::size_t> stride; ///< column-major strides, for indexing
};

/**
 * @param Phi balance function of the existing OI station, column-major over
 *            the lattice 0 <= n <= N (length prod(N+1))
 * @param N   (R) closed population vector
 * @param f   target queue-dependent function with f(0) = 0; empty for sum(n)
 */
template <class T>
OiFncResult<T> pfqn_oi_fnc(const std::vector<T>& Phi, const std::vector<int>& N,
                           const std::function<T(const std::vector<int>&)>& f) {
    const std::size_t R = N.size();
    if (R == 0) throw InputError("pfqn_oi_fnc: empty population vector");
    std::vector<std::size_t> shp(R), stride(R, 1);
    std::size_t total = 1;
    for (std::size_t d = 0; d < R; ++d) {
        if (N[d] < 0) throw InputError("pfqn_oi_fnc: negative population");
        shp[d] = static_cast<std::size_t>(N[d]) + 1;
    }
    for (std::size_t d = 1; d < R; ++d) stride[d] = stride[d - 1] * shp[d - 1];
    for (std::size_t d = 0; d < R; ++d) total *= shp[d];
    if (Phi.size() != total) throw InputError("pfqn_oi_fnc: numel(Phi) must equal prod(N+1)");

    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    std::vector<std::vector<int>> subs(total, std::vector<int>(R, 0));
    for (std::size_t i = 0; i < total; ++i) {
        std::size_t li = i;
        for (std::size_t d = 0; d < R; ++d) {
            subs[i][d] = static_cast<int>(li % shp[d]);
            li /= shp[d];
        }
    }

    // Step 1: deconvolve (Psi * Phi)(n) = (1 + f(n)) Phi(n).
    std::vector<T> Psi(total, zero);
    for (std::size_t i = 0; i < total; ++i) {
        const std::vector<int>& n = subs[i];
        const T fv = f ? f(n) : [&]() {
            T s = zero;
            for (int v : n) s += num_traits<T>::from_int(v);
            return s;
        }();
        T acc = T(T(one + fv) * Phi[i]);
        for (std::size_t j = 0; j < i; ++j) {
            const std::vector<int>& k = subs[j];
            bool le = true;
            for (std::size_t d = 0; d < R && le; ++d)
                if (k[d] > n[d]) le = false;
            if (!le) continue;
            std::size_t idx = 0;
            for (std::size_t d = 0; d < R; ++d)
                idx += static_cast<std::size_t>(n[d] - k[d]) * stride[d];
            acc -= Psi[j] * Phi[idx];
        }
        Psi[i] = acc;
    }

    // Step 2: balanced-fairness inversion to the rate.
    const T inf = detail::num_inf_marker<T>();
    std::vector<T> mu(total, inf);
    for (std::size_t i = 0; i < total; ++i) {
        const std::vector<int>& n = subs[i];
        int tot = 0;
        for (int v : n) tot += v;
        if (tot == 0) {
            mu[i] = zero;
            continue;
        }
        if (Psi[i] == zero) continue;  // non-physical / undefined rate
        T num = zero;
        for (std::size_t r = 0; r < R; ++r)
            if (n[r] > 0) num += Psi[i - stride[r]];
        mu[i] = T(num / Psi[i]);
    }

    OiFncResult<T> res;
    res.Psi = Psi;
    res.mu = mu;
    res.stride = stride;
    return res;
}

template <class T>
OiFncResult<T> pfqn_oi_fnc(const std::vector<T>& Phi, const std::vector<int>& N) {
    return pfqn_oi_fnc(Phi, N, std::function<T(const std::vector<int>&)>());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PFQN_OI_FNC_H
