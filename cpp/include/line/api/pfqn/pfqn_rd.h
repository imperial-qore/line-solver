/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_RD_H
#define LINE_API_PFQN_RD_H

/**
 * Reduction heuristic (RD) for the normalizing constant of a closed
 * LOAD-DEPENDENT product-form network.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_rd.m. MATLAB is the ONLY usable
 * reference here: jar/src/main/java/jline/api/pfqn/nc/Pfqn_rd.java carries four
 * recorded defects (the load-independent flag is not reset per station, the
 * rate comparison runs over the wrong axis, the demand division skips class 0,
 * and the rate matrix is not reset), so it was not used to adjudicate anything
 * in this port.
 *
 * Method. Each station's rate profile mu_i(k) is split into a constant part
 * and a residual. Let s_i be the first population at which the rate reaches its
 * terminal value mu_i(sum N), so that mu_i(k) = mu_i(s_i) for k >= s_i. The
 * demands are rescaled by that terminal rate, y = L / mu_i(s_i), which turns
 * the station load independent above s_i, and the residual profile
 *
 *   gamma_i(k) = mu_i(k) / mu_i(s_i),
 *   beta_i(1)  = gamma_i(1) / (1 - gamma_i(1)),
 *   beta_i(j)  = (1 - gamma_i(j-1)) gamma_i(j) / (1 - gamma_i(j))
 *
 * carries what the rescaling threw away. The heuristic then writes
 *
 *   G(N) = G_LI(y, N, Z) * Cgamma,
 *   Cgamma = sum_{v=0}^{vmax} (sum(N) - max(0, v-1))/sum(N) * E_v,
 *
 * with E_v the single-class load-dependent constant pfqn_lldsingle(rho, v,
 * beta) evaluated at the single-class utilizations rho = y X, X the exact MVA
 * throughput of the rescaled model, and vmax = min(sum_i (s_i - 1), sum(N))
 * over the stations that are genuinely load dependent. A station whose rate is
 * already constant is folded into the demands up front and contributes
 * nothing.
 *
 * Reference behaviour preserved verbatim: a not-a-number rate becomes
 * infinite, an infinite terminal rate pushes s_i back to the last finite
 * column, a not-a-number beta becomes infinite (this is the 0 * Inf that every
 * column past s_i produces), and an all-infinite beta means the residual is
 * empty and the plain load-independent constant is returned unchanged. The
 * reference also relies on MATLAB auto-growing lEN so that its first entry is
 * zero, i.e. E_0 = 1; that is written explicitly here. THAT RELIANCE HOLDS ONLY
 * FOR vmax >= 1, where the first assignment lands at index 2 and index 1 is
 * filled with zero as a side effect. At vmax = 0 the assigning loop never runs,
 * so nothing auto-grows and the read raises "Unrecognized function or variable
 * 'lEN'"; the explicit zeros here are what make that case defined, and the
 * reference has been corrected to preallocate the same way.
 *
 * pfqn_lldsingle is called on beta, which can be negative, so its linear
 * (non-logarithmic) recursion is the one that applies; the reference wraps the
 * call in real() for exactly that reason. The port has only the linear
 * recursion, so no wrapper is needed.
 *
 * Arithmetic: INEXACT BY CONSTRUCTION. This is a heuristic reduction, not an
 * identity: Cgamma is a truncated and reweighted correction series, the
 * comparison that locates s_i is against a tolerance, and the result is
 * reported as a log. Gated on has_transcendental accordingly.
 */

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "line/api/pfqn/pfqn_lldsingle.h"
#include "line/api/pfqn/pfqn_mva.h"
#include "line/api/pfqn/pfqn_nc.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Return value of pfqn_rd, mirroring [lGN, Cgamma]. */
template <class T>
struct RdResult {
    double lGN;  ///< log of the normalizing constant
    T Cgamma;    ///< the correction factor, one when the residual is empty
};

/**
 * @param L0     (M x R) service demands
 * @param N      (R) population per class
 * @param Z      (K x R) think times, summed over rows; empty for none
 * @param mu0    (M x >= sum N) load-dependent rates
 * @param tol    tolerance used to locate the terminal rate (reference 1e-6)
 * @param method the load-independent constant algorithm to reduce to
 */
template <class T>
RdResult<T> pfqn_rd(const Matrix<T>& L0, const std::vector<int>& N, const Matrix<T>& Z,
                    const Matrix<T>& mu0, double tol, NcMethod method) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_rd requires transcendental arithmetic: it is a heuristic reduction whose "
                  "correction series is truncated and reweighted, it locates the terminal rate by "
                  "a tolerance comparison, and it reports a logarithm");

    const std::size_t M = L0.rows();
    const std::size_t R = N.size();
    if (!L0.empty() && L0.cols() != R)
        throw InputError("pfqn_rd: L and N disagree on the class count");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T inf = num_traits<T>::from_double(std::numeric_limits<double>::infinity());
    const T tolT = num_traits<T>::from_double(tol);
    const auto is_nan = [](const T& v) { return !(v == v); };
    const auto is_finite = [&](const T& v) { return !is_nan(v) && v < inf && v > T(-inf); };

    long Ntot = 0;
    for (int n : N) Ntot += n;

    RdResult<T> res;
    res.Cgamma = one;
    if (Ntot < 0) {
        res.lGN = -std::numeric_limits<double>::infinity();
        return res;
    }

    Matrix<T> L = L0;
    Matrix<T> mu = mu0;
    // ---- fold the load-independent stations into the demands ---------------
    for (std::size_t i = 0; i < M; ++i) {
        bool constant = true;
        for (std::size_t k = 1; k < mu.cols(); ++k)
            if (mu(i, k) != mu(i, 0)) constant = false;
        if (!constant) continue;
        for (std::size_t r = 0; r < R; ++r) L(i, r) = L(i, r) / mu(i, 0);
        for (std::size_t k = 0; k < mu.cols(); ++k) mu(i, k) = one;
    }
    if (Ntot == 0) {
        res.lGN = 0.0;
        return res;
    }
    const std::size_t Nt = static_cast<std::size_t>(Ntot);
    if (mu.rows() != M) throw InputError("pfqn_rd: mu has the wrong station count");
    if (mu.cols() < Nt) throw InputError("pfqn_rd: mu has fewer rate columns than the population");

    // ---- truncate the rate table and normalize the missing entries ---------
    Matrix<T> muT(M, Nt, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < Nt; ++k) muT(i, k) = is_nan(mu(i, k)) ? inf : mu(i, k);

    // ---- s_i: first column at which the terminal rate is reached -----------
    std::vector<std::size_t> s(M, Nt);
    for (std::size_t i = 0; i < M; ++i) {
        if (!is_finite(muT(i, Nt - 1))) {
            s[i] = Nt;  // 1-based sum(N)
            continue;
        }
        std::size_t found = Nt;
        for (std::size_t k = 0; k < Nt; ++k)
            if (num_abs(T(muT(i, k) - muT(i, Nt - 1))) < tolT) {
                found = k + 1;  // 1-based, as in the reference
                break;
            }
        s[i] = found;
    }

    // ---- rescale the demands by the terminal rate --------------------------
    Matrix<T> y = L;
    for (std::size_t i = 0; i < M; ++i) {
        if (!is_finite(muT(i, s[i] - 1))) {
            std::size_t lastfinite = 0;
            for (std::size_t k = 0; k < Nt; ++k)
                if (is_finite(muT(i, k))) lastfinite = k + 1;
            if (lastfinite == 0)
                throw NumericError("pfqn_rd: a station has no finite load-dependent rate");
            s[i] = lastfinite;
        }
        for (std::size_t r = 0; r < R; ++r) y(i, r) = y(i, r) / muT(i, s[i] - 1);
    }

    // ---- residual profile and its beta transform ---------------------------
    Matrix<T> gamma(M, Nt, one);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < Nt; ++k) gamma(i, k) = muT(i, k) / muT(i, s[i] - 1);

    Matrix<T> beta(M, Nt, one);
    for (std::size_t i = 0; i < M; ++i) {
        beta(i, 0) = gamma(i, 0) / (one - gamma(i, 0));
        for (std::size_t j = 1; j < Nt; ++j)
            beta(i, j) = (one - gamma(i, j - 1)) * (gamma(i, j) / (one - gamma(i, j)));
    }
    bool allInf = true;
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t j = 0; j < Nt; ++j) {
            if (is_nan(beta(i, j))) beta(i, j) = inf;
            if (beta(i, j) != inf) allInf = false;
        }

    const std::vector<T> lambda;
    if (allInf) {
        // The residual is empty: the rescaled model IS the model.
        res.lGN = pfqn_nc(lambda, L, N, Z, method, zero).lG;
        return res;
    }

    // ---- the correction series --------------------------------------------
    long vmax_l = 0;
    for (std::size_t i = 0; i < M; ++i)
        if (s[i] > 1) vmax_l += static_cast<long>(s[i]) - 1;
    if (vmax_l > Ntot) vmax_l = Ntot;
    const std::size_t vmax = static_cast<std::size_t>(vmax_l < 0 ? 0 : vmax_l);

    const MvaResult<T> Y = pfqn_mva(y, N, Matrix<T>());
    Matrix<T> rhoN(M, 1, zero);
    for (std::size_t i = 0; i < M; ++i) {
        T acc = zero;
        for (std::size_t r = 0; r < R; ++r) acc += y(i, r) * Y.XN[r];
        rhoN(i, 0) = acc;
    }

    // lEN[0] = 0, i.e. E_0 = 1: the reference gets this from MATLAB's
    // auto-growth of an unassigned array slot.
    std::vector<double> lEN(vmax + 1, 0.0);
    for (std::size_t v = 1; v <= vmax; ++v) {
        const NcResult<T> e = pfqn_lldsingle(rhoN, static_cast<int>(v), beta);
        // negative-beta real() rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
        lEN[v] = e.G < zero ? num_traits<T>::log_as_double(T(-e.G)) : e.lG;
    }

    T Cgamma = zero;
    for (std::size_t v = 0; v <= vmax; ++v) {
        const T EN = num_traits<T>::from_double(std::exp(lEN[v]));
        const long shift = v >= 1 ? static_cast<long>(v) - 1 : 0;
        Cgamma += num_traits<T>::from_int(Ntot - (shift > 0 ? shift : 0)) /
                  num_traits<T>::from_int(Ntot) * EN;
    }
    res.Cgamma = Cgamma;
    res.lGN = pfqn_nc(lambda, y, N, Z, method, zero).lG + num_traits<T>::log_as_double(Cgamma);
    return res;
}

/**
 * Reference defaults: tol 1e-6, and the exact convolution for the reduced
 * load-independent constant. The reference sets options.method = 'default',
 * whose multi-station branch in this tree dispatches to the cub / le family
 * that is not ported; 'ca' is what MATLAB's 'default' itself selects for
 * models of the size this heuristic targets, and it is exact.
 */
template <class T>
RdResult<T> pfqn_rd(const Matrix<T>& L, const std::vector<int>& N, const Matrix<T>& Z,
                    const Matrix<T>& mu) {
    return pfqn_rd(L, N, Z, mu, 1e-6, NcMethod::Ca);
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_RD_H
