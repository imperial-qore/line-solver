/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_CLWOI_H
#define LINE_API_PFQN_CLWOI_H

/**
 * Normalizing constant of a closed network of ORDER-INDEPENDENT (OI) stations
 * plus one aggregated delay, by numerical inversion of the multichain
 * generating function (Choudhury-Leung-Whitt, J. ACM 42(5):935-970, 1995).
 *
 * Templated port of matlab/src/api/pfqn/pfqn_clwoi.m. It is the OI counterpart
 * of pfqn_clw_lld and the transform counterpart of the convolution routine
 * pfqn_ncoi; both return the same G(N) and differ only in cost.
 *
 *   G(z) = exp( sum_r Z_r z_r ) prod_i F_i(z),   F_i(z) = sum_n Phi_i(n) z^n,
 *
 * with Phi_i the v-weighted balanced-fairness balance function of station i,
 * mu_i(n) Phi_i(n) = sum_{r: n_r>0} v_{i,r} Phi_i(n - e_r).
 *
 * Unlike a load-dependent station, whose factor collapses to a function of the
 * single argument sum_r rho_{ri} z_r, an OI station factor depends on the whole
 * vector z, because mu_i(n) depends on the occupancy only through its SUPPORT
 * supp(n) = {r : n_r > 0}. It is nevertheless rational and available in closed
 * form. Splitting the count lattice by support, on which mu_i(n) = mu_{i,S} is
 * constant, and writing F_{i,S} for the part of F_i carried by the states of
 * support S, the balance recursion gives
 *
 *   ( mu_{i,S} - sum_{r in S} v_{i,r} z_r ) F_{i,S}(z)
 *        = sum_{r in S} v_{i,r} z_r F_{i,S minus r}(z),   F_{i,{}} = 1,
 *   F_i(z) = sum_S F_{i,S}(z),
 *
 * since removing a class-r job from a state of support S lands on support S
 * when n_r >= 2 and on S minus r when n_r = 1. The singularities are the hyperplanes
 * sum_{r in S} v_{i,r} z_r = mu_{i,S}, one per support, in place of the single
 * pole x = c_i of the load-dependent case; a load-independent single-server
 * queue (mu_{i,S} = 1) gives back 1/(1 - sum_r v_{i,r} z_r).
 *
 * The restrictive static scaling of eqs. 5.41-5.46 is reused verbatim on the
 * EXPANDED constraint matrix that lists one row per (station, nonempty support)
 * pair with unit-pole intensities v_{i,r}/mu_{i,S}, after dropping the rows
 * dominated by a superset of no larger rate. Each surviving row is a binding
 * singular hyperplane, so the contour stays inside the domain of analyticity
 * exactly as the single-pole normalization does for pfqn_clw_lld.
 *
 * SCOPE. The rates must be support-only, mu_i(n) = mu_i(supp(n)), which is the
 * defining property of an OI station and what makes the transform a finite
 * rational function. Every rate handle is verified EXHAUSTIVELY on the count
 * lattice 0 < n <= N before the inversion, at prod_r (N_r+1) evaluations per
 * station (below the contour points spent afterwards), and a state whose rate
 * differs from that of its support is an error naming that state: a rate that
 * varies inside a support is a general balanced-fairness station and belongs to
 * pfqn_ncoi. It is refused rather than warned-and-inverted because the
 * inversion would otherwise return a plausible but wrong G(N). A non-empty swap
 * graph breaks the closure of Phi on the count vector altogether and requires
 * the microstate routine pfqn_pas_nc.
 *
 * Arithmetic: TRANSCENDENTAL, double and Real only, for the same reason as the
 * rest of the CLW family (the contour radius 10^{-gamma/(2 l K)} is not in the
 * field of the inputs). Accuracy against the exact pfqn_ncoi is ~1e-9 for two
 * chains and ~1e-8 for three, i.e. the accuracy the method itself has.
 *
 * COST. prod_r 2 l_r N_r contour points, each costing O(M R 2^R), against
 * O(M prod_r (N_r+1)(N_r+2)/2) for pfqn_ncoi: linear rather than quadratic in
 * each population, so it wins on large populations with few chains and loses as
 * the chain count grows. Unlike pfqn_ncoi it returns G at the single population
 * N; throughputs need the R additional inversions at N - e_r.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "line/api/pfqn/pfqn_clw.h"
#include "line/api/pfqn/pfqn_ncoi.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

namespace detail {

/**
 * The rate of one support, read at its 0/1 indicator. The indicator is itself a
 * lattice point of that support, every retained chain having N_r >= 1;
 * constancy over the support is verified by clwoi_check_support.
 */
template <class T>
T clwoi_support_rate(const OiRate<T>& murate, const std::vector<int>& chi) {
    const T zero = num_traits<T>::from_int(0);
    const T rate = murate(chi);
    if (!(rate > zero))
        throw InputError("pfqn_clwoi: a station has a non-positive rate on a reachable support");
    return rate;
}

/** Formats a count vector as "[n1 n2 ...]" for the diagnostics below. */
inline std::string clwoi_state_str(const std::vector<int>& n) {
    std::string s = "[";
    for (std::size_t r = 0; r < n.size(); ++r) {
        if (r) s += " ";
        s += std::to_string(n[r]);
    }
    return s + "]";
}

/**
 * Exhaustive support-only check: every state 0 < n <= N is compared against the
 * rate of its own support. The transform is exact only if mu_i is constant on
 * each support, and a rate that violates that yields a wrong G with no other
 * symptom, so it is refused rather than inverted. The scan costs
 * prod_r (N_r+1) rate evaluations per station, below the prod_r 2 l_r N_r
 * contour points the inversion itself spends (l_r >= 1 gives 2 l_r N_r >= N_r+1).
 */
template <class T>
void clwoi_check_support(const std::vector<OiRate<T>>& mu, const Matrix<T>& muS,
                          const std::vector<std::size_t>& keep, const std::vector<int>& Nk,
                          std::size_t R) {
    const std::size_t M = mu.size(), p = keep.size();
    if (M == 0) return;
    std::size_t L = 1;
    for (std::size_t j = 0; j < p; ++j) L *= static_cast<std::size_t>(Nk[j]) + 1;
    for (std::size_t i = 0; i < M; ++i) {
        std::vector<int> n(R, 0);
        for (std::size_t idx = 1; idx < L; ++idx) {  // idx 0 is the empty support, unused
            std::size_t rem = idx, mask = 0;
            for (std::size_t j = 0; j < p; ++j) {
                const std::size_t base = static_cast<std::size_t>(Nk[j]) + 1;
                const std::size_t nj = rem % base;
                rem /= base;
                n[keep[j]] = static_cast<int>(nj);
                if (nj > 0) mask |= static_cast<std::size_t>(1) << j;
            }
            const T rate = mu[i](n);
            const T ref = muS(i, mask);
            T diff = T(rate - ref);
            if (diff < num_traits<T>::from_int(0)) diff = T(-diff);
            const T scale = (ref > num_traits<T>::from_int(1)) ? ref : num_traits<T>::from_int(1);
            if (diff > num_traits<T>::from_double(1e-9) * scale) {
                std::vector<int> chi(R, 0);
                for (std::size_t j = 0; j < p; ++j)
                    if (mask & (static_cast<std::size_t>(1) << j)) chi[keep[j]] = 1;
                throw InputError(
                    "pfqn_clwoi: station " + std::to_string(i) +
                    " has a rate that varies within a support (state " + clwoi_state_str(n) +
                    " against the indicator " + clwoi_state_str(chi) +
                    "); pfqn_clwoi requires order-independent (support-only) rates, "
                    "mu(n)=mu(supp(n)), use pfqn_ncoi for a general balanced-fairness station");
            }
        }
    }
}

}  // namespace detail

/**
 * @param Z      (R) think-time demand of the aggregated delay node
 * @param N      (R) closed population, finite
 * @param mu     one rate handle per OI station, mapping a per-class count vector
 *               to the total service rate; must depend on the count vector only
 *               through its support. Empty for a pure delay network
 * @param visits (M x R) per-station class visit ratios weighting the balance
 *               recursion; empty for unit visits
 * @param opt    lattice and aliasing parameters
 */
template <class T>
ClwResult<T> pfqn_clwoi(const std::vector<T>& Z, const std::vector<int>& N,
                         const std::vector<OiRate<T>>& mu, const Matrix<T>& visits,
                         const ClwOptions& opt) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_clwoi requires transcendental arithmetic (contour integration of a "
                  "generating function)");
    using std::exp;
    using std::log;
    const std::size_t R = N.size();
    const std::size_t M = mu.size();
    if (!Z.empty() && Z.size() != R)
        throw InputError("pfqn_clwoi: Z and N disagree on the chain count");
    if (!visits.empty() && (visits.rows() != M || visits.cols() != R))
        throw InputError("pfqn_clwoi: visits must be M x R");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);

    ClwResult<T> res;
    long Ntot = 0;
    for (std::size_t r = 0; r < R; ++r) {
        if (N[r] < 0) {
            res.G = zero;
            res.lG = T(-std::numeric_limits<T>::infinity());
            return res;
        }
        Ntot += N[r];
    }
    if (Ntot == 0) {
        res.G = one;
        res.lG = zero;
        return res;
    }

    std::vector<int> lfull;
    std::vector<double> gfull;
    detail::clw_defaults(R, opt, lfull, gfull);

    // drop zero-population chains: the coefficient of z_r^0 is the generating
    // function at z_r = 0, which kills every F_{i,S} with r in S
    std::vector<std::size_t> keep;
    for (std::size_t r = 0; r < R; ++r)
        if (N[r] > 0) keep.push_back(r);
    const std::size_t p = keep.size();
    std::vector<int> Nk(p, 0), l(p, 1);
    std::vector<T> Zk(p, zero);
    std::vector<double> gam(p, 0.0);
    for (std::size_t j = 0; j < p; ++j) {
        Nk[j] = N[keep[j]];
        Zk[j] = Z.empty() ? zero : Z[keep[j]];
        l[j] = lfull[keep[j]];
        gam[j] = gfull[keep[j]];
    }
    const std::size_t nmask = static_cast<std::size_t>(1) << p;

    // support rate table mu_{i,S}, S a bitmask over the retained chains
    Matrix<T> muS(M ? M : 1, nmask, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t mask = 1; mask < nmask; ++mask) {
            std::vector<int> chi(R, 0);
            for (std::size_t b = 0; b < p; ++b)
                if (mask & (static_cast<std::size_t>(1) << b)) chi[keep[b]] = 1;
            muS(i, mask) = detail::clwoi_support_rate<T>(mu[i], chi);
        }
    detail::clwoi_check_support<T>(mu, muS, keep, Nk, R);

    // per-mask chain lists and the S minus r column indices of the recursion
    std::vector<std::vector<std::size_t>> bits(nmask), subcol(nmask);
    for (std::size_t mask = 1; mask < nmask; ++mask)
        for (std::size_t b = 0; b < p; ++b)
            if (mask & (static_cast<std::size_t>(1) << b)) {
                bits[mask].push_back(b);
                subcol[mask].push_back(mask & ~(static_cast<std::size_t>(1) << b));
            }

    // per-station visit vectors restricted to the retained chains
    Matrix<T> V(M ? M : 1, p, one);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t b = 0; b < p; ++b) V(i, b) = visits.empty() ? one : visits(i, keep[b]);

    std::vector<T> r(p, one);
    for (std::size_t j = 0; j < p; ++j)
        r[j] = detail::clw_pow_real(
            num_traits<T>::from_int(10),
            T(num_traits<T>::from_double(-gam[j]) /
              num_traits<T>::from_int(2 * static_cast<long>(l[j]) * Nk[j])));

    // One constraint row per (station, nonempty support), holding the unit-pole
    // intensities v_{i,r}/mu_{i,S} of that singular hyperplane. Support S is
    // dominated by a superset S' with mu_{i,S'} <= mu_{i,S}, since then
    // v/mu_{S'} >= v/mu_S on all of S; keeping such slack rows would perturb the
    // group averages of eq. 5.44 (for load-independent stations only the
    // full-support pole survives, reproducing the pfqn_clw_lld scaling).
    std::vector<std::vector<T>> rows;
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t mask = 1; mask < nmask; ++mask) {
            bool dominated = false;
            for (std::size_t mask2 = 1; mask2 < nmask && !dominated; ++mask2)
                if (mask2 != mask && (mask & mask2) == mask &&
                    muS(i, mask2) <= muS(i, mask) * num_traits<T>::from_double(1 + 1e-12))
                    dominated = true;
            if (dominated) continue;
            std::vector<T> row(p, zero);
            for (std::size_t t = 0; t < bits[mask].size(); ++t)
                row[bits[mask][t]] = V(i, bits[mask][t]) / muS(i, mask);
            rows.push_back(row);
        }
    const std::size_t nrow = rows.size();
    Matrix<T> Lt(nrow ? nrow : 1, p, zero);
    for (std::size_t i = 0; i < nrow; ++i)
        for (std::size_t j = 0; j < p; ++j) Lt(i, j) = rows[i][j];
    if (nrow == 0)
        for (std::size_t j = 0; j < p; ++j) Lt(0, j) = zero;  // pure delay network

    const std::vector<long> mult(nrow ? nrow : 1, 1);
    const std::vector<T> alpha = detail::clw_scaling(Lt, Lt, Nk, Zk, l, r, mult);

    std::vector<T> arho0(p, zero);
    Matrix<T> vs(M ? M : 1, p, zero);
    for (std::size_t j = 0; j < p; ++j) {
        arho0[j] = alpha[j] * Zk[j];
        for (std::size_t i = 0; i < M; ++i) vs(i, j) = V(i, j) * alpha[j];
    }

    // Gbar(w) = exp(sum_r alpha_r Z_r (w_r - 1)) prod_i F_i(alpha_r v_{i,r} w_r).
    // Summing logs is legitimate for the principal complex log because
    // exp(log a + log b) = a b.
    std::vector<detail::Cx<T>> FS(nmask, detail::Cx<T>(zero, zero));
    const auto gbar = [&](const std::vector<detail::Cx<T>>& w) {
        detail::Cx<T> expo(zero, zero);
        for (std::size_t j = 0; j < p; ++j)
            expo = detail::cx_add(
                expo, detail::cx_scale(detail::Cx<T>(T(w[j].re - one), w[j].im), arho0[j]));
        detail::Cx<T> logF(zero, zero);
        for (std::size_t i = 0; i < M; ++i) {
            FS[0] = detail::Cx<T>(one, zero);
            detail::Cx<T> tot(one, zero);
            for (std::size_t mask = 1; mask < nmask; ++mask) {
                detail::Cx<T> num(zero, zero);
                detail::Cx<T> den(muS(i, mask), zero);
                for (std::size_t t = 0; t < bits[mask].size(); ++t) {
                    const std::size_t b = bits[mask][t];
                    const detail::Cx<T> x = detail::cx_scale(w[b], vs(i, b));
                    num = detail::cx_add(num, detail::cx_mul(x, FS[subcol[mask][t]]));
                    den = detail::cx_sub(den, x);
                }
                FS[mask] = detail::cx_div(num, den);
                tot = detail::cx_add(tot, FS[mask]);
            }
            logF = detail::cx_add(logF, detail::cx_log(tot));
        }
        return detail::cx_exp(detail::cx_add(expo, logF));
    };

    std::vector<detail::Cx<T>> w(p);
    const detail::Cx<T> gv = detail::clw_invert(0, w, Nk, l, r, p, gbar);
    if (!(gv.re > zero))
        throw NumericError("pfqn_clwoi: the inverted generating function is not positive");

    T lG = log(gv.re);
    for (std::size_t j = 0; j < p; ++j)
        lG += arho0[j] - num_traits<T>::from_int(Nk[j]) * log(alpha[j]);
    res.lG = lG;
    res.G = (lG > num_traits<T>::from_int(709)) ? T(std::numeric_limits<T>::infinity()) : T(exp(lG));
    return res;
}

/** Overload with the CLW default parameters. */
template <class T>
ClwResult<T> pfqn_clwoi(const std::vector<T>& Z, const std::vector<int>& N,
                         const std::vector<OiRate<T>>& mu, const Matrix<T>& visits) {
    return pfqn_clwoi(Z, N, mu, visits, ClwOptions());
}

/** Overload with unit visits. */
template <class T>
ClwResult<T> pfqn_clwoi(const std::vector<T>& Z, const std::vector<int>& N,
                         const std::vector<OiRate<T>>& mu) {
    return pfqn_clwoi(Z, N, mu, Matrix<T>(), ClwOptions());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_CLWOI_H
