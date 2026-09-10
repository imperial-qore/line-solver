/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_CLWJD_H
#define LINE_API_PFQN_CLWJD_H

/**
 * Normalizing constant of a closed network of LIMITED JOINT-DEPENDENT (LJD)
 * stations plus one aggregated delay, by numerical inversion of the multichain
 * generating function (Choudhury-Leung-Whitt, J. ACM 42(5):935-970, 1995).
 *
 * Templated port of matlab/src/api/pfqn/pfqn_clwjd.m. This stands to
 * pfqn_clwoi as pfqn_clw_lld stands to pfqn_clw: a per-station cutoff beyond
 * which the rate stops changing turns an infinite series into a rational
 * function of the same denominators.
 *
 * Station i has a rate mu_i(n) that reads the whole per-class occupancy but
 * saturates coordinatewise: with a cutoff vector l_i,
 *
 *   mu_i(n) = c_{i,t},  t = ( min(n_1,l_{i,1}), ..., min(n_R,l_{i,R}) ),
 *
 * so past l_{i,r} further class-r jobs no longer change the rate. Order
 * independence is l_i = 1 (t is the support indicator); a multiserver station
 * with c servers is l_i = c, min(sum n, c) being a function of the clipped
 * vector once every l_{i,r} >= c.
 *
 * Splitting the count lattice by clipped region, on which mu_i is constant, and
 * writing F_{i,t} for the part of F_i carried by the states with t_i(n) = t,
 *
 *   ( mu_{i,t} - sum_{r: t_r = l_{i,r}} v_{i,r} z_r ) F_{i,t}(z)
 *        = sum_{r: t_r >= 1} v_{i,r} z_r F_{i,t-e_r}(z),   F_{i,0} = 1,
 *   F_i(z) = sum_t F_{i,t}(z),
 *
 * the two sides differing because removing a class-r job leaves the region only
 * on an UNsaturated coordinate: for t_r < l_{i,r} the region pins n_r = t_r so
 * n - e_r lands in t - e_r, while for t_r = l_{i,r} the region is n_r >= l and
 * n - e_r lands in t or in t - e_r. The singularities are therefore the
 * hyperplanes sum_{r in S} v_{i,r} z_r = mu_{i,t} over the SATURATED sets
 * S = {r : t_r = l_{i,r}}: at most 2^R per station, however large the cutoffs
 * are. Setting l_i = 1 reduces this to the support recursion of pfqn_clwoi and
 * R = 1 reduces it to Bertozzi-McKenna eq. 2.19.
 *
 * The restrictive static scaling of eqs. 5.41-5.46 runs on one row per
 * (station, saturated set), carrying v_{i,r}/min{mu_{i,t} : saturated set of t
 * is S}, the smallest rate over regions sharing a saturated set being the
 * binding one; rows dominated by a superset of no larger rate are dropped.
 *
 * SCOPE. The rate must be constant on each clipped region, which is checked on
 * probe states. Any rate is admissible with an all-N cutoff (the default), the
 * clipping being vacuous on the reachable lattice.
 *
 * Arithmetic: TRANSCENDENTAL, double and Real only, as for the rest of the CLW
 * family.
 *
 * COST. prod_r 2 l_r N_r contour points, each costing
 * O(M R prod_r (l_{i,r}+1)), against O(M prod_r (N_r+1)(N_r+2)/2) for
 * pfqn_ncjd: the inversion is linear rather than quadratic in each population,
 * but the per-point region box grows with the cutoff, so it pays off exactly
 * when the joint dependence saturates early and loses outright at cutoff N.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
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
 * The rate of one clipped region, with the constancy check. The region pins
 * every unsaturated coordinate and leaves the saturated ones free above the
 * cutoff, so the rate is probed at the region representative and above it.
 */
template <class T>
T clwjd_region_rate(const OiRate<T>& murate, const std::vector<int>& t,
                    const std::vector<std::size_t>& keep, std::size_t R,
                    const std::vector<int>& N, const std::vector<int>& lrow) {
    const T zero = num_traits<T>::from_int(0);
    std::vector<int> nrep(R, 0);
    for (std::size_t b = 0; b < t.size(); ++b) nrep[keep[b]] = t[b];
    const T rate = murate(nrep);
    if (!(rate > zero))
        throw InputError("pfqn_clwjd: a station has a non-positive rate on a reachable region");
    for (int pass = 0; pass < 2; ++pass) {
        std::vector<int> probe = nrep;
        bool differs = false;
        for (std::size_t b = 0; b < t.size(); ++b) {
            if (t[b] != lrow[b]) continue;
            const int nr = N[keep[b]];
            probe[keep[b]] = (pass == 0) ? nr : std::max(t[b], (t[b] + nr) / 2);
            if (probe[keep[b]] != nrep[keep[b]]) differs = true;
        }
        if (!differs) continue;
        const T rt = murate(probe);
        T diff = T(rt - rate);
        if (diff < zero) diff = T(-diff);
        const T scale = (rate > num_traits<T>::from_int(1)) ? rate : num_traits<T>::from_int(1);
        if (diff > num_traits<T>::from_double(1e-9) * scale)
            throw InputError(
                "pfqn_clwjd: a station has a rate that varies within a clipped region; raise the "
                "cutoff or use pfqn_ncjd");
    }
    return rate;
}

}  // namespace detail

/**
 * @param Z      (R) think-time demand of the aggregated delay node
 * @param N      (R) closed population, finite
 * @param mu     one rate handle per joint-dependent station; empty for a pure delay
 * @param visits (M x R) per-station class visit ratios; empty for unit visits
 * @param lcut   (M x R) per-station per-class saturation cutoffs l_{i,r} >= 1,
 *               clipped to N_r (exact: a rate difference at n_r > N_r can only
 *               move coefficients with n_r > N_r); empty for the all-N cutoff
 * @param opt    lattice and aliasing parameters
 */
template <class T>
ClwResult<T> pfqn_clwjd(const std::vector<T>& Z, const std::vector<int>& N,
                        const std::vector<OiRate<T>>& mu, const Matrix<T>& visits,
                        const Matrix<int>& lcut, const ClwOptions& opt) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_clwjd requires transcendental arithmetic (contour integration of a "
                  "generating function)");
    using std::exp;
    using std::log;
    const std::size_t R = N.size();
    const std::size_t M = mu.size();
    if (!Z.empty() && Z.size() != R)
        throw InputError("pfqn_clwjd: Z and N disagree on the chain count");
    if (!visits.empty() && (visits.rows() != M || visits.cols() != R))
        throw InputError("pfqn_clwjd: visits must be M x R");
    if (!lcut.empty() && (lcut.rows() != M || lcut.cols() != R))
        throw InputError("pfqn_clwjd: lcut must be M x R");
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
    // function at z_r = 0, which kills every F_{i,t} with t_r >= 1
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

    // saturation cutoffs, broadcast and clipped to the reachable lattice
    std::vector<std::vector<int>> Lk(M, std::vector<int>(p, 1));
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t b = 0; b < p; ++b) {
            int li = lcut.empty() ? Nk[b] : lcut(i, keep[b]);
            if (li < 1) li = 1;
            if (li > Nk[b]) li = Nk[b];
            Lk[i][b] = li;
        }

    // per-station region tables over the clipped box prod_r {0,...,l_{i,r}}, in
    // mixed radix so that t - e_r always precedes t
    std::vector<std::size_t> ntreg(M, 1);
    std::vector<std::vector<T>> muT(M);
    std::vector<std::vector<std::vector<std::size_t>>> regDec(M), regSat(M);
    std::vector<std::vector<std::size_t>> satMask(M);
    for (std::size_t i = 0; i < M; ++i) {
        std::vector<std::size_t> st(p, 1);
        std::size_t nt = 1;
        for (std::size_t b = 0; b < p; ++b) {
            st[b] = nt;
            nt *= static_cast<std::size_t>(Lk[i][b] + 1);
        }
        ntreg[i] = nt;
        muT[i].assign(nt, one);
        regDec[i].assign(nt, std::vector<std::size_t>());
        regSat[i].assign(nt, std::vector<std::size_t>());
        satMask[i].assign(nt, 0);
        std::vector<int> t(p, 0);
        for (std::size_t tl = 0; tl < nt; ++tl) {
            for (std::size_t b = 0; b < p; ++b)
                t[b] = static_cast<int>((tl / st[b]) % static_cast<std::size_t>(Lk[i][b] + 1));
            for (std::size_t b = 0; b < p; ++b) {
                if (t[b] >= 1) {
                    regDec[i][tl].push_back(b);
                    regDec[i][tl].push_back(tl - st[b]);
                }
                if (t[b] == Lk[i][b]) {
                    regSat[i][tl].push_back(b);
                    satMask[i][tl] |= (static_cast<std::size_t>(1) << b);
                }
            }
            if (tl > 0) muT[i][tl] = detail::clwjd_region_rate<T>(mu[i], t, keep, R, N, Lk[i]);
        }
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

    // binding rate of each saturated set: regions sharing a saturated set share the
    // hyperplane sum_{r in S} v_r z_r = mu, so the smallest rate constrains
    const std::size_t nmask = static_cast<std::size_t>(1) << p;
    std::vector<std::vector<T>> muS(M, std::vector<T>(nmask, zero));
    std::vector<std::vector<bool>> muSet(M, std::vector<bool>(nmask, false));
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t tl = 0; tl < ntreg[i]; ++tl) {
            const std::size_t sm = satMask[i][tl];
            if (sm == 0) continue;
            if (!muSet[i][sm] || muT[i][tl] < muS[i][sm]) {
                muS[i][sm] = muT[i][tl];
                muSet[i][sm] = true;
            }
        }

    // one constraint row per (station, saturated set), dominated rows dropped:
    // set S is implied by a superset S' with mu_{i,S'} <= mu_{i,S}
    std::vector<std::vector<T>> rows;
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t mask = 1; mask < nmask; ++mask) {
            if (!muSet[i][mask]) continue;
            bool dominated = false;
            for (std::size_t mask2 = 1; mask2 < nmask && !dominated; ++mask2)
                if (mask2 != mask && (mask & mask2) == mask && muSet[i][mask2] &&
                    muS[i][mask2] <= muS[i][mask] * num_traits<T>::from_double(1 + 1e-12))
                    dominated = true;
            if (dominated) continue;
            std::vector<T> row(p, zero);
            for (std::size_t b = 0; b < p; ++b)
                if (mask & (static_cast<std::size_t>(1) << b)) row[b] = V(i, b) / muS[i][mask];
            rows.push_back(row);
        }
    const std::size_t nrow = rows.size();
    Matrix<T> Lt(nrow ? nrow : 1, p, zero);
    for (std::size_t i = 0; i < nrow; ++i)
        for (std::size_t j = 0; j < p; ++j) Lt(i, j) = rows[i][j];

    const std::vector<long> mult(nrow ? nrow : 1, 1);
    const std::vector<T> alpha = detail::clw_scaling(Lt, Lt, Nk, Zk, l, r, mult);

    std::vector<T> arho0(p, zero);
    Matrix<T> vs(M ? M : 1, p, zero);
    for (std::size_t j = 0; j < p; ++j) {
        arho0[j] = alpha[j] * Zk[j];
        for (std::size_t i = 0; i < M; ++i) vs(i, j) = V(i, j) * alpha[j];
    }

    std::size_t ntmax = 1;
    for (std::size_t i = 0; i < M; ++i) ntmax = std::max(ntmax, ntreg[i]);
    std::vector<detail::Cx<T>> FT(ntmax, detail::Cx<T>(zero, zero));
    const auto gbar = [&](const std::vector<detail::Cx<T>>& w) {
        detail::Cx<T> expo(zero, zero);
        for (std::size_t j = 0; j < p; ++j)
            expo = detail::cx_add(
                expo, detail::cx_scale(detail::Cx<T>(T(w[j].re - one), w[j].im), arho0[j]));
        detail::Cx<T> logF(zero, zero);
        for (std::size_t i = 0; i < M; ++i) {
            FT[0] = detail::Cx<T>(one, zero);
            detail::Cx<T> tot(one, zero);
            for (std::size_t tl = 1; tl < ntreg[i]; ++tl) {
                detail::Cx<T> num(zero, zero);
                const std::vector<std::size_t>& dec = regDec[i][tl];
                for (std::size_t k = 0; k < dec.size(); k += 2)
                    num = detail::cx_add(
                        num, detail::cx_mul(detail::cx_scale(w[dec[k]], vs(i, dec[k])),
                                            FT[dec[k + 1]]));
                detail::Cx<T> den(muT[i][tl], zero);
                const std::vector<std::size_t>& sat = regSat[i][tl];
                for (std::size_t k = 0; k < sat.size(); ++k)
                    den = detail::cx_sub(den, detail::cx_scale(w[sat[k]], vs(i, sat[k])));
                FT[tl] = detail::cx_div(num, den);
                tot = detail::cx_add(tot, FT[tl]);
            }
            logF = detail::cx_add(logF, detail::cx_log(tot));
        }
        return detail::cx_exp(detail::cx_add(expo, logF));
    };

    std::vector<detail::Cx<T>> w(p);
    const detail::Cx<T> gv = detail::clw_invert(0, w, Nk, l, r, p, gbar);
    if (!(gv.re > zero))
        throw NumericError("pfqn_clwjd: the inverted generating function is not positive");

    T lG = log(gv.re);
    for (std::size_t j = 0; j < p; ++j)
        lG += arho0[j] - num_traits<T>::from_int(Nk[j]) * log(alpha[j]);
    res.lG = lG;
    res.G = (lG > num_traits<T>::from_int(709)) ? T(std::numeric_limits<T>::infinity()) : T(exp(lG));
    return res;
}

/** Overload with the CLW default parameters. */
template <class T>
ClwResult<T> pfqn_clwjd(const std::vector<T>& Z, const std::vector<int>& N,
                        const std::vector<OiRate<T>>& mu, const Matrix<T>& visits,
                        const Matrix<int>& lcut) {
    return pfqn_clwjd(Z, N, mu, visits, lcut, ClwOptions());
}

/** Overload with the all-N cutoff, i.e. no truncation of the joint dependence. */
template <class T>
ClwResult<T> pfqn_clwjd(const std::vector<T>& Z, const std::vector<int>& N,
                        const std::vector<OiRate<T>>& mu, const Matrix<T>& visits) {
    return pfqn_clwjd(Z, N, mu, visits, Matrix<int>(), ClwOptions());
}

/** Overload with unit visits and the all-N cutoff. */
template <class T>
ClwResult<T> pfqn_clwjd(const std::vector<T>& Z, const std::vector<int>& N,
                        const std::vector<OiRate<T>>& mu) {
    return pfqn_clwjd(Z, N, mu, Matrix<T>(), Matrix<int>(), ClwOptions());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_CLWJD_H
