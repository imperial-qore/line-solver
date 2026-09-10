/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_CACHE_CACHE_MISS_POS_RMF_H
#define LINE_API_CACHE_CACHE_MISS_POS_RMF_H

/**
 * Position-resolved mean-field miss rates for FIFO(m) and strict FIFO(m).
 *
 * Templated port of `matlab/src/api/cache/cache_miss_fifo_rmf.m`,
 * `cache_miss_sfifo_rmf.m` and the drift they share,
 * `cache_pos_drift_graph.m`.
 *
 * WHY THESE EXIST AT ALL, GIVEN `cache_miss_rmf`. Gast and Van Houdt
 * (SIGMETRICS 2015, Thm 1) prove pi_FIFO(m) = pi_RAND(m) EXACTLY, so on the
 * linear access graph the FIFO STEADY STATE is served from the cheaper refined
 * mean field of `cache_miss_rmf` and this file is not reached. Two things break
 * that equality:
 *
 *   THE TRANSIENT. FIFO evicts the deterministic tail -- residence is exactly
 *   m insertions -- while RANDOM evicts a uniformly drawn victim, so residence
 *   is geometric. H(inf) agrees; H(t) from a cold cache does not, and a
 *   trajectory read off the RANDOM drift would ramp at the wrong rate.
 *
 *   A NON-LINEAR ACCESS GRAPH. The equality is proved for the linear chain. Once
 *   admission or promotion is item-dependent, the per-item per-list occupancy
 *   that RANDOM(m) tracks no longer determines the dynamics, and FIFO needs its
 *   own position-resolved state.
 *
 * Strict FIFO(m) is a THIRD policy, not a spelling of FIFO(m). Gast and Van
 * Houdt show it differs from RANDOM(m) and give it no mean-field model. The
 * difference is the within-list age ordering: on a hit at position j of list
 * i < h the demoted tail of list i+1 is reinserted at position 1 of list i and
 * positions 1..j-1 shift back (strict), whereas FIFO(m) drops it into the
 * VACATED position j with no shift. That single choice is the `reinsert`
 * parameter of the shared drift, and it is why strict FIFO(m) is never served
 * from `cache_miss_rmf` even on the linear chain -- it degenerates to FIFO(m)
 * only when m_1 = ... = m_{h-1} = 1, where there is no within-list order to
 * disagree about.
 *
 * THE STATE. x[k,i,j] = P(item k occupies position j of list i), over the
 * `sum(m)` in-cache slots only; the out-of-cache mass is the complement
 * 1 - sum_{i,j} x[k,i,j], which is what `pos_out` returns and what the miss
 * probability is read from. That complement is CLIPPED to [0,1] rather than
 * asserted, exactly as the reference does: the mean-field trajectory can leave
 * the simplex by an integration tolerance without the fixed point being wrong.
 *
 * THE INITIAL CONDITION DIFFERS BY PATH, and deliberately. The linear drift
 * starts POPULARITY-ORDERED (the S most requested items pre-loaded, one per
 * slot), which is near its own fixed point and integrates quickly. The general
 * graph drift starts COLD (an empty cache), because a graph may make some items
 * non-admissible and a pre-loaded non-admissible item has no outflow term that
 * can drain it -- it would sit in the cache forever and report a hit rate the
 * policy never produces.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <string>
#include <vector>

#include "line/api/cache/cache_miss_rmf.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/util/ode.h"

namespace line {
namespace cache {

/** Return value of the position-resolved routines, as CacheMissRmfResult. */
template <class T>
using CacheMissPosRmfResult = CacheMissRmfResult<T>;

namespace pos_detail {

/** The (list, position) slot map: `slots[s] = (i,j)` with i, j 1-based. */
struct SlotMap {
    std::vector<std::pair<std::size_t, std::size_t> > slots;  ///< (list, position), 1-based
    std::vector<std::vector<std::size_t> > sidx;              ///< sidx[i-1][j-1] = s, 0-based
    std::size_t S = 0;
};

inline SlotMap build_slots(const std::vector<int>& m) {
    SlotMap sm;
    const std::size_t h = m.size();
    std::size_t mmax = 0;
    for (std::size_t i = 0; i < h; ++i) mmax = std::max(mmax, static_cast<std::size_t>(m[i]));
    sm.sidx.assign(h, std::vector<std::size_t>(mmax, 0));
    for (std::size_t i = 1; i <= h; ++i)
        for (int j = 1; j <= m[i - 1]; ++j) {
            sm.slots.push_back(std::make_pair(i, static_cast<std::size_t>(j)));
            sm.sidx[i - 1][static_cast<std::size_t>(j) - 1] = sm.slots.size() - 1;
        }
    sm.S = sm.slots.size();
    return sm;
}

/** `(k-1)*S + sidx(i,j)` in 0-based form. */
inline std::size_t kidx(std::size_t k, std::size_t i, std::size_t j, const SlotMap& sm) {
    return k * sm.S + sm.sidx[i - 1][j - 1];
}

/** `fifo_out` / `sfifo_out`: out-of-cache occupancy of item k, clipped. */
template <class T>
T pos_out(const std::vector<T>& x, std::size_t k, const SlotMap& sm) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    T acc = zero;
    for (std::size_t s = 0; s < sm.S; ++s) acc += x[k * sm.S + s];
    T o = one - acc;
    if (o < zero) o = zero;
    if (o > one) o = one;
    return o;
}

/** Total in-cache occupancy of item k over the whole of list i. */
template <class T>
T pos_list_occ(const std::vector<T>& x, std::size_t k, std::size_t i, const std::vector<int>& m,
               const SlotMap& sm) {
    T acc = num_traits<T>::from_int(0);
    for (int j = 1; j <= m[i - 1]; ++j) acc += x[kidx(k, i, static_cast<std::size_t>(j), sm)];
    return acc;
}

/** `sfifo_gi` / `pos_gg`: the aggregate rate at positions DEEPER than jp. */
template <class T>
T pos_deeper(const Matrix<T>& H, std::size_t i, std::size_t jp, const std::vector<int>& m,
             std::size_t h) {
    const T zero = num_traits<T>::from_int(0);
    if (i == h) return zero;  // the top list never shifts on a hit
    T g = zero;
    for (int jj = static_cast<int>(jp) + 1; jj <= m[i - 1]; ++jj)
        g += H(i - 1, static_cast<std::size_t>(jj) - 1);
    return g;
}

/**
 * `fifo_drift` and `sfifo_drift`, which differ only in the two `strict` arms.
 *
 * FIFO(m) shifts a whole list on every insertion into it, so its outflow term
 * is `Sfull(i) * x`. Strict FIFO(m) ALSO shifts the prefix of a list when a hit
 * lands deeper in that same list, which is the `pos_deeper` term; and the tail
 * demoted from list i+1 arrives at position 1 rather than at the position the
 * promoted job vacated, which is the `j == 1` arm. Everything else is common,
 * and writing it once is what keeps the two policies' shared terms from
 * drifting apart.
 */
template <class T>
std::vector<T> pos_drift_linear(const std::vector<T>& x_in, const std::vector<T>& p,
                                const std::vector<int>& m, std::size_t n, std::size_t h,
                                const SlotMap& sm, bool strict) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    std::vector<T> x = x_in;
    for (std::size_t a = 0; a < x.size(); ++a) {
        if (x[a] < zero) x[a] = zero;
        if (x[a] > one) x[a] = one;
    }
    std::size_t mmax = 0;
    for (std::size_t i = 0; i < h; ++i) mmax = std::max(mmax, static_cast<std::size_t>(m[i]));

    Matrix<T> Hpos(h, mmax, zero);
    std::vector<T> Hi(h + 1, zero);
    for (std::size_t s = 0; s < sm.S; ++s) {
        const std::size_t i = sm.slots[s].first, j = sm.slots[s].second;
        T acc = zero;
        for (std::size_t k = 0; k < n; ++k) acc += p[k] * x[kidx(k, i, j, sm)];
        Hpos(i - 1, j - 1) = acc;
        Hi[i - 1] += acc;
    }
    T Mrate = zero;
    for (std::size_t k = 0; k < n; ++k) Mrate += p[k] * pos_out(x, k, sm);

    // `Sfull(i)`: the rate at which list i shifts as a whole. List 1 shifts on
    // every miss; list i > 1 shifts on every hit in list i-1, i.e. on every
    // promotion INTO it.
    std::vector<T> Sfull(h + 1, zero);
    Sfull[0] = Mrate;
    for (std::size_t i = 2; i <= h; ++i) Sfull[i - 1] = Hi[i - 2];

    std::vector<T> dX(n * sm.S, zero);
    for (std::size_t k = 0; k < n; ++k)
        for (std::size_t s = 0; s < sm.S; ++s) {
            const std::size_t i = sm.slots[s].first, j = sm.slots[s].second;
            const std::size_t at = kidx(k, i, j, sm);
            const T xk = x[at];
            T o = Sfull[i - 1] * xk;
            if (strict) o += pos_deeper(Hpos, i, j, m, h) * xk;
            if (i < h) o += p[k] * xk;
            dX[at] -= o;
            if (j >= 2) {
                T rate = Sfull[i - 1];
                if (strict) rate += pos_deeper(Hpos, i, j - 1, m, h);
                dX[at] += rate * x[kidx(k, i, j - 1, sm)];
            } else {
                if (i == 1)
                    dX[kidx(k, 1, 1, sm)] += p[k] * pos_out(x, k, sm);
                else
                    dX[kidx(k, i, 1, sm)] += p[k] * pos_list_occ(x, k, i - 1, m, sm);
                // Strict reinsertion: the demoted tail of list i+1 lands at
                // position 1 and the prefix shifts back.
                if (strict && i < h)
                    dX[kidx(k, i, 1, sm)] +=
                        Hi[i - 1] * x[kidx(k, i + 1, static_cast<std::size_t>(m[i]), sm)];
            }
            // FIFO reinsertion: the demoted tail lands IN PLACE, at the position
            // the promoted item just vacated, so the flow is per-position.
            if (!strict && i < h)
                dX[at] += Hpos(i - 1, j - 1) *
                          x[kidx(k, i + 1, static_cast<std::size_t>(m[i]), sm)];
        }
    return dX;
}

/**
 * Port of `cache_pos_drift_graph.m`: the same two policies under a per-item
 * access graph.
 *
 * `strict` is the reference's `reinsert` argument, 'head' (strict FIFO) against
 * 'pos' (FIFO). Row 0 of G[k] is miss admission (column 0 rejects, column 1+l
 * admits to list l); row 1+i is a hit in list i, where column 1+i means STAY IN
 * PLACE -- that diagonal entry is why the outflow carries `1 - G[k](i,i)` and
 * not `1`, and is the FIFO/SFIFO convention rather than the RANDOM(m) one.
 */
template <class T>
std::vector<T> pos_drift_graph(const std::vector<T>& x_in, const std::vector<T>& p,
                               const std::vector<Matrix<T> >& G, const std::vector<int>& m,
                               std::size_t n, std::size_t h, const SlotMap& sm, bool strict) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    std::vector<T> x = x_in;
    for (std::size_t a = 0; a < x.size(); ++a) {
        if (x[a] < zero) x[a] = zero;
        if (x[a] > one) x[a] = one;
    }
    std::size_t mmax = 0;
    for (std::size_t i = 0; i < h; ++i) mmax = std::max(mmax, static_cast<std::size_t>(m[i]));

    std::vector<T> MI(h + 1, zero);      // MI[l]: miss admission into list l
    Matrix<T> HP(h + 1, h + 1, zero);    // HP(i,b): promotion i -> b, b > i
    for (std::size_t k = 0; k < n; ++k) {
        const T ok = pos_out(x, k, sm);
        for (std::size_t l = 1; l <= h; ++l) MI[l] += p[k] * ok * G[k](0, l);
        for (std::size_t i = 1; i <= h; ++i) {
            const T oc = pos_list_occ(x, k, i, m, sm);
            for (std::size_t b = i + 1; b <= h; ++b) HP(i, b) += p[k] * oc * G[k](i, b);
        }
    }
    std::vector<T> Sin(h + 1, zero);
    for (std::size_t l = 1; l <= h; ++l) {
        Sin[l] = MI[l];
        for (std::size_t s = 1; s + 1 <= l; ++s) Sin[l] += HP(s, l);
    }
    // POp(i,j): the rate at which the occupant of (i,j) LEAVES its position on a
    // hit, i.e. weighted by 1 - G(i,i) rather than by 1.
    Matrix<T> POp(h, mmax, zero);
    for (std::size_t s = 0; s < sm.S; ++s) {
        const std::size_t i = sm.slots[s].first, j = sm.slots[s].second;
        T acc = zero;
        for (std::size_t k = 0; k < n; ++k)
            acc += p[k] * x[kidx(k, i, j, sm)] * T(one - G[k](i, i));
        POp(i - 1, j - 1) = acc;
    }

    std::vector<T> dX(n * sm.S, zero);
    for (std::size_t k = 0; k < n; ++k) {
        const T ok = pos_out(x, k, sm);
        for (std::size_t s = 0; s < sm.S; ++s) {
            const std::size_t i = sm.slots[s].first, j = sm.slots[s].second;
            const std::size_t at = kidx(k, i, j, sm);
            const T xk = x[at];
            T o = p[k] * xk * T(one - G[k](i, i));
            o += (strict ? T(Sin[i] + pos_deeper(POp, i, j, m, h)) : Sin[i]) * xk;
            dX[at] -= o;
            if (j >= 2) {
                const T rate = strict ? T(Sin[i] + pos_deeper(POp, i, j - 1, m, h)) : Sin[i];
                dX[at] += rate * x[kidx(k, i, j - 1, sm)];
            } else {
                dX[kidx(k, i, 1, sm)] += p[k] * ok * G[k](0, i);
                for (std::size_t ss = 1; ss + 1 <= i; ++ss)
                    dX[kidx(k, i, 1, sm)] +=
                        p[k] * pos_list_occ(x, k, ss, m, sm) * G[k](ss, i);
            }
            for (std::size_t b = i + 1; b <= h; ++b) {
                if (strict) {
                    if (j == 1)
                        dX[kidx(k, i, 1, sm)] +=
                            HP(i, b) * x[kidx(k, b, static_cast<std::size_t>(m[b - 1]), sm)];
                } else {
                    T poj = zero;
                    for (std::size_t kk = 0; kk < n; ++kk)
                        poj += p[kk] * x[kidx(kk, i, j, sm)] * G[kk](i, b);
                    dX[at] += poj * x[kidx(k, b, static_cast<std::size_t>(m[b - 1]), sm)];
                }
            }
        }
    }
    return dX;
}

/** The shared body of the two entry points; `strict` selects the policy. */
template <class T>
CacheMissPosRmfResult<T> pos_rmf(const std::vector<T>& gamma, const std::vector<int>& m,
                                 const Matrix<T>& lambda,
                                 const std::vector<std::vector<Matrix<T> > >& accost,
                                 bool strict, bool want_transient, const T& t0, const T& t1,
                                 const std::vector<T>& x0init) {
    static_assert(num_traits<T>::has_transcendental,
                  "cache_miss_fifo_rmf / cache_miss_sfifo_rmf require transcendental arithmetic: "
                  "the fixed point is reached by a tolerance-driven integration");
    (void)gamma;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t u = lambda.rows(), n = lambda.cols(), h = m.size();
    const std::string who = strict ? "cache_miss_sfifo_rmf" : "cache_miss_fifo_rmf";
    if (u == 0 || n == 0) throw InputError(who + ": empty request-rate matrix");
    if (h == 0) throw InputError(who + ": at least one cache list is required");
    for (std::size_t i = 0; i < h; ++i)
        if (m[i] <= 0) throw InputError(who + ": a list has non-positive capacity");

    // non-finite rates are dropped, as the reference's `row(~isfinite(row)) = 0`
    std::vector<T> lam_i(n, zero);
    T tot = zero;
    for (std::size_t v = 0; v < u; ++v)
        for (std::size_t k = 0; k < n; ++k) {
            if (!std::isfinite(num_traits<T>::to_double(lambda(v, k)))) continue;
            lam_i[k] += lambda(v, k);
            tot += lambda(v, k);
        }
    std::vector<T> p(n, zero);
    if (tot > zero)
        for (std::size_t k = 0; k < n; ++k) p[k] = lam_i[k] / tot;
    else
        for (std::size_t k = 0; k < n; ++k) p[k] = one / num_traits<T>::from_int(static_cast<long>(n));

    const SlotMap sm = build_slots(m);
    const std::size_t dim = n * sm.S;

    // Popularity-ordered warm start: the S most requested items, one per slot.
    std::vector<std::size_t> order(n);
    for (std::size_t k = 0; k < n; ++k) order[k] = k;
    std::stable_sort(order.begin(), order.end(),
                     [&](std::size_t a, std::size_t b) { return p[a] > p[b]; });
    std::vector<T> x0(dim, zero);
    for (std::size_t s = 0; s < sm.S && s < n; ++s) x0[order[s] * sm.S + s] = one;

    const std::vector<Matrix<T> > G = rmf_detail::build_item_graphs(accost, lambda, n, h);
    const bool graph = !G.empty();
    // A cold start on the graph path: see the header comment.
    const std::vector<T> x0s = graph ? std::vector<T>(dim, zero) : x0;
    const auto f = [&](const T& t, const std::vector<T>& xx) {
        (void)t;
        return graph ? pos_drift_graph(xx, p, G, m, n, h, sm, strict)
                     : pos_drift_linear(xx, p, m, n, h, sm, strict);
    };

    OdeOptions<T> opt;
    opt.rtol = num_traits<T>::from_double(1e-8);
    opt.atol = num_traits<T>::from_double(1e-10);
    opt.store_trajectory = false;

    CacheMissPosRmfResult<T> res;
    const std::vector<T> xss =
        ode_rosenbrock4(f, zero, T(num_traits<T>::from_int(20000)), x0s, opt).final_state();
    res.xss = xss;
    res.pi0.assign(n, zero);
    for (std::size_t k = 0; k < n; ++k) res.pi0[k] = pos_out(xss, k, sm);
    res.MI.assign(n, zero);
    res.M = zero;
    for (std::size_t k = 0; k < n; ++k) {
        res.MI[k] = lam_i[k] * res.pi0[k];
        res.M += res.MI[k];
    }
    res.MU.assign(u, zero);
    for (std::size_t v = 0; v < u; ++v) {
        T s = zero;
        for (std::size_t k = 0; k < n; ++k) {
            if (!std::isfinite(num_traits<T>::to_double(lambda(v, k)))) continue;
            s += lambda(v, k) * res.pi0[k];
        }
        res.MU[v] = s;
    }
    if (!want_transient) return res;

    OdeOptions<T> topt = opt;
    topt.store_trajectory = true;
    const std::vector<T> xt0 = x0init.empty() ? x0s : x0init;
    if (xt0.size() != dim)
        throw InputError(who + ": the initial occupancy has the wrong dimension");
    const OdeSolution<T> tr = ode_rosenbrock4(f, t0, t1, xt0, topt);
    const std::size_t nt = tr.t.size();
    res.tout = tr.t;
    res.xtraj = Matrix<T>(dim, nt, zero);
    for (std::size_t c = 0; c < nt; ++c)
        for (std::size_t a = 0; a < dim; ++a) res.xtraj(a, c) = tr.y[c][a];
    res.pi0_t = Matrix<T>(n, nt, zero);
    for (std::size_t c = 0; c < nt; ++c)
        for (std::size_t k = 0; k < n; ++k) res.pi0_t(k, c) = pos_out(tr.y[c], k, sm);
    res.MU_t = Matrix<T>(u, nt, zero);
    for (std::size_t v = 0; v < u; ++v)
        for (std::size_t c = 0; c < nt; ++c) {
            T s = zero;
            for (std::size_t k = 0; k < n; ++k) {
                if (!std::isfinite(num_traits<T>::to_double(lambda(v, k)))) continue;
                s += lambda(v, k) * res.pi0_t(k, c);
            }
            res.MU_t(v, c) = s;
        }
    return res;
}

}  // namespace pos_detail

/**
 * Port of `cache_miss_fifo_rmf.m`: the FIFO(m) position-resolved mean field.
 *
 * @param gamma  present for signature compatibility; the reference marks it
 *               unused and so is it here
 * @param m      (h) list capacities
 * @param lambda (u x n) per-user per-item request rates, the reference's
 *               `lambda(:,:,1)` page
 * @param accost per-(user,item) access graph; empty is the linear chain
 */
template <class T>
CacheMissPosRmfResult<T> cache_miss_fifo_rmf(
    const std::vector<T>& gamma, const std::vector<int>& m, const Matrix<T>& lambda,
    const std::vector<std::vector<Matrix<T> > >& accost = std::vector<std::vector<Matrix<T> > >()) {
    return pos_detail::pos_rmf(gamma, m, lambda, accost, false, false,
                               num_traits<T>::from_int(0), num_traits<T>::from_int(0),
                               std::vector<T>());
}

/** `cache_miss_fifo_rmf` with the optional TSPAN/X0INIT transient. */
template <class T>
CacheMissPosRmfResult<T> cache_miss_fifo_rmf_transient(
    const std::vector<T>& gamma, const std::vector<int>& m, const Matrix<T>& lambda, const T& t0,
    const T& t1, const std::vector<T>& x0init,
    const std::vector<std::vector<Matrix<T> > >& accost = std::vector<std::vector<Matrix<T> > >()) {
    return pos_detail::pos_rmf(gamma, m, lambda, accost, false, true, t0, t1, x0init);
}

/**
 * Port of `cache_miss_sfifo_rmf.m`: the strict FIFO(m) position-resolved mean
 * field. Same arguments as `cache_miss_fifo_rmf`; the policies differ only in
 * where a demoted tail is reinserted (see the header comment).
 */
template <class T>
CacheMissPosRmfResult<T> cache_miss_sfifo_rmf(
    const std::vector<T>& gamma, const std::vector<int>& m, const Matrix<T>& lambda,
    const std::vector<std::vector<Matrix<T> > >& accost = std::vector<std::vector<Matrix<T> > >()) {
    return pos_detail::pos_rmf(gamma, m, lambda, accost, true, false,
                               num_traits<T>::from_int(0), num_traits<T>::from_int(0),
                               std::vector<T>());
}

/** `cache_miss_sfifo_rmf` with the optional TSPAN/X0INIT transient. */
template <class T>
CacheMissPosRmfResult<T> cache_miss_sfifo_rmf_transient(
    const std::vector<T>& gamma, const std::vector<int>& m, const Matrix<T>& lambda, const T& t0,
    const T& t1, const std::vector<T>& x0init,
    const std::vector<std::vector<Matrix<T> > >& accost = std::vector<std::vector<Matrix<T> > >()) {
    return pos_detail::pos_rmf(gamma, m, lambda, accost, true, true, t0, t1, x0init);
}

}  // namespace cache
}  // namespace line

#endif  // LINE_API_CACHE_CACHE_MISS_POS_RMF_H
