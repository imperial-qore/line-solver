/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_EXPLICIT_LD_H
#define LINE_API_PFQN_PFQN_EXPLICIT_LD_H

/**
 * Explicit closed-form normalizing constant of a multiclass LIMITED LOAD-DEPENDENT
 * network.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_explicit_ld.m. Load-dependent
 * counterpart of pfqn_explicit: it evaluates the same divided-difference form of
 * G. Casale, "Accelerating Performance Inference over Closed Systems by
 * Asymptotic Methods", ACM SIGMETRICS 2017, Corollary 3.2,
 *
 *   G(N) = sum_{0<=t<=N} (-1)^(|N|-|t|)/(N_1!...N_R!) prod_r C(N_r,t_r) h_t(|N|)
 *
 * but substitutes for the single-class constant h_t(|N|) the limited
 * load-dependent closed form of G. Casale, P. G. Harrison, W. H. Ong,
 * "Facilitating Load-Dependent Queueing Analysis Through Factorization",
 * Perform. Eval. 2021, Theorem 1, Eq. (8),
 *
 *   h_theta(N) = sum_{0<=v<s} g_sigma(N-|v|) prod_k phi_k(v_k)
 *   phi_k(v_k) = theta_k^v_k / prod_{t=1..v_k} alpha_k(t) * (1 - alpha_k(v_k)/alpha_k(s_k))
 *
 * at the induced demands theta_k(t) = sum_r t_r L(k,r). Here alpha_k(.) = mu(k,.)
 * is the load-dependent scaling of station k, s_k the population past which it
 * stays constant, sigma_k = theta_k/alpha_k(s_k) the SCALED demands, and g_sigma
 * the FIXED-RATE single-class constant at those scaled demands, which is exactly
 * what pfqn_explicit evaluates in closed form (Eqs. 15 and 16). The result is
 * explicit throughout, with no recursion over population; pfqn_gldsingle is the
 * same constant by an O(M|N|^2) recursion instead.
 *
 * TWO CONVENTIONS OF THEOREM 1 ARE NOT THOSE OF THE EQUILIBRIUM DISTRIBUTION.
 * alpha_k(0) is taken as ZERO inside the bracket of phi_k, so that phi_k(0) = 1,
 * even though the state probabilities use alpha_k(0) = 1; and g_sigma(n) = 0 for
 * n < 0, which caps the outer sum at |v| <= |N|. With alpha_k(n) = min(n,s_k) the
 * expression collapses to Gordon's multi-server formula, Oper. Res. 38(5), 1990,
 * Eq. (29), but unlike that one it needs neither a multi-server shape nor
 * distinct scaled demands.
 *
 * LIMITED LOAD DEPENDENCE. Theorem 1 holds for any s_k with
 * alpha_k(n) = alpha_k(s_k) for all n >= s_k, and a LARGER s_k is always
 * admissible, so s_k is detected here as the smallest index whose value the tail
 * of mu(k,:) repeats to within tol. A station whose rates never settle (an
 * infinite server, mu(k,n) = n) gets s_k = |N|, which is still exact: populations
 * above |N| do not occur, so redefining alpha_k there changes nothing. It is
 * merely expensive, since the inner sum costs prod_k s_k terms, capped by
 * |v| <= |N|. Think time is not admissible: a delay would have to enter g_sigma,
 * whose closed form covers queues only.
 *
 * ARITHMETIC. Both sums alternate in sign with terms far larger than the result,
 * so they are evaluated as SIGNED log-sum-exps, and the routine is gated on
 * num_traits<T>::has_transcendental exactly as pfqn_explicit is; a caller that
 * wants the same constant in exact arithmetic wants pfqn_gld. phi_k is
 * sign-definite when alpha_k increases, as a multi-server station does, and
 * changes sign where alpha_k decreases, so a decreasing rate function costs
 * digits in the inner sum too.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "line/api/pfqn/pfqn_explicit.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

namespace detail {

/**
 * Theorem 1 of Casale-Harrison-Ong (2021), Eq. (8): the single-class limited
 * load-dependent constant at induced demands th and total population Nt, as the
 * finite sum over 0 <= v < s of the fixed-rate constant at the scaled demands
 * th/alpha(s), one population level lower for every job held back by v.
 */
template <class T>
SignedLse<T> explicit_hlld(const std::vector<T>& th, std::size_t M, long Nt,
                           const std::vector<T>& alphaS, const std::vector<int>& vcap,
                           const std::vector<std::vector<T> >& lcum,
                           const std::vector<std::vector<T> >& lbr,
                           const std::vector<std::vector<double> >& sbr, const std::string& expr,
                           double tol) {
    using std::log;
    const T zero = num_traits<T>::from_int(0);
    const double dinf = std::numeric_limits<double>::infinity();
    std::vector<T> sigma(M), lth(M);
    for (std::size_t i = 0; i < M; ++i) {
        sigma[i] = T(th[i] / alphaS[i]);
        lth[i] = (th[i] > zero) ? T(log(th[i])) : num_traits<T>::from_double(-dinf);
    }
    std::vector<T> lterm;
    std::vector<double> sterm;
    double lossDigits = 0.0;
    std::vector<int> v(M, 0);
    while (true) {
        long nv = 0;
        for (std::size_t i = 0; i < M; ++i) nv += v[i];
        if (nv <= Nt) {
            T lval = zero;
            double sval = 1.0;
            bool dead = false;
            for (std::size_t i = 0; i < M; ++i) {
                const int vi = v[i];
                if (vi > 0) {
                    if (!std::isfinite(num_traits<T>::to_double(lth[i]))) {
                        // theta_k = 0 kills every v_k>0, and 0^0=1 keeps v_k=0
                        dead = true;
                        break;
                    }
                    // kept inside the guard because 0*(-Inf) is NaN, not 0
                    lval = T(lval + num_traits<T>::from_int(vi) * lth[i]);
                }
                lval = T(lval - lcum[i][vi] + lbr[i][vi]);
                sval *= sbr[i][vi];
            }
            if (!dead && sval != 0.0 && std::isfinite(num_traits<T>::to_double(lval))) {
                T lg = zero;
                int sg = 1;
                double dl = 0.0;
                if (Nt - nv > 0) {
                    // g_sigma(0) = 1 by definition, so the Nt == nv case is left
                    // exact: reading it off the partial fraction instead would
                    // spend digits on an alternating sum of known value.
                    const SignedLse<T> g =
                        (expr == "distinct")
                            ? explicit_gdistinct(sigma, num_traits<T>::from_int(Nt - nv), M)
                            : explicit_grepeated(sigma, num_traits<T>::from_int(Nt - nv), M, tol);
                    lg = g.lS;
                    sg = g.sgn;
                    dl = g.lossDigits;
                }
                lossDigits = std::max(lossDigits, dl);
                if (sg != 0) {
                    lterm.push_back(T(lval + lg));
                    sterm.push_back(sval * static_cast<double>(sg));
                }
            }
        }
        std::size_t i = M;
        while (i > 0 && v[i - 1] == vcap[i - 1]) v[--i] = 0;
        if (i == 0) break;
        ++v[i - 1];
    }
    SignedLse<T> h = explicit_signed_logsumexp(lterm, sterm);
    h.lossDigits = std::max(lossDigits, h.lossDigits);
    return h;
}

}  // namespace detail

/**
 * @param L       (M x R) service demands
 * @param N       population per class
 * @param mu      (M x >= sum(N)) load-dependent rate lattice, alpha_i(j) = mu(i,j-1);
 *                an empty matrix means all ones
 * @param tol     relative tolerance declaring two scaled demands redundant, and the
 *                rate tail constant
 * @param method  "auto", "distinct" (force Eq. 15) or "repeated" (force Eq. 16)
 * @param maxloss cancellation budget in decimal digits; a finite value turns the
 *                overrun into a silent REFUSAL (valid = false) for callers that
 *                hold a fallback, the default keeps the result whatever it costs
 */
template <class T>
ExplicitResult<T> pfqn_explicit_ld(const Matrix<T>& L, const std::vector<int>& N,
                                   const Matrix<T>& mu,
                                   double tol = std::numeric_limits<double>::epsilon(),
                                   const std::string& method = "auto",
                                   double maxloss = std::numeric_limits<double>::infinity()) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_explicit_ld requires transcendental arithmetic: the alternating sums are "
                  "carried as signed log-sum-exps so that no intermediate can overflow. Use "
                  "pfqn_gld for the same constant in exact arithmetic");
    using std::exp;
    using std::log;
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const double dinf = std::numeric_limits<double>::infinity();
    const std::size_t R = N.size();

    ExplicitResult<T> res;
    res.method = "distinct";
    res.lG = num_traits<T>::from_double(-dinf);
    res.G = zero;

    if (method != "auto" && method != "distinct" && method != "repeated")
        throw InputError(
            "pfqn_explicit_ld: unrecognized method, use 'auto', 'distinct' (Eq. 15) or 'repeated' "
            "(Eq. 16)");
    long Nsum = 0;
    for (int v : N) Nsum += v;
    if (Nsum < 0) return res;
    if (Nsum == 0) {
        res.lG = zero;
        res.G = one;
        return res;
    }
    if (L.rows() == 0 || L.cols() == 0) return res;
    if (static_cast<std::size_t>(L.cols()) != R)
        throw InputError("pfqn_explicit_ld: the demand matrix must have one column per class of N");
    const std::size_t M = static_cast<std::size_t>(L.rows());
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r)
            if (L(i, r) < zero)
                throw InputError("pfqn_explicit_ld: the demand matrix must be nonnegative");
    const std::size_t Nt = static_cast<std::size_t>(Nsum);

    // ---- the rate lattice, defaulting to a fixed-rate model ----
    std::vector<std::vector<T> > alpha(M, std::vector<T>(Nt, one));
    if (mu.rows() != 0 && mu.cols() != 0) {
        if (static_cast<std::size_t>(mu.rows()) != M)
            throw InputError(
                "pfqn_explicit_ld: the load-dependent rate matrix must have one row per station "
                "of L");
        if (static_cast<std::size_t>(mu.cols()) < Nt)
            throw InputError(
                "pfqn_explicit_ld: the load-dependent rate matrix must have at least sum(N) "
                "columns");
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t k = 0; k < Nt; ++k) {
                alpha[i][k] = mu(i, k);
                if (!(alpha[i][k] > zero))
                    throw InputError(
                        "pfqn_explicit_ld: the load-dependent rates must be strictly positive");
            }
    }

    // ---- s_k: the smallest index whose value the tail of the rate row repeats ----
    // Any larger s_k also satisfies alpha_k(n)=alpha_k(s_k) for n>=s_k, so a missed
    // tie only adds terms; a false tie would be a wrong answer, hence the strict tol.
    std::vector<std::size_t> s(M, Nt);
    std::vector<T> alphaS(M, one);
    for (std::size_t i = 0; i < M; ++i) {
        const T tail = alpha[i][Nt - 1];
        const double atail = std::fabs(num_traits<T>::to_double(tail));
        for (std::size_t n = Nt; n > 1; --n) {
            const T d = T(alpha[i][n - 2] - tail);
            if (std::fabs(num_traits<T>::to_double(d)) <= tol * std::max(atail, 1.0))
                s[i] = n - 1;
            else
                break;
        }
        alphaS[i] = alpha[i][s[i] - 1];
    }

    // ---- per-station phi tables, in the log domain, indexed by v_k = 0..s_k-1 ----
    std::vector<std::vector<T> > lcum(M), lbr(M);
    std::vector<std::vector<double> > sbr(M);
    std::vector<int> vcap(M, 0);
    for (std::size_t i = 0; i < M; ++i) {
        lcum[i].resize(s[i]);
        lbr[i].resize(s[i]);
        sbr[i].resize(s[i]);
        for (std::size_t v = 0; v < s[i]; ++v) {
            lcum[i][v] = (v == 0) ? zero : T(lcum[i][v - 1] + log(alpha[i][v - 1]));
            const T br = T(one - (((v == 0) ? zero : alpha[i][v - 1]) / alphaS[i]));
            if (br == zero) {
                lbr[i][v] = num_traits<T>::from_double(-dinf);
                sbr[i][v] = 0.0;
            } else if (br > zero) {
                lbr[i][v] = log(br);
                sbr[i][v] = 1.0;
            } else {
                lbr[i][v] = log(T(zero - br));
                sbr[i][v] = -1.0;
            }
        }
        // g_sigma vanishes below zero population, Eq. (8) caps |v| <= |N|
        vcap[i] = static_cast<int>(std::min<std::size_t>(s[i] - 1, Nt));
    }

    // ---- redundancy scan: are the SCALED induced demands pairwise distinct? ----
    // The scan MUST form sigma exactly as explicit_hlld does, (L*t)/alphaS and
    // not (L/alphaS)*t: the two orderings differ in the last ulp, so an exact
    // tie can clear an eps-relative gap under one and not the other, and
    // Eq. (15) would then divide by that ulp.
    const auto induced = [&](const std::vector<int>& t) {
        std::vector<T> th(M, zero);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < R; ++r)
                th[i] = T(th[i] + L(i, r) * num_traits<T>::from_int(t[r]));
        return th;
    };
    const auto scaled = [&](const std::vector<int>& t) {
        std::vector<T> th = induced(t);
        for (std::size_t i = 0; i < M; ++i) th[i] = T(th[i] / alphaS[i]);
        return th;
    };
    const auto redundant_at = [&](std::vector<T> th) {
        std::sort(th.begin(), th.end(), [](const T& a, const T& b) { return a < b; });
        const T scale = th.back();
        // scale == 0 leaves every scaled demand at zero, so the term takes no part
        // in the sum
        if (!(scale > zero)) return false;
        const T gap = T(num_traits<T>::from_double(tol) * scale);
        for (std::size_t i = 1; i < th.size(); ++i)
            if (!(T(th[i] - th[i - 1]) > gap)) return true;
        return false;
    };
    bool isRedundant = false;
    if (R == 1) {
        // the scaled demands at t are t*sigma, so both the tie structure and the
        // relative tolerance are those of sigma itself, at every t at once
        std::vector<T> th(M);
        for (std::size_t i = 0; i < M; ++i) th[i] = T(L(i, 0) / alphaS[i]);
        isRedundant = redundant_at(th);
    } else {
        std::vector<int> t(R, 0);
        while (true) {
            long ts = 0;
            for (int val : t) ts += val;
            if (ts > 0 && redundant_at(scaled(t))) {
                isRedundant = true;
                break;
            }
            std::size_t r = R;
            while (r > 0 && t[r - 1] == N[r - 1]) t[--r] = 0;
            if (r == 0) break;
            ++t[r - 1];
        }
    }
    std::string expr = method;
    if (expr == "auto") {
        expr = isRedundant ? "repeated" : "distinct";
    } else if (expr == "distinct" && isRedundant) {
        throw InputError(
            "pfqn_explicit_ld: Eq. (15) requires pairwise distinct scaled demands, but two of them "
            "agree to within tol. Use 'auto' or 'repeated'");
    }
    res.method = expr;

    detail::SignedLse<T> total;
    if (R == 1) {
        // ---- single class: the divided difference is the identity ----
        std::vector<T> th(M);
        for (std::size_t i = 0; i < M; ++i) th[i] = L(i, 0);
        total = detail::explicit_hlld(th, M, Nsum, alphaS, vcap, lcum, lbr, sbr, expr, tol);
    } else {
        // ---- outer divided-difference sum over 0 <= t <= N ----
        std::vector<T> lterm;
        std::vector<double> sterm;
        double innerLoss = 0.0;
        std::vector<int> t(R, 0);
        while (true) {
            long ts = 0;
            for (int val : t) ts += val;
            if (ts > 0) {
                std::vector<T> th = induced(t);
                T thmax = th[0];
                for (const T& val : th)
                    if (val > thmax) thmax = val;
                if (thmax > zero) {
                    const detail::SignedLse<T> h = detail::explicit_hlld(
                        th, M, Nsum, alphaS, vcap, lcum, lbr, sbr, expr, tol);
                    innerLoss = std::max(innerLoss, h.lossDigits);
                    if (h.sgn != 0) {
                        T l = h.lS;
                        for (std::size_t r = 0; r < R; ++r) {
                            l = T(l - detail::num_factln<T>(num_traits<T>::from_int(t[r])));
                            l = T(l -
                                  detail::num_factln<T>(num_traits<T>::from_int(N[r] - t[r])));
                        }
                        lterm.push_back(l);
                        sterm.push_back(static_cast<double>(h.sgn) *
                                        (((Nsum - ts) % 2 == 0) ? 1.0 : -1.0));
                    }
                }
            }
            std::size_t r = R;
            while (r > 0 && t[r - 1] == N[r - 1]) t[--r] = 0;
            if (r == 0) break;
            ++t[r - 1];
        }
        total = detail::explicit_signed_logsumexp(lterm, sterm);
        total.lossDigits = std::max(total.lossDigits, innerLoss);
    }
    res.lossDigits = total.lossDigits;

    // A caller that named a cancellation budget has a fallback and wants a verdict,
    // not a warning: refuse quietly. lossDigits is infinite when the sum vanished
    // identically, which is a total loss rather than a legitimate G = 0.
    if (std::isfinite(maxloss) && (total.sgn < 0 || total.lossDigits > maxloss)) {
        res.valid = false;
        res.lG = num_traits<T>::from_double(std::numeric_limits<double>::quiet_NaN());
        res.G = res.lG;
        return res;
    }
    if (total.sgn == 0) {
        res.lG = num_traits<T>::from_double(-dinf);
        res.G = zero;
        return res;
    }
    if (total.sgn < 0) {
        // Double precision is exhausted by cancellation; the sign itself is wrong,
        // so there is no result to hand back.
        res.valid = false;
        res.lG = num_traits<T>::from_double(std::numeric_limits<double>::quiet_NaN());
        res.G = res.lG;
        return res;
    }
    res.lG = total.lS;
    res.G = exp(res.lG);
    return res;
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PFQN_EXPLICIT_LD_H
