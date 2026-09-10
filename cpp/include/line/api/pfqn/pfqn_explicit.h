/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_EXPLICIT_H
#define LINE_API_PFQN_PFQN_EXPLICIT_H

/**
 * Explicit closed-form normalizing constant of a multiclass closed network.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_explicit.m, i.e. Eqs. (15) and
 * (16) of G. Casale, "Accelerating Performance Inference over Closed Systems by
 * Asymptotic Methods", ACM SIGMETRICS 2017. Both instantiate the
 * divided-difference form of Corollary 3.2,
 *
 *   G(N) = sum_{0<=t<=N} (-1)^(|N|-|t|)/(N_1!...N_R!) prod_r C(N_r,t_r) g_t(|N|)
 *
 * by substituting a closed form for the single-class constant g_t(|N|) at the
 * induced demands theta_k(t) = sum_r t_r L(k,r). Eq. (15) is Gordon's partial
 * fraction and needs the induced demands PAIRWISE DISTINCT; Eq. (16) is the
 * general partial-fraction expansion over the distinct values and their
 * multiplicities, and reduces term by term to Eq. (15) when every multiplicity
 * is one. The choice is automatic.
 *
 * SINGLE CLASS. At R=1 the multiclass constant IS the single-class constant at
 * demands L, so the outer sum is skipped: g_t(N) = t^N g_1(N) and
 * sum_t (-1)^(N-t) t^N/(t!(N-t)!) = S(N,N) = 1. Running the difference anyway
 * would add N alternating terms, and their cancellation, to a closed form that
 * carries none of them. What is left is O(K^2) work at any population, which is
 * why the single-class route is the cheap one on large populations.
 *
 * ARITHMETIC. Both expressions alternate in sign with terms far larger than the
 * result, so they are evaluated as SIGNED log-sum-exps: that removes the
 * floating-point RANGE problem but not the cancellation, which is what makes
 * multiprecision arithmetic necessary on all but small models. The routine is
 * therefore gated on num_traits<T>::has_transcendental; a caller that wants the
 * same constant in exact arithmetic wants pfqn_ca.
 *
 * ADMISSIBILITY. Single-server load-independent queues only: infinite servers
 * need the integral form of Corollary 3.4 and load-dependent rates need
 * pfqn_explicit_ld, which keeps this closed form as its inner kernel.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "line/api/pfqn/pfqn_asympt_common.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Return value of pfqn_explicit, mirroring [lG, G, method, lossDigits]. */
template <class T>
struct ExplicitResult {
    T lG;                    ///< logarithm of the normalizing constant
    T G;                     ///< the normalizing constant
    std::string method;      ///< expression used, "distinct" (Eq. 15) or "repeated" (Eq. 16)
    double lossDigits = 0.0; ///< decimal digits lost to cancellation
    /// False when a caller's cancellation budget was exceeded: lG and G are then
    /// meaningless and the caller is expected to fall back.
    bool valid = true;
};

namespace detail {

/** Logarithm of the binomial coefficient C(n,m); MATLAB nchoosekln.m. */
template <class T>
T explicit_nchoosekln(const T& n, const T& m) {
    const T one = num_traits<T>::from_int(1);
    return T(num_lgamma<T>(T(one + n)) - num_lgamma<T>(T(one + n - m)) -
             num_lgamma<T>(T(one + m)));
}

/** Signed log-sum-exp of S = sum_i s_i exp(l_i): log|S|, sign(S), digits lost. */
template <class T>
struct SignedLse {
    T lS;
    int sgn;
    double lossDigits;
};

template <class T>
SignedLse<T> explicit_signed_logsumexp(const std::vector<T>& lterm,
                                       const std::vector<double>& sterm) {
    using std::exp;
    using std::log;
    const double dinf = std::numeric_limits<double>::infinity();
    SignedLse<T> out;
    out.lS = num_traits<T>::from_double(-dinf);
    out.sgn = 0;
    out.lossDigits = 0.0;
    std::vector<std::size_t> keep;
    for (std::size_t i = 0; i < lterm.size(); ++i) {
        const double li = num_traits<T>::to_double(lterm[i]);
        if (std::isfinite(li) && sterm[i] != 0.0) keep.push_back(i);
    }
    if (keep.empty()) return out;
    T a = lterm[keep[0]];
    for (std::size_t i : keep)
        if (lterm[i] > a) a = lterm[i];
    T s = num_traits<T>::from_int(0);
    for (std::size_t i : keep)
        s = T(s + num_traits<T>::from_double(sterm[i]) * exp(T(lterm[i] - a)));
    const T zero = num_traits<T>::from_int(0);
    if (s == zero) {
        out.lossDigits = dinf;
        return out;
    }
    out.sgn = (s > zero) ? 1 : -1;
    const T abss = (s > zero) ? s : T(zero - s);
    out.lS = T(a + log(abss));
    // max(exp(l_i - a)) is 1, so -log10|s| is the shortfall of the sum against
    // its largest term. Every one of the n terms carries a rounding error of
    // order eps*max_term, so the digits actually lost are that shortfall PLUS
    // log10(n); dropping the count understates the loss and lets a wrong answer
    // past the guard.
    out.lossDigits =
        std::max(0.0, std::log10(static_cast<double>(keep.size()) /
                                 std::fabs(num_traits<T>::to_double(abss))));
    return out;
}

/**
 * Eq. (14): single-class constant at pairwise distinct demands th, population
 * Nt over K queues. A zero demand contributes nothing, which also realizes the
 * 0/0 = 0 convention of Eq. (15) when the zero is repeated.
 */
template <class T>
SignedLse<T> explicit_gdistinct(const std::vector<T>& th, const T& Nt, std::size_t K) {
    using std::log;
    const T zero = num_traits<T>::from_int(0);
    const double dinf = std::numeric_limits<double>::infinity();
    std::vector<T> lin(K, num_traits<T>::from_double(-dinf));
    std::vector<double> sgv(K, 0.0);
    for (std::size_t k = 0; k < K; ++k) {
        if (!(th[k] > zero)) continue;
        T acc = T(T(Nt + num_traits<T>::from_int(static_cast<long>(K) - 1)) * log(th[k]));
        double sign = 1.0;
        for (std::size_t i = 0; i < K; ++i) {
            if (i == k) continue;
            const T d = T(th[k] - th[i]);
            const T ad = (d > zero) ? d : T(zero - d);
            acc = T(acc - log(ad));
            sign *= (d > zero) ? 1.0 : ((d < zero) ? -1.0 : 0.0);
        }
        lin[k] = acc;
        sgv[k] = sign;
    }
    return explicit_signed_logsumexp(lin, sgv);
}

/** Every Kp-vector r >= 0 with sum(r) = k, i.e. MATLAB multichoose(Kp,k). */
inline void explicit_multichoose(std::size_t Kp, int k, std::vector<int>& current,
                                 std::size_t idx, std::vector<std::vector<int> >& out) {
    if (idx + 1 == Kp) {
        current[idx] = k;
        out.push_back(current);
        return;
    }
    for (int i = 0; i <= k; ++i) {
        current[idx] = i;
        explicit_multichoose(Kp, k - i, current, idx + 1, out);
    }
}

inline std::vector<std::vector<int> > explicit_multichoose(std::size_t Kp, int k) {
    std::vector<std::vector<int> > out;
    if (Kp == 0 || k < 0) return out;
    std::vector<int> current(Kp, 0);
    explicit_multichoose(Kp, k, current, 0, out);
    return out;
}

/**
 * Eq. (16): single-class constant at demands th of arbitrary multiplicity,
 * population Nt over K queues. Demands within tol of each other, relatively to
 * the largest one, are merged into one distinct value carrying their count.
 */
template <class T>
SignedLse<T> explicit_grepeated(const std::vector<T>& th, const T& Nt, std::size_t K,
                                double tol) {
    using std::log;
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const double dinf = std::numeric_limits<double>::infinity();
    std::vector<T> ths(th);
    std::sort(ths.begin(), ths.end(), [](const T& a, const T& b) { return a < b; });
    T scale = ths.back();
    if (!(scale > zero)) scale = one;
    const T gap = T(num_traits<T>::from_double(tol) * scale);
    // clusters of consecutive near-ties, represented by their centroid
    std::vector<T> thd;
    std::vector<int> m;
    for (std::size_t i = 0; i < ths.size();) {
        std::size_t j = i + 1;
        while (j < ths.size() && !(T(ths[j] - ths[j - 1]) > gap)) ++j;
        T sum = zero;
        for (std::size_t q = i; q < j; ++q) sum = T(sum + ths[q]);
        thd.push_back(T(sum / num_traits<T>::from_int(static_cast<long>(j - i))));
        m.push_back(static_cast<int>(j - i));
        i = j;
    }
    const std::size_t Kp = thd.size();
    std::vector<T> lin;
    std::vector<double> sgv;
    for (std::size_t j = 0; j < Kp; ++j) {
        if (!(thd[j] > zero)) {
            // the exponent Nt+K-m_j is at least Nt>=1, so a zero cluster
            // contributes nothing
            continue;
        }
        const T louter =
            T(T(Nt + num_traits<T>::from_int(static_cast<long>(K) - m[j])) * log(thd[j]));
        const double souter = ((m[j] - 1) % 2 == 0) ? 1.0 : -1.0;
        const std::vector<std::vector<int> > rs = explicit_multichoose(Kp, m[j] - 1);
        for (std::size_t i = 0; i < rs.size(); ++i) {
            const std::vector<int>& r = rs[i];
            T lval = T(louter + explicit_nchoosekln<T>(T(Nt + num_traits<T>::from_int(r[j])),
                                                       num_traits<T>::from_int(r[j])));
            double sval = souter * ((r[j] % 2 == 0) ? 1.0 : -1.0);
            bool vanished = false;
            for (std::size_t k = 0; k < Kp; ++k) {
                if (k == j) continue;
                lval = T(lval + explicit_nchoosekln<T>(
                                    num_traits<T>::from_int(m[k] + r[k] - 1),
                                    num_traits<T>::from_int(r[k])));
                if (r[k] > 0) {
                    if (!(thd[k] > zero)) {
                        // theta_k^r_k vanishes; 0^0 = 1 is the r_k = 0 case
                        vanished = true;
                        break;
                    }
                    lval = T(lval + num_traits<T>::from_int(r[k]) * log(thd[k]));
                }
                const T dd = T(thd[j] - thd[k]);
                const T add = (dd > zero) ? dd : T(zero - dd);
                lval = T(lval - num_traits<T>::from_int(m[k] + r[k]) * log(add));
                if (!(dd > zero) && ((m[k] + r[k]) % 2 != 0)) sval = -sval;
            }
            if (vanished) {
                lin.push_back(num_traits<T>::from_double(-dinf));
                sgv.push_back(0.0);
            } else {
                lin.push_back(lval);
                sgv.push_back(sval);
            }
        }
    }
    return explicit_signed_logsumexp(lin, sgv);
}

}  // namespace detail

/**
 * @param L       (K x R) service demands of single-server load-independent queues
 * @param N       population per class
 * @param tol     relative tolerance declaring two induced demands redundant
 * @param method  "auto", "distinct" (force Eq. 15) or "repeated" (force Eq. 16)
 * @param maxloss cancellation budget in decimal digits; a finite value turns the
 *                overrun into a silent REFUSAL (valid = false) for callers that
 *                hold a fallback, the default keeps the result whatever it costs
 */
template <class T>
ExplicitResult<T> pfqn_explicit(const Matrix<T>& L, const std::vector<int>& N,
                                double tol = std::numeric_limits<double>::epsilon(),
                                const std::string& method = "auto",
                                double maxloss = std::numeric_limits<double>::infinity()) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_explicit requires transcendental arithmetic: the alternating sums are "
                  "carried as signed log-sum-exps so that no intermediate can overflow. Use "
                  "pfqn_ca for the same constant in exact arithmetic");
    using std::exp;
    using std::log;
    const T zero = num_traits<T>::from_int(0);
    const double dinf = std::numeric_limits<double>::infinity();
    const std::size_t R = N.size();

    ExplicitResult<T> res;
    res.method = "distinct";
    res.lG = num_traits<T>::from_double(-dinf);
    res.G = zero;

    if (method != "auto" && method != "distinct" && method != "repeated")
        throw InputError(
            "pfqn_explicit: unrecognized method, use 'auto', 'distinct' (Eq. 15) or 'repeated' "
            "(Eq. 16)");
    long Nsum = 0;
    for (int v : N) Nsum += v;
    if (Nsum < 0) return res;
    if (Nsum == 0) {
        res.lG = zero;
        res.G = num_traits<T>::from_int(1);
        return res;
    }
    if (L.rows() == 0 || L.cols() == 0) return res;
    if (static_cast<std::size_t>(L.cols()) != R)
        throw InputError("pfqn_explicit: the demand matrix must have one column per class of N");
    const std::size_t K = static_cast<std::size_t>(L.rows());
    for (std::size_t i = 0; i < K; ++i)
        for (std::size_t r = 0; r < R; ++r)
            if (L(i, r) < zero)
                throw InputError("pfqn_explicit: the demand matrix must be nonnegative");
    const T Nt = num_traits<T>::from_int(Nsum);

    // ---- redundancy scan: are the induced demands pairwise distinct at every t? ----
    const auto induced = [&](const std::vector<int>& t) {
        std::vector<T> th(K, zero);
        for (std::size_t i = 0; i < K; ++i)
            for (std::size_t r = 0; r < R; ++r)
                th[i] = T(th[i] + L(i, r) * num_traits<T>::from_int(t[r]));
        return th;
    };
    const auto redundant_at = [&](std::vector<T> th) {
        std::sort(th.begin(), th.end(), [](const T& a, const T& b) { return a < b; });
        const T scale = th.back();
        // scale == 0 leaves every induced demand at zero, so g_t(|N|) = 0 at
        // sum(N) > 0 and the term takes no part in the sum
        if (!(scale > zero)) return false;
        const T gap = T(num_traits<T>::from_double(tol) * scale);
        for (std::size_t i = 1; i < th.size(); ++i)
            if (!(T(th[i] - th[i - 1]) > gap)) return true;
        return false;
    };
    bool isRedundant = false;
    if (R == 1) {
        // the induced demands at t are t*L, so both the tie structure and the
        // relative tolerance are those of L itself, at every t at once
        std::vector<T> th(K);
        for (std::size_t i = 0; i < K; ++i) th[i] = L(i, 0);
        isRedundant = redundant_at(th);
    } else {
        std::vector<int> t(R, 0);
        while (true) {
            long ts = 0;
            for (int v : t) ts += v;
            if (ts > 0 && redundant_at(induced(t))) {
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
            "pfqn_explicit: Eq. (15) requires pairwise distinct induced demands, but two of them "
            "agree to within tol. Use 'auto' or 'repeated'");
    }
    res.method = expr;

    detail::SignedLse<T> total;
    if (R == 1) {
        // ---- single class: the divided difference is the identity, evaluate g ----
        std::vector<T> th(K);
        for (std::size_t i = 0; i < K; ++i) th[i] = L(i, 0);
        total = (expr == "distinct") ? detail::explicit_gdistinct(th, Nt, K)
                                     : detail::explicit_grepeated(th, Nt, K, tol);
    } else {
        // ---- outer divided-difference sum over 0 <= t <= N ----
        std::vector<T> lterm;
        std::vector<double> sterm;
        double innerLoss = 0.0;
        std::vector<int> t(R, 0);
        while (true) {
            long ts = 0;
            for (int v : t) ts += v;
            if (ts > 0) {
                std::vector<T> th = induced(t);
                T thmax = th[0];
                for (const T& v : th)
                    if (v > thmax) thmax = v;
                if (thmax > zero) {
                    const detail::SignedLse<T> g =
                        (expr == "distinct") ? detail::explicit_gdistinct(th, Nt, K)
                                             : detail::explicit_grepeated(th, Nt, K, tol);
                    innerLoss = std::max(innerLoss, g.lossDigits);
                    if (g.sgn != 0) {
                        T l = g.lS;
                        for (std::size_t r = 0; r < R; ++r) {
                            l = T(l - detail::num_factln<T>(num_traits<T>::from_int(t[r])));
                            l = T(l - detail::num_factln<T>(
                                          num_traits<T>::from_int(N[r] - t[r])));
                        }
                        lterm.push_back(l);
                        sterm.push_back(static_cast<double>(g.sgn) *
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

    // A caller that named a cancellation budget has a fallback and wants a
    // verdict, not a warning: refuse quietly. lossDigits is infinite when the
    // sum vanished identically, which is a total loss rather than a legitimate
    // G = 0.
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
        // Double precision is exhausted by cancellation; the sign itself is
        // wrong, so there is no result to hand back.
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

#endif  // LINE_API_PFQN_PFQN_EXPLICIT_H
