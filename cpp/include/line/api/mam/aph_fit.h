/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_APH_FIT_H
#define LINE_API_MAM_APH_FIT_H

/**
 * Minimal-order acyclic phase-type fit of the first three moments
 * (matlab/lib/kpctoolbox/aph/aph_fit.m).
 *
 * Implements A. Bobbio, A. Horvath, M. Telek, "Matching three moments with
 * minimal acyclic phase type distributions", Stochastic Models 21:303-326,
 * 2005. The routine searches the smallest order n <= nmax whose normalized
 * moment region contains (n2, n3) = (e2/e1^2, e3/(e1 e2)) and then evaluates
 * one of the paper's two closed forms for the canonical APH(n).
 *
 * Gated on transcendental arithmetic: both cases take square roots of moment
 * discriminants, the order search takes square roots of the region bounds,
 * and the second case additionally takes cube roots (of complex arguments).
 * None of these has an exact rational counterpart.
 *
 * The second case follows MATLAB and the JAR (Aph_fit.java) in evaluating the
 * chain K9..K22 in complex arithmetic: several of those radicands are negative
 * for feasible moment sets and the imaginary parts cancel in f, so real
 * arithmetic would produce NaN. Only the real part of f is used, as in both
 * references.
 *
 * Reference defect carried over deliberately: the third branch of the f
 * selection reads `n3 == 2*n2/2`, i.e. n3 == n2, where the surrounding
 * branches make `3*n2/2` the intended boundary. MATLAB and the JAR agree on
 * the text, so the port reproduces it rather than silently "fixing" the
 * reference; the branch is unreachable in practice because exact equality of
 * two computed doubles is required to enter it.
 */

#include <cstddef>
#include <limits>
#include <vector>

#include "line/api/mam/map_fit_detail.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_transform.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/** Result of aph_fit. */
template <class T>
struct AphFitResult {
    Map<T> aph;      ///< the fitted APH as a MAP pair (D0, D1)
    bool isexact;    ///< false when the moment set had to be relaxed
    unsigned order;  ///< number of phases actually used
};

namespace fitdetail {

/** Canonical bidiagonal APH assembled from (alpha, T) as in aph_fit.m. */
template <class T>
Map<T> aph_canonical(const Matrix<T>& Tm, const std::vector<T>& alpha, const T& e1) {
    const std::size_t n = Tm.rows();
    Map<T> m;
    m.D0 = Tm;
    m.D1 = Matrix<T>(n, n, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < n; ++i) {
        T rs = num_traits<T>::from_int(0);
        for (std::size_t k = 0; k < n; ++k) rs -= Tm(i, k);
        for (std::size_t j = 0; j < n; ++j) m.D1(i, j) = rs * alpha[j];
    }
    return map_scale(map_normalize(m), e1);
}

/**
 * Slack for the order search's boundary tests, derived from the number type.
 *
 * WHY THIS EXISTS. The Bobbio-Horvath-Telek feasibility conditions compare the
 * normalized third moment against bounds ln(n) and un(n) that a moment set can
 * ATTAIN exactly: the moments of an Erlang, and of anything else sitting on the
 * APH(n) boundary, make n3 equal to a bound rather than merely close to it.
 * Both sides of such a comparison are computed with rounding, so their
 * difference is pure noise of the order of the working epsilon, with an
 * arbitrary sign, and an exact >= or <= then decides feasibility by coin toss.
 *
 * Measured on (e1, e2, e3) = (1/3, 8/45, 2/15), where n3 = 9/4 and un(2) = 9/4
 * exactly: un - n3 came out as +8.88e-16 at double, -1.07e-50 at Real50 and
 * exactly 0 at Real100, so the fitted order was 2, then 3, then 2 again. The
 * flip is not monotone in the precision because the residual is noise, not a
 * bias, which is precisely why a frozen constant cannot fix it and a
 * precision-relative one can.
 *
 * The slack is zero in an exact field, where the comparison is decidable and
 * the boundary is attained exactly.
 */
template <class T>
T aph_boundary_slack() {
    if constexpr (num_traits<T>::is_exact) {
        return num_traits<T>::from_int(0);
    } else {
        return T(std::numeric_limits<T>::epsilon() * num_traits<T>::from_int(16));
    }
}

/** a >= b, admitting an attained boundary within the slack. */
template <class T>
bool ge_slack(const T& a, const T& b, const T& slack) {
    const T one = num_traits<T>::from_int(1);
    const T ab = num_abs(b);
    const T scale = ab > one ? ab : one;
    return a >= T(b - slack * scale);
}

/** a <= b, admitting an attained boundary within the slack. */
template <class T>
bool le_slack(const T& a, const T& b, const T& slack) {
    const T one = num_traits<T>::from_int(1);
    const T ab = num_abs(b);
    const T scale = ab > one ? ab : one;
    return a <= T(b + slack * scale);
}

/**
 * Read a radicand that is provably nonnegative where it is consumed.
 *
 * SAME PHENOMENON AS aph_boundary_slack, ONE LEVEL DOWN. The bound formulas
 * below do not only COMPARE against an attained boundary, they take square
 * roots of quantities that VANISH on it: at n2 = (n+1)/n the radicand of un is
 * exactly 1 + n(n2-2)/(n-1) = 0, and case 1's discriminant is 0 on the same
 * moment sets. Each is nonnegative on the whole region where its value is used
 * -- that is a property of the Bobbio-Horvath-Telek feasibility region, not an
 * empirical observation -- so a negative value within rounding noise of zero is
 * noise and nothing else, and reading it as zero is the only reading consistent
 * with the mathematics at every precision.
 *
 * Measured on the Erlang(2) moment set (1/3, 1/6, 1/9), where n2 = 3/2 sits ON
 * the n2 >= (n+1)/n boundary: the un radicand is exactly 0 at double, -1.07e-50
 * at Real50 and exactly 0 at Real100. Untreated, Real50 replaced the true bound
 * un = 2 with the placeholder 0, the n3 <= un test failed and the fit returned
 * order 3 for a distribution that IS an APH(2). The sign of that residual is
 * noise, so, as with the comparisons, only a precision-derived rule fixes it.
 *
 * `scale` is the magnitude of the terms whose cancellation produced x. A
 * negative value beyond the slack is NOT noise and is returned unchanged, so
 * the caller's own infeasibility handling still sees it.
 */
template <class T>
T read_nonneg_radicand(const T& x, const T& slack, const T& scale) {
    const T zero = num_traits<T>::from_int(0);
    if (x >= zero) return x;
    const T one = num_traits<T>::from_int(1);
    const T as = num_abs(scale);
    const T s = as > one ? as : one;
    return x >= T(-slack * s) ? zero : x;
}

}  // namespace fitdetail

/**
 * Fit an APH(n) with n <= nmax to the raw moments e1, e2, e3.
 *
 * tol is the tolerance used only for the exponential degeneracy screen
 * (scv == 1 with the matching third moment), where the general APH(2)
 * formulas divide by zero.
 *
 * The order search's own boundary tests carry a separate, precision-derived
 * slack (fitdetail::aph_boundary_slack) so that a moment set lying ON an
 * APH(n) bound is accepted at every arithmetic rather than at whichever ones
 * happen to round the residual the right way. See that function for the
 * measurement that motivates it.
 */
template <class T>
AphFitResult<T> aph_fit(const T& e1, const T& e2, const T& e3, unsigned nmax, const T& tol) {
    static_assert(num_traits<T>::has_transcendental, "aph_fit requires transcendental arithmetic");
    using fitdetail::Cplx;
    using fitdetail::cplx_pow_real;
    using fitdetail::cplx_sqrt;
    using fitdetail::num_sqrt;
    using fitdetail::pw;

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T three = num_traits<T>::from_int(3);
    const T six = num_traits<T>::from_int(6);

    if (e1 <= zero) throw InputError("aph_fit: the first moment must be positive");
    if (nmax < 2) throw InputError("aph_fit: nmax must be at least 2");

    AphFitResult<T> res;
    res.isexact = true;
    const T slack = fitdetail::aph_boundary_slack<T>();

    const T scv = (e2 - e1 * e1) / (e1 * e1);
    if (num_abs(T(scv - one)) < tol && num_abs(T(e3 - six * pw(e1, 3))) < tol) {
        res.aph = map_exponential_mean(e1);
        res.order = 1;
        return res;
    }

    T n2 = e2 / (e1 * e1);
    T n3 = e3 / e1 / e2;

    bool n2_feas = false, n3_ubfeas = false, n3_lbfeas = false;
    unsigned n = 1;
    T un = zero;
    T un_1 = zero;
    while ((!n2_feas || !n3_lbfeas || !n3_ubfeas) && n < nmax) {
        ++n;
        const T nn = num_traits<T>::from_int(static_cast<long>(n));
        un_1 = un;
        // bound radicand boundary rationale: see _kb/03-api-layer.md (cpp port notes: mam)
        const T uarg = fitdetail::read_nonneg_radicand(
            T(one + nn * (n2 - two) / (nn - one)), slack, T(nn * (n2 - two) / (nn - one)));
        un = uarg < zero ? zero
                         : (one / (nn * nn * n2)) *
                               (two * (nn - two) * (nn * n2 - nn - one) * num_sqrt(uarg) +
                                (nn + two) * (three * nn * n2 - two * nn - two));

        if (fitdetail::ge_slack(n2, T((nn + one) / nn), slack) &&
            fitdetail::le_slack(n2, T((nn + num_traits<T>::from_int(4)) / (nn + one)), slack)) {
            n2_feas = true;
            // On this interval 4(n+1) - 3 n n2 >= (n-2)^2/(n+1) >= 0, so both
            // radicands below are real.
            const T pn = ((nn + one) * (n2 - two) / (three * n2 * (nn - one))) *
                         (-two * num_sqrt(T(nn + one)) /
                              num_sqrt(T(num_traits<T>::from_int(4) * (nn + one) - three * nn * n2)) -
                          one);
            // boundary rationale: see _kb/03-api-layer.md (cpp port notes: mam)
            const T arad = fitdetail::read_nonneg_radicand(
                T(pn * pn + pn * nn * (n2 - two) / (nn - one)), slack, T(pn * pn));
            const T an = (n2 - two) / (pn * (one - n2) + num_sqrt(arad));
            const T ln = ((three + an) * (nn - one) + two * an) / ((nn - one) * (one + an * pn)) -
                         (two * an * (nn + one)) /
                             (two * (nn - one) + an * pn * (nn * an + two * nn - two));
            if (fitdetail::ge_slack(n3, ln, slack)) n3_lbfeas = true;
        } else if (fitdetail::ge_slack(n2, T((nn + num_traits<T>::from_int(4)) / (nn + one)),
                                       slack)) {
            n2_feas = true;
            if (fitdetail::ge_slack(n3, T(n2 * (nn + one) / nn), slack)) n3_lbfeas = true;
        }
        if (fitdetail::ge_slack(n2, T((nn + one) / nn), slack) &&
            fitdetail::le_slack(n2, T(nn / (nn - one)), slack)) {
            n2_feas = true;
            if (fitdetail::le_slack(n3, un, slack)) n3_ubfeas = true;
        } else if (fitdetail::ge_slack(n2, T(nn / (nn - one)), slack)) {
            n2_feas = true;
            n3_ubfeas = true;  // MATLAB: n3 < Inf, always true for a finite n3
        }
    }

    if (!n2_feas || !n3_lbfeas || !n3_ubfeas || n == nmax) {
        const T nn = num_traits<T>::from_int(static_cast<long>(n));
        n2 = (nn + one) / nn;
        n3 = two * n2 - one;
        res.isexact = false;
    }

    const T nn = num_traits<T>::from_int(static_cast<long>(n));
    res.order = n;

    std::vector<T> alpha(n, zero);
    Matrix<T> Tm(n, n, zero);
    const T lambda = one;

    if (n2 <= nn / (nn - one) || n3 <= two * n2 - one) {
        // discriminant boundary rationale: see _kb/03-api-layer.md (cpp port notes: mam)
        const T rad = fitdetail::read_nonneg_radicand(
            T(num_traits<T>::from_int(12) * n2 * n2 * (nn + one) +
              num_traits<T>::from_int(16) * n3 * (nn + one) +
              n2 * (nn * (n3 - num_traits<T>::from_int(15)) * (n3 + one) -
                    num_traits<T>::from_int(8) * (n3 + three))),
            slack, T(num_traits<T>::from_int(12) * n2 * n2 * (nn + one)));
        if (rad < zero) throw NumericError("aph_fit: negative discriminant in case 1");
        const T b = two * (num_traits<T>::from_int(4) - nn * (three * n2 - num_traits<T>::from_int(4))) /
                    (n2 * (num_traits<T>::from_int(4) + nn - nn * n3) + num_sqrt(T(nn * n2)) * num_sqrt(rad));
        const T a = (b * n2 - two) * (nn - one) * b / ((b - one) * nn);
        if (a == zero) throw NumericError("aph_fit: degenerate case-1 parameter");
        const T p = (b - one) / a;
        const T mu = lambda * (nn - one) / a;
        alpha[0] = p;
        alpha[n - 1] = one - p;
        for (std::size_t i = 0; i < n; ++i) Tm(i, i) = -mu;
        for (std::size_t i = 0; i + 1 < n; ++i) Tm(i, i + 1) = mu;
        Tm(n - 1, n - 1) = -lambda;
    } else if (n2 > nn / (nn - one) && n3 > un_1) {
        // case 2 of 2: the Bobbio-Horvath-Telek chain, evaluated in C
        const T K1 = nn - one;
        const T K2 = nn - two;
        const T K3 = three * n2 - two * n3;
        const T K4 = n3 - three;
        const T K5 = nn - n2;
        const T K6 = one + n2 - n3;
        const T K7 = nn + n2 - nn * n2;
        const T K8 = three + three * n2 * n2 + n3 - three * n2 * n3;
        if (K3 == zero || K4 == zero) throw NumericError("aph_fit: degenerate case-2 parameter");

        const T k9a = num_traits<T>::from_int(4) * K1 * pw(K5, 3) +
                      K1 * K1 * K2 * K4 * K4 * nn * n2 * n2 +
                      num_traits<T>::from_int(4) * K2 * nn * n2 *
                          (K4 * nn * nn - three * K6 * n2 + K8 * nn);
        const Cplx<T> K9inner =
            cplx_sqrt(Cplx<T>(T(-num_traits<T>::from_int(16) * K1 * K1 * pw(K7, 6) + k9a * k9a)));
        const Cplx<T> K9 =
            (K9inner + Cplx<T>(T(num_traits<T>::from_int(4) * K2 * K2 * K3 * nn * nn * n2 +
                                 K1 * K1 * K2 * K4 * K4 * nn * n2 * n2 +
                                 num_traits<T>::from_int(4) * K1 * K5 *
                                     (K5 * K5 - three * K2 * K6 * nn * n2)))) *
            T(num_traits<T>::from_int(108) * K1 * K1);
        const T K10 = K4 * K4 / (num_traits<T>::from_int(4) * K3 * K3) - K5 / (K1 * K3 * n2);
        const T third = one / three;
        const Cplx<T> K9cbrt = cplx_pow_real(K9, third);
        const T c2 = fitdetail::num_exp(T(third * fitdetail::num_log(two)));         // 2^(1/3)
        const T c7 = fitdetail::num_exp(T((num_traits<T>::from_int(7) / three) * fitdetail::num_log(two)));
        const Cplx<T> K11 =
            fitdetail::cplx_inv(K9cbrt) * T(c2 * (three * K5 * K5 + K2 * (K3 + two * K4) * nn * n2) / (K3 * n2));
        const Cplx<T> K12 = K9cbrt / T(three * c7 * K1 * K1 * K3 * n2);
        const Cplx<T> K13 = cplx_sqrt(K11 + K12 + Cplx<T>(K10));
        const T K20 = six * K1 * K3 * K4 * K5 + num_traits<T>::from_int(4) * K2 * K3 * K3 * nn -
                      K1 * K1 * pw(K4, 3) * n2;
        const Cplx<T> K14 =
            fitdetail::cplx_inv(K13) * T(K20 / (num_traits<T>::from_int(4) * K1 * K1 * pw(K3, 3) * n2));
        const T K15 = -K4 / (two * K3);
        const Cplx<T> twoK10(T(two * K10));
        const Cplx<T> K16 = cplx_sqrt(twoK10 - K11 - K12 - K14);
        const Cplx<T> K17 = cplx_sqrt(twoK10 - K11 - K12 + K14);
        const T k18a = num_traits<T>::from_int(4) * pw(K5, 3) +
                       num_traits<T>::from_int(4) * K2 * K4 * K5 * nn * n2 +
                       K1 * K2 * K4 * K4 * nn * n2 * n2;
        const Cplx<T> K18 =
            -cplx_sqrt(Cplx<T>(T(num_traits<T>::from_int(81) * k18a * k18a -
                                 num_traits<T>::from_int(48) *
                                     pw(T(three * K5 * K5 + two * K2 * K4 * nn * n2), 3)))) +
            Cplx<T>(T(num_traits<T>::from_int(36) * pw(K5, 3) +
                      num_traits<T>::from_int(36) * K2 * K4 * K5 * nn * n2 +
                      num_traits<T>::from_int(9) * K1 * K2 * K4 * K4 * nn * n2 * n2));
        const Cplx<T> K18cbrt = cplx_pow_real(K18, third);
        const T c23 = fitdetail::num_exp(T((two / three) * fitdetail::num_log(two)));   // 2^(2/3)
        const T c31 = fitdetail::num_exp(T(third * fitdetail::num_log(three)));         // 3^(1/3)
        const T c62 = fitdetail::num_exp(T((two / three) * fitdetail::num_log(six)));   // 6^(2/3)
        const Cplx<T> K19 =
            Cplx<T>(T(-K5 / (K1 * K4 * n2))) -
            fitdetail::cplx_inv(K18cbrt) *
                T(c23 * (three * K5 * K5 + two * K2 * K4 * nn * n2) / (c31 * K1 * K4 * n2)) -
            K18cbrt / T(c62 * K1 * K4 * n2);
        const Cplx<T> K21 = K11 + K12 + Cplx<T>(T(K5 / (two * nn * K1 * K3)));
        const Cplx<T> K22 = cplx_sqrt(
            cplx_sqrt(K21 * K21 * num_traits<T>::from_int(4) - Cplx<T>(T(nn * K2 / (n2 * K1 * K1 * K3)))) +
            Cplx<T>(T(three * K4 * K4 / (num_traits<T>::from_int(4) * K3 * K3) - three * K5 / (K1 * K3 * n2))));

        Cplx<T> fC;
        bool picked = false;
        if (n3 > un_1 && n3 < three * n2 / two) {
            fC = K13 + Cplx<T>(K15) - K17;
            picked = true;
        } else if (n3 == two * n2 / two) {  // reference text; see header comment
            fC = K19;
            picked = true;
        } else if (n3 > three * n2 / two && K20 > zero) {
            fC = -K13 + Cplx<T>(K15) + K16;
            picked = true;
        } else if (K20 == zero) {
            fC = K22 + Cplx<T>(K15);
            picked = true;
        } else if (K20 < zero) {
            fC = K13 + Cplx<T>(K15) + K17;
            picked = true;
        }
        if (!picked) throw NumericError("aph_fit: no branch of the case-2 selection applies");
        const T f = fC.re;
        const T denom = (nn - one) * (n2 * f * f - two * f + two) - nn;
        if (denom == zero) throw NumericError("aph_fit: degenerate case-2 denominator");
        const T a = two * (f - one) * (nn - one) / denom;
        if (a == zero) throw NumericError("aph_fit: degenerate case-2 parameter");
        const T p = (f - one) * a;
        const T mu = lambda * (nn - one) / a;
        alpha[0] = p;
        alpha[1] = one - p;
        for (std::size_t i = 0; i < n; ++i) Tm(i, i) = -mu;
        for (std::size_t i = 0; i + 1 < n; ++i) Tm(i, i + 1) = mu;
        Tm(0, 0) = -lambda;
        Tm(0, 1) = lambda;
    } else {
        throw NumericError("aph_fit: the moment set cannot be matched with an APH distribution");
    }

    res.aph = fitdetail::aph_canonical(Tm, alpha, e1);
    return res;
}

/** aph_fit with the MATLAB defaults nmax = 10 and a 1e-12 degeneracy tolerance. */
template <class T>
AphFitResult<T> aph_fit(const T& e1, const T& e2, const T& e3) {
    return aph_fit(e1, e2, e3, 10u, T(num_traits<T>::from_double(1e-12)));
}

/** aph_fit with an explicit order cap and the default degeneracy tolerance. */
template <class T>
AphFitResult<T> aph_fit(const T& e1, const T& e2, const T& e3, unsigned nmax) {
    return aph_fit(e1, e2, e3, nmax, T(num_traits<T>::from_double(1e-12)));
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_APH_FIT_H
