/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_NPFQN_DPS_MORRISON_H
#define LINE_API_NPFQN_DPS_MORRISON_H

/**
 * Two-term heavy-usage asymptotic approximation for a closed queueing network with one
 * infinite-server (think) station and one discriminatory processor-sharing (DPS) station.
 *
 * Templated port of matlab/src/api/npfqn/npfqn_dps_morrison.m, cross-checked against
 * jar/src/main/java/jline/api/npfqn/Npfqn_dps_morrison.java (identical term for term, including
 * the sigma solve and the W_m recursion).
 *
 * Reference: J.A. Morrison, "Asymptotic analysis of a large closed queueing network with
 * discriminatory processor sharing", Queueing Systems 9 (1991) 191-214.
 *
 * The network is NOT product-form, so nothing here computes a normalizing constant: the method
 * expands the GENERATING FUNCTION of the balance equations. The substitution P(n) = <w,n> f(n)
 * clears the DPS denominator and turns the balance recursion into a linear PDE with affine
 * coefficients (eq. 2.5); rescaling z = 1 - xi/sqrt(N) and expanding in powers of N^(-1/2) leaves
 * a degenerate leading operator whose kernel is the functions of the similarity variable eta, and
 * the solvability condition along its characteristic gives an ODE for the amplitude (eq. 2.20).
 * RESULT 1 (eq. 4.11) and RESULT 2 (eq. 4.17) are the two-term approximations returned here.
 *
 * Scaling. Morrison writes K_j = N b_j and lambda_j = N r_j g_j with usage
 * rho = sum_j b_j/g_j = 1 - a/sqrt(N). N is bookkeeping only and the approximation is invariant to
 * it, so this routine fixes N = 1: b = N_pop, g = Z/S, r = 1/Z, a = 1 - rho. Accuracy is governed
 * by the PHYSICAL regime -- large populations with rho near 1. rho > 1 is admissible, being the
 * saturated regime of appendix A.
 *
 * Arithmetic. The W_m of eq. (3.23) need an erfc and an exp, so this requires transcendental
 * arithmetic and cannot be instantiated at T = Rational.
 */

#include <cstddef>
#include <vector>

#include "line/api/npfqn/npfqn_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace npfqn {

/** Mean queue lengths, sojourn times and throughputs, with Morrison's intermediate constants. */
template <class T>
struct DpsMorrisonResult {
    std::vector<T> Q;      ///< mean number of class-k jobs at the DPS station
    std::vector<T> R;      ///< mean class-k sojourn time per visit to the DPS station
    std::vector<T> X;      ///< per-class throughput
    std::vector<T> Qlead;  ///< leading-order (one-term) queue lengths
    std::vector<T> Rlead;  ///< leading-order (one-term) sojourn times
    std::vector<T> sigma;  ///< the vector sigma of eq. (4.12)
    std::vector<T> W;      ///< W_0..W_4 of eq. (3.23)
    T rho;                 ///< usage sum_j b_j/g_j
    T a;                   ///< heavy-usage parameter a = 1 - rho at the N = 1 scale
    T cB, cC, cD, cH, cI, cJ, cK, cL, cM, cQ, delta, cR, cS, cU, cA, cV;  ///< Morrison's constants
};

namespace detail {

/**
 * Scaled complementary error function exp(x^2) erfc(x), for either sign. The direct product
 * overflows past x ~ 26, where the asymptotic series is already exact to double precision; for
 * x < 0 the reflection erfcx(x) = 2 exp(x^2) - erfcx(-x) is used, which overflows below x ~ -26
 * and is reported as such by the caller.
 */
template <class T>
T num_erfcx(const T& x) {
    const T zero = num_traits<T>::from_int(0);
    if (x < zero) {
        return num_traits<T>::from_int(2) * num_exp(T(x * x)) - num_erfcx(T(-x));
    }
    if (x < num_traits<T>::from_int(25)) {
        return num_exp(T(x * x)) * num_erfc(x);
    }
    const T y = num_traits<T>::from_int(1) / (num_traits<T>::from_int(2) * x * x);
    T term = num_traits<T>::from_int(1);
    T sum = num_traits<T>::from_int(1);
    for (int k = 1; k <= 12; ++k) {
        term *= -num_traits<T>::from_int(2 * k - 1) * y;
        sum += term;
    }
    return sum / (x * num_sqrt(boost::math::constants::pi<T>()));
}

/** I_m(yh) = int_0^inf s^m exp(-s^2/2 - yh s) ds by composite Simpson, the recursion's fallback. */
template <class T>
T dps_quad_I(std::size_t m, const T& yh) {
    const T zero = num_traits<T>::from_int(0);
    const T twelve = num_traits<T>::from_int(12);
    T hi = twelve;
    if (yh >= num_traits<T>::from_int(1)) {
        const T shrunk = num_traits<T>::from_int(40) / yh;
        hi = (shrunk < twelve) ? shrunk : twelve;
    } else if (yh < zero) {
        hi = twelve - yh;
    }
    const int n = 4096;
    const T h = hi / num_traits<T>::from_int(n);
    T sum = zero;
    for (int k = 0; k <= n; ++k) {
        const T s = num_traits<T>::from_int(k) * h;
        T f = num_exp(T(-s * s / num_traits<T>::from_int(2) - yh * s));
        for (std::size_t j = 0; j < m; ++j) f *= s;
        const T wgt = (k == 0 || k == n) ? num_traits<T>::from_int(1)
                                         : num_traits<T>::from_int((k % 2 == 1) ? 4 : 2);
        sum += wgt * f;
    }
    return sum * h / num_traits<T>::from_int(3);
}

/** Gaussian elimination with partial pivoting on a small dense square system. */
template <class T>
std::vector<T> dps_solve_square(std::vector<std::vector<T> > A, std::vector<T> rhs) {
    const std::size_t n = rhs.size();
    for (std::size_t i = 0; i < n; ++i) A[i].push_back(rhs[i]);
    for (std::size_t c = 0; c < n; ++c) {
        std::size_t piv = c;
        for (std::size_t i = c + 1; i < n; ++i) {
            T ai = A[i][c] < num_traits<T>::from_int(0) ? T(-A[i][c]) : A[i][c];
            T ap = A[piv][c] < num_traits<T>::from_int(0) ? T(-A[piv][c]) : A[piv][c];
            if (ai > ap) piv = i;
        }
        T ap = A[piv][c] < num_traits<T>::from_int(0) ? T(-A[piv][c]) : A[piv][c];
        if (!(ap > num_traits<T>::from_double(1e-300))) {
            throw InputError("npfqn_dps_morrison: singular system while solving Morrison's "
                             "sigma equations (4.12).");
        }
        A[c].swap(A[piv]);
        for (std::size_t i = c + 1; i < n; ++i) {
            const T f = A[i][c] / A[c][c];
            for (std::size_t j = c; j <= n; ++j) A[i][j] -= f * A[c][j];
        }
    }
    std::vector<T> x(n, num_traits<T>::from_int(0));
    for (std::size_t ii = n; ii-- > 0;) {
        T acc = A[ii][n];
        for (std::size_t j = ii + 1; j < n; ++j) acc -= A[ii][j] * x[j];
        x[ii] = acc / A[ii][ii];
    }
    return x;
}

/** W_m of eq. (3.23), m = 0..4; see the MATLAB local_W for the normalization used. */
template <class T>
std::vector<T> dps_W(const T& cB, const T& cC, const T& cD, const T& y) {
    const std::size_t mmax = 4;
    const T sig = num_sqrt(T(cD / (cB * cC)));
    const T yh = y * num_sqrt(T(cB / (cC * cD)));
    const T two = num_traits<T>::from_int(2);

    std::vector<T> Iv(mmax + 1, num_traits<T>::from_int(0));
    Iv[0] = num_sqrt(T(boost::math::constants::pi<T>() / two)) * num_erfcx(T(yh / num_sqrt(two)));
    if (!num_isfinite(Iv[0])) {
        throw InputError("npfqn_dps_morrison: the usage is so far above saturation that the "
                         "expansion overflows. This model is outside the moderately-heavy regime the "
                         "approximation is derived for; use SolverFLD, SolverMVA or SolverCTMC.");
    }
    Iv[1] = num_traits<T>::from_int(1) - yh * Iv[0];
    for (std::size_t m = 2; m <= mmax; ++m) {
        Iv[m] = num_traits<T>::from_int(static_cast<int>(m) - 1) * Iv[m - 2] - yh * Iv[m - 1];
    }
    bool positive = true;
    for (std::size_t m = 0; m <= mmax; ++m) {
        if (!(Iv[m] > num_traits<T>::from_int(0))) positive = false;
    }
    if (!positive) {
        for (std::size_t m = 0; m <= mmax; ++m) Iv[m] = dps_quad_I(m, yh);
    }
    std::vector<T> W(mmax + 1, num_traits<T>::from_int(0));
    T sp = sig;
    for (std::size_t m = 0; m <= mmax; ++m) {
        W[m] = (cB / cD) * (cB / cD) * sp * Iv[m];
        sp *= sig;
    }
    return W;
}

}  // namespace detail

/**
 * Evaluates Morrison's two-term approximation.
 *
 * @param N per-class populations, finite and positive
 * @param Z per-class mean think times, finite and positive
 * @param S per-class mean DPS service times, finite and positive
 * @param w per-class DPS weights, finite and positive
 * @return the mean performance measures and the intermediate constants
 */
template <class T>
DpsMorrisonResult<T> npfqn_dps_morrison(const std::vector<T>& N, const std::vector<T>& Z,
                                        const std::vector<T>& S, const std::vector<T>& w) {
    static_assert(num_traits<T>::has_transcendental,
                  "npfqn_dps_morrison requires transcendental arithmetic");
    const std::size_t p = N.size();
    if (Z.size() != p || S.size() != p || w.size() != p) {
        throw InputError("npfqn_dps_morrison: N, Z, S and w must have the same number of classes.");
    }
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T three = num_traits<T>::from_int(3);

    for (std::size_t i = 0; i < p; ++i) {
        if (!detail::num_isfinite(N[i]) || !(N[i] > zero)) {
            throw InputError("npfqn_dps_morrison: the approximation requires finite positive class "
                             "populations.");
        }
        if (!detail::num_isfinite(Z[i]) || !(Z[i] > zero) || !detail::num_isfinite(S[i]) ||
            !(S[i] > zero)) {
            throw InputError("npfqn_dps_morrison: think times Z and DPS service times S must be "
                             "finite and positive.");
        }
        if (!detail::num_isfinite(w[i]) || !(w[i] > zero)) {
            throw InputError("npfqn_dps_morrison: DPS weights must be finite and positive.");
        }
    }

    // Morrison's parameters at the bookkeeping scale N = 1
    std::vector<T> b(N), r(p, zero), g(p, zero);
    T rho = zero;
    for (std::size_t i = 0; i < p; ++i) {
        r[i] = one / Z[i];
        g[i] = Z[i] / S[i];
        rho += b[i] / g[i];
    }
    const T a = one - rho;

    // constants, eqs. (2.18), (2.19), (3.11), (3.13)-(3.15)
    T cB = zero, cC = zero, cD = zero, cH = zero, cI = zero, cJ = zero;
    T cK = zero, cL = zero, cM = zero, cQ = zero;
    for (std::size_t i = 0; i < p; ++i) {
        const T g2 = g[i] * g[i], g3 = g2 * g[i], w2 = w[i] * w[i], r2 = r[i] * r[i];
        cB += b[i] / (r[i] * g2 * w[i]);
        cC += b[i] / (g2 * w[i]);
        cD += b[i] / (r[i] * g2);
        cH += b[i] / (r2 * g3 * w[i]);
        cI += b[i] / (r2 * g3 * w2);
        cJ += b[i] / (r[i] * g3 * w2);
        cK += b[i] / (g3 * w2);
        cL += b[i] / (r[i] * g3 * w[i]);
        cM += b[i] / (r2 * g3);
        cQ += b[i] / g2;
    }

    // sigma: eq. (4.12) with the normalization (4.13). The p equations have rank p-1, so the last
    // one -- implied by the others -- is REPLACED by (4.13), giving a square nonsingular system.
    // All four codebases use this same scheme so their sigma agree.
    std::vector<std::vector<T> > A(p, std::vector<T>(p, zero));
    std::vector<T> rhs(p, zero);
    for (std::size_t i = 0; i < p; ++i) {
        A[i][i] += rho;
        for (std::size_t j = 0; j < p; ++j) {
            const T den = r[i] * g[i] * w[i] + r[j] * g[j] * w[j];
            A[i][i] -= w[j] * b[j] * r[j] / den;
            A[i][j] -= w[j] * b[i] * r[i] / den;
        }
        rhs[i] = rho * (b[i] / g[i]) * (cD / (cB * w[i]) - one);
    }
    for (std::size_t j = 0; j < p; ++j) A[p - 1][j] = one / (r[j] * g[j]);
    rhs[p - 1] = zero;
    const std::vector<T> sigma = detail::dps_solve_square(A, rhs);

    // alpha from eq. (4.9), then delta of eq. (3.15)
    T delta = zero;
    for (std::size_t i = 0; i < p; ++i) {
        const T alpha_i = sigma[i] - (b[i] / g[i]) * (cD / (cB * w[i]) - one);
        delta += alpha_i / g[i];
    }

    // eqs. (3.19)-(3.21)
    const T cR = three * (cB * cL - cD * cJ) / (cB * cD);
    const T cS = (two * cB * (cD * cH - cB * cM) - cD * (cD * cI - cB * cH)) / (two * cB * cB * cD * cD);
    const T cU = (cQ - cC * cD / cB - delta) / rho - cD * cR / cB + (a * a - cC * cD / cB) * cS;
    const T cA = cS * cC * cC + cR * cC - cK;
    const T cV = cR + two * cS * cC;

    const std::vector<T> W = detail::dps_W(cB, cC, cD, a);

    // RESULT 1 (4.11) and RESULT 2 (4.17), at sqrt(N) = 1. NOTE the numerator bracket carries
    // U*W2: eq. (4.10) of the paper misprints it as U*W1, but (4.7), (4.11), (A6) and (B2) all
    // agree on U*W2, and it is what the derivation from (4.4)-(4.9) gives.
    const T eps = cB / cD;
    const T num = W[1] - eps * (cA / three * W[4] + a / two * cV * W[3] + cU * W[2]);
    const T den = W[0] - eps * (cA / three * W[3] + a / two * cV * W[2] + cU * W[1] + cS);
    if (!detail::num_isfinite(den) || den == zero) {
        throw InputError("npfqn_dps_morrison: the expansion is degenerate for this model (vanishing "
                         "denominator); the usage is too far from the moderately-heavy regime.");
    }

    DpsMorrisonResult<T> res;
    res.Q.resize(p); res.R.resize(p); res.X.resize(p);
    res.Qlead.resize(p); res.Rlead.resize(p); res.sigma = sigma; res.W = W;
    for (std::size_t j = 0; j < p; ++j) {
        const T gw = g[j] * w[j];
        res.Qlead[j] = b[j] * W[1] / (gw * W[0]);
        res.Q[j] = b[j] * num / (gw * den) - b[j] * W[2] / (g[j] * g[j] * w[j] * w[j] * W[0]) -
                   sigma[j] / rho;
        res.Rlead[j] = W[1] / (r[j] * gw * W[0]);
        res.R[j] = num / (r[j] * gw * den) +
                   ((W[1] / W[0]) * (W[1] / W[0]) - W[2] / W[0]) /
                       (r[j] * g[j] * g[j] * w[j] * w[j]) -
                   sigma[j] / (rho * r[j] * b[j]);
        res.X[j] = r[j] * (b[j] - res.Q[j]);
    }
    res.rho = rho; res.a = a;
    res.cB = cB; res.cC = cC; res.cD = cD; res.cH = cH; res.cI = cI; res.cJ = cJ;
    res.cK = cK; res.cL = cL; res.cM = cM; res.cQ = cQ; res.delta = delta;
    res.cR = cR; res.cS = cS; res.cU = cU; res.cA = cA; res.cV = cV;
    return res;
}

}  // namespace npfqn
}  // namespace line

#endif  // LINE_API_NPFQN_DPS_MORRISON_H
