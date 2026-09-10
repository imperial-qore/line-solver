/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_NRE_H
#define LINE_API_PFQN_PFQN_NRE_H

/**
 * Norlund-Rice inversion of the normalizing constant on a SADDLE-TILTED
 * contour, with a second-order Edgeworth correction.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_nre.m. Two corrections separate
 * it from pfqn_nrl and pfqn_nrp, which Laplace-approximate the same integral
 * on the untilted contour X = 1:
 *
 *   1. The integrand is invariant under t -> t + c*1, since h is homogeneous
 *      of degree sum(N) in the class variables and that degree cancels
 *      against exp(-i N t). The redundant direction is quotiented out, so the
 *      integral is (R-1)-dimensional; pfqn_nrl and pfqn_nrp integrate over R
 *      dimensions and let the substitution Jacobian supply curvature along
 *      the null direction, which is an artifact of the change of variables.
 *   2. The contour radii are tilted per class to the saddle point, the X
 *      solving X_r dlog(h)/dX_r = N_r, so the origin is a stationary point of
 *      the phase. On X = 1 it is not, which is the leading bias of nrl / nrp.
 *
 * NO COMPLEX ARITHMETIC. Every integrand evaluation sits at real positive
 * demands, so unlike pfqn_nrl this routine only needs pfqn_lldsingle. It is
 * still gated on num_traits<T>::has_transcendental: the cumulant generating
 * function, the Gaussian curvature term and the tilt all take logs and exps.
 *
 * OVERFLOW CEILING. The cumulant generating function is log(G) of the tilted
 * single-class model, and this port's pfqn_lldsingle accumulates G itself
 * rather than its logarithm (the MATLAB reference switches to the log domain).
 * A double instantiation therefore loses the saddle search once log G passes
 * ~709; a Real<D> instantiation does not.
 *
 * WARNINGS. This layer has no warning channel, so the two conditions the
 * reference warns about -- a saddle search that runs out of iterations and a
 * non-positive Edgeworth correction -- take the same fallback silently: the
 * current estimate, and the bare saddlepoint term respectively.
 *
 * COST is O(I R^2 + R^4) evaluations of a single-class LLD constant, hence
 * polynomial in the class count. The fourth-cumulant tensor caps the port at
 * 8 classes, as in the reference.
 */

#include <cmath>
#include <cstddef>
#include <map>
#include <vector>

#include "line/api/pfqn/pfqn_asympt_common.h"
#include "line/api/pfqn/pfqn_gld.h"
#include "line/api/pfqn/pfqn_lldsingle.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

namespace detail {

/** Finite-difference step; the reference is flat over [5e-3,5e-2]. */
template <class T>
inline T nre_hstep() {
    return num_traits<T>::from_rational(2, 100);
}

/** Beyond 8 classes the fourth-cumulant tensor is no longer affordable. */
inline constexpr std::size_t nre_max_dim() { return 7; }

/**
 * Cumulant generating function of the tilted single-class model, memoised on
 * the finite-difference stencil around the current expansion point.
 */
template <class T>
class NreCgf {
public:
    NreCgf(const Matrix<T>& L, int Nt, const Matrix<T>& alpha, std::size_t d)
        : L_(L), Nt_(Nt), alpha_(alpha), d_(d), vbase_(d, num_traits<T>::from_int(0)) {}

    /** Move the expansion point, which invalidates every memoised value. */
    void reset(const std::vector<T>& v) {
        vbase_ = v;
        cache_.clear();
    }

    /** Value at vbase + off*hstep. */
    T at(const std::vector<int>& off) {
        const typename std::map<std::vector<int>, T>::const_iterator it = cache_.find(off);
        if (it != cache_.end()) return it->second;
        const T h = nre_hstep<T>();
        std::vector<T> v(d_);
        for (std::size_t r = 0; r < d_; ++r)
            v[r] = vbase_[r] + num_traits<T>::from_int(off[r]) * h;
        const T y = at_point(v);
        cache_.insert(std::make_pair(off, y));
        return y;
    }

    /** log of the single-class LLD constant at the class tilt X = [exp(v),1]. */
    T at_point(const std::vector<T>& v) const {
        using std::exp;
        using std::log;
        const std::size_t M = L_.rows(), R = L_.cols();
        Matrix<T> Lx(M, 1, num_traits<T>::from_int(0));
        for (std::size_t i = 0; i < M; ++i) {
            T acc = num_traits<T>::from_int(0);
            for (std::size_t r = 0; r < R; ++r)
                acc += L_(i, r) * (r < d_ ? exp(v[r]) : num_traits<T>::from_int(1));
            Lx(i, 0) = acc;
        }
        const NcResult<T> res = pfqn_lldsingle(Lx, Nt_, alpha_);
        return log(res.G);
    }

private:
    const Matrix<T>& L_;
    int Nt_;
    const Matrix<T>& alpha_;
    std::size_t d_;
    std::vector<T> vbase_;
    std::map<std::vector<int>, T> cache_;
};

/** Offset vector with sgn at coordinate a. */
inline std::vector<int> nre_unitoff(std::size_t d, std::size_t a, int sgn) {
    std::vector<int> off(d, 0);
    off[a] = sgn;
    return off;
}

/** Mixed second difference of the cumulant generating function. */
template <class T>
T nre_second_diff(NreCgf<T>& cgf, std::size_t d, std::size_t a, std::size_t b) {
    std::vector<int> pp(d, 0), pm(d, 0), mp(d, 0), mm(d, 0);
    pp[a] += 1;
    pp[b] += 1;
    pm[a] += 1;
    pm[b] -= 1;
    mp[a] -= 1;
    mp[b] += 1;
    mm[a] -= 1;
    mm[b] -= 1;
    const T h = nre_hstep<T>();
    return (cgf.at(pp) - cgf.at(pm) - cgf.at(mp) + cgf.at(mm)) /
           (num_traits<T>::from_int(4) * h * h);
}

/**
 * Cholesky factor of a symmetric matrix, or `false` if it is not positive
 * definite. Stands in for the reference's min(eig(Sigma)) <= 0 test: a
 * symmetric matrix admits a Cholesky factorization exactly when its smallest
 * eigenvalue is positive, and the factor also gives log det as 2 sum log l_ii
 * without forming the determinant.
 */
template <class T>
bool nre_chol(const Matrix<T>& A, Matrix<T>& Lo) {
    using std::sqrt;
    const std::size_t n = A.rows();
    const T zero = num_traits<T>::from_int(0);
    Lo = Matrix<T>(n, n, zero);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            T acc = A(i, j);
            for (std::size_t k = 0; k < j; ++k) acc -= Lo(i, k) * Lo(j, k);
            if (i == j) {
                if (!(acc > zero)) return false;
                Lo(i, i) = sqrt(acc);
            } else {
                Lo(i, j) = acc / Lo(j, j);
            }
        }
    }
    return true;
}

}  // namespace detail

/**
 * The reference's `[lG,G,lGs,vsad]`.
 *
 * `lG - lGs` is the Edgeworth correction, so a caller wanting the plain
 * saddlepoint estimate reads `lGs` rather than re-deriving it. `vsad` is empty
 * on the shortcut arms that never solve a saddle point, matching the
 * reference's empty `[]`.
 */
template <class T>
struct PfqnNreResult {
    T lG;                    ///< log G, Edgeworth correction included
    T lGs;                   ///< log of the saddlepoint term alone
    std::vector<T> vsad;     ///< the tilt actually used, empty when none was solved
};

/**
 * Saddle-tilted Edgeworth approximation of log G for a limited load-dependent
 * model: the full form of the reference's outputs, named alike in the JAR and
 * the native python port.
 *
 * @param L     (M x R) demands
 * @param N     (R) population
 * @param Z     (R) think times, empty for zero
 * @param alpha (M x Ntot) load-dependent rates, empty for all ones
 * @param vfix  tilt to use instead of solving the saddle-point equation, empty
 *              for the standard estimator. Supplying the tilt obtained at a
 *              nearby population makes numerator and denominator of a ratio
 *              share one expansion point, the Tierney-Kadane arrangement.
 */
template <class T>
PfqnNreResult<T> pfqn_nre_full(const Matrix<T>& L0, const std::vector<T>& N,
                               const std::vector<T>& Z, const Matrix<T>& alpha0,
                               const std::vector<T>& vfix) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_nre requires transcendental arithmetic (a saddlepoint expansion of a "
                  "coefficient-extraction integral)");
    using std::log;
    using std::sqrt;
    const std::size_t M0 = L0.rows(), R = L0.cols();
    if (N.size() != R) throw InputError("pfqn_nre: L and N disagree on the class count");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);

    T Ntsum = zero, Zsum = zero;
    for (std::size_t r = 0; r < R; ++r) Ntsum += N[r];
    for (std::size_t r = 0; r < Z.size(); ++r) Zsum += Z[r];
    if (Ntsum < zero) throw InputError("pfqn_nre: negative population");
    if (Ntsum == zero) return PfqnNreResult<T>{zero, zero, std::vector<T>()};
    const int Nt = static_cast<int>(num_traits<T>::to_double(Ntsum) + 0.5);
    const std::size_t Ntot = static_cast<std::size_t>(Nt);

    // Append the delay as an infinite-server station, and trim the rate matrix
    // so that every rate used downstream is positive.
    const std::size_t M = M0 + (Zsum > zero ? 1 : 0);
    Matrix<T> L(M, R, zero);
    for (std::size_t i = 0; i < M0; ++i)
        for (std::size_t r = 0; r < R; ++r) L(i, r) = L0(i, r);
    Matrix<T> alpha(M, Ntot, one);
    for (std::size_t i = 0; i < M0 && i < alpha0.rows(); ++i)
        for (std::size_t k = 0; k < Ntot; ++k)
            alpha(i, k) = k < alpha0.cols() ? alpha0(i, k) : one;
    if (Zsum > zero) {
        if (Z.size() != R) throw InputError("pfqn_nre: Z has the wrong length");
        for (std::size_t r = 0; r < R; ++r) L(M0, r) = Z[r];
        for (std::size_t k = 0; k < Ntot; ++k)
            alpha(M0, k) = num_traits<T>::from_int(static_cast<long>(k) + 1);
    }

    if (M == 1) {
        std::vector<int> Ni(R, 0);
        for (std::size_t r = 0; r < R; ++r)
            Ni[r] = static_cast<int>(num_traits<T>::to_double(N[r]) + 0.5);
        const T lGone = num_traits<T>::from_double(pfqn_gld(L, Ni, alpha).lG);
        return PfqnNreResult<T>{lGone, lGone, std::vector<T>()};
    }

    // Scale demands into [0,1] per class; the residual factor is exact by
    // homogeneity of the integrand.
    T lGscale = zero;
    for (std::size_t r = 0; r < R; ++r) {
        T m = zero;
        for (std::size_t i = 0; i < M; ++i)
            if (L(i, r) > m) m = L(i, r);
        if (!(m > zero)) m = one;
        for (std::size_t i = 0; i < M; ++i) L(i, r) = L(i, r) / m;
        lGscale += N[r] * log(m);
    }

    if (R == 1) {
        // coefficient extraction is the identity in a single class
        const T lGone = num_traits<T>::from_double(pfqn_lldsingle(L, Nt, alpha).lG) + lGscale;
        return PfqnNreResult<T>{lGone, lGone, std::vector<T>()};
    }

    const std::size_t d = R - 1;  // dimension of the quotient torus
    if (d > detail::nre_max_dim())
        throw UnsupportedError(
            "pfqn_nre: pfqn_nre is limited to 8 classes, use nrl or clw beyond that");

    std::vector<T> Nd(d);
    for (std::size_t a = 0; a < d; ++a) Nd[a] = N[a];
    detail::NreCgf<T> cgf(L, Nt, alpha, d);
    const T h = detail::nre_hstep<T>();
    const std::vector<int> origin(d, 0);

    // ---- saddle point: minimise the convex F(v) = K(v) - Nd*v ----
    // A tilt supplied by the caller is used as given, so that a ratio of two
    // constants can be expanded about one common point rather than two.
    std::vector<T> vbase(d, zero);
    const T tol = num_traits<T>::from_double(1e-10);
    bool converged = false;
    const bool tilt_given = !vfix.empty();
    if (tilt_given) {
        if (vfix.size() < d)
            throw InputError(
                "pfqn_nre: the supplied tilt must have one entry per quotient dimension (R-1)");
        for (std::size_t a = 0; a < d; ++a) vbase[a] = vfix[a];
        converged = true;
    }
    for (int it = 0; !tilt_given && it < 100; ++it) {
        cgf.reset(vbase);
        std::vector<T> grad(d);
        Matrix<T> hess(d, d, zero);
        for (std::size_t a = 0; a < d; ++a)
            grad[a] = (cgf.at(detail::nre_unitoff(d, a, 1)) - cgf.at(detail::nre_unitoff(d, a, -1))) /
                          (num_traits<T>::from_int(2) * h) -
                      Nd[a];
        for (std::size_t a = 0; a < d; ++a)
            for (std::size_t b = 0; b < d; ++b) hess(a, b) = detail::nre_second_diff(cgf, d, a, b);
        std::vector<T> step = solve(hess, grad);
        for (std::size_t a = 0; a < d; ++a) step[a] = -step[a];

        T F0 = cgf.at(origin);
        for (std::size_t a = 0; a < d; ++a) F0 -= Nd[a] * vbase[a];
        T tau = one;
        std::vector<T> vtry(d);
        while (tau > tol) {
            T obj = zero;
            for (std::size_t a = 0; a < d; ++a) vtry[a] = vbase[a] + tau * step[a];
            obj = cgf.at_point(vtry);
            for (std::size_t a = 0; a < d; ++a) obj -= Nd[a] * vtry[a];
            if (obj <= F0) break;
            tau = tau / num_traits<T>::from_int(2);
        }
        T stepNorm = zero;
        for (std::size_t a = 0; a < d; ++a) {
            const T delta = tau * step[a];
            vbase[a] += delta;
            stepNorm += delta * delta;
        }
        // Newton converges to the root of the DIFFERENCED gradient, whose own
        // O(hstep^2) bias puts any absolute gradient target out of reach.
        if (sqrt(stepNorm) < tol) {
            converged = true;
            break;
        }
    }
    // NO WARNING CHANNEL in this layer: the reference warns here and returns
    // its current estimate anyway, so the estimate is what the port returns.
    (void)converged;

    // ---- cumulants of the tilted distribution at the saddle ----
    cgf.reset(vbase);
    const T K0 = cgf.at(origin);
    Matrix<T> Sigma(d, d, zero);
    for (std::size_t a = 0; a < d; ++a)
        for (std::size_t b = 0; b < d; ++b) Sigma(a, b) = detail::nre_second_diff(cgf, d, a, b);
    for (std::size_t a = 0; a < d; ++a)
        for (std::size_t b = a + 1; b < d; ++b) {
            const T sym = (Sigma(a, b) + Sigma(b, a)) / num_traits<T>::from_int(2);
            Sigma(a, b) = sym;
            Sigma(b, a) = sym;
        }
    Matrix<T> chol;
    if (!detail::nre_chol(Sigma, chol))
        throw NumericError(
            "pfqn_nre: the tilted covariance is singular, a class has no demand at any station");
    T logdet = zero;
    for (std::size_t a = 0; a < d; ++a) logdet += num_traits<T>::from_int(2) * log(chol(a, a));

    std::vector<T> k3(d * d * d, zero);
    for (std::size_t a = 0; a < d; ++a)
        for (std::size_t b = 0; b < d; ++b)
            for (std::size_t c = 0; c < d; ++c) {
                T acc = zero;
                for (int s = 0; s < 8; ++s) {
                    const int s1 = 1 - 2 * ((s >> 0) & 1);
                    const int s2 = 1 - 2 * ((s >> 1) & 1);
                    const int s3 = 1 - 2 * ((s >> 2) & 1);
                    std::vector<int> off(d, 0);
                    off[a] += s1;
                    off[b] += s2;
                    off[c] += s3;
                    acc += num_traits<T>::from_int(s1 * s2 * s3) * cgf.at(off);
                }
                k3[(a * d + b) * d + c] = acc / (num_traits<T>::from_int(8) * h * h * h);
            }

    std::vector<T> k4(d * d * d * d, zero);
    for (std::size_t a = 0; a < d; ++a)
        for (std::size_t b = 0; b < d; ++b)
            for (std::size_t c = 0; c < d; ++c)
                for (std::size_t e = 0; e < d; ++e) {
                    T acc = zero;
                    for (int s = 0; s < 16; ++s) {
                        const int s1 = 1 - 2 * ((s >> 0) & 1);
                        const int s2 = 1 - 2 * ((s >> 1) & 1);
                        const int s3 = 1 - 2 * ((s >> 2) & 1);
                        const int s4 = 1 - 2 * ((s >> 3) & 1);
                        std::vector<int> off(d, 0);
                        off[a] += s1;
                        off[b] += s2;
                        off[c] += s3;
                        off[e] += s4;
                        acc += num_traits<T>::from_int(s1 * s2 * s3 * s4) * cgf.at(off);
                    }
                    k4[((a * d + b) * d + c) * d + e] =
                        acc / (num_traits<T>::from_int(16) * h * h * h * h);
                }

    // ---- second-order Edgeworth factor, see the header ----
    const Matrix<T> S = inverse(Sigma);
    T rho4 = zero;
    for (std::size_t a = 0; a < d; ++a)
        for (std::size_t b = 0; b < d; ++b)
            for (std::size_t c = 0; c < d; ++c)
                for (std::size_t e = 0; e < d; ++e)
                    rho4 += k4[((a * d + b) * d + c) * d + e] * S(a, b) * S(c, e);
    std::vector<T> u(d, zero);
    for (std::size_t c = 0; c < d; ++c)
        for (std::size_t a = 0; a < d; ++a)
            for (std::size_t b = 0; b < d; ++b) u[c] += S(a, b) * k3[(a * d + b) * d + c];
    T rhoA = zero;
    for (std::size_t c = 0; c < d; ++c)
        for (std::size_t e = 0; e < d; ++e) rhoA += u[c] * S(c, e) * u[e];
    // staged contraction of k3 against three copies of Sigma^{-1}
    std::vector<T> T1(d * d * d, zero), T2(d * d * d, zero);
    for (std::size_t i = 0; i < d; ++i)
        for (std::size_t b = 0; b < d; ++b)
            for (std::size_t c = 0; c < d; ++c) {
                T acc = zero;
                for (std::size_t a = 0; a < d; ++a) acc += S(i, a) * k3[(a * d + b) * d + c];
                T1[(i * d + b) * d + c] = acc;
            }
    for (std::size_t i = 0; i < d; ++i)
        for (std::size_t j = 0; j < d; ++j)
            for (std::size_t c = 0; c < d; ++c) {
                T acc = zero;
                for (std::size_t b = 0; b < d; ++b) acc += S(j, b) * T1[(i * d + b) * d + c];
                T2[(i * d + j) * d + c] = acc;
            }
    T rhoB = zero;
    for (std::size_t i = 0; i < d; ++i)
        for (std::size_t j = 0; j < d; ++j)
            for (std::size_t k = 0; k < d; ++k) {
                T acc = zero;
                for (std::size_t c = 0; c < d; ++c) acc += S(k, c) * T2[(i * d + j) * d + c];
                rhoB += k3[(i * d + j) * d + k] * acc;
            }
    T corr = one + rho4 / num_traits<T>::from_int(8) -
             (num_traits<T>::from_int(3) * rhoA + num_traits<T>::from_int(2) * rhoB) /
                 num_traits<T>::from_int(24);
    // Same as above: a non-positive correction falls back on the bare
    // saddlepoint term, which is the reference's behaviour after its warning.
    if (!(corr > zero)) corr = one;

    T lGs = K0 - num_traits<T>::from_double(0.5 * static_cast<double>(d) *
                                            std::log(2.0 * 3.14159265358979323846)) -
            logdet / num_traits<T>::from_int(2) + lGscale;
    for (std::size_t a = 0; a < d; ++a) lGs -= Nd[a] * vbase[a];
    return PfqnNreResult<T>{lGs + log(corr), lGs, vbase};
}

/**
 * Saddle-tilted Edgeworth approximation of log G for a limited load-dependent
 * model.
 *
 * @param L     (M x R) demands
 * @param N     (R) population
 * @param Z     (R) think times, empty for zero
 * @param alpha (M x Ntot) load-dependent rates, empty for all ones
 */
template <class T>
T pfqn_nre(const Matrix<T>& L0, const std::vector<T>& N, const std::vector<T>& Z,
           const Matrix<T>& alpha0) {
    return pfqn_nre_full(L0, N, Z, alpha0, std::vector<T>()).lG;
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PFQN_NRE_H
