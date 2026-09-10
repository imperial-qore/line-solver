/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MC_CTMC_SADDLEPOINT_H
#define LINE_API_MC_CTMC_SADDLEPOINT_H

/**
 * Saddlepoint approximation of Pr{N(t)=k} for the counting process of a MAP.
 *
 * Templated port of matlab/src/api/mc/ctmc_saddlepoint.m. The probability that
 * the Markovian arrival process (D0,D1) records exactly k events in (0,t],
 * obtained by steepest-descent inversion of the counting generating function
 * instead of by forming the k-th superdiagonal block of expm(t*X).
 *
 * The counting generating function is the matrix exponential
 *
 *   sum_k P(k,t) z^k = expm(t*(D0 + z*D1)),
 *
 * so the cumulant generating function of N(t) is eta(theta) = spectral abscissa
 * of A(theta) = D0 + exp(theta)*D1, the Perron root of an irreducible Metzler
 * matrix: real, simple, strictly convex in theta, with eta(0)=0 and
 * eta'(0)=lambda. Inverting by steepest descent gives Daniels (1954),
 *
 *   Pr{N(t)=k} ~ g(theta*) exp(t eta(theta*) - k theta*)
 *                / sqrt(2 pi t eta''(theta*)),
 *
 * with the saddle theta* solving eta'(theta*) = k/t and g the amplitude of the
 * Perron projection, g(theta) = (pi0 v)(u 1), u and v the left and right Perron
 * vectors normalised by u v = 1.
 *
 * THE EXPANSION PARAMETER IS K2 = t*eta''(theta*), THE VARIANCE OF THE COUNT,
 * not its mean and not t. Measured error laws, with the constants flat to two
 * digits over Erlang orders 1..8 and horizons 10..160:
 *
 *   err(daniels) = 0.083 / K2        err(daniels2) = 0.017 / K2^2
 *
 * For a renewal Erlang(r) the count variance rate is lambda/r, so
 * K2 = lambda*t/r and an Erlang-4 at t=50 is as accurate as a Poisson at
 * t=12.5: low variability shrinks the parameter, it does not break the method.
 * Below K2 = 5 the expansion is out of its regime and the result carries a
 * flag saying so.
 *
 * ATTRIBUTION. The first-order form is Daniels (1954). The amplitude g and the
 * whole 'daniels2' bracket are NOT a rederivation: they are Jensen, "Saddlepoint
 * Expansions for Sums of Markov Dependent Variables on a Continuous State Space",
 * Probab. Th. Rel. Fields 89, 1991, Eq. (4.4) with the coefficients on p.191. His
 * gamma_0(s) = (sum_i c_i)(sum_i r_i P(Y_0=i)) is exactly g under his own
 * normalisation sum_i r_i c_i = 1, and expanding his
 * alpha_0 + (1/n){-alpha_3/2 + alpha_4/8 - 5*alpha_5/24} reproduces
 * g*(1 + lam4/8 - 5*lam3^2/24) - g''/(2*K2) + g'*K3/(2*K2^2) term for term; his
 * Theorem 4.1 gives the O(n^-2) error measured here as 0.017/K2^2. Jensen works
 * with discrete-n sums over a Markov chain, so the continuous-time MAP counting
 * process is that result transcribed, n -> t and the kernel eigenvalue -> the
 * Perron root of D0+exp(theta)*D1.
 *
 * GATED ON TRANSCENDENTAL ARITHMETIC (exp, log, sqrt) and, for the Perron root,
 * on LAPACK through eig_values, as several other api/mam headers already are.
 *
 * This is an asymptotic method, not a quadrature: use it for rare-event and
 * large-deviation coefficients, where k/t is away from lambda or where the
 * probability underflows. For the bulk of the transient distribution, i.e.
 * every block k=0..N-1 at once at moderate t, uniformization
 * (ctmc_uniformization, ctmc_foxglynn) is both exact and faster.
 */

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <cstddef>
#include <numeric>
#include <string>
#include <vector>

#include "line/num/number.h"
#include "line/util/eig.h"
#include "line/util/error.h"
#include "line/util/expm.h"
#include "line/api/mc/ctmc_solve.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace mc {

/**
 * Below this value of K2 = t*eta''(theta*) the expansion is out of its regime.
 * Do NOT threshold on lambda*t: for Erlang(r) the count variance rate is
 * lambda/r, so K2 = lambda*t/r, and lambda*t over-warns on Poisson-like
 * processes while under-warning on low-variability ones.
 */
static const double CTMC_SADDLEPOINT_K2_MIN = 5.0;

/** Which term of the steepest-descent expansion to stop at. */
enum SaddlepointMethod {
    SADDLEPOINT_DANIELS2 = 0,  ///< second order, error O(1/K2^2) -- the DEFAULT
    SADDLEPOINT_DANIELS = 1,   ///< first order with the Perron amplitude, O(1/K2)
    SADDLEPOINT_PLAIN = 2      ///< bare first order, amplitude set to 1
};

namespace detail {

/**
 * exp/log/isfinite through ADL, the idiom the rest of api/mc uses: `using
 * std::exp` then an unqualified call, so a Real<D> or Rational picks up its own
 * overload. num_traits carries no transcendental entry points.
 */
template <class T>
inline T tx_exp(const T& x) {
    using std::exp;
    return exp(x);
}

template <class T>
inline T tx_log(const T& x) {
    using std::log;
    return log(x);
}

template <class T>
inline bool tx_finite(const T& x) {
    const double d = num_traits<T>::to_double(x);
    return d == d && d != std::numeric_limits<double>::infinity() &&
           d != -std::numeric_limits<double>::infinity();
}

}  // namespace detail

using detail::tx_exp;
using detail::tx_log;
using detail::tx_finite;

/** Perron root of A(theta) with its first two derivatives and the amplitude. */
template <class T>
struct PerronState {
    T eta;
    T deta;
    T d2eta;
    T ampl;
};

/** One entry per (t,k) pair. */
template <class T>
struct SaddlepointResult {
    std::vector<T> p;         ///< approximation of Pr{N(t)=k}
    std::vector<T> logp;      ///< its natural logarithm, accurate below the floor
    std::vector<T> theta;     ///< the saddle theta*, -inf where k=0
    std::vector<T> eta;       ///< eta(theta*)
    std::vector<T> deta;      ///< eta'(theta*), equal to k/t at convergence
    std::vector<T> d2eta;     ///< eta''(theta*)
    std::vector<T> ampl;      ///< the Perron amplitude g(theta*)
    std::vector<T> corr;      ///< the bracket multiplying the leading term
    std::vector<T> k2;        ///< K2 = t*eta''(theta*), the expansion parameter
    std::vector<int> iter;    ///< Newton steps taken
    std::vector<bool> exact;  ///< true where the value is exact, not approximated
    T lambda;                 ///< the stationary event rate eta'(0)
    /// true when some point fell below CTMC_SADDLEPOINT_K2_MIN
    bool out_of_regime;
    T worst_k2;               ///< the smallest K2 met
    double worst_t;           ///< horizon at which it was met
    long worst_k;             ///< count at which it was met
};

namespace detail {

/** Dense solve with partial pivoting, orders K and K+1 only. */
template <class T>
inline std::vector<T> saddle_solve(Matrix<T> A, std::vector<T> b) {
    std::vector<std::size_t> piv = lu_factor(A);
    lu_solve(A, piv, b);
    return b;
}

}  // namespace detail

/**
 * Perron root of A(th) = D0 + exp(th)*D1 with deta, d2eta and the amplitude.
 *
 * Only EIGENVALUES are taken from the eigensolver; the Perron vectors come from
 * bordered solves, the idiom ctmc_solve already uses. That keeps all four
 * codebases on ONE algorithm: eig.h exposes values only, and the JAR's
 * commons-math hands back Schur blocks rather than eigenvectors as soon as a
 * complex pair appears, so neither can supply a left eigenvector.
 */
template <class T>
inline PerronState<T> ctmc_saddlepoint_perron(const Matrix<T>& D0, const Matrix<T>& D1,
                                              const std::vector<T>& pi0, const T& th) {
    const std::size_t n = D0.rows();
    const T ex = tx_exp(th);
    Matrix<T> W(n, n), A(n, n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            W(i, j) = ex * D1(i, j);  // A'(th) = A''(th) = exp(th)*D1
            A(i, j) = D0(i, j) + W(i, j);
        }

    Matrix<double> Ad(n, n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) Ad(i, j) = num_traits<T>::to_double(A(i, j));
    const std::vector<std::complex<double> > ev = eig_values(Ad);
    double etad = -std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < ev.size(); ++i)
        if (ev[i].real() > etad) etad = ev[i].real();
    const T eta = num_traits<T>::from_double(etad);

    Matrix<T> Ashift(n, n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            Ashift(i, j) = A(i, j) - (i == j ? eta : num_traits<T>::from_int(0));

    std::vector<T> rhs(n, num_traits<T>::from_int(0));
    rhs[n - 1] = num_traits<T>::from_int(1);
    // (A-eta*I)v = 0 with the last row replaced by sum(v)=1. A row may be
    // dropped because A-eta*I is a singular irreducible M-matrix, every proper
    // principal submatrix of which is nonsingular.
    Matrix<T> M = Ashift;
    for (std::size_t j = 0; j < n; ++j) M(n - 1, j) = num_traits<T>::from_int(1);
    const std::vector<T> v = detail::saddle_solve(M, rhs);
    // u(A-eta*I) = 0 by the same construction on the transpose
    Matrix<T> Mt(n, n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            Mt(i, j) = (i == n - 1) ? num_traits<T>::from_int(1) : Ashift(j, i);
    std::vector<T> u = detail::saddle_solve(Mt, rhs);
    T uv = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < n; ++i) uv += u[i] * v[i];
    for (std::size_t i = 0; i < n; ++i) u[i] /= uv;  // u v = 1 fixes the scale

    PerronState<T> st;
    st.eta = eta;
    st.deta = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < n; ++i) {
        T inner = num_traits<T>::from_int(0);
        for (std::size_t j = 0; j < n; ++j) inner += W(i, j) * v[j];
        st.deta += u[i] * inner;
    }

    // First-order eigenvector perturbation (A-eta*I)v' = (eta'*I-W)v taken with
    // u v' = 0; the bordered system is nonsingular because the Perron root of an
    // irreducible Metzler matrix is simple.
    Matrix<T> B(n + 1, n + 1, num_traits<T>::from_int(0));
    std::vector<T> r(n + 1, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) B(i, j) = Ashift(i, j);
        B(i, n) = v[i];
        B(n, i) = u[i];
        T acc = st.deta * v[i];
        for (std::size_t j = 0; j < n; ++j) acc -= W(i, j) * v[j];
        r[i] = acc;
    }
    const std::vector<T> sol = detail::saddle_solve(B, r);
    T corr2 = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < n; ++i) {
        T inner = num_traits<T>::from_int(0);
        for (std::size_t j = 0; j < n; ++j) inner += W(i, j) * sol[j];
        corr2 += u[i] * inner;
    }
    st.d2eta = st.deta + num_traits<T>::from_int(2) * corr2;

    T pv = num_traits<T>::from_int(0), u1 = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < n; ++i) {
        pv += pi0[i] * v[i];
        u1 += u[i];
    }
    st.ampl = pv * u1;
    return st;
}

namespace detail {

/**
 * Saddle of the counting cumulant generating function at rate r, the root of
 * eta'(th) = r. eta' is continuous and strictly increasing from 0 to +inf, so
 * the root exists and is unique for every r>0; it is bracketed by geometric
 * expansion from th0 and refined by Newton on log(eta'), safeguarded by
 * bisection.
 */
template <class T>
inline T saddle_root(const Matrix<T>& D0, const Matrix<T>& D1, const std::vector<T>& pi0,
                     const T& r, const T& th0, const T& thmin, const T& thmax, int& iters) {
    const T TOL = num_traits<T>::from_double(1e-13);
    const int MAXIT = 200;
    T th = std::min(std::max(th0, thmin), thmax);
    T d1 = ctmc_saddlepoint_perron(D0, D1, pi0, th).deta;
    T lo = th, hi = th, dlo = d1, dhi = d1;
    T step = num_traits<T>::from_int(1);
    while (dlo > r) {
        hi = lo;
        dhi = dlo;
        lo = lo - step;
        if (lo <= thmin) {
            lo = thmin;
            dlo = ctmc_saddlepoint_perron(D0, D1, pi0, lo).deta;
            if (dlo > r)
                throw InputError("ctmc_saddlepoint: the rate k/t is below the representable "
                                 "range of eta'");
            break;
        }
        dlo = ctmc_saddlepoint_perron(D0, D1, pi0, lo).deta;
        step += step;
    }
    step = num_traits<T>::from_int(1);
    while (dhi < r) {
        lo = hi;
        dlo = dhi;
        hi = hi + step;
        if (hi >= thmax) {
            hi = thmax;
            dhi = ctmc_saddlepoint_perron(D0, D1, pi0, hi).deta;
            if (dhi < r)
                throw InputError("ctmc_saddlepoint: the rate k/t is above the representable "
                                 "range of eta'");
            break;
        }
        dhi = ctmc_saddlepoint_perron(D0, D1, pi0, hi).deta;
        step += step;
    }
    th = std::min(std::max(th, lo), hi);
    const T logr = tx_log(r);
    iters = 0;
    for (int it = 1; it <= MAXIT; ++it) {
        iters = it;
        const PerronState<T> si = ctmc_saddlepoint_perron(D0, D1, pi0, th);
        const T f = tx_log(si.deta) - logr;
        if (num_abs(f) <= TOL) break;
        if (f > num_traits<T>::from_int(0))
            hi = th;
        else
            lo = th;
        T thn = th - f * si.deta / si.d2eta;
        if (!tx_finite(thn) || thn <= lo || thn >= hi)
            thn = (lo + hi) / num_traits<T>::from_int(2);
        if (num_abs(thn - th) <=
            TOL * std::max(num_traits<T>::from_int(1), num_abs(th))) {
            th = thn;
            break;
        }
        th = thn;
    }
    return th;
}

}  // namespace detail

/**
 * Pr{N(t)=k} over arrays of horizons and counts.
 *
 * @param D0     generator of the phase process with the counted transitions removed
 * @param D1     rates of the counted transitions; D0+D1 must be an irreducible generator
 * @param t      time horizons; length 1 broadcasts against k
 * @param k      event counts; length 1 broadcasts against t
 * @param method SADDLEPOINT_DANIELS2 (the default), SADDLEPOINT_DANIELS or SADDLEPOINT_PLAIN
 * @param pi0    initial phase distribution; empty selects the stationary distribution of D0+D1
 */
template <class T>
inline SaddlepointResult<T> ctmc_saddlepoint(const Matrix<T>& D0, const Matrix<T>& D1,
                                             const std::vector<T>& t, const std::vector<long>& k,
                                             SaddlepointMethod method = SADDLEPOINT_DANIELS2,
                                             const std::vector<T>& pi0 = std::vector<T>()) {
    static_assert(num_traits<T>::has_transcendental,
                  "ctmc_saddlepoint is unavailable over the rational field: the saddlepoint "
                  "involves exp, log and sqrt of the Perron root, none a rational function "
                  "of the rates");
    const std::size_t nph = D0.rows();
    if (D0.cols() != nph || D1.rows() != nph || D1.cols() != nph)
        throw InputError("ctmc_saddlepoint: D0 and D1 must be square matrices of the same order");
    const T zero = num_traits<T>::from_int(0);
    T maxrate = zero, maxq = zero;
    for (std::size_t i = 0; i < nph; ++i)
        for (std::size_t j = 0; j < nph; ++j) {
            if (D1(i, j) < zero) throw InputError("ctmc_saddlepoint: D1 must be nonnegative");
            maxrate = std::max(maxrate, D1(i, j));
            maxq = std::max(maxq, num_abs(D0(i, j) + D1(i, j)));
        }
    if (!(maxrate > zero))
        throw InputError("ctmc_saddlepoint: D1 has no counted transitions, the counting process "
                         "is identically zero");
    Matrix<T> Q(nph, nph);
    for (std::size_t i = 0; i < nph; ++i) {
        T rowsum = zero;
        for (std::size_t j = 0; j < nph; ++j) {
            Q(i, j) = D0(i, j) + D1(i, j);
            rowsum += Q(i, j);
        }
        if (num_abs(rowsum) >
            num_traits<T>::from_double(1e-8) * std::max(num_traits<T>::from_int(1), maxq))
            throw InputError("ctmc_saddlepoint: D0+D1 must be an infinitesimal generator "
                             "(zero row sums)");
    }

    std::vector<T> pi = pi0;
    if (pi.empty()) pi = ctmc_solve(Q);
    if (pi.size() != nph)
        throw InputError("ctmc_saddlepoint: pi0 must have one entry per phase");
    T pisum = zero;
    for (std::size_t i = 0; i < nph; ++i) pisum += pi[i];
    if (num_abs(pisum - num_traits<T>::from_int(1)) > num_traits<T>::from_double(1e-8))
        throw InputError("ctmc_saddlepoint: pi0 must sum to one");

    // Broadcast the horizons against the counts
    const std::size_t n = std::max(t.size(), k.size());
    if ((t.size() != n && t.size() != 1) || (k.size() != n && k.size() != 1))
        throw InputError("ctmc_saddlepoint: t and k must be scalars or arrays of the same size");
    std::vector<T> tv(n);
    std::vector<long> kv(n);
    for (std::size_t i = 0; i < n; ++i) {
        tv[i] = t.size() == 1 ? t[0] : t[i];
        kv[i] = k.size() == 1 ? k[0] : k[i];
        if (tv[i] < zero) throw InputError("ctmc_saddlepoint: the horizon t must be nonnegative");
        if (kv[i] < 0) throw InputError("ctmc_saddlepoint: the count k must be nonnegative");
    }

    SaddlepointResult<T> res;
    const T ninf = num_traits<T>::from_double(-std::numeric_limits<double>::infinity());
    const T nan = num_traits<T>::from_double(std::numeric_limits<double>::quiet_NaN());
    res.p.assign(n, zero);
    res.logp.assign(n, ninf);
    res.theta.assign(n, ninf);
    res.eta.assign(n, nan);
    res.deta.assign(n, nan);
    res.d2eta.assign(n, nan);
    res.ampl.assign(n, nan);
    res.corr.assign(n, nan);
    res.k2.assign(n, nan);
    res.iter.assign(n, 0);
    res.exact.assign(n, false);
    res.lambda = ctmc_saddlepoint_perron(D0, D1, pi, zero).deta;
    res.out_of_regime = false;
    res.worst_k2 = num_traits<T>::from_double(std::numeric_limits<double>::infinity());
    res.worst_t = 0.0;
    res.worst_k = 0;

    // exp(theta) multiplies D1, so the saddle is confined to the range over
    // which A(theta) is representable; never active for a feasible k/t
    const T thmax = num_traits<T>::from_double(std::log(std::numeric_limits<double>::max() / 1e6)) -
                    tx_log(maxrate);
    const T thmin = num_traits<T>::from_double(std::log(std::numeric_limits<double>::min() * 1e6)) -
                    tx_log(maxrate);

    // Sorting by the rate k/t lets each Newton solve warm-start from the
    // previous saddle, the saddle being a monotone function of that rate alone
    std::vector<std::size_t> ord(n);
    for (std::size_t i = 0; i < n; ++i) ord[i] = i;
    std::vector<double> rate(n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
        if (tv[i] > zero) rate[i] = double(kv[i]) / num_traits<T>::to_double(tv[i]);
    std::stable_sort(ord.begin(), ord.end(),
                     [&rate](std::size_t a, std::size_t b) { return rate[a] < rate[b]; });

    T thprev = zero;
    for (std::size_t idx = 0; idx < n; ++idx) {
        const std::size_t i = ord[idx];
        const T ti = tv[i];
        const long ki = kv[i];
        if (ti == zero) {
            // No time has elapsed, so the count is zero with probability one
            res.exact[i] = true;
            if (ki == 0) {
                res.p[i] = num_traits<T>::from_int(1);
                res.logp[i] = zero;
            }
            continue;
        }
        if (ki == 0) {
            // The saddle runs off to -inf; the exact value is one matrix
            // exponential of the taboo generator and costs no more than a step
            // of the approximation itself
            res.exact[i] = true;
            const Matrix<T> E = expm(D0, ti);
            T acc = zero;
            for (std::size_t a = 0; a < nph; ++a)
                for (std::size_t b = 0; b < nph; ++b) acc += pi[a] * E(a, b);
            res.p[i] = acc;
            res.logp[i] = acc > zero ? tx_log(acc) : ninf;
            continue;
        }

        int iters = 0;
        const T th = detail::saddle_root(D0, D1, pi, num_traits<T>::from_double(rate[i]), thprev,
                                         thmin, thmax, iters);
        thprev = th;
        res.theta[i] = th;
        res.iter[i] = iters;

        const PerronState<T> s = ctmc_saddlepoint_perron(D0, D1, pi, th);
        res.eta[i] = s.eta;
        res.deta[i] = s.deta;
        res.d2eta[i] = s.d2eta;

        const T K2 = ti * s.d2eta;
        res.k2[i] = K2;
        if (K2 < res.worst_k2) {
            res.worst_k2 = K2;
            res.worst_t = num_traits<T>::to_double(ti);
            res.worst_k = ki;
        }
        if (!(K2 > zero))
            throw NumericError("ctmc_saddlepoint: the cumulant generating function is not "
                               "strictly convex at the saddle; D0+D1 is probably reducible");
        const T kiT = num_traits<T>::from_double(double(ki));
        const T base = ti * s.eta - kiT * th -
                       num_traits<T>::from_double(0.5) *
                           tx_log(num_traits<T>::from_double(2.0 * 3.14159265358979323846) * K2);
        const T ampl = (method == SADDLEPOINT_PLAIN) ? num_traits<T>::from_int(1) : s.ampl;
        res.ampl[i] = s.ampl;

        T corr;
        if (method != SADDLEPOINT_DANIELS2) {
            corr = ampl;
        } else {
            // The higher cumulants and the derivatives of the amplitude come
            // from central differences of the analytic eta'' and g, both of
            // which carry full precision at each evaluation point
            const T h = num_traits<T>::from_double(1e-3) *
                        std::max(num_traits<T>::from_int(1), num_abs(th));
            const PerronState<T> sp = ctmc_saddlepoint_perron(D0, D1, pi, th + h);
            const PerronState<T> sm = ctmc_saddlepoint_perron(D0, D1, pi, th - h);
            const T two = num_traits<T>::from_int(2);
            const T d3 = (sp.d2eta - sm.d2eta) / (two * h);
            const T d4 = (sp.d2eta - two * s.d2eta + sm.d2eta) / (h * h);
            const T K3 = ti * d3, K4 = ti * d4;
            const T lam3sq = K3 * K3 / (K2 * K2 * K2);
            const T lam4 = K4 / (K2 * K2);
            const T gp = (sp.ampl - sm.ampl) / (two * h);
            const T gpp = (sp.ampl - two * s.ampl + sm.ampl) / (h * h);
            // Steepest descent to O(1/K2), Jensen (1991) Eq. (4.4): the Daniels
            // bracket on the amplitude, plus the two terms the amplitude
            // contributes through its own curvature along the contour
            corr = ampl * (num_traits<T>::from_int(1) + lam4 / num_traits<T>::from_int(8) -
                           num_traits<T>::from_int(5) * lam3sq / num_traits<T>::from_int(24)) -
                   gpp / (two * K2) + gp * K3 / (two * K2 * K2);
            if (corr <= zero) corr = ampl;  // expansion broken down, fall back
        }
        res.corr[i] = corr;
        res.logp[i] = base + tx_log(corr);
        res.p[i] = tx_exp(res.logp[i]);
    }

    res.out_of_regime =
        num_traits<T>::to_double(res.worst_k2) < CTMC_SADDLEPOINT_K2_MIN;
    return res;
}

/** Pr{N(t)=k} at a single (t,k), with the default method. */
template <class T>
inline T ctmc_saddlepoint(const Matrix<T>& D0, const Matrix<T>& D1, const T& t, long k,
                          SaddlepointMethod method = SADDLEPOINT_DANIELS2) {
    return ctmc_saddlepoint(D0, D1, std::vector<T>(1, t), std::vector<long>(1, k), method).p[0];
}

}  // namespace mc
}  // namespace line

#endif  // LINE_API_MC_CTMC_SADDLEPOINT_H
