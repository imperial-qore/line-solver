/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_INFER_INFER_VARIATIONAL_H
#define LINE_API_INFER_INFER_VARIATIONAL_H

/**
 * Variational inference for Markovian queueing networks.
 *
 * Port of matlab/src/api/infer/infer_variational.m and of the JAR twin
 * jline.inference.api.Infer_variational, following I. Perez, G. Casale,
 * "Variational Inference for Markovian Queueing Networks", Advances in Applied
 * Probability 53(3), 2021.
 *
 * The network trajectory is reparameterised by the transition counts Y^eta,
 * eta = (i,j,c), so that the station marginals decouple:
 *
 *   x_{i,c}(t) = x_{i,c}(0) + sum_{eta in In(i,c)} Y^eta(t)
 *                           - sum_{eta in Out(i,c)} Y^eta(t)
 *
 * The variational family is a product of inhomogeneous pure-birth processes,
 * one per transition, with rate nu^eta(t,y), times a product of Gamma
 * densities over the unknown service rates. The state space is expanded by
 * adding DELTA to every feasible rate, so that queue lengths may go negative
 * and the approximating measure stays mutually absolutely continuous with the
 * target; the original model is recovered as DELTA -> 0.
 *
 * Each iteration performs, per transition, a backward pass for the Lagrange
 * multipliers r^eta with multiplicative jumps at the observation epochs, the
 * rate update nu^eta(t,y) = exp(E log Xi^eta(t,y)) r^eta(t,y+1)/r^eta(t,y),
 * and a forward pass of the master equation for the marginal. The conjugate
 * Gamma posteriors are then refreshed from the expected number of firings and
 * the expected exposure time of each station-class pair.
 *
 * ARITHMETIC: the backward and forward passes are uniformizations, i.e. convex
 * combinations of sub-stochastic matrix actions, so they stay positive and
 * bounded whatever the rate scale. The transcendental content is exp, log and
 * the digamma of the Gamma posteriors; there is no random-number stream, since
 * the expectations over the other transitions are taken on a Halton lattice
 * mapped through the inverse marginal c.d.f. The estimator therefore
 * reproduces the MATLAB, Java and Python implementations digit for digit.
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace infer {

/** Service discipline of a station, as seen by the load factor Upsilon. */
enum class VariationalSched { INF = 0, SHARED = 1, EXTERNAL = 2 };

/**
 * Inference problem handed to infer_variational.
 *
 * Station-class pairs are flattened column-major, so that pair (m,r) sits at
 * index r*M+m, matching the MATLAB, Java and Python specifications. Station
 * and class indices inside `arcs` are one-based; index 0 marks the external
 * source or sink.
 */
template <class T>
struct VariationalSpec {
    /** (T x 3) transitions [i j c]; i==0 external source, j==0 sink. */
    std::vector<std::array<std::size_t, 3>> arcs;
    /** (M x R) initial queue lengths. */
    Matrix<T> x0;
    /** (M) discipline codes. */
    std::vector<int> sched;
    /** (M) number of servers. */
    std::vector<T> nservers;
    /** (T) routing probability of each transition. */
    std::vector<T> routeprob;
    /** (T) index in 1..P of the rate governing the transition, 0 when known. */
    std::vector<std::size_t> arcparam;
    /** (T) known rate for transitions with arcparam == 0. */
    std::vector<T> arcrate;
    /** (P) Gamma prior shapes and rates. */
    std::vector<T> alpha0, beta0;
    /** (K) observation epochs. */
    std::vector<T> obsTimes;
    /** (K x M*R) observed queue lengths; an unobserved entry is `unobserved`. */
    Matrix<T> obsData;
    /** (M*R) support size of the uniform contamination. */
    std::vector<T> obsRange;
    /** probability that a reading is faulty. */
    T epsilon;
    /**
     * (M*R) upper bound on the queue length, empty for none. In a closed
     * network this is the chain population, and clamping the load there keeps
     * the expanded state space from crediting a station with more jobs than
     * the network holds.
     */
    std::vector<T> capacity;

    /** Sentinel marking an unobserved entry of obsData. */
    static T unobserved() { return num_traits<T>::from_double(-1.0); }
};

/** Options of infer_variational; a negative box means "derive a default". */
template <class T>
struct VariationalOptions {
    int verbose = 0;
    std::size_t iter_max = 20;
    double tol = 1e-3;
    std::size_t nsamples = 200;
    /** rate added to every feasible transition by the space expansion. */
    double delta = 1e-3;
    double floor = 1e-4;
    /** cap on the variational rates; derived from ymax and tmax when negative. */
    double rate_max = -1.0;
    double rate_cap_factor = 10.0;
    double unifmax = 30.0;
    double unif_tol = 1e-12;
    std::size_t unif_max_terms = 2000;
    double tmax = -1.0;
    double dt = -1.0;
    long ngrid = -1;
    long ymax = -1;
};

/** Outcome of infer_variational. */
template <class T>
struct VariationalResult {
    std::vector<T> alpha, beta, rates, mean_service_time;
    std::vector<T> bound;
    Matrix<T> alpha_trace, beta_trace;  // (P x iter)
    std::vector<Matrix<T>> Y, nu;       // one (G x ymax+1) block per transition
    std::vector<T> tgrid;
    Matrix<T> qlen;  // (G x M*R)
    std::size_t iter = 0;
    bool converged = false;
    T tailmass;
};

namespace detail {

/** Recurrence threshold: the Stirling tails are below 1e-17 from here on. */
inline int iv_gamma_shift() { return 20; }

/** psi(x) for x > 0, upward recurrence to x >= 20 then the Stirling series. */
template <class T>
T iv_digamma(const T& x0) {
    const T one = num_traits<T>::from_int(1);
    const T lim = num_traits<T>::from_int(iv_gamma_shift());
    T x = x0, acc = num_traits<T>::from_int(0);
    while (x < lim) {
        acc -= one / x;
        x += one;
    }
    const T inv = one / x, inv2 = inv * inv;
    using std::log;
    T s = log(x) - inv / num_traits<T>::from_int(2);
    T p = inv2;
    s -= p / num_traits<T>::from_int(12);
    p *= inv2;
    s += p / num_traits<T>::from_int(120);
    p *= inv2;
    s -= p / num_traits<T>::from_int(252);
    p *= inv2;
    s += p / num_traits<T>::from_int(240);
    p *= inv2;
    s -= p / num_traits<T>::from_int(132);
    return acc + s;
}

/** log Gamma(x) for x > 0, the same construction on log Gamma(x)=log Gamma(x+1)-log x. */
template <class T>
T iv_lgamma(const T& x0) {
    using std::log;
    const T one = num_traits<T>::from_int(1);
    const T lim = num_traits<T>::from_int(iv_gamma_shift());
    T x = x0, acc = num_traits<T>::from_int(0);
    while (x < lim) {
        acc -= log(x);
        x += one;
    }
    const T half = one / num_traits<T>::from_int(2);
    const T log2pi = num_traits<T>::from_double(std::log(2.0 * 3.14159265358979323846));
    const T inv = one / x, inv2 = inv * inv;
    T s = (x - half) * log(x) - x + log2pi / num_traits<T>::from_int(2);
    T p = inv;
    s += p / num_traits<T>::from_int(12);
    p *= inv2;
    s -= p / num_traits<T>::from_int(360);
    p *= inv2;
    s += p / num_traits<T>::from_int(1260);
    p *= inv2;
    s -= p / num_traits<T>::from_int(1680);
    p *= inv2;
    s += p / num_traits<T>::from_int(1188);
    return acc + s;
}

/** Load factor Upsilon of a transition leaving a station-class pair. */
template <class T>
T iv_ups(const T& xic, const T& xis, const T& nservers, int sched, const T& cap,
         const T& capstat) {
    const T zero = num_traits<T>::from_int(0);
    if (sched == static_cast<int>(VariationalSched::EXTERNAL)) return num_traits<T>::from_int(1);
    T a = xic < zero ? zero : xic;
    if (a > cap) a = cap;
    if (sched == static_cast<int>(VariationalSched::INF)) return a;
    T b = xis < zero ? zero : xis;
    if (b > capstat) b = capstat;
    if (!(b > zero)) return zero;
    const T srv = nservers < b ? nservers : b;
    return a / b * srv;
}

/** k-th prime, k >= 1. */
inline std::size_t iv_prime(std::size_t k) {
    std::size_t n = 0, c = 1, p = 2;
    while (n < k) {
        ++c;
        bool isp = true;
        for (std::size_t d = 2; d * d <= c; ++d) {
            if (c % d == 0) {
                isp = false;
                break;
            }
        }
        if (isp) {
            ++n;
            p = c;
        }
    }
    return p;
}

/** Van der Corput radical inverse of i in the given base. */
inline double iv_radical_inverse(std::size_t i, std::size_t base) {
    double r = 0.0, f = 1.0 / static_cast<double>(base);
    while (i > 0) {
        r += f * static_cast<double>(i % base);
        i /= base;
        f /= static_cast<double>(base);
    }
    return r;
}

/**
 * Inverse-c.d.f. samples of a marginal on a Halton lattice. Each transition
 * uses its own prime base, so the samples of distinct transitions are jointly
 * equidistributed rather than comonotone.
 */
template <class T>
Matrix<T> iv_sample(const Matrix<T>& q, std::size_t S, std::size_t e) {
    const std::size_t G = q.rows(), ny = q.cols();
    const std::size_t base = iv_prime(e + 1);
    std::vector<double> u(S);
    std::vector<std::size_t> ord(S);
    for (std::size_t s = 0; s < S; ++s) {
        u[s] = iv_radical_inverse(s + 1, base);
        ord[s] = s;
    }
    std::stable_sort(ord.begin(), ord.end(),
                     [&u](std::size_t a, std::size_t b) { return u[a] < u[b]; });
    std::vector<double> us(S);
    for (std::size_t s = 0; s < S; ++s) us[s] = u[ord[s]];
    Matrix<T> ys(G, S, num_traits<T>::from_int(0));
    std::vector<double> c(ny);
    for (std::size_t g = 0; g < G; ++g) {
        double acc = 0.0;
        for (std::size_t y = 0; y < ny; ++y) {
            acc += num_traits<T>::to_double(q(g, y));
            c[y] = acc;
        }
        if (c[ny - 1] > 0) {
            for (std::size_t y = 0; y < ny; ++y) c[y] /= c[ny - 1];
        }
        c[ny - 1] = 1.0;
        std::size_t j = 0;
        for (std::size_t s = 0; s < S; ++s) {
            while (j + 1 < ny && c[j] < us[s]) ++j;
            ys(g, ord[s]) = num_traits<T>::from_int(static_cast<int>(j));
        }
    }
    return ys;
}

/** Slack multiplier of the rate cap; unity when the cap is inactive. */
template <class T>
T iv_damp(const T& sl, const T& ye, double floor) {
    const T zero = num_traits<T>::from_int(0);
    if (!(sl > zero)) return num_traits<T>::from_int(1);
    using std::exp;
    const T fl = num_traits<T>::from_double(floor);
    const T z = sl / (ye > fl ? ye : fl);
    const T d = (num_traits<T>::from_int(1) + z) / exp(z);
    if (!(d > zero)) return zero;
    return d;
}

/** Normalise by the largest entry; only the ratios of the multipliers matter. */
template <class T>
void iv_rescale(std::vector<T>& v) {
    const T zero = num_traits<T>::from_int(0);
    T m = zero;
    for (std::size_t i = 0; i < v.size(); ++i)
        if (v[i] > m) m = v[i];
    if (m > zero) {
        for (std::size_t i = 0; i < v.size(); ++i) v[i] = v[i] / m;
    }
}

/** One uniformization step of the backward sub-generator. */
template <class T>
std::vector<T> iv_back_uniformize(const std::vector<T>& v0, const std::vector<T>& pd,
                                  const std::vector<T>& pu, const T& lt,
                                  const VariationalOptions<T>& opt) {
    using std::exp;
    const std::size_t ny = v0.size();
    const T one = num_traits<T>::from_int(1);
    T w = exp(-lt);
    std::vector<T> v(ny), u(v0), un(ny);
    for (std::size_t i = 0; i < ny; ++i) v[i] = w * v0[i];
    T cum = w;
    std::size_t n = 1;
    const T tol = num_traits<T>::from_double(opt.unif_tol);
    while ((one - cum) > tol && n < opt.unif_max_terms) {
        for (std::size_t i = 0; i < ny; ++i) un[i] = u[i] * (one - pd[i]);
        for (std::size_t i = 0; i + 1 < ny; ++i) un[i] += u[i + 1] * pu[i];
        u = un;
        w = w * lt / num_traits<T>::from_int(static_cast<int>(n));
        for (std::size_t i = 0; i < ny; ++i) v[i] += w * u[i];
        cum += w;
        ++n;
    }
    return v;
}

/** One uniformization step of the pure-birth chain. */
template <class T>
std::vector<T> iv_uniformize(const std::vector<T>& v0, const std::vector<T>& p, const T& lt,
                             const VariationalOptions<T>& opt) {
    using std::exp;
    const std::size_t ny = v0.size();
    const T one = num_traits<T>::from_int(1);
    T w = exp(-lt);
    std::vector<T> v(ny), u(v0), un(ny);
    for (std::size_t i = 0; i < ny; ++i) v[i] = w * v0[i];
    T cum = w;
    std::size_t n = 1;
    const T tol = num_traits<T>::from_double(opt.unif_tol);
    while ((one - cum) > tol && n < opt.unif_max_terms) {
        for (std::size_t i = 0; i < ny; ++i) un[i] = u[i] * (one - p[i]);
        for (std::size_t i = ny; i-- > 1;) un[i] += u[i - 1] * p[i - 1];
        u = un;
        w = w * lt / num_traits<T>::from_int(static_cast<int>(n));
        for (std::size_t i = 0; i < ny; ++i) v[i] += w * u[i];
        cum += w;
        ++n;
    }
    return v;
}

/**
 * Backward pass for the Lagrange multipliers. The equation is linear in r, so
 * on a grid cell with frozen coefficients it is the action of a matrix
 * exponential. The generator has non-positive row sums by Jensen, so
 * uniformization evaluates it without the stiffness an explicit rule suffers
 * when exp(E log Xi) falls orders of magnitude below E[Xi].
 */
template <class T>
Matrix<T> iv_backward(const Matrix<T>& ge, const Matrix<T>& he, const Matrix<T>& sl,
                      const Matrix<T>& Ye, const std::vector<std::size_t>& obsIdx,
                      const Matrix<T>& obsw, const T& dt, const VariationalOptions<T>& opt) {
    const std::size_t G = ge.rows(), ny = ge.cols();
    const T zero = num_traits<T>::from_int(0);
    Matrix<T> r(G, ny, zero);
    std::vector<T> v(ny, num_traits<T>::from_int(1));
    for (std::size_t q = 0; q < obsIdx.size(); ++q) {
        if (obsIdx[q] + 1 == G) {
            for (std::size_t y = 0; y < ny; ++y)
                v[y] = v[y] * (obsw(q, y) > zero ? obsw(q, y) : zero);
        }
    }
    iv_rescale(v);
    for (std::size_t y = 0; y < ny; ++y) r(G - 1, y) = v[y];
    std::vector<T> pd(ny), pu(ny);
    for (std::size_t gi = G - 1; gi-- > 0;) {
        T lam = zero;
        for (std::size_t y = 0; y < ny; ++y) {
            const T gv = ge(gi, y) > zero ? ge(gi, y) : zero;
            T hv = he(gi, y) * iv_damp(sl(gi, y), Ye(gi, y), opt.floor);
            if (!(hv > zero)) hv = zero;
            if (hv > gv) hv = gv;
            pd[y] = gv;
            pu[y] = hv;
            if (gv > lam) lam = gv;
        }
        if (lam > zero) {
            const double lamd = num_traits<T>::to_double(lam);
            const double dtd = num_traits<T>::to_double(dt);
            const std::size_t ncell =
                std::max<std::size_t>(1, static_cast<std::size_t>(std::ceil(lamd * dtd / opt.unifmax)));
            const T h = dt / num_traits<T>::from_int(static_cast<int>(ncell));
            std::vector<T> pdn(ny), pun(ny);
            for (std::size_t y = 0; y < ny; ++y) {
                pdn[y] = pd[y] / lam;
                pun[y] = pu[y] / lam;
            }
            const T lt = lam * h;
            for (std::size_t c = 0; c < ncell; ++c) v = iv_back_uniformize(v, pdn, pun, lt, opt);
            iv_rescale(v);
        }
        for (std::size_t q = 0; q < obsIdx.size(); ++q) {
            if (obsIdx[q] == gi) {
                for (std::size_t y = 0; y < ny; ++y)
                    v[y] = v[y] * (obsw(q, y) > zero ? obsw(q, y) : zero);
                iv_rescale(v);
            }
        }
        for (std::size_t y = 0; y < ny; ++y) r(gi, y) = v[y];
    }
    return r;
}

/** Forward master equation of an inhomogeneous pure-birth process. */
template <class T>
Matrix<T> iv_forward(const Matrix<T>& nue, const T& dt, const VariationalOptions<T>& opt) {
    const std::size_t G = nue.rows(), ny = nue.cols();
    const T zero = num_traits<T>::from_int(0);
    Matrix<T> q(G, ny, zero);
    std::vector<T> v(ny, zero);
    v[0] = num_traits<T>::from_int(1);
    for (std::size_t y = 0; y < ny; ++y) q(0, y) = v[y];
    std::vector<T> p(ny);
    for (std::size_t g = 0; g + 1 < G; ++g) {
        T lam = zero;
        for (std::size_t y = 0; y < ny; ++y) {
            p[y] = nue(g, y) > zero ? nue(g, y) : zero;
            if (p[y] > lam) lam = p[y];
        }
        if (!(lam > zero)) {
            for (std::size_t y = 0; y < ny; ++y) q(g + 1, y) = v[y];
            continue;
        }
        const double lamd = num_traits<T>::to_double(lam);
        const double dtd = num_traits<T>::to_double(dt);
        const std::size_t ncell =
            std::max<std::size_t>(1, static_cast<std::size_t>(std::ceil(lamd * dtd / opt.unifmax)));
        const T h = dt / num_traits<T>::from_int(static_cast<int>(ncell));
        std::vector<T> pn(ny);
        for (std::size_t y = 0; y < ny; ++y) pn[y] = p[y] / lam;
        const T lt = lam * h;
        for (std::size_t c = 0; c < ncell; ++c) v = iv_uniformize(v, pn, lt, opt);
        T s = zero;
        for (std::size_t y = 0; y < ny; ++y) {
            if (v[y] < zero) v[y] = zero;
            s += v[y];
        }
        if (s > zero) {
            for (std::size_t y = 0; y < ny; ++y) v[y] = v[y] / s;
        }
        for (std::size_t y = 0; y < ny; ++y) q(g + 1, y) = v[y];
    }
    return q;
}

/** Trapezoidal integral of a grid function. */
template <class T>
T iv_trapz(const std::vector<T>& f, const T& dt) {
    if (f.size() < 2) return num_traits<T>::from_int(0);
    T s = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < f.size(); ++i) s += f[i];
    const T half = num_traits<T>::from_int(1) / num_traits<T>::from_int(2);
    return dt * (s - half * f.front() - half * f.back());
}

/** KL(Gamma(a,b) || Gamma(a0,b0)) with rate parameterisation. */
template <class T>
T iv_kl_gamma(const T& a, const T& b, const T& a0, const T& b0) {
    using std::log;
    return (a - a0) * iv_digamma(a) - iv_lgamma(a) + iv_lgamma(a0) + a0 * (log(b) - log(b0)) +
           a * (b0 - b) / b;
}

}  // namespace detail

/**
 * Run the variational inference procedure.
 *
 * @param spec    the inference problem
 * @param options solver options; defaults are derived from the specification
 * @return        posterior Gamma parameters, the bound trace and the marginals
 */
template <class T>
VariationalResult<T> infer_variational(VariationalSpec<T> spec,
                                       VariationalOptions<T> opt = VariationalOptions<T>()) {
    using std::exp;
    using std::log;
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const std::size_t M = spec.x0.rows(), R = spec.x0.cols();
    const std::size_t narcs = spec.arcs.size(), P = spec.alpha0.size();
    const std::size_t MR = M * R;
    const std::size_t K = spec.obsTimes.size();

    if (spec.sched.size() != M || spec.nservers.size() != M)
        throw InputError("infer_variational: sched and nservers must have one entry per station");
    if (spec.routeprob.size() != narcs || spec.arcparam.size() != narcs ||
        spec.arcrate.size() != narcs)
        throw InputError(
            "infer_variational: routeprob, arcparam and arcrate must have one entry per "
            "transition");
    if (spec.alpha0.size() != spec.beta0.size())
        throw InputError("infer_variational: alpha0 and beta0 must have the same length");
    if (spec.obsData.rows() != K || (K > 0 && spec.obsData.cols() != MR))
        throw InputError("infer_variational: obsData must be K x M*R");
    if (spec.obsRange.size() != MR)
        throw InputError("infer_variational: obsRange must have M*R entries");
    for (std::size_t e = 0; e < narcs; ++e) {
        if (spec.arcs[e][0] == 0 && spec.arcs[e][1] == 0)
            throw InputError("infer_variational: a transition cannot be external at both ends");
        if (spec.arcparam[e] == 0 && !(spec.arcrate[e] > zero))
            throw InputError(
                "infer_variational: a transition without a parameter needs a positive known rate");
    }
    if (spec.capacity.empty())
        spec.capacity.assign(MR, num_traits<T>::from_double(std::numeric_limits<double>::infinity()));
    if (spec.capacity.size() != MR)
        throw InputError("infer_variational: capacity must have M*R entries");

    const T unobs = VariationalSpec<T>::unobserved();
    std::vector<std::size_t> src(narcs), dst(narcs), cls(narcs);
    for (std::size_t e = 0; e < narcs; ++e) {
        src[e] = spec.arcs[e][0];
        dst[e] = spec.arcs[e][1];
        cls[e] = spec.arcs[e][2];
    }

    Matrix<T> sgnClass(narcs, MR, zero), sgnStat(narcs, M, zero);
    for (std::size_t e = 0; e < narcs; ++e) {
        if (dst[e] > 0) {
            sgnClass(e, (cls[e] - 1) * M + dst[e] - 1) += one;
            sgnStat(e, dst[e] - 1) += one;
        }
        if (src[e] > 0) {
            sgnClass(e, (cls[e] - 1) * M + src[e] - 1) -= one;
            sgnStat(e, src[e] - 1) -= one;
        }
    }

    std::vector<T> x0v(MR, zero), x0s(M, zero), capStat(M, zero);
    for (std::size_t m = 0; m < M; ++m) {
        for (std::size_t r = 0; r < R; ++r) {
            x0v[r * M + m] = spec.x0(m, r);
            x0s[m] += spec.x0(m, r);
            capStat[m] += spec.capacity[r * M + m];
        }
    }

    // mean occupancy used to size the truncation and the initial rates
    std::vector<T> xbar(x0v), xbars(M, zero);
    for (std::size_t k = 0; k < MR; ++k) {
        T s = zero;
        std::size_t n = 0;
        for (std::size_t q = 0; q < K; ++q) {
            if (!(spec.obsData(q, k) == unobs)) {
                s += spec.obsData(q, k);
                ++n;
            }
        }
        if (n > 0) xbar[k] = s / num_traits<T>::from_int(static_cast<int>(n));
    }
    for (std::size_t m = 0; m < M; ++m)
        for (std::size_t r = 0; r < R; ++r) xbars[m] += xbar[r * M + m];

    if (opt.tmax < 0) {
        if (K == 0) throw InputError("infer_variational: tmax is required without observations");
        double t = 0.0;
        for (std::size_t k = 0; k < K; ++k)
            t = std::max(t, num_traits<T>::to_double(spec.obsTimes[k]));
        opt.tmax = t;
    }
    if (!(opt.tmax > 0)) throw InputError("infer_variational: tmax must be positive");
    if (opt.ngrid < 0 && opt.dt < 0) opt.ngrid = 201;
    if (opt.ngrid < 0) opt.ngrid = static_cast<long>(std::llround(opt.tmax / opt.dt)) + 1;
    opt.ngrid = std::max<long>(2, opt.ngrid);
    opt.dt = opt.tmax / static_cast<double>(opt.ngrid - 1);

    if (opt.ymax < 0) {
        double fmax = 0.0;
        for (std::size_t e = 0; e < narcs; ++e) {
            double lam;
            if (spec.arcparam[e] > 0) {
                const std::size_t p = spec.arcparam[e] - 1;
                lam = num_traits<T>::to_double(spec.routeprob[e] * spec.alpha0[p] / spec.beta0[p]);
            } else {
                lam = num_traits<T>::to_double(spec.routeprob[e] * spec.arcrate[e]);
            }
            double u = 1.0;
            if (src[e] > 0) {
                const std::size_t kc = (cls[e] - 1) * M + src[e] - 1;
                u = num_traits<T>::to_double(detail::iv_ups(xbar[kc], xbars[src[e] - 1],
                                                            spec.nservers[src[e] - 1],
                                                            spec.sched[src[e] - 1],
                                                            spec.capacity[kc], capStat[src[e] - 1]));
            }
            fmax = std::max(fmax, lam * u * opt.tmax);
        }
        opt.ymax = std::max<long>(
            20, static_cast<long>(std::ceil(2 * fmax + 5 * std::sqrt(std::max(1.0, fmax)))));
    }
    opt.ymax = std::max<long>(2, opt.ymax);
    if (opt.rate_max < 0)
        opt.rate_max = opt.rate_cap_factor * static_cast<double>(opt.ymax) / opt.tmax;

    const std::size_t G = static_cast<std::size_t>(opt.ngrid);
    const std::size_t ymax = static_cast<std::size_t>(opt.ymax);
    const std::size_t ny = ymax + 1;
    const std::size_t S = opt.nsamples;
    const T dt = num_traits<T>::from_double(opt.dt);
    const T deltaT = num_traits<T>::from_double(opt.delta);
    const T floorT = num_traits<T>::from_double(opt.floor);
    const T rateMaxT = num_traits<T>::from_double(opt.rate_max);

    std::vector<T> yvec(ny), tgrid(G);
    for (std::size_t y = 0; y < ny; ++y) yvec[y] = num_traits<T>::from_int(static_cast<int>(y));
    for (std::size_t g = 0; g < G; ++g) tgrid[g] = num_traits<T>::from_int(static_cast<int>(g)) * dt;

    std::vector<std::size_t> obsIdx(K, 0);
    for (std::size_t k = 0; k < K; ++k) {
        const long idx = std::lround(num_traits<T>::to_double(spec.obsTimes[k]) / opt.dt);
        obsIdx[k] = static_cast<std::size_t>(std::min<long>(std::max<long>(idx, 0),
                                                            static_cast<long>(G) - 1));
    }

    std::vector<int> arcSched(narcs, static_cast<int>(VariationalSched::EXTERNAL));
    std::vector<T> arcServers(narcs, one);
    for (std::size_t e = 0; e < narcs; ++e) {
        if (src[e] > 0) {
            arcSched[e] = spec.sched[src[e] - 1];
            arcServers[e] = spec.nservers[src[e] - 1];
        }
    }

    std::vector<T> alpha(spec.alpha0), beta(spec.beta0);

    std::vector<Matrix<T>> Y, nu, slack, gexp, hexp;
    Y.reserve(narcs);
    nu.reserve(narcs);
    slack.reserve(narcs);
    gexp.reserve(narcs);
    hexp.reserve(narcs);

    // E[lambda_eta] and E[log lambda_eta] under the current Gamma posterior
    auto rate_mean = [&](std::size_t e) {
        const std::size_t p = spec.arcparam[e];
        if (p == 0) return spec.routeprob[e] * spec.arcrate[e];
        return spec.routeprob[e] * alpha[p - 1] / beta[p - 1];
    };
    auto rate_log_mean = [&](std::size_t e) {
        const std::size_t p = spec.arcparam[e];
        if (p == 0) return log(spec.routeprob[e] * spec.arcrate[e]);
        return log(spec.routeprob[e]) + detail::iv_digamma(alpha[p - 1]) - log(beta[p - 1]);
    };

    for (std::size_t e = 0; e < narcs; ++e) {
        const T lam = rate_mean(e);
        T u0 = one;
        if (src[e] > 0) {
            const std::size_t kc = (cls[e] - 1) * M + src[e] - 1;
            u0 = detail::iv_ups(xbar[kc], xbars[src[e] - 1], arcServers[e], arcSched[e],
                                spec.capacity[kc], capStat[src[e] - 1]);
        }
        T nu0 = lam * u0;
        if (nu0 < deltaT) nu0 = deltaT;
        Matrix<T> nue(G, ny, zero);
        for (std::size_t g = 0; g < G; ++g)
            for (std::size_t y = 0; y < ymax; ++y) nue(g, y) = nu0;
        nu.push_back(nue);
        Y.push_back(detail::iv_forward(nue, dt, opt));
        slack.push_back(Matrix<T>(G, ny, zero));
        gexp.push_back(Matrix<T>(G, ny, zero));
        hexp.push_back(Matrix<T>(G, ny, zero));
    }

    /**
     * Conditional rate moments of one transition, and its observation jumps.
     * Fills E[Xi | Y^eta=y] and exp(E[log Xi | Y^eta=y]) on the time grid,
     * both under Q with the transition's own contribution removed.
     */
    auto rate_moments = [&](std::size_t e, const std::vector<Matrix<T>>& Ys, Matrix<T>& ge,
                            Matrix<T>& he, Matrix<T>& obsw, bool want_obs) {
        const T lam = rate_mean(e);
        const T loglam = rate_log_mean(e);
        const bool has_origin = src[e] > 0;
        const std::size_t kclass = has_origin ? (cls[e] - 1) * M + src[e] - 1 : 0;
        const T sgnOwnClass = has_origin ? sgnClass(e, kclass) : zero;
        const T sgnOwnStat = has_origin ? sgnStat(e, src[e] - 1) : zero;
        const T d0 = deltaT / lam;
        for (std::size_t q = 0; q < K; ++q)
            for (std::size_t y = 0; y < ny; ++y) obsw(q, y) = one;
        std::vector<T> accg(ny), acch(ny);
        Matrix<T> aStore(MR, S, zero);
        for (std::size_t g = 0; g < G; ++g) {
            bool hasObs = false;
            if (want_obs) {
                for (std::size_t q = 0; q < K; ++q)
                    if (obsIdx[q] == g) hasObs = true;
            }
            std::fill(accg.begin(), accg.end(), zero);
            std::fill(acch.begin(), acch.end(), zero);
            for (std::size_t s = 0; s < S; ++s) {
                T a = zero, b = zero;
                if (has_origin) {
                    a = x0v[kclass];
                    b = x0s[src[e] - 1];
                    for (std::size_t f = 0; f < narcs; ++f) {
                        if (f == e) continue;
                        a += sgnClass(f, kclass) * Ys[f](g, s);
                        b += sgnStat(f, src[e] - 1) * Ys[f](g, s);
                    }
                }
                if (hasObs) {
                    for (std::size_t k = 0; k < MR; ++k) {
                        T acc = x0v[k];
                        for (std::size_t f = 0; f < narcs; ++f)
                            if (f != e) acc += sgnClass(f, k) * Ys[f](g, s);
                        aStore(k, s) = acc;
                    }
                }
                for (std::size_t y = 0; y < ny; ++y) {
                    T u = one;
                    if (has_origin)
                        u = detail::iv_ups(a + sgnOwnClass * yvec[y], b + sgnOwnStat * yvec[y],
                                           arcServers[e], arcSched[e], spec.capacity[kclass],
                                           capStat[src[e] - 1]);
                    accg[y] += u;
                    acch[y] += log(u + d0);
                }
            }
            const T Sn = num_traits<T>::from_int(static_cast<int>(S));
            // E[Xi] and exp(E[log Xi]) of the SAME rate Xi = delta + lam*Ups;
            // writing the second as exp(E[log lam]) exp(E[log(Ups + delta/E[lam])])
            // keeps the two consistent wherever Ups is deterministic, which is
            // what stops the backward equation from developing a gradient away
            // from the empty-station boundary
            for (std::size_t y = 0; y < ny; ++y) {
                ge(g, y) = deltaT + lam * (accg[y] / Sn);
                he(g, y) = exp(loglam + acch[y] / Sn);
            }
            if (hasObs) {
                for (std::size_t q = 0; q < K; ++q) {
                    if (obsIdx[q] != g) continue;
                    for (std::size_t y = 0; y < ny; ++y) {
                        T acc = zero;
                        for (std::size_t k = 0; k < MR; ++k) {
                            if (spec.obsData(q, k) == unobs) continue;
                            const T range =
                                spec.obsRange[k] > one ? spec.obsRange[k] : one;
                            T ak = zero;
                            for (std::size_t s = 0; s < S; ++s) {
                                const T x = aStore(k, s) + sgnClass(e, k) * yvec[y];
                                T p = zero;
                                if (x == spec.obsData(q, k)) {
                                    p = one - spec.epsilon;
                                } else if (!(x < zero) && !(x > spec.obsRange[k])) {
                                    p = spec.epsilon / range;
                                }
                                ak += log(floorT + p);
                            }
                            acc += ak / Sn;
                        }
                        obsw(q, y) = exp(acc);
                    }
                }
            }
        }
    };

    std::vector<T> bound;
    Matrix<T> alphaTrace(P, opt.iter_max, zero), betaTrace(P, opt.iter_max, zero);
    bool converged = false;
    std::size_t iter = 0;
    std::vector<Matrix<T>> Ys(narcs, Matrix<T>(G, S, zero));
    Matrix<T> ge(G, ny, zero), he(G, ny, zero), obsw(std::max<std::size_t>(K, 1), ny, zero);

    for (std::size_t it = 1; it <= opt.iter_max; ++it) {
        iter = it;
        for (std::size_t e = 0; e < narcs; ++e) {
            for (std::size_t f = 0; f < narcs; ++f) Ys[f] = detail::iv_sample(Y[f], S, f);
            rate_moments(e, Ys, ge, he, obsw, true);
            const Matrix<T> r = detail::iv_backward(ge, he, slack[e], Y[e], obsIdx, obsw, dt, opt);

            // Eq. (15). A vanishing multiplier marks a count the future
            // observations rule out; the rate there is zero, which is what
            // keeps the forward pass from placing mass on it.
            Matrix<T> nue(G, ny, zero), sl(G, ny, zero);
            for (std::size_t g = 0; g < G; ++g) {
                for (std::size_t y = 0; y < ymax; ++y) {
                    T val = zero;
                    if (r(g, y) > zero) val = he(g, y) * r(g, y + 1) / r(g, y);
                    if (!(val > zero)) val = zero;
                    if (val > rateMaxT) {
                        const T w = Y[e](g, y) > floorT ? Y[e](g, y) : floorT;
                        sl(g, y) = w * log(val / rateMaxT);
                        val = rateMaxT;
                    }
                    nue(g, y) = val;
                }
            }
            nu[e] = nue;
            slack[e] = sl;
            Y[e] = detail::iv_forward(nue, dt, opt);
        }

        // conjugate Gamma updates: the shape gains the expected number of
        // firings, the rate the expected exposure time of the station-class
        // pair that the parameter governs
        for (std::size_t f = 0; f < narcs; ++f) Ys[f] = detail::iv_sample(Y[f], S, f);
        std::vector<T> firings(P, zero), exposure(P, zero);
        std::vector<char> seen(P * MR, 0);
        for (std::size_t e = 0; e < narcs; ++e) {
            const std::size_t p = spec.arcparam[e];
            if (p == 0) continue;
            // expected number of firings over the horizon, taken from the
            // marginal itself, which is exact, rather than by quadrature of
            // the intensity, which a near-deterministic marginal makes
            // inaccurate
            T m1 = zero, m0 = zero;
            for (std::size_t y = 0; y < ny; ++y) {
                m1 += Y[e](G - 1, y) * yvec[y];
                m0 += Y[e](0, y) * yvec[y];
            }
            firings[p - 1] += m1 - m0;
            const std::size_t kclass = (cls[e] - 1) * M + src[e] - 1;
            if (!seen[(p - 1) * MR + kclass]) {
                seen[(p - 1) * MR + kclass] = 1;
                std::vector<T> ue(G, zero);
                const T Sn = num_traits<T>::from_int(static_cast<int>(S));
                for (std::size_t g = 0; g < G; ++g) {
                    T acc = zero;
                    for (std::size_t s = 0; s < S; ++s) {
                        T a = x0v[kclass], b = x0s[src[e] - 1];
                        for (std::size_t f = 0; f < narcs; ++f) {
                            a += sgnClass(f, kclass) * Ys[f](g, s);
                            b += sgnStat(f, src[e] - 1) * Ys[f](g, s);
                        }
                        acc += detail::iv_ups(a, b, spec.nservers[src[e] - 1],
                                              spec.sched[src[e] - 1], spec.capacity[kclass],
                                              capStat[src[e] - 1]);
                    }
                    ue[g] = acc / Sn;
                }
                exposure[p - 1] += detail::iv_trapz(ue, dt);
            }
        }
        for (std::size_t p = 0; p < P; ++p) {
            alpha[p] = spec.alpha0[p] + firings[p];
            beta[p] = spec.beta0[p] + exposure[p];
            alphaTrace(p, it - 1) = alpha[p];
            betaTrace(p, it - 1) = beta[p];
        }

        // the bound is evaluated at the state the iteration ended in, so the
        // rate moments are recomputed against the updated marginals rather
        // than reused from the sweep that produced them
        for (std::size_t e = 0; e < narcs; ++e) {
            rate_moments(e, Ys, ge, he, obsw, false);
            gexp[e] = ge;
            hexp[e] = he;
        }

        // evidence lower bound: path term, observation term and the divergence
        // of the rate posteriors from their priors
        T b = zero;
        {
            std::vector<T> acc(G, zero);
            for (std::size_t e = 0; e < narcs; ++e) {
                for (std::size_t g = 0; g < G; ++g) {
                    T s = zero;
                    for (std::size_t y = 0; y < ny; ++y) {
                        const T n = nu[e](g, y);
                        T term = n - gexp[e](g, y);
                        if (n > zero) {
                            const T hh = hexp[e](g, y) > floorT ? hexp[e](g, y) : floorT;
                            term -= n * log(n / hh);
                        }
                        s += Y[e](g, y) * term;
                    }
                    acc[g] = s;
                }
                b += detail::iv_trapz(acc, dt);
            }
            const T Sn = num_traits<T>::from_int(static_cast<int>(S));
            for (std::size_t k = 0; k < K; ++k) {
                const std::size_t g = obsIdx[k];
                T tot = zero;
                for (std::size_t s = 0; s < S; ++s) {
                    T a = zero;
                    for (std::size_t j = 0; j < MR; ++j) {
                        if (spec.obsData(k, j) == unobs) continue;
                        T x = x0v[j];
                        for (std::size_t f = 0; f < narcs; ++f) x += sgnClass(f, j) * Ys[f](g, s);
                        T p = zero;
                        if (x == spec.obsData(k, j)) {
                            p = one - spec.epsilon;
                        } else if (!(x < zero) && !(x > spec.obsRange[j])) {
                            p = spec.epsilon / (spec.obsRange[j] > one ? spec.obsRange[j] : one);
                        }
                        a += log(floorT + p);
                    }
                    tot += a;
                }
                b += tot / Sn;
            }
            for (std::size_t p = 0; p < P; ++p)
                b -= detail::iv_kl_gamma(alpha[p], beta[p], spec.alpha0[p], spec.beta0[p]);
        }
        bound.push_back(b);

        // The rate update solves a stationarity condition rather than
        // maximising the bound in a block, so the bound need not ascend;
        // convergence is judged on the bound AND on the rate posteriors.
        if (it > 1) {
            const double prevb = num_traits<T>::to_double(bound[it - 2]);
            double crit = std::abs(num_traits<T>::to_double(b) - prevb) / std::max(1.0, std::abs(prevb));
            for (std::size_t p = 0; p < P; ++p) {
                const double prev = num_traits<T>::to_double(alphaTrace(p, it - 2) /
                                                             betaTrace(p, it - 2));
                const double cur = num_traits<T>::to_double(alpha[p] / beta[p]);
                crit = std::max(crit, std::abs(cur - prev) / std::max(1e-12, prev));
            }
            if (crit <= opt.tol) {
                converged = true;
                break;
            }
        }
    }

    VariationalResult<T> out;
    out.alpha = alpha;
    out.beta = beta;
    out.rates.resize(P);
    out.mean_service_time.resize(P);
    for (std::size_t p = 0; p < P; ++p) {
        out.rates[p] = alpha[p] / beta[p];
        out.mean_service_time[p] = beta[p] / alpha[p];
    }
    out.bound = bound;
    out.alpha_trace = Matrix<T>(P, iter, zero);
    out.beta_trace = Matrix<T>(P, iter, zero);
    for (std::size_t p = 0; p < P; ++p) {
        for (std::size_t i = 0; i < iter; ++i) {
            out.alpha_trace(p, i) = alphaTrace(p, i);
            out.beta_trace(p, i) = betaTrace(p, i);
        }
    }
    out.Y = Y;
    out.nu = nu;
    out.tgrid = tgrid;
    out.iter = iter;
    out.converged = converged;

    T tailmass = zero;
    for (std::size_t e = 0; e < narcs; ++e)
        for (std::size_t g = 0; g < G; ++g)
            if (Y[e](g, ny - 1) > tailmass) tailmass = Y[e](g, ny - 1);
    out.tailmass = tailmass;

    out.qlen = Matrix<T>(G, MR, zero);
    for (std::size_t g = 0; g < G; ++g)
        for (std::size_t k = 0; k < MR; ++k) out.qlen(g, k) = x0v[k];
    for (std::size_t e = 0; e < narcs; ++e) {
        for (std::size_t g = 0; g < G; ++g) {
            T my = zero;
            for (std::size_t y = 0; y < ny; ++y) my += Y[e](g, y) * yvec[y];
            for (std::size_t k = 0; k < MR; ++k) out.qlen(g, k) += my * sgnClass(e, k);
        }
    }
    return out;
}

}  // namespace infer
}  // namespace line

#endif  // LINE_API_INFER_INFER_VARIATIONAL_H
