/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_MAPG1K_H
#define LINE_API_QSYS_QSYS_MAPG1K_H

/**
 * The MAP/G/1/K queue with tail drop: Markovian arrivals, an arbitrary service
 * law F, and a buffer of K packets counting the one in transmission.
 *
 * Port of matlab/src/api/qsys/qsys_mapg1k.m, which is self-contained (it calls
 * no Q-MAM and no BUTools), so this is a faithful transcription rather than a
 * reconstruction. The service law is NOT fitted to a phase-type distribution:
 * F enters exactly through the functionals A_m and Q_m, evaluated by
 * uniformizing the arrival MAP at theta = max_i (-D0(i,i)).
 *
 * METHOD. The chain embedded at departure epochs has state (n, j) with
 * n = 0..K-1 the packets left behind and j the MAP phase. With A_m the matrix
 * of "m arrivals during one service, phase i -> phase j",
 *   n >= 1: n' = n - 1 + min(m, K - n), the overflow being sum_{m >= K-n} A_m
 *   n == 0: the phase first jumps through Psi = (-D0)^-1 D1, because the idle
 *           period ends at an arrival, and the service then proceeds as from 1.
 * Its stationary law sigma gives, by Markov renewal reward over one
 * inter-departure cycle,
 *   E[cycle] = S + sigma_0 (-D0)^-1 e,  T = 1/E[cycle],  p0 = 1 - T S,
 * and the level-holding times come from Q_m, the expected time within a service
 * during which exactly m arrivals have occurred. No PASTA argument is used
 * anywhere: the MAP phase resolution does that work instead, which is what lets
 * qsys_mmapg1k read exact PER-CLASS loss ratios off pKvec.
 *
 * WHERE THE PORT DIFFERS FROM THE REFERENCE, AND WHY.
 *  - The uniformization coefficients c_n = E[e^{-theta S}(theta S)^n/n!] are
 *    built by RECURSION rather than from log-gamma. For the gamma law they are
 *    the negative binomial pmf, c_0 = (1+theta th)^-al and
 *    c_n = c_{n-1} p (al+n-1)/n with p = theta th/(1+theta th); for the
 *    deterministic law they are the Poisson pmf, c_0 = e^{-theta d} and
 *    c_n = c_{n-1} theta d/n; for a PH law c_n = theta^n alpha M^{n+1} t with
 *    M = (theta I - T)^-1, accumulated as a running row vector. The reference
 *    evaluates each term independently through gammaln, which is the same
 *    number to rounding but costs a special function the port would otherwise
 *    not need at Real50. The recursions are also monotone in n and cannot lose
 *    the leading digits to cancellation.
 *  - The reference WARNS and continues when the c_n series is truncated with a
 *    residual above 1e-6; the port throws. A residual that large means the
 *    uniformization has not converged and every downstream quantity is wrong by
 *    an unknown amount, and a warning nobody reads is the worse failure mode.
 *    The residual is reported in the result either way.
 *  - The density path integrates over an EXPANDING finite window in
 *    u = log x rather than over (-Inf, log tmax] in one call, because the
 *    port's adaptive Gauss-Kronrod rule takes finite endpoints. The window is
 *    widened geometrically until a whole new panel contributes less than the
 *    tolerance, so the truncation is measured rather than assumed. The
 *    substitution itself is the reference's: it turns an integrable
 *    singularity x^(al-1) at the origin into e^(al u), which decays smoothly,
 *    so the singularity disappears instead of being resolved.
 *
 * ARITHMETIC. Gated on num_traits<T>::has_transcendental: c_0 needs exp or a
 * non-integer power for every service law, and the series truncation is a
 * tolerance. Everything after the c_n -- the A_m and Q_m sums, the embedded
 * chain, its stationary vector and the reward averaging -- is finite exact
 * matrix algebra and adds no error of its own.
 *
 * MEASURED AGREEMENT: see cpp/tests/test_qsys_mapg1k.cpp. The M/M/1/K collapse
 * is checked against the closed form (1-rho)rho^K/(1-rho^(K+1)), the M/G/1/K
 * and MAP/G/1/K instances against MATLAB, and three identities are asserted on
 * every instance: sum_l plevel(l) = 1, p0 = 1 - T S, and
 * lossProbability = 1 - T/lambda = pK-weighted arrival loss.
 *
 * References:
 * [1] Chydzinski, A. Per-Flow Throughput of a FIFO Buffer. Applied System
 *     Innovation 2026, 9, 112.
 * [2] Niu, Z.; Cooper, R.B. Transform-Free Analysis of M/G/1/K and Related
 *     Queues. Mathematics of Operations Research 1993, 18, 486-510.
 */

#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <string>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/api/mc/dtmc_solve.h"
#include "line/api/qsys/qsys_quadrature.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace qsys {

/** Which family the service law belongs to. */
enum class ServiceKind { Gamma, Deterministic, PhaseType, Density };

/**
 * Service-time descriptor, the C++ form of the MATLAB svc struct.
 *
 * Build one with the static factories; the fields that a given kind does not
 * use are left at their default and are never read.
 */
template <class T>
struct ServiceLaw {
    ServiceKind kind = ServiceKind::Deterministic;
    T shape = num_traits<T>::from_int(1);   ///< gamma shape alpha
    T scale = num_traits<T>::from_int(1);   ///< gamma scale theta
    T det = num_traits<T>::from_int(1);     ///< deterministic service time
    std::vector<T> ph_alpha;                ///< PH initial probability row
    Matrix<T> ph_T;                         ///< PH subgenerator
    std::function<T(const T&)> pdf;         ///< density of the service time
    T tmax = num_traits<T>::from_int(0);    ///< upper support limit of the density
    bool tmax_finite = false;               ///< whether tmax is used

    /** Gamma(shape, scale). shape = 1 is the exponential, integer shape Erlang. */
    static ServiceLaw gamma(const T& shape, const T& scale) {
        ServiceLaw s;
        s.kind = ServiceKind::Gamma;
        s.shape = shape;
        s.scale = scale;
        return s;
    }
    /** Constant service time d. */
    static ServiceLaw deterministic(const T& d) {
        ServiceLaw s;
        s.kind = ServiceKind::Deterministic;
        s.det = d;
        return s;
    }
    /** Phase type (alpha, Tmat). */
    static ServiceLaw phase_type(const std::vector<T>& alpha, const Matrix<T>& Tmat) {
        ServiceLaw s;
        s.kind = ServiceKind::PhaseType;
        s.ph_alpha = alpha;
        s.ph_T = Tmat;
        return s;
    }
    /** Arbitrary density on (0, inf), or on (0, tmax] when tmax is supplied. */
    static ServiceLaw density(std::function<T(const T&)> f) {
        ServiceLaw s;
        s.kind = ServiceKind::Density;
        s.pdf = f;
        return s;
    }
    static ServiceLaw density(std::function<T(const T&)> f, const T& tmax) {
        ServiceLaw s = density(f);
        s.tmax = tmax;
        s.tmax_finite = true;
        return s;
    }
};

/** Return value of qsys_mapg1k, mirroring the MATLAB result struct. */
template <class T>
struct MapG1kResult {
    T p0;                     ///< P(buffer empty)
    T pK;                     ///< P(buffer full)
    T throughput;             ///< aggregate departure rate
    T lossProbability;        ///< 1 - throughput/lambda
    T lambda;                 ///< aggregate MAP arrival rate
    T meanServiceTime;        ///< S
    T utilization;            ///< 1 - p0
    T rho;                    ///< offered load lambda S
    T meanQueueLength;        ///< E[number in system]
    std::size_t nmax;         ///< uniformization order used
    T countingResidual;       ///< |1 - sum_n c_n| at the truncation
    std::vector<T> sigma;     ///< stationary law of the embedded chain, K M entries
    std::vector<T> pKvec;     ///< P(level = K, phase j), summing to pK
    std::vector<T> p0vec;     ///< P(level = 0, phase j), summing to p0
    std::vector<T> plevel;    ///< P(level = l), l = 0..K
};

namespace mapg1k_detail {

/** Machine epsilon of T, the "a whole block added nothing" threshold. */
template <class T>
T eps_of() {
    return std::numeric_limits<T>::epsilon();
}

/** Mean of the service law. */
template <class T>
T service_mean(const ServiceLaw<T>& svc) {
    switch (svc.kind) {
        case ServiceKind::Gamma:
            return svc.shape * svc.scale;
        case ServiceKind::Deterministic:
            return svc.det;
        case ServiceKind::PhaseType: {
            const std::size_t p = svc.ph_T.rows();
            const std::vector<T> e = ones<T>(p);
            const std::vector<T> x = solve(svc.ph_T, e);  // T^-1 e
            T s = num_traits<T>::from_int(0);
            for (std::size_t i = 0; i < p; ++i) s -= svc.ph_alpha[i] * x[i];
            return s;
        }
        case ServiceKind::Density:
            break;
    }
    // Density: E[S] under the same log substitution used for the c_n.
    return num_traits<T>::from_int(0);  // replaced by the caller, see density_moment
}

/**
 * E[w(S)] for a density-specified service law, under the substitution
 * x = exp(u): the integrand becomes w(e^u) f(e^u) e^u, which tends to zero at
 * both ends because w is bounded and integrability of f forces x f(x) -> 0.
 * The window is widened geometrically to the left (and to the right when the
 * support is unbounded) until a whole new panel adds less than the tolerance.
 */
template <class T, class W>
T density_moment(const ServiceLaw<T>& svc, W&& w, const T& reltol) {
    const T zero = num_traits<T>::from_int(0);
    using std::exp;
    using std::log;
    std::function<T(const T&)> g = [&](const T& u) -> T {
        const T x = exp(u);
        const T v = w(x) * svc.pdf(x) * x;
        return (v == v && num_abs(v) < std::numeric_limits<T>::infinity()) ? v : zero;
    };
    // Seed window around log of a representative scale: start at [-1, 1] and
    // widen by 4 in u each round, so the covered range of x doubles in decades.
    const T hi0 = svc.tmax_finite ? T(log(svc.tmax)) : num_traits<T>::from_int(1);
    T lo = hi0 - num_traits<T>::from_int(2);
    T hi = hi0;
    const T abstol = num_traits<T>::from_double(1e-300);
    T total = detail::num_integral<T>(g, lo, hi, reltol, abstol, 40u);
    const unsigned rounds = 60u;
    for (unsigned k = 0; k < rounds; ++k) {
        const T lonew = lo - num_traits<T>::from_int(4);
        const T add_lo = detail::num_integral<T>(g, lonew, lo, reltol, abstol, 40u);
        lo = lonew;
        T add_hi = zero;
        if (!svc.tmax_finite) {
            const T hinew = hi + num_traits<T>::from_int(4);
            add_hi = detail::num_integral<T>(g, hi, hinew, reltol, abstol, 40u);
            hi = hinew;
        }
        total += add_lo + add_hi;
        const T added = num_abs(T(add_lo)) + num_abs(T(add_hi));
        if (added <= reltol * num_abs(total)) break;
    }
    return total;
}

/**
 * c_n = E[e^{-theta S}(theta S)^n/n!] for n = 0..N, together with the mean.
 * sum_n c_n = E[e^{-theta S} e^{theta S}] = 1 exactly, which both sets the
 * truncation order and certifies it.
 */
template <class T>
struct ServiceCoefficients {
    std::vector<T> cn;
    T mean;
    T residual;  ///< |1 - sum_n c_n|
};

template <class T>
ServiceCoefficients<T> service_coefficients(const ServiceLaw<T>& svc, const T& theta, const T& tol,
                                            std::size_t nmaxCap) {
    using std::exp;
    using std::log;
    using std::pow;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    ServiceCoefficients<T> out;
    const T qreltol = num_traits<T>::from_double(1e-13);
    out.mean = (svc.kind == ServiceKind::Density)
                   ? density_moment(svc, [](const T& x) { return x; }, qreltol)
                   : service_mean(svc);
    if (out.mean <= zero) throw InputError("qsys_mapg1k: service law has non-positive mean");

    // Initial order: the mean of the Poisson-like count plus a deviation
    // allowance, exactly the reference's i_guess.
    std::size_t n0;
    {
        const double md = num_traits<T>::to_double(T(theta * out.mean));
        const double g = md + 10.0 * std::sqrt(md > 1.0 ? md : 1.0) + 32.0;
        n0 = static_cast<std::size_t>(g < 32.0 ? 32.0 : std::ceil(g));
        if (n0 + 1 > nmaxCap) n0 = nmaxCap > 0 ? nmaxCap - 1 : 0;
    }

    std::vector<T>& cn = out.cn;
    // Per-family generator of c_n from c_{n-1}, plus the carried state.
    switch (svc.kind) {
        case ServiceKind::Gamma: {
            const T th = svc.scale, al = svc.shape;
            if (al <= zero || th <= zero)
                throw InputError("qsys_mapg1k: gamma shape and scale must be positive");
            const T q = one + th * theta;
            const T p = th * theta / q;
            cn.push_back(exp(-al * log(q)));
            for (std::size_t n = 1; n <= n0; ++n) {
                const T nT = num_traits<T>::from_int(static_cast<long>(n));
                cn.push_back(cn[n - 1] * p * (al + nT - one) / nT);
            }
            break;
        }
        case ServiceKind::Deterministic: {
            if (svc.det <= zero)
                throw InputError("qsys_mapg1k: deterministic service time must be positive");
            const T m = theta * svc.det;
            cn.push_back(exp(-m));
            for (std::size_t n = 1; n <= n0; ++n)
                cn.push_back(cn[n - 1] * m / num_traits<T>::from_int(static_cast<long>(n)));
            break;
        }
        case ServiceKind::PhaseType: {
            const std::size_t p = svc.ph_T.rows();
            if (svc.ph_alpha.size() != p || svc.ph_T.cols() != p)
                throw InputError("qsys_mapg1k: PH alpha and T are inconsistent");
            Matrix<T> ThI(p, p);
            for (std::size_t i = 0; i < p; ++i)
                for (std::size_t j = 0; j < p; ++j)
                    ThI(i, j) = (i == j ? theta : zero) - svc.ph_T(i, j);
            const Matrix<T> Minv = inverse(ThI);
            const std::vector<T> e = ones<T>(p);
            std::vector<T> t = mulvec(svc.ph_T, e);
            for (T& v : t) v = -v;
            std::vector<T> row = vecmul(svc.ph_alpha, Minv);  // alpha M
            for (std::size_t n = 0; n <= n0; ++n) {
                T v = zero;
                for (std::size_t i = 0; i < p; ++i) v += row[i] * t[i];
                cn.push_back(v);
                row = vecmul(row, Minv);
                for (T& x : row) x *= theta;
            }
            break;
        }
        case ServiceKind::Density: {
            if (!svc.pdf) throw InputError("qsys_mapg1k: density service law has no pdf");
            for (std::size_t n = 0; n <= n0; ++n) {
                const T nT = num_traits<T>::from_int(static_cast<long>(n));
                const T lfact = log(num_factorial<T>(static_cast<unsigned>(n)));
                cn.push_back(density_moment(
                    svc,
                    [&](const T& x) {
                        return exp(-theta * x + nT * log(theta * x) - lfact);
                    },
                    qreltol));
            }
            break;
        }
    }
    // series growth termination rationale: see _kb/03-api-layer.md (cpp port notes: qsys)
    T total = zero;
    for (const T& v : cn) total += v;
    while (cn.size() < nmaxCap) {
        if (num_abs(T(one - total)) <= tol) break;
        const std::size_t first = cn.size();
        const std::size_t last = (first + 63 < nmaxCap) ? first + 63 : nmaxCap - 1;
        T added = zero;
        switch (svc.kind) {
            case ServiceKind::Gamma: {
                const T q = one + svc.scale * theta;
                const T p = svc.scale * theta / q;
                for (std::size_t n = first; n <= last; ++n) {
                    const T nT = num_traits<T>::from_int(static_cast<long>(n));
                    cn.push_back(cn[n - 1] * p * (svc.shape + nT - one) / nT);
                    added += cn.back();
                }
                break;
            }
            case ServiceKind::Deterministic: {
                const T m = theta * svc.det;
                for (std::size_t n = first; n <= last; ++n) {
                    cn.push_back(cn[n - 1] * m / num_traits<T>::from_int(static_cast<long>(n)));
                    added += cn.back();
                }
                break;
            }
            case ServiceKind::PhaseType: {
                // block-restart rationale: see _kb/03-api-layer.md (cpp port notes: qsys)
                const std::size_t p = svc.ph_T.rows();
                Matrix<T> ThI(p, p);
                for (std::size_t i = 0; i < p; ++i)
                    for (std::size_t j = 0; j < p; ++j)
                        ThI(i, j) = (i == j ? theta : zero) - svc.ph_T(i, j);
                const Matrix<T> Minv = inverse(ThI);
                const std::vector<T> e = ones<T>(p);
                std::vector<T> t = mulvec(svc.ph_T, e);
                for (T& v : t) v = -v;
                std::vector<T> row = vecmul(svc.ph_alpha, Minv);
                for (std::size_t n = 0; n < first; ++n) {
                    row = vecmul(row, Minv);
                    for (T& x : row) x *= theta;
                }
                for (std::size_t n = first; n <= last; ++n) {
                    T v = zero;
                    for (std::size_t i = 0; i < p; ++i) v += row[i] * t[i];
                    cn.push_back(v);
                    added += v;
                    row = vecmul(row, Minv);
                    for (T& x : row) x *= theta;
                }
                break;
            }
            case ServiceKind::Density: {
                for (std::size_t n = first; n <= last; ++n) {
                    const T nT = num_traits<T>::from_int(static_cast<long>(n));
                    const T lfact = log(num_factorial<T>(static_cast<unsigned>(n)));
                    cn.push_back(density_moment(
                        svc,
                        [&](const T& x) { return exp(-theta * x + nT * log(theta * x) - lfact); },
                        qreltol));
                    added += cn.back();
                }
                break;
            }
        }
        total += added;
        if (added <= eps_of<T>() * total) break;
    }
    out.residual = num_abs(T(one - total));
    return out;
}

}  // namespace mapg1k_detail

/**
 * MAP/G/1/K with tail drop.
 *
 * @param arrival arrival MAP (D0, D1)
 * @param svc     service law
 * @param K       buffer size in packets, K >= 1, the one in service included
 * @param tol     uniformization truncation tolerance
 * @param nmaxCap cap on the uniformization order
 */
template <class T>
MapG1kResult<T> qsys_mapg1k(const mam::Map<T>& arrival, const ServiceLaw<T>& svc, std::size_t K,
                            const T& tol, std::size_t nmaxCap) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_mapg1k requires transcendental arithmetic");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t M = arrival.D0.rows();
    if (arrival.D0.cols() != M || arrival.D1.rows() != M || arrival.D1.cols() != M)
        throw InputError("qsys_mapg1k: D0 and D1 must be square matrices of equal size");
    if (K < 1) throw InputError("qsys_mapg1k: buffer size K must be a positive integer");
    const Matrix<T>& D0 = arrival.D0;
    const Matrix<T>& D1 = arrival.D1;

    T theta = zero;
    for (std::size_t i = 0; i < M; ++i) {
        const T b = -D0(i, i);
        if (b <= zero)
            throw InputError("qsys_mapg1k: D0 must have strictly negative diagonal entries");
        if (b > theta) theta = b;
    }

    const mapg1k_detail::ServiceCoefficients<T> sc =
        mapg1k_detail::service_coefficients(svc, theta, tol, nmaxCap);
    const std::vector<T>& cn = sc.cn;
    const T Smean = sc.mean;
    const std::size_t nmax = cn.size() - 1;
    if (sc.residual > num_traits<T>::from_double(1e-6))
        throw NumericError(
            "qsys_mapg1k: the uniformization series for c_n was truncated with residual " +
            num_traits<T>::to_string(sc.residual) + "; raise nmax");

    // telescoping-form rationale: see _kb/03-api-layer.md (cpp port notes: qsys)
    std::vector<T> dn(nmax + 1, zero);
    {
        T tail = zero;
        for (std::size_t n = nmax + 1; n-- > 0;) {
            dn[n] = tail / theta;  // tail is sum_{k > n} c_k at this point
            tail += cn[n];
        }
    }

    // A_m and Q_m for m = 0..K-1, plus B0 = sum_m A_m and Qtot = sum_m Q_m.
    const std::size_t mmax = K - 1;
    std::vector<Matrix<T>> A(mmax + 1, Matrix<T>(M, M, zero));
    std::vector<Matrix<T>> Q(mmax + 1, Matrix<T>(M, M, zero));
    std::vector<Matrix<T>> Sn(mmax + 1, Matrix<T>(M, M, zero));
    Sn[0] = eye<T>(M);
    Matrix<T> B0(M, M, zero), Qtot(M, M, zero);
    Matrix<T> Pn = eye<T>(M);
    Matrix<T> Pt0(M, M), Pt1(M, M), PD(M, M);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t j = 0; j < M; ++j) {
            Pt0(i, j) = (i == j ? one : zero) + D0(i, j) / theta;
            Pt1(i, j) = D1(i, j) / theta;
            PD(i, j) = (i == j ? one : zero) + (D0(i, j) + D1(i, j)) / theta;
        }
    for (std::size_t n = 0; n <= nmax; ++n) {
        const std::size_t mtop = (n < mmax) ? n : mmax;
        for (std::size_t m = 0; m <= mtop; ++m)
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t j = 0; j < M; ++j) {
                    A[m](i, j) += Sn[m](i, j) * cn[n];
                    Q[m](i, j) += Sn[m](i, j) * dn[n];
                }
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t j = 0; j < M; ++j) {
                B0(i, j) += Pn(i, j) * cn[n];
                Qtot(i, j) += Pn(i, j) * dn[n];
            }
        if (n < nmax) {
            std::vector<Matrix<T>> Snew(mmax + 1, Matrix<T>(M, M, zero));
            const std::size_t mt = (n + 1 < mmax) ? n + 1 : mmax;
            for (std::size_t m = 0; m <= mt; ++m) {
                Matrix<T> acc(M, M, zero);
                if (m <= n) acc = matmul(Sn[m], Pt0);
                if (m >= 1 && m - 1 <= n) {
                    const Matrix<T> add = matmul(Sn[m - 1], Pt1);
                    for (std::size_t i = 0; i < M; ++i)
                        for (std::size_t j = 0; j < M; ++j) acc(i, j) += add(i, j);
                }
                Snew[m] = acc;
            }
            Sn = Snew;
            Pn = matmul(Pn, PD);
        }
    }

    const std::vector<T> e = ones<T>(M);
    Matrix<T> negD0(M, M);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t j = 0; j < M; ++j) negD0(i, j) = -D0(i, j);
    const Matrix<T> negD0inv = inverse(negD0);
    const Matrix<T> Psi = matmul(negD0inv, D1);  // phase at the arrival ending an idle period
    const std::vector<T> idle = mulvec(negD0inv, e);

    // Embedded chain at departure epochs, state (n,j) -> index n M + j.
    Matrix<T> P(K * M, K * M, zero);
    const std::size_t lastblk = (K - 1) * M;
    for (std::size_t n = 1; n + 1 <= K; ++n) {
        Matrix<T> Bacc = B0;
        for (std::size_t m = 0; m + n + 1 <= K; ++m) {
            const std::size_t col = (n - 1 + m) * M;
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t j = 0; j < M; ++j) {
                    P(n * M + i, col + j) += A[m](i, j);
                    Bacc(i, j) -= A[m](i, j);
                }
        }
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t j = 0; j < M; ++j) P(n * M + i, lastblk + j) += Bacc(i, j);
    }
    {
        Matrix<T> Bacc = B0;
        for (std::size_t m = 0; m + 2 <= K; ++m) {
            const Matrix<T> blk = matmul(Psi, A[m]);
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t j = 0; j < M; ++j) {
                    P(i, m * M + j) += blk(i, j);
                    Bacc(i, j) -= A[m](i, j);
                }
        }
        const Matrix<T> tailblk = matmul(Psi, Bacc);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t j = 0; j < M; ++j) P(i, lastblk + j) += tailblk(i, j);
    }

    {
        T rowdev = zero;
        for (std::size_t i = 0; i < K * M; ++i) {
            T s = zero;
            for (std::size_t j = 0; j < K * M; ++j) s += P(i, j);
            const T d = num_abs(T(s - one));
            if (d > rowdev) rowdev = d;
        }
        if (rowdev > num_traits<T>::from_double(1e-8))
            throw NumericError(
                "qsys_mapg1k: embedded chain rows deviate from 1 by " +
                num_traits<T>::to_string(rowdev) +
                "; the uniformization series for A_m has not converged, raise nmax");
    }

    const std::vector<T> sigma = mc::dtmc_solve(P);
    std::vector<T> sigma0(sigma.begin(), sigma.begin() + M);

    // Markov renewal reward over one inter-departure cycle.
    T idleTime = zero;
    for (std::size_t j = 0; j < M; ++j) idleTime += sigma0[j] * idle[j];
    const T Ecyc = Smean + idleTime;
    const T Tput = one / Ecyc;
    const T p0 = idleTime / Ecyc;

    // Expected time with the buffer full during a cycle, resolved by phase.
    std::vector<Matrix<T>> Qcum(mmax + 1, Matrix<T>(M, M, zero));
    {
        Matrix<T> acc(M, M, zero);
        for (std::size_t m = 0; m <= mmax; ++m) {
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t j = 0; j < M; ++j) acc(i, j) += Q[m](i, j);
            Qcum[m] = acc;
        }
    }
    std::vector<T> timeKvec(M, zero);
    for (std::size_t n = 1; n + 1 <= K; ++n) {
        const std::size_t rr = K - n - 1;
        std::vector<T> row(sigma.begin() + n * M, sigma.begin() + (n + 1) * M);
        Matrix<T> D = Qtot;
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t j = 0; j < M; ++j) D(i, j) -= Qcum[rr](i, j);
        const std::vector<T> t = vecmul(row, D);
        for (std::size_t j = 0; j < M; ++j) timeKvec[j] += t[j];
    }
    {
        const std::vector<T> s0Psi = vecmul(sigma0, Psi);
        Matrix<T> D = Qtot;
        if (K >= 2)
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t j = 0; j < M; ++j) D(i, j) -= Qcum[K - 2](i, j);
        const std::vector<T> t = vecmul(s0Psi, D);
        for (std::size_t j = 0; j < M; ++j) timeKvec[j] += t[j];
    }
    std::vector<T> pKvec(M, zero);
    T pK = zero;
    for (std::size_t j = 0; j < M; ++j) {
        pKvec[j] = timeKvec[j] / Ecyc;
        pK += pKvec[j];
    }

    // Expected time at every level, from the same Q_m.
    std::vector<T> timeL(K + 1, zero);
    timeL[0] = idleTime;
    for (std::size_t n = 1; n + 1 <= K; ++n) {
        std::vector<T> row(sigma.begin() + n * M, sigma.begin() + (n + 1) * M);
        for (std::size_t l = n; l + 1 <= K; ++l) {
            const std::vector<T> t = vecmul(row, Q[l - n]);
            for (std::size_t j = 0; j < M; ++j) timeL[l] += t[j];
        }
    }
    {
        const std::vector<T> s0Psi = vecmul(sigma0, Psi);
        for (std::size_t l = 1; l + 1 <= K; ++l) {
            const std::vector<T> t = vecmul(s0Psi, Q[l - 1]);
            for (std::size_t j = 0; j < M; ++j) timeL[l] += t[j];
        }
    }
    {
        T s = zero;
        for (std::size_t j = 0; j < M; ++j) s += timeKvec[j];
        timeL[K] = s;
    }
    std::vector<T> plevel(K + 1, zero);
    T mass = zero;
    for (std::size_t l = 0; l <= K; ++l) {
        plevel[l] = timeL[l] / Ecyc;
        mass += plevel[l];
    }
    if (num_abs(T(mass - one)) > num_traits<T>::from_double(1e-8))
        throw NumericError("qsys_mapg1k: the level distribution has mass " +
                           num_traits<T>::to_string(mass) +
                           "; the Q_m series has not converged, raise nmax");
    T meanQ = zero;
    for (std::size_t l = 0; l <= K; ++l)
        meanQ += num_traits<T>::from_int(static_cast<long>(l)) * plevel[l];

    std::vector<T> p0vec = vecmul(sigma0, negD0inv);
    for (T& v : p0vec) v /= Ecyc;

    const T lambda = mam::map_lambda(arrival);

    MapG1kResult<T> r;
    r.p0 = p0;
    r.pK = pK;
    r.throughput = Tput;
    r.lossProbability = one - Tput / lambda;
    r.lambda = lambda;
    r.meanServiceTime = Smean;
    r.utilization = one - p0;
    r.rho = lambda * Smean;
    r.meanQueueLength = meanQ;
    r.nmax = nmax;
    r.countingResidual = sc.residual;
    r.sigma = sigma;
    r.pKvec = pKvec;
    r.p0vec = p0vec;
    r.plevel = plevel;
    return r;
}

/** qsys_mapg1k with the reference defaults tol = 1e-12, nmax = 200000. */
template <class T>
MapG1kResult<T> qsys_mapg1k(const mam::Map<T>& arrival, const ServiceLaw<T>& svc, std::size_t K) {
    return qsys_mapg1k(arrival, svc, K, T(num_traits<T>::from_double(1e-12)),
                       static_cast<std::size_t>(200000));
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_MAPG1K_H
