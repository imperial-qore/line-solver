/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MAP_ANFIT_H
#define LINE_API_MAM_MAP_ANFIT_H

/**
 * Fit a superposition of interrupted Poisson processes to a Hurst parameter.
 *
 * Templated port of matlab/lib/kpctoolbox/map/map_anfit.m (Andersen-Nielsen).
 * The construction builds d two-state IPPs whose switching rates form a
 * geometric ladder k(2,i) = a^(1-i) k(2,1) with a = 10^(n/(d-1)), and whose
 * weights phi(i) are chosen so that the aggregate variance-time curve follows
 * t^beta with beta = 2 - 2H over n decades. Superposing them with an optional
 * Poisson remainder gives a MAP whose autocorrelation decays like a fractional
 * process over that range, which is what makes it a long-range-dependence
 * surrogate rather than a moment fit.
 *
 * The phi ladder is built by the reference's backward recursion
 *   D = a^(i beta) - sum_{j<i} phi(d-j)^2 exp(1 - a^(i-j)),
 *   phi(d-i) = sqrt(D) when D >= 0, else 0 and the depth d is grown by one,
 * which is a greedy fill: a level that cannot carry positive variance is set to
 * zero and the ladder is extended instead of failing.
 *
 * REFERENCE INCONSISTENCY, resolved in favour of the least-squares branch.
 * `map_anfit.m` builds the Poisson remainder as `map_exponential(lP)`, and the
 * global `map_exponential` takes a MEAN, whereas the file's own `objfun` builds
 * the SAME component as `{[-lP], [lP]}`, i.e. with lP as a RATE. The two cannot
 * both be right. The rate reading is taken here, for two reasons: `objfun` is
 * unambiguous, and lP = 0 -- which the `ls < L` branch produces deliberately to
 * mean "no Poisson remainder" -- is the null stream under the rate reading and a
 * division by zero under the mean reading.
 *
 * THE LEAST-SQUARES BRANCH replaces MATLAB's `fmincon` with
 * `line/util/auglag.h`, whose header carries the acceptance contract for that
 * substitution. It minimizes the 2-norm between the fitted autocorrelation and
 * the supplied one over the per-IPP ratios r(i), subject to the reference's
 * feasibility constraint ls(i) >= sqrt(k1(i) k2(i) / r(i)), which is what keeps
 * the second IPP rate non-negative.
 *
 * ARITHMETIC: transcendental.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/api/mam/map_transform.h"
#include "line/num/number.h"
#include "line/util/auglag.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

namespace anfitdetail {

/** Superposition of two MAPs: the Kronecker sum of each block. */
template <class T>
Map<T> super2(const Map<T>& a, const Map<T>& b) {
    Map<T> m;
    m.D0 = krons(a.D0, b.D0);
    m.D1 = krons(a.D1, b.D1);
    return map_normalize(m);
}

/** One interrupted Poisson process: on/off switching c1, c2 and rate l. */
template <class T>
Map<T> ipp(const T& c1, const T& c2, const T& l) {
    const T zero = num_traits<T>::from_int(0);
    Map<T> m;
    m.D0 = Matrix<T>(2, 2, zero);
    m.D1 = Matrix<T>(2, 2, zero);
    m.D0(0, 1) = c1;
    m.D0(1, 0) = c2;
    m.D1(0, 0) = l;
    return map_normalize(m);
}

/** The ladder weights phi(1..d), growing d when a level cannot be filled. */
template <class T>
std::vector<double> phi_ladder(double beta, double n, std::size_t ds, std::size_t* dOut,
                               double* aOut) {
    std::size_t d = ds;
    for (std::size_t guard = 0; guard < 64; ++guard) {
        const double a = std::pow(10.0, n / static_cast<double>(d - 1));
        std::vector<double> phi(d + 1, 0.0);  // 1-based, as the reference indexes it
        phi[d] = 1.0;
        bool grew = false;
        for (std::size_t i = 1; i < d; ++i) {
            double S = 0.0;
            for (std::size_t j = 0; j < i; ++j)
                S += phi[d - j] * phi[d - j] *
                     std::exp(1.0 - std::pow(a, static_cast<double>(i - j)));
            const double D = std::pow(a, static_cast<double>(i) * beta) - S;
            if (D < 0.0) {
                phi[d - i] = 0.0;
                // The reference grows the ladder while it still has room to.
                if (ds > d - 1) {
                    ++d;
                    grew = true;
                    break;
                }
            } else {
                phi[d - i] = std::sqrt(D);
            }
        }
        if (grew) continue;
        *dOut = d;
        *aOut = a;
        std::vector<double> out(d + 1, 0.0);
        for (std::size_t i = 1; i <= d; ++i) out[i] = phi[i];
        return out;
    }
    throw NumericError("map_anfit: the phi ladder did not terminate");
}

}  // namespace anfitdetail

/** What map_anfit returns: the fitted MAP and the ladder it was built on. */
template <class T>
struct MapAnfitResult {
    Map<T> map;
    std::size_t d = 0;           ///< number of IPPs actually used
    std::vector<T> switching;    ///< k(2,i), the switching-rate ladder
    std::vector<T> rates;        ///< l(i), the on-state arrival rates
    T poisson_rate = num_traits<T>::from_int(0);  ///< lP, the Poisson remainder
};

/**
 * @param ls  target arrival rate
 * @param rho target burstiness parameter, below 1/2 for the construction to hold
 * @param H   Hurst parameter; beta = 2 - 2H
 * @param n   number of decades the variance-time curve should follow t^beta
 * @param ds  starting number of IPPs; grown when a ladder level cannot be filled
 */
template <class T>
MapAnfitResult<T> map_anfit(const T& ls, const T& rho, const T& H, double n, std::size_t ds) {
    static_assert(num_traits<T>::has_transcendental, "map_anfit builds a geometric rate ladder");
    if (ds < 2) throw InputError("map_anfit: at least two IPPs are required");
    const double lsd = num_traits<T>::to_double(ls), rhod = num_traits<T>::to_double(rho);
    if (!(lsd > 0.0)) throw InputError("map_anfit: the arrival rate must be positive");
    const double beta = 2.0 - 2.0 * num_traits<T>::to_double(H);

    std::size_t d = ds;
    double a = 0.0;
    const std::vector<double> phi = anfitdetail::phi_ladder<T>(beta, n, ds, &d, &a);

    const double k21 = 0.8;
    std::vector<double> k2(d + 1, 0.0);
    for (std::size_t i = 1; i <= d; ++i) k2[i] = std::pow(a, 1.0 - static_cast<double>(i)) * k21;

    double S = 0.0, phisum = 0.0;
    for (std::size_t i = 1; i <= d; ++i) {
        const double kappa = k2[i], e = std::exp(-kappa);
        S += phi[i] * phi[i] / (kappa * kappa) *
             ((1.0 - e) * (1.0 - e) - 2.0 * rhod * (kappa - (1.0 - e)));
        phisum += phi[i];
    }
    if (!(S > 0.0))
        throw NumericError(
            "map_anfit: the ladder carries no variance at this (rho, H, n); the construction has "
            "no interrupted-Poisson superposition for these targets");
    const double eta = std::sqrt(4.0 * rhod * lsd) / std::sqrt(S);
    const double L = eta * phisum / 2.0;

    std::vector<double> c1(d + 1, 0.0), c2v(d + 1, 0.0), l(d + 1, 0.0);
    double lP = 0.0;
    if (lsd < L) {
        for (std::size_t i = 1; i <= d; ++i) {
            c1[i] = L * L / (lsd * lsd + L * L) * k2[i];
            c2v[i] = k2[i] - c1[i];
            l[i] = phi[i] * (lsd * lsd + L * L) / (lsd * phisum);
        }
    } else {
        lP = lsd - L;
        for (std::size_t i = 1; i <= d; ++i) {
            c2v[i] = 0.5 * k2[i];
            c1[i] = c2v[i];
            l[i] = eta * phi[i];
        }
    }

    // The Poisson remainder, as a RATE; see the header note. lP = 0 gives the
    // null stream, which is the identity for superposition.
    Map<T> acc;
    acc.D0 = Matrix<T>(1, 1, num_traits<T>::from_double(-lP));
    acc.D1 = Matrix<T>(1, 1, num_traits<T>::from_double(lP));
    for (std::size_t i = 1; i <= d; ++i)
        acc = anfitdetail::super2(acc, anfitdetail::ipp(num_traits<T>::from_double(c1[i]),
                                                        num_traits<T>::from_double(c2v[i]),
                                                        num_traits<T>::from_double(l[i])));

    MapAnfitResult<T> out;
    out.map = map_normalize(acc);
    out.d = d;
    out.poisson_rate = num_traits<T>::from_double(lP);
    for (std::size_t i = 1; i <= d; ++i) {
        out.switching.push_back(num_traits<T>::from_double(k2[i]));
        out.rates.push_back(num_traits<T>::from_double(l[i]));
    }
    return out;
}

/**
 * The least-squares variant: after the deterministic construction, the per-IPP
 * ratios are tuned so the fitted autocorrelation matches a supplied one.
 *
 * @param SA     target autocorrelation values
 * @param SAlags the lags they were measured at
 */
template <class T>
MapAnfitResult<T> map_anfit_lsq(const T& ls, const T& rho, const T& H, double n, std::size_t ds,
                                const std::vector<T>& SA, const std::vector<unsigned>& SAlags,
                                unsigned iter_max = 100) {
    static_assert(num_traits<T>::has_transcendental, "map_anfit_lsq minimizes over the ladder");
    if (SA.size() != SAlags.size())
        throw InputError("map_anfit_lsq: the autocorrelation values and lags must agree in length");
    if (SA.empty()) throw InputError("map_anfit_lsq: no autocorrelation targets given");

    MapAnfitResult<T> base = map_anfit(ls, rho, H, n, ds);
    const std::size_t d = base.d;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);

    // The reference's per-IPP characteristics, from the deterministic fit.
    std::vector<double> k1(d, 0.0), k2(d, 0.0), lsi(d, 0.0);
    for (std::size_t i = 0; i < d; ++i) {
        const double li = num_traits<T>::to_double(base.rates[i]);
        const double kk = num_traits<T>::to_double(base.switching[i]);
        // c1 = c2 = kk/2 in the branch the reference derives these from.
        const double c1 = 0.5 * kk, c2 = 0.5 * kk;
        k1[i] = li * li * (c1 * c2) / std::pow(c1 + c2, 3.0);
        k2[i] = kk;
        lsi[i] = (c2 * li) / (c1 + c2);
    }

    auto build = [&](const std::vector<T>& r) {
        Map<T> acc;
        acc.D0 = Matrix<T>(1, 1, T(-base.poisson_rate));
        acc.D1 = Matrix<T>(1, 1, base.poisson_rate);
        for (std::size_t i = 0; i < d; ++i) {
            const double ri = num_traits<T>::to_double(r[i]);
            if (!(ri > 0.0)) throw NumericError("map_anfit_lsq: a ratio left the feasible region");
            Map<T> m;
            m.D0 = Matrix<T>(2, 2, zero);
            m.D1 = Matrix<T>(2, 2, zero);
            m.D0(0, 1) = num_traits<T>::from_double(k2[i] * ri / (1.0 + ri));
            m.D0(1, 0) = num_traits<T>::from_double(k2[i] / (1.0 + ri));
            m.D1(0, 0) = num_traits<T>::from_double(lsi[i] + std::sqrt(k1[i] * k2[i] * ri));
            m.D1(1, 1) = num_traits<T>::from_double(lsi[i] - std::sqrt(k1[i] * k2[i] / ri));
            acc = anfitdetail::super2(acc, map_normalize(m));
        }
        return map_normalize(acc);
    };

    auto fobj = [&](const std::vector<T>& r) {
        T s = zero;
        try {
            const std::vector<T> acf = map_acf(build(r), SAlags);
            for (std::size_t i = 0; i < SA.size(); ++i) {
                const T diff = T(acf[i] - SA[i]);
                s += diff * diff;
            }
        } catch (const Error&) {
            return num_traits<T>::from_double(1e30);  // outside the feasible region
        }
        return s;
    };
    auto heq = [&](const std::vector<T>&) { return std::vector<T>(); };
    auto gineq = [&](const std::vector<T>& r) {
        std::vector<T> g(d, zero);
        for (std::size_t i = 0; i < d; ++i) {
            const double ri = num_traits<T>::to_double(r[i]);
            const double lim = ri > 0.0 ? std::sqrt(k1[i] * k2[i] / ri) : 1e30;
            g[i] = num_traits<T>::from_double(lim - lsi[i]);  // <= 0
        }
        return g;
    };

    std::vector<T> r0(d, one), bestx;
    std::vector<Bound<T>> bounds(d);
    for (std::size_t i = 0; i < d; ++i) {
        r0[i] = num_traits<T>::from_double(
            1.0 + (lsi[i] > 0.0 ? k1[i] * k2[i] / (lsi[i] * lsi[i]) : 1.0));
        bounds[i].lo = num_traits<T>::from_double(1e-7);
        bounds[i].hi = num_traits<T>::from_double(1e7);
    }
    AugLagOptions<T> opt = auglag_defaults<T>();
    opt.max_outer = iter_max;
    const AugLagResult<T> res = auglag(fobj, heq, gineq, r0, bounds, opt);

    MapAnfitResult<T> out = base;
    out.map = build(res.x);
    return out;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MAP_ANFIT_H
