/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_CYCLET_OFREE_H
#define LINE_API_PFQN_PFQN_CYCLET_OFREE_H

/**
 * Exact passage-time law along an OVERTAKE-FREE PATH of a closed single-chain
 * tree-like product-form network.
 *
 * Port of matlab/src/api/pfqn/pfqn_cyclet_ofree.m, which is the reference.
 * Source: P. G. Harrison and W. J. Knottenbelt, "Passage Time Distributions in
 * Large Markov Chains", 2002, Sec. 7.1, Theorems 1 and 2, after P. G. Harrison,
 * J. Appl. Prob. 27, 1990 and H. Duduna, Adv. Appl. Prob. 14, 1982. The
 * underlying sojourn-time result for overtake-free paths is F. Kelly and
 * P. Pollett, Adv. Appl. Prob. 15, 1983.
 *
 * THE ONE FACT THAT MAKES ALL THREE ROUTES WORK. Conditional on the path,
 *
 *     T | z  =  sum_{j in z} Erlang(u_{z_j} + 1, mu_{z_j})
 *
 * with u distributed as the network's equilibrium population vector AT N-1 (the
 * arrival theorem). Hence the transform of Theorem 1 collapses to
 *
 *     L(s|z) = prod_{j in z} mu_j/(s+mu_j) * G(y(s), N-1) / G(x, N-1)
 *
 * where x_i = v_i/mu_i and y_i(s) = x_i mu_i/(s+mu_i) on the path, x_i off it.
 * One Buzen convolution per value of s.
 *
 * MOMENTS ARE NEVER TAKEN FROM THE DENSITY. They come from running the same
 * Buzen convolution in the ring of truncated power series in s, so they are
 * exact to machine precision, are unaffected by the time grid, and stay valid
 * when the rates coincide and Theorem 2 does not apply.
 *
 * NOTE ON THE PAPER. The inner sum of Theorem 2 reads (v_j t)^(c-i)/(c-i)! and
 * that is CORRECT as printed, however odd the visit ratio looks against a time:
 * substituting the service rate instead returns negative densities. Verified
 * against a direct mixture-of-Erlangs oracle to 1e-15, and at the paper's own
 * N = 18 example against the transform route to 1e-11.
 *
 * ARITHMETIC: double. The density needs exp and an incomplete gamma.
 */

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <string>
#include <vector>

#include "line/api/lti/laplace_invert.h"
#include "line/util/error.h"

namespace line {
namespace pfqn {

/** One path's outcome: which route ran, and the network constant it used. */
struct CycletPathInfo {
    std::string method;
    double lG = 0.0;
    std::vector<std::size_t> path;
};

/** Density, distribution and moments of the passage time along a path. */
struct CycletResult {
    std::vector<double> f;
    std::vector<double> F;
    std::vector<double> mom;
    std::vector<CycletPathInfo> info;
};

namespace cyclet_detail {

/**
 * Buzen's convolution: g[k] = G at population k for the node set y, k = 0..n.
 * This is the k(y,a,b) recursion of Sec. 7.1 with the node index rolled up,
 * k(y,a,b) = k(y,a-1,b) + y_a k(y,a,b-1), k(y,a,0) = 1, k(y,0,b>0) = 0.
 */
template <class S>
std::vector<S> buzen(const std::vector<S>& y, std::size_t n) {
    std::vector<S> g(n + 1, S(0.0));
    g[0] = S(1.0);
    for (std::size_t i = 0; i < y.size(); ++i)
        for (std::size_t k = 1; k <= n; ++k) g[k] = g[k] + y[i] * g[k - 1];
    return g;
}

inline std::vector<double> series_mul(const std::vector<double>& a, const std::vector<double>& b,
                                      std::size_t K) {
    std::vector<double> c(K + 1, 0.0);
    for (std::size_t i = 0; i <= K; ++i) {
        if (a[i] == 0.0) continue;
        for (std::size_t j = 0; j + i <= K; ++j) c[i + j] += a[i] * b[j];
    }
    return c;
}

/**
 * The regularized lower incomplete gamma P(a, x) for INTEGER a, which is all
 * this file needs: P(k+1, x) = 1 - exp(-x) sum_{j<k+1} x^j/j!. The finite sum is
 * exact for the integer shapes the Erlang terms produce, so no series/continued
 * fraction switch is needed.
 */
inline double gammainc_int(std::size_t k, double x) {
    if (!(x > 0.0)) return 0.0;
    double term = std::exp(-x);
    double acc = term;
    for (std::size_t j = 1; j <= k; ++j) {
        term *= x / double(j);
        acc += term;
    }
    return std::min(1.0, std::max(0.0, 1.0 - acc));
}

inline double factorial_d(std::size_t k) {
    double r = 1.0;
    for (std::size_t i = 2; i <= k; ++i) r *= double(i);
    return r;
}

}  // namespace cyclet_detail

/**
 * Exact passage-time density, CDF and moments along the overtake-free paths
 * `paths`, mixed by `pathprob`.
 *
 * @param method "auto" (default) uses "exact" when the path rates are separated
 *               and "lt" otherwise; "exact" is Theorem 2 in closed form and
 *               REQUIRES DISTINCT RATES on the path, since its partial fractions
 *               divide by prod_{i!=j}(mu_i - mu_j)
 */
inline CycletResult pfqn_cyclet_ofree(const std::vector<double>& v, const std::vector<double>& mu,
                                      std::size_t N,
                                      const std::vector<std::vector<std::size_t>>& paths,
                                      const std::vector<double>& tset,
                                      const std::string& method = "auto", std::size_t nmom = 3,
                                      const std::vector<double>& pathprob = {},
                                      const std::string& lti_method = "euler",
                                      double tol = 1e-8) {
    const std::size_t M = v.size();
    if (mu.size() != M)
        throw InputError("pfqn_cyclet_ofree: v and mu must name the same number of nodes");
    for (std::size_t i = 0; i < M; ++i)
        if (!(mu[i] > 0.0))
            throw InputError("pfqn_cyclet_ofree: every service rate must be positive");
    if (N < 1) throw InputError("pfqn_cyclet_ofree: the population N must be positive");
    if (paths.empty()) throw InputError("pfqn_cyclet_ofree: no path was given");

    std::vector<double> pp = pathprob;
    if (pp.empty()) pp.assign(paths.size(), 1.0 / double(paths.size()));
    if (pp.size() != paths.size())
        throw InputError("pfqn_cyclet_ofree: pathprob must carry one probability per path");

    std::vector<double> x(M);
    for (std::size_t i = 0; i < M; ++i) x[i] = v[i] / mu[i];
    const double Gn1 = cyclet_detail::buzen(x, N - 1).back();
    if (!(Gn1 > 0.0))
        throw NumericError(
            "pfqn_cyclet_ofree: the network normalizing constant at population N-1 vanished; "
            "check v and mu");

    CycletResult out;
    out.f.assign(tset.size(), 0.0);
    out.F.assign(tset.size(), 0.0);
    out.mom.assign(nmom, 0.0);

    for (std::size_t ip = 0; ip < paths.size(); ++ip) {
        const std::vector<std::size_t>& z = paths[ip];
        if (z.empty())
            throw InputError(
                "pfqn_cyclet_ofree: an overtake-free path must contain at least the root node");
        for (std::size_t a = 0; a < z.size(); ++a) {
            if (z[a] >= M)
                throw InputError("pfqn_cyclet_ofree: a path node is outside the network");
            for (std::size_t b = a + 1; b < z.size(); ++b)
                if (z[a] == z[b])
                    throw InputError("pfqn_cyclet_ofree: a path must have distinct nodes");
        }
        const std::size_t m = z.size();
        std::vector<bool> onpath(M, false);
        for (std::size_t j : z) onpath[j] = true;

        std::string mth = method;
        if (mth == "auto") {
            if (m == 1) {
                mth = "exact";
            } else {
                double sep = std::numeric_limits<double>::infinity();
                double mx = 0.0;
                for (std::size_t a = 0; a < m; ++a) {
                    mx = std::max(mx, mu[z[a]]);
                    for (std::size_t b = a + 1; b < m; ++b)
                        sep = std::min(sep, std::abs(mu[z[a]] - mu[z[b]]));
                }
                mth = (sep > tol * mx) ? "exact" : "lt";
            }
        }

        std::vector<double> fi(tset.size(), 0.0), Fi(tset.size(), 0.0);
        if (mth == "exact") {
            std::vector<double> xoff;
            for (std::size_t i = 0; i < M; ++i)
                if (!onpath[i]) xoff.push_back(x[i]);
            const std::vector<double> Gm = cyclet_detail::buzen(xoff, N - 1);

            // coef[j][k] multiplies t^k exp(-mu_j t)
            std::vector<std::vector<double>> coef(m, std::vector<double>(N, 0.0));
            for (std::size_t j = 0; j < m; ++j) {
                double den = 1.0;
                for (std::size_t i = 0; i < m; ++i)
                    if (i != j) den *= (mu[z[i]] - mu[z[j]]);
                if (den == 0.0)
                    throw InputError(
                        "pfqn_cyclet_ofree: Theorem 2 needs distinct service rates on the path; "
                        "two coincide. Use method 'lt'");
                std::vector<double> w;
                for (std::size_t i = 0; i < m; ++i)
                    if (i != j) w.push_back((v[z[i]] - v[z[j]]) / (mu[z[i]] - mu[z[j]]));
                const std::vector<double> K = cyclet_detail::buzen(w, N - 1);
                for (std::size_t c = 0; c < N; ++c) {
                    const double Gmc = Gm[N - 1 - c];
                    if (Gmc == 0.0) continue;
                    for (std::size_t i = 0; i <= c; ++i) coef[j][c - i] += Gmc * K[i] / den;
                }
            }

            double pref = 1.0 / Gn1;
            for (std::size_t j = 0; j < m; ++j) pref *= mu[z[j]];
            for (std::size_t j = 0; j < m; ++j) {
                const double mj = mu[z[j]], vj = v[z[j]];
                for (std::size_t k = 0; k < N; ++k) {
                    const double c = coef[j][k];
                    if (c == 0.0) continue;
                    const double vk = std::pow(vj, double(k));
                    const double kf = cyclet_detail::factorial_d(k);
                    for (std::size_t it = 0; it < tset.size(); ++it) {
                        const double t = tset[it];
                        if (t < 0.0) continue;
                        fi[it] += pref * c * vk * std::pow(t, double(k)) / kf * std::exp(-mj * t);
                        // int_0^t s^k exp(-mu s) ds = k!/mu^(k+1) P(k+1, mu t)
                        Fi[it] += pref * c * vk / std::pow(mj, double(k + 1)) *
                                  cyclet_detail::gammainc_int(k, mj * t);
                    }
                }
            }
        } else if (mth == "lt") {
            using C = std::complex<double>;
            const lti::LaplaceFn L = [&](C s) {
                std::vector<C> y(M);
                for (std::size_t i = 0; i < M; ++i) y[i] = C(x[i], 0.0);
                for (std::size_t j : z) y[j] = y[j] * (C(mu[j], 0.0) / (s + C(mu[j], 0.0)));
                C acc = cyclet_detail::buzen(y, N - 1).back() / C(Gn1, 0.0);
                for (std::size_t j : z) acc = acc * (C(mu[j], 0.0) / (s + C(mu[j], 0.0)));
                return acc;
            };
            const lti::LaplaceMethod lm = lti::laplace_method(lti_method);
            fi = lti::laplace_invert_pdf(L, tset, lm);
            Fi = lti::laplace_invert_cdf(L, tset, lm);
        } else {
            throw InputError("pfqn_cyclet_ofree: unknown method '" + method +
                             "', expected auto, exact or lt");
        }

        // Moments by the series-ring Buzen convolution.
        const std::size_t K = nmom;
        std::vector<std::vector<double>> Y(M, std::vector<double>(K + 1, 0.0));
        for (std::size_t i = 0; i < M; ++i) {
            if (onpath[i]) {
                double p = 1.0;
                for (std::size_t k = 0; k <= K; ++k) {
                    Y[i][k] = x[i] * p;
                    p *= (-1.0 / mu[i]);
                }
            } else {
                Y[i][0] = x[i];
            }
        }
        std::vector<std::vector<double>> G(N, std::vector<double>(K + 1, 0.0));
        G[0][0] = 1.0;
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t nn = 1; nn < N; ++nn) {
                const std::vector<double> add = cyclet_detail::series_mul(Y[i], G[nn - 1], K);
                for (std::size_t k = 0; k <= K; ++k) G[nn][k] += add[k];
            }
        std::vector<double> L(K + 1, 0.0);
        for (std::size_t k = 0; k <= K; ++k) L[k] = G[N - 1][k] / Gn1;
        for (std::size_t j : z) {
            std::vector<double> e(K + 1, 0.0);
            double p = 1.0;
            for (std::size_t k = 0; k <= K; ++k) {
                e[k] = p;
                p *= (-1.0 / mu[j]);
            }
            L = cyclet_detail::series_mul(L, e, K);
        }
        std::vector<double> momi(nmom, 0.0);
        for (std::size_t q = 1; q <= nmom; ++q)
            momi[q - 1] = ((q % 2) ? -1.0 : 1.0) * cyclet_detail::factorial_d(q) * L[q];

        for (std::size_t it = 0; it < tset.size(); ++it) {
            out.f[it] += pp[ip] * fi[it];
            out.F[it] += pp[ip] * Fi[it];
        }
        for (std::size_t q = 0; q < nmom; ++q) out.mom[q] += pp[ip] * momi[q];

        CycletPathInfo pi;
        pi.method = mth;
        pi.lG = std::log(Gn1);
        pi.path = z;
        out.info.push_back(pi);
    }

    for (std::size_t it = 0; it < tset.size(); ++it) {
        out.f[it] = std::max(0.0, out.f[it]);
        out.F[it] = std::min(1.0, std::max(0.0, out.F[it]));
    }
    return out;
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PFQN_CYCLET_OFREE_H
