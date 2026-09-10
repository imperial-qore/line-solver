/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MC_CTMC_FAU_H
#define LINE_API_MC_CTMC_FAU_H

/**
 * Transient distribution of a CTMC by fast adaptive uniformization.
 *
 * Templated port of matlab/src/api/mc/ctmc_fau.m and
 * jar/src/main/java/jline/api/mc/Ctmc_fau.java.
 *
 * Ordinary uniformization fixes one rate q >= max_i |q_ii| over the WHOLE state
 * space and mixes the powers of P = I + Q/q against a Poisson(q t) law, so its
 * cost is set by the fastest state anywhere, including states carrying no
 * probability at time t. Adaptive uniformization (van Moorsel and Sanders,
 * 1994) picks a rate per step from the states the iterate occupies,
 *
 *   Lambda_n >= max{|q_ii| : i in supp(u^(n))},  u^(n+1) = u^(n)(I + Q/Lambda_n),
 *
 * which keeps every entry of u^(n+1) nonnegative. The subordinating process is
 * then the pure birth process N(t) with rates Lambda_0, Lambda_1, ... and
 * pi(t) = sum_n P{N(t) = n} u^(n). The fast variant (Mateescu, Wolf, Didier and
 * Henzinger, 2010) drops an entry below delta rather than propagating it, so
 * the support tracks the states of non-negligible occupancy instead of the
 * reachable set.
 *
 * NOTHING IS RENORMALIZED, so the error is measured rather than estimated: the
 * birth index truncated at K, the Poisson window of the weight computation and
 * the delta threshold each remove mass and none puts any back, whence
 * 0 <= pi(t) - pit componentwise and |pi(t) - pit|_1 = sum(pi0) - sum(pit),
 * which is what errorBound reports.
 *
 * The birth weights are exact rather than quadratured: the rates generate a
 * bidiagonal generator on the birth index plus one absorbing overflow index,
 * and its transient distribution is obtained by uniformizing that scalar chain
 * at Lstar = max_n Lambda_n and mixing the shipped Fox-Glynn weights. The sweep
 * runs twice because b_n(t) needs the rates up to n, which are not known before
 * the sweep ends, while u^(n) is needed after them, and storing every iterate
 * would cost K times the support.
 *
 * GATED ON TRANSCENDENTAL ARITHMETIC, exactly as ctmc_foxglynn is: the stopping
 * rule is a logarithmic tail bound and the mixture is an approximation of
 * exp(Qt) controlled by epsilon, so the rate bookkeeping and the truncation
 * decisions are taken in double while the iterate and the weighted sum run in
 * T, which is what a high-precision instantiation buys.
 *
 * This is a transient method: it produces no stationary distribution.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <vector>

#include "line/api/mc/ctmc_foxglynn.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mc {

/** Default cap on birth steps, so a pathological horizon reports truncation. */
constexpr long FAU_MAX_STEPS = 1000000;

template <class T>
struct FauResult {
    std::vector<T> pit;      ///< defective distribution at t, a lower bound on pi(t)
    long steps = 0;          ///< number of birth steps K+1 actually taken
    double lambdaMin = 0.0;  ///< smallest adaptive rate used
    double lambdaMax = 0.0;  ///< largest adaptive rate used, the Lstar of the weights
    double uniformRate = 0.0;///< max_i |q_ii|, the rate ordinary uniformization would use
    T weightTail;            ///< mass reaching the overflow index, i.e. P{N(t) > K}
    T weightWindow;          ///< Poisson mass outside the Fox-Glynn window
    T droppedMass;           ///< probability removed by the occupancy threshold
    T errorBound;            ///< sum(pi0) - sum(pit), which IS the L1 error
    std::size_t supportMax = 0;   ///< largest occupied support over the sweep
    std::size_t supportFinal = 0; ///< support at the last step
    bool truncated = false;  ///< maxsteps stopped the sweep
    bool absorbed = false;   ///< the support emptied or became absorbing
};

namespace detail {

/**
 * Upper bound on P{N(t) >= k} for the birth process, through the stochastic
 * domination of its k-th jump epoch by an Erlang(k, lstar): the bound is the
 * Poisson(lstar t) upper tail P{X >= k} at its Chernoff exponent. That exponent
 * bounds the upper tail only above the mean, so below it the bound is vacuous.
 */
inline double fau_tailbound(double lstar, double t, long k) {
    const double lambda = lstar * t;
    const double kd = static_cast<double>(k);
    if (lambda <= 0.0 || kd <= lambda) return 1.0;
    return std::exp(-(lambda - kd + kd * std::log(kd / lambda)));
}

/** Indices holding positive mass, ascending. */
template <class T>
std::vector<std::size_t> fau_support(const std::vector<T>& u) {
    std::vector<std::size_t> act;
    for (std::size_t i = 0; i < u.size(); ++i)
        if (u[i] > num_traits<T>::from_int(0)) act.push_back(i);
    return act;
}

/** Largest exit rate over the occupied states. */
inline double fau_max_exit(const std::vector<double>& d, const std::vector<std::size_t>& act) {
    double L = 0.0;
    for (std::size_t k = 0; k < act.size(); ++k)
        if (d[act[k]] > L) L = d[act[k]];
    return L;
}

/**
 * One adaptive uniformization step u <- u(I + Q/L), reading only the rows of Q
 * in the current support, followed by the drop rule. A state with a zero exit
 * rate is absorbing: its row of Q is empty, so it holds its mass and stays in
 * the support. DROPPED, when not null, accumulates the mass the threshold
 * removes.
 */
template <class T>
std::vector<std::size_t> fau_step(std::vector<T>& u, const std::vector<std::size_t>& act,
                                  const Matrix<T>& Q, const T& L, const T& delta,
                                  std::vector<T>& scratch, std::vector<char>& touchedFlag,
                                  T* dropped) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t n = u.size();
    std::vector<std::size_t> touched;
    for (std::size_t k = 0; k < act.size(); ++k) {
        const std::size_t i = act[k];
        const T& ui = u[i];
        for (std::size_t j = 0; j < n; ++j) {
            const T& qij = Q(i, j);
            if (qij == zero) continue;
            if (!touchedFlag[j]) {
                touchedFlag[j] = 1;
                touched.push_back(j);
            }
            scratch[j] += ui * qij;
        }
    }
    if (touched.empty()) return act;
    std::sort(touched.begin(), touched.end());

    std::vector<std::size_t> written;
    written.reserve(touched.size());
    for (std::size_t k = 0; k < touched.size(); ++k) {
        const std::size_t j = touched[k];
        const T contrib = scratch[j];
        scratch[j] = zero;
        touchedFlag[j] = 0;
        // A row that cancels exactly leaves its state untouched; the union
        // below carries the surviving part of the old support anyway.
        if (contrib == zero) continue;
        T v = u[j] + contrib / L;
        if (v < delta) {
            if (dropped != nullptr && v > zero) *dropped += v;
            v = zero;
        }
        u[j] = v;
        if (v > zero) written.push_back(j);
    }
    std::vector<std::size_t> survivors;
    survivors.reserve(act.size());
    for (std::size_t k = 0; k < act.size(); ++k)
        if (u[act[k]] > zero) survivors.push_back(act[k]);

    std::vector<std::size_t> out;
    out.reserve(survivors.size() + written.size());
    std::set_union(survivors.begin(), survivors.end(), written.begin(), written.end(),
                   std::back_inserter(out));
    return out;
}

/**
 * Transient distribution of the pure birth process with rates LAMBDA at time T,
 * that is b_n = P{N(t) = n} for n = 0..K, plus the mass that reached the
 * absorbing overflow index and therefore measures P{N(t) > K}, and the Poisson
 * mass left outside the Fox-Glynn window. The chain is uniformized at
 * Lstar = max(LAMBDA), so the kernel entries 1 - Lambda_n/Lstar and
 * Lambda_n/Lstar are probabilities and nothing cancels; the weights are taken
 * unnormalized so that the window loss stays visible as missing mass.
 */
template <class T>
void fau_weights(const std::vector<double>& lambda, const T& t, double tDouble, double tol,
                 std::vector<T>& b, T& tail, T& window) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t k1 = lambda.size();
    b.assign(k1, zero);
    tail = zero;
    window = zero;
    if (k1 == 0) return;
    double lstar = 0.0;
    for (std::size_t m = 0; m < k1; ++m) lstar = std::max(lstar, lambda[m]);
    if (lstar <= 0.0 || tDouble <= 0.0) {
        b[0] = num_traits<T>::from_int(1);
        return;
    }
    const double lambdaDouble = lstar * tDouble;
    const long left = detail::foxglynn_left(lambdaDouble, tol);
    const long right = detail::foxglynn_right(lambdaDouble, tol);
    const T lstarT = num_traits<T>::from_double(lstar);
    const std::vector<T> w =
        detail::foxglynn_poisson(lstarT * t, left, right, lambdaDouble, false);

    T wsum = zero;
    for (std::size_t i = 0; i < w.size(); ++i) wsum += w[i];
    window = (num_traits<T>::from_int(1) > wsum) ? num_traits<T>::from_int(1) - wsum : zero;

    std::vector<T> v(k1 + 1, zero);
    std::vector<T> acc(k1 + 1, zero);
    v[0] = num_traits<T>::from_int(1);
    std::vector<T> c(k1, zero);
    std::vector<T> a(k1, zero);
    for (std::size_t m = 0; m < k1; ++m) {
        c[m] = num_traits<T>::from_double(lambda[m]) / lstarT;
        a[m] = num_traits<T>::from_int(1) - c[m];
    }
    for (long k = 0; k <= right; ++k) {
        if (k >= left) {
            const T& wk = w[static_cast<std::size_t>(k - left)];
            for (std::size_t m = 0; m <= k1; ++m) acc[m] += wk * v[m];
        }
        if (k < right) {
            for (std::size_t m = k1; m >= 1; --m) {
                const T forward = v[m - 1] * c[m - 1];
                const T stay = (m < k1) ? v[m] * a[m] : v[m];
                v[m] = stay + forward;
            }
            v[0] = v[0] * a[0];
        }
    }
    for (std::size_t m = 0; m < k1; ++m) b[m] = acc[m];
    tail = acc[k1];
}

}  // namespace detail

/**
 * @param pi0      initial distribution (row vector)
 * @param Q        generator
 * @param t        time horizon, t >= 0
 * @param epsilon  birth-process truncation tolerance (MATLAB default 1e-6)
 * @param delta    occupancy threshold below which a state is dropped (default 1e-12)
 * @param maxsteps cap on birth steps; <= 0 for the default cap
 */
template <class T>
FauResult<T> ctmc_fau(const std::vector<T>& pi0, const Matrix<T>& Q, const T& t,
                      double epsilon = 1e-6, double delta = 1e-12, long maxsteps = -1) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t n = Q.rows();
    if (Q.cols() != n) throw InputError("ctmc_fau: Q must be square");
    if (pi0.size() != n) throw InputError("ctmc_fau: pi0 and Q have inconsistent sizes");
    const double tDouble = num_traits<T>::to_double(t);
    if (tDouble < 0.0) throw InputError("ctmc_fau: t must be nonnegative");
    const double eps = (epsilon > 0.0) ? epsilon : 1e-6;
    const T dropThreshold = num_traits<T>::from_double((delta > 0.0) ? delta : 0.0);
    const long cap = (maxsteps > 0) ? maxsteps : FAU_MAX_STEPS;

    FauResult<T> r;
    r.weightTail = zero;
    r.weightWindow = zero;
    r.droppedMass = zero;
    r.errorBound = zero;

    std::vector<double> d(n, 0.0);
    std::vector<T> exitRate(n, zero);
    for (std::size_t i = 0; i < n; ++i) {
        exitRate[i] = zero - Q(i, i);
        d[i] = num_traits<T>::to_double(exitRate[i]);
        r.uniformRate = std::max(r.uniformRate, d[i]);
    }
    if (tDouble == 0.0 || n == 0) {
        r.pit = pi0;
        r.steps = 1;
        r.supportMax = detail::fau_support(pi0).size();
        r.supportFinal = r.supportMax;
        return r;
    }

    // Pass one: the adaptive rate sequence, and where it stops.
    std::vector<double> lambda;
    {
        std::vector<T> u = pi0;
        std::vector<T> scratch(n, zero);
        std::vector<char> touchedFlag(n, 0);
        std::vector<std::size_t> act = detail::fau_support(u);
        double lstar = 0.0;
        while (true) {
            if (act.empty()) {
                r.absorbed = true;
                break;
            }
            const double L = detail::fau_max_exit(d, act);
            lambda.push_back(L);
            if (L <= 0.0) {
                // Every occupied state is absorbing: the birth process stops
                // here and the remaining weight falls on this iterate.
                r.absorbed = true;
                break;
            }
            lstar = std::max(lstar, L);
            if (detail::fau_tailbound(lstar, tDouble, static_cast<long>(lambda.size())) <= eps) break;
            if (static_cast<long>(lambda.size()) >= cap) {
                r.truncated = true;
                break;
            }
            act = detail::fau_step(u, act, Q, num_traits<T>::from_double(L), dropThreshold, scratch,
                                   touchedFlag, static_cast<T*>(nullptr));
        }
    }

    // The birth-process weights of that rate sequence, exactly.
    std::vector<T> b;
    detail::fau_weights(lambda, t, tDouble, eps, b, r.weightTail, r.weightWindow);

    // Pass two: the same sweep again, accumulating sum_n b_n u^(n).
    {
        std::vector<T> u = pi0;
        std::vector<T> scratch(n, zero);
        std::vector<char> touchedFlag(n, 0);
        std::vector<std::size_t> act = detail::fau_support(u);
        r.pit.assign(n, zero);
        r.supportMax = act.size();
        r.supportFinal = act.size();
        for (std::size_t m = 0; m < lambda.size(); ++m) {
            if (act.empty()) break;
            r.supportMax = std::max(r.supportMax, act.size());
            r.supportFinal = act.size();
            for (std::size_t k = 0; k < act.size(); ++k) r.pit[act[k]] += b[m] * u[act[k]];
            if (m + 1 < lambda.size()) {
                const double L = detail::fau_max_exit(d, act);
                if (L <= 0.0) break;
                act = detail::fau_step(u, act, Q, num_traits<T>::from_double(L), dropThreshold,
                                       scratch, touchedFlag, &r.droppedMass);
            }
        }
    }

    r.steps = static_cast<long>(lambda.size());
    if (!lambda.empty()) {
        r.lambdaMin = *std::min_element(lambda.begin(), lambda.end());
        r.lambdaMax = *std::max_element(lambda.begin(), lambda.end());
    }
    T mass0 = zero;
    T massT = zero;
    for (std::size_t i = 0; i < n; ++i) {
        mass0 += pi0[i];
        massT += r.pit[i];
    }
    r.errorBound = mass0 - massT;
    return r;
}

}  // namespace mc
}  // namespace line

#endif  // LINE_API_MC_CTMC_FAU_H
