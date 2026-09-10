/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MAP_SAMPLE_H
#define LINE_API_MAM_MAP_SAMPLE_H

/**
 * Sample the inter-arrival times of a MAP, a RAP or a matrix exponential.
 *
 * Templated port of matlab/lib/kpctoolbox/map/map_sample.m, rap_sample.m and
 * me_sample.m, together with m3a's randp.m.
 *
 * A MAP is simulated on its own state space: from the current phase the process
 * takes hidden transitions (D0's off-diagonal) until one of the arrival
 * transitions (D1) fires, and the inter-arrival time is the sum of the
 * exponential holding times spent on the way. The reference draws the whole path
 * and then adds the holding times; this port adds them as it goes, which is the
 * same sum with no path buffer and no growing `visits` matrix.
 *
 * A RAP IS NOT SIMULATED THAT WAY. Its D0 has negative off-diagonals, so there
 * is no embedded jump chain to walk: the "phase" is a signed vector, not a
 * state. The reference samples it by INVERSE TRANSFORM on the conditional
 * distribution, advancing the row vector
 *   a <- a exp(D0 x) D1 / (a exp(D0 x) D1 e)
 * after each arrival at x. That is what `rap_sample` does here, and it is why it
 * costs a matrix exponential per sample where `map_sample` costs a handful of
 * deviates. An ME is a RENEWAL process, so it does not pay that price: it lives
 * in me_sample.h, which tabulates its CDF once.
 *
 * RANDOMNESS. The generator is `line::pfqn::McRng` passed by reference and
 * advanced by the call, the convention every Monte Carlo entry point in this
 * tree uses. The stream is NOT comparable with MATLAB's -- different generator,
 * different mapping from bits to deviates -- so the oracle for these functions
 * is distributional, never path-for-path.
 *
 * ARITHMETIC: transcendental. Exponential deviates and, for the RAP and ME
 * paths, a matrix exponential.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/api/pfqn/pfqn_mc_common.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/**
 * Draw an index from a discrete law, m3a's `randp`.
 *
 * The weights need not be normalized; they must be non-negative and not all
 * zero, which the reference reports rather than assumes.
 */
template <class T>
std::size_t randp(const std::vector<T>& P, pfqn::McRng& rng) {
    const T zero = num_traits<T>::from_int(0);
    T total = zero;
    for (std::size_t i = 0; i < P.size(); ++i) {
        if (P[i] < zero) throw InputError("randp: all probabilities should be 0 or larger");
        total += P[i];
    }
    if (P.empty() || !(total > zero)) throw InputError("randp: all zero probabilities");
    const double u = pfqn::mc_uniform01(rng);
    T acc = zero;
    for (std::size_t i = 0; i < P.size(); ++i) {
        acc += P[i];
        if (num_traits<T>::to_double(acc / total) >= u) return i;
    }
    return P.size() - 1;
}

/** The state a sampled inter-arrival began and ended in. */
struct SampleTrace {
    std::vector<std::size_t> first;  ///< phase at the start of each interval, 0-based
    std::vector<std::size_t> last;   ///< phase entered on each arrival, 0-based
};

/**
 * @param m        the MAP
 * @param n        number of inter-arrival times to draw
 * @param rng      generator, advanced by the call
 * @param pie0     initial phase law; empty selects map_pie, the reference's
 *                 interval-stationary initialization
 * @param trace    optional per-sample start and end phases
 */
template <class T>
std::vector<T> map_sample(const Map<T>& m, std::size_t n, pfqn::McRng& rng,
                          const std::vector<T>& pie0 = std::vector<T>(),
                          SampleTrace* trace = 0) {
    static_assert(num_traits<T>::has_transcendental,
                  "map_sample draws exponential deviates");
    const std::size_t K = m.order();
    if (K == 0 || m.D1.rows() != K) throw InputError("map_sample: D0 and D1 disagree");
    const T zero = num_traits<T>::from_int(0);

    std::vector<T> out;
    out.reserve(n);
    if (trace != 0) {
        trace->first.assign(n, 0);
        trace->last.assign(n, 0);
    }

    // The exponential case has no phase to walk.
    if (K == 1) {
        const T mean = map_mean(m);
        for (std::size_t i = 0; i < n; ++i)
            out.push_back(T(mean * num_traits<T>::from_double(-std::log(pfqn::mc_uniform01(rng)))));
        return out;
    }

    const std::vector<T> start = pie0.empty() ? map_pie(m) : pie0;
    if (start.size() != K) throw InputError("map_sample: the initial law has the wrong length");
    std::size_t cur = randp(start, rng);

    // Row i of the jump law over the 2K destinations: hidden moves first, then
    // the arrival moves, each divided by the total rate out of i.
    std::vector<std::vector<T>> jump(K, std::vector<T>(2 * K, zero));
    std::vector<T> hold(K, zero);
    for (std::size_t i = 0; i < K; ++i) {
        const T rate = -m.D0(i, i);
        if (!(num_traits<T>::to_double(rate) > 0.0))
            throw InputError("map_sample: a phase has no exit rate");
        hold[i] = T(num_traits<T>::from_int(1) / rate);
        for (std::size_t j = 0; j < K; ++j) {
            if (i != j) jump[i][j] = T(m.D0(i, j) / rate);
            jump[i][K + j] = T(m.D1(i, j) / rate);
        }
    }

    for (std::size_t s = 0; s < n; ++s) {
        if (trace != 0) trace->first[s] = cur;
        T acc = zero;
        for (;;) {
            // One exponential holding time in the current phase.
            acc += T(hold[cur] * num_traits<T>::from_double(-std::log(pfqn::mc_uniform01(rng))));
            const std::size_t d = randp(jump[cur], rng);
            if (d >= K) {  // an arrival: the interval ends here
                cur = d - K;
                break;
            }
            cur = d;
        }
        if (trace != 0) trace->last[s] = cur;
        out.push_back(acc);
    }
    return out;
}

namespace sampledetail {

/** exp(A t) by scaling and squaring around a truncated Taylor series. */
template <class T>
Matrix<T> expm(const Matrix<T>& A, const T& t) {
    const std::size_t n = A.rows();
    double nrm = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        double r = 0.0;
        for (std::size_t j = 0; j < n; ++j)
            r += std::fabs(num_traits<T>::to_double(A(i, j)) * num_traits<T>::to_double(t));
        nrm = std::max(nrm, r);
    }
    int s = 0;
    while (nrm > 0.5) {
        nrm /= 2.0;
        ++s;
    }
    const T h = T(t / num_traits<T>::from_double(std::pow(2.0, s)));
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    Matrix<T> R(n, n, zero), term(n, n, zero);
    for (std::size_t i = 0; i < n; ++i) {
        R(i, i) = one;
        term(i, i) = one;
    }
    for (int k = 1; k <= 40; ++k) {
        Matrix<T> nx(n, n, zero);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) {
                T v = zero;
                for (std::size_t q = 0; q < n; ++q) v += term(i, q) * A(q, j) * h;
                nx(i, j) = T(v / num_traits<T>::from_int(k));
            }
        term = nx;
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) R(i, j) += term(i, j);
    }
    for (int k = 0; k < s; ++k) {
        Matrix<T> sq(n, n, zero);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) {
                T v = zero;
                for (std::size_t q = 0; q < n; ++q) v += R(i, q) * R(q, j);
                sq(i, j) = v;
            }
        R = sq;
    }
    return R;
}

/** P(X > x) = a exp(D0 x) e, the conditional survival from entry law a. */
template <class T>
T survival(const Matrix<T>& D0, const std::vector<T>& a, const T& x) {
    // QUALIFIED deliberately. Unqualified, ordinary lookup finds this
    // namespace's `expm` while ADL on `line::Matrix<T>` also drags in
    // `line::expm` from util/expm.h, and the call is AMBIGUOUS in any
    // translation unit that includes both headers. Nothing did until the
    // native LDES engine, which is why it compiled for so long.
    const Matrix<T> E = sampledetail::expm(D0, x);
    T s = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t j = 0; j < a.size(); ++j) s += a[i] * E(i, j);
    return s;
}

}  // namespace sampledetail

/**
 * Sample a RAP or a matrix exponential by inverse transform.
 *
 * There is no embedded jump chain to walk -- D0 may carry negative
 * off-diagonals -- so each sample is the root of a exp(D0 x) e = u, found by
 * bracketing and bisection, after which the entry law is advanced to
 * a exp(D0 x) D1, renormalized.
 *
 * @param m   the RAP or ME as a (D0, D1) pair
 * @param n   number of inter-arrival times to draw
 * @param rng generator, advanced by the call
 * @param a0  initial entry law; empty selects map_pie
 * @param a_out when non-null, receives the entry law AFTER the last sample, so
 *              a caller drawing one variate at a time can chain the calls and
 *              keep the process correlated. A RAP's state is this real-valued
 *              vector and not a discrete phase, so it cannot be recovered from
 *              the sampled times the way a MAP's can from `SampleTrace`:
 *              without it, repeated n=1 calls silently restart the process
 *              from its stationary entry law every time and deliver a RENEWAL
 *              stream with the right marginal and no autocorrelation.
 */
template <class T>
std::vector<T> rap_sample(const Map<T>& m, std::size_t n, pfqn::McRng& rng,
                          const std::vector<T>& a0 = std::vector<T>(),
                          std::vector<T>* a_out = 0) {
    static_assert(num_traits<T>::has_transcendental,
                  "rap_sample inverts a matrix-exponential survival function");
    const std::size_t K = m.order();
    if (K == 0 || m.D1.rows() != K) throw InputError("rap_sample: D0 and D1 disagree");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);

    std::vector<T> a = a0.empty() ? map_pie(m) : a0;
    if (a.size() != K) throw InputError("rap_sample: the initial law has the wrong length");
    const double mean = num_traits<T>::to_double(map_mean(m));
    if (!(mean > 0.0)) throw InputError("rap_sample: the process has no positive mean");

    std::vector<T> out;
    out.reserve(n);
    for (std::size_t s = 0; s < n; ++s) {
        const double u = pfqn::mc_uniform01(rng);
        // Bracket the root of survival(x) = u, then bisect.
        double lo = 0.0, hi = mean;
        for (int k = 0; k < 200; ++k) {
            if (num_traits<T>::to_double(
                    sampledetail::survival(m.D0, a, num_traits<T>::from_double(hi))) <= u)
                break;
            lo = hi;
            hi *= 2.0;
        }
        for (int k = 0; k < 200; ++k) {
            const double mid = 0.5 * (lo + hi);
            const double sv = num_traits<T>::to_double(
                sampledetail::survival(m.D0, a, num_traits<T>::from_double(mid)));
            if (sv > u)
                lo = mid;
            else
                hi = mid;
            if (hi - lo < 1e-14 * (1.0 + hi)) break;
        }
        const T x = num_traits<T>::from_double(0.5 * (lo + hi));
        out.push_back(x);

        // Advance the entry law: a <- a exp(D0 x) D1, renormalized.
        const Matrix<T> E = sampledetail::expm(m.D0, x);
        std::vector<T> b(K, zero), c(K, zero);
        for (std::size_t j = 0; j < K; ++j)
            for (std::size_t i = 0; i < K; ++i) b[j] += a[i] * E(i, j);
        T tot = zero;
        for (std::size_t j = 0; j < K; ++j) {
            for (std::size_t i = 0; i < K; ++i) c[j] += b[i] * m.D1(i, j);
            tot += c[j];
        }
        if (!(num_traits<T>::to_double(tot) > 0.0))
            throw NumericError("rap_sample: the entry law lost all its mass");
        for (std::size_t j = 0; j < K; ++j) a[j] = T(c[j] / tot);
    }
    (void)one;
    if (a_out != 0) *a_out = a;
    return out;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MAP_SAMPLE_H
