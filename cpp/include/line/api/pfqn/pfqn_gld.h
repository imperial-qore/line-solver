/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_GLD_H
#define LINE_API_PFQN_GLD_H

/**
 * Exact normalizing constant of a closed product-form network whose stations
 * may be load dependent (generalized Buzen, Reiser-Kobayashi 1975).
 *
 * Templated port of matlab/src/api/pfqn/pfqn_gld.m (and of the single-class
 * specialization matlab/src/api/pfqn/pfqn_gldsingle.m, which the recursion
 * below subsumes), cross-checked term for term against mp_pfqn's
 * gld/gld_multi.c, the exact GMP reference.
 *
 * Model. Station i serves at rate mu(i,k) when it holds k jobs, k = 1 ... Nt
 * with Nt = sum_r N_r. The single-station balance function for a class vector
 * k with j = sum_r k_r jobs is
 *
 *   Y_i(k) = j! / prod_r k_r! * prod_r L(i,r)^{k_r} / prod_{a=1}^{j} mu(i,a)
 *
 * and the constant is the convolution of the M station factors,
 *
 *   G_0(n) = [n == 0],   G_i(n) = sum_{0 <= k <= n} Y_i(k) G_{i-1}(n - k),
 *   G(N)   = G_M(N).
 *
 * Algorithm. MATLAB writes this as a recursion on (M, N, mu) with no
 * memoization, whose cost is exponential in Nt; mp_pfqn replaced it by the
 * station-by-station convolution above, which is what this port implements.
 * The value is identical, the cost is O(M P^2 R) with P = prod_r (N_r + 1).
 * A station whose rates are all 1 is load independent and is folded in by the
 * classical in-place Buzen update
 *
 *   G_i(n) = G_{i-1}(n) + sum_r L(i,r) G_i(n - e_r)
 *
 * in O(P R) instead, so a model with no load-dependent station reduces
 * operation for operation to pfqn_ca on the same demands.
 *
 * NO pfqn_lld HERE, deliberately. MATLAB, python and the JAR carry a pfqn_lld
 * alongside pfqn_gld: there the recursion is the unmemoised one described
 * above, and saturating the rate shift at the LLD threshold s_k makes its
 * state repeat, so a memo turns an exponential tree into a bounded one. That
 * is a cure for an algorithm this port does not use. The convolution here is
 * already polynomial and visits each station once, so there is nothing for the
 * threshold to collapse; the LLD structure would have to be exploited by a
 * different device, splitting a station's balance function into the s_k terms
 * below the threshold and a geometric tail folded in by the Buzen update, and
 * that is a change of algorithm rather than a port. pfqn_lldsingle IS ported,
 * because the single-class kernel it accelerates is the same recursion in
 * every language.
 *
 * Delay stations. MATLAB's pfqn_gld takes no think-time argument: a delay is
 * an ordinary row of L whose rates are mu(i,k) = k, for which the factorials
 * cancel and Y_i(k) collapses to prod_r Z_r^{k_r} / prod_r k_r!. The port
 * keeps that convention, so an infinite-server station is expressed by giving
 * it the rate row 1, 2, ..., Nt and needs no special case anywhere.
 *
 * Arithmetic. Every operation is an addition, a multiplication or a division
 * in the field of the inputs, so the algorithm is exact in rational arithmetic
 * with no reformulation; nothing here needs a transcendental function. As in
 * pfqn_ca, IEEE double is the only arithmetic that can overflow, and it gets
 * the same power-of-two rescaling of the demands: dividing every demand by
 * 2^k divides every Y_i(k) of total degree j by 2^{jk}, hence divides G(N) by
 * exactly 2^{Nt k}, and ldexp moves the exponent without touching a mantissa
 * bit. The estimate that picks k accounts for the rates as well as the
 * demands, since a delay station depresses G by Nt!.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <type_traits>
#include <vector>

#include "line/api/pfqn/pfqn_ca.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/util/population.h"

namespace line {
namespace pfqn {

namespace detail {

/**
 * Power-of-two scale exponent for the load-dependent recursion, from the
 * largest single-station term log Y_i(N). Returns 0 for every arithmetic other
 * than double, and also whenever a demand or a rate is non-positive, since the
 * estimate is then unavailable and the recursion is run unscaled.
 */
template <class T>
int gld_scale_exponent(const Matrix<T>& L, const std::vector<int>& N, const Matrix<T>& mu,
                       long Nt) {
    if (!std::is_same<T, double>::value) return 0;
    if (Nt <= 0) return 0;
    const std::size_t M = L.rows(), R = L.cols();

    double lmulti = std::lgamma(static_cast<double>(Nt) + 1.0);
    for (std::size_t r = 0; r < R; ++r) lmulti -= std::lgamma(static_cast<double>(N[r]) + 1.0);

    double lGest = -std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < M; ++i) {
        double t = lmulti;
        bool ok = true;
        for (std::size_t r = 0; r < R && ok; ++r) {
            if (N[r] > 0) {
                const double lir = num_traits<T>::to_double(L(i, r));
                if (lir > 0)
                    t += N[r] * std::log(lir);
                else
                    ok = false;
            }
        }
        for (long a = 0; a < Nt && ok; ++a) {
            const double m = num_traits<T>::to_double(mu(i, static_cast<std::size_t>(a)));
            if (m > 0)
                t -= std::log(m);
            else
                ok = false;
        }
        if (ok && t > lGest) lGest = t;
    }
    if (!std::isfinite(lGest)) return 0;
    return static_cast<int>(std::lround(lGest / (static_cast<double>(Nt) * std::log(2.0))));
}

}  // namespace detail

/**
 * @param L  (M x R) service demands, M stations and R closed classes
 * @param N  (R) population per class
 * @param mu (M x Nt') load-dependent service rates, Nt' >= sum(N); mu(i,k-1)
 *           is the rate of station i while it holds k jobs. A row of all ones
 *           is a single server, the row 1, 2, ..., Nt is an infinite server.
 *
 * @throws InputError on a dimension mismatch, on a rate matrix with fewer
 *         columns than the total population, or on a zero rate (which would
 *         make the balance function undefined rather than infinite).
 */
template <class T>
NcResult<T> pfqn_gld(const Matrix<T>& L, const std::vector<int>& N, const Matrix<T>& mu) {
    const std::size_t M = L.rows();
    const std::size_t R = N.size();
    if (!L.empty() && L.cols() != R)
        throw InputError("pfqn_gld: demand matrix and population vector disagree on the class count");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    long Nt = 0;
    bool negative = false;
    for (int v : N) {
        if (v < 0) negative = true;
        Nt += v;
    }
    // Same contract as pfqn_ca: an unreachable population has no states.
    if (negative) return {zero, -std::numeric_limits<double>::infinity()};
    if (Nt == 0) return {one, 0.0};
    // MATLAB returns G = 0 for an empty demand matrix and a positive population.
    if (M == 0) return {zero, -std::numeric_limits<double>::infinity()};

    if (mu.rows() != M)
        throw InputError("pfqn_gld: rate matrix and demand matrix disagree on the station count");
    if (static_cast<long>(mu.cols()) < Nt)
        throw InputError("pfqn_gld: rate matrix needs one column per job in the total population");
    for (std::size_t i = 0; i < M; ++i)
        for (long a = 0; a < Nt; ++a)
            if (mu(i, static_cast<std::size_t>(a)) == zero)
                throw InputError("pfqn_gld: load-dependent service rate must be nonzero");

    const int kscale = detail::gld_scale_exponent(L, N, mu, Nt);
    Matrix<T> Ls = L;
    if constexpr (std::is_same<T, double>::value) {
        if (kscale != 0) {
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t r = 0; r < R; ++r) Ls(i, r) = std::ldexp(Ls(i, r), -kscale);
        }
    }

    const std::vector<std::size_t> prods = plane_sizes(N);
    const std::size_t total = population_count(N);

    // Running constant G_i over the population lattice, starting at G_0 = delta_0.
    std::vector<T> Gv(total, zero);
    Gv[0] = one;

    std::vector<T> Y, Gs;
    std::vector<T> muprod;
    std::vector<int> n(R, 0), k(R, 0);

    for (std::size_t m = 0; m < M; ++m) {
        bool loadIndependent = true;
        for (long a = 0; a < Nt && loadIndependent; ++a)
            if (!(mu(m, static_cast<std::size_t>(a)) == one)) loadIndependent = false;

        if (loadIndependent) {
            // In-place Buzen, ascending: n - e_r is already updated for this
            // station by the time n is reached, which is what the recursion asks.
            std::fill(n.begin(), n.end(), 0);
            while (next_pop(n, N)) {
                const std::size_t idx = pop_index(n, prods);
                T acc = Gv[idx];
                for (std::size_t r = 0; r < R; ++r)
                    if (n[r] > 0) acc += Ls(m, r) * Gv[idx - prods[r]];
                Gv[idx] = acc;
            }
            continue;
        }

        // Load-dependent station: build its balance function over the lattice,
        // then convolve. muprod[j] = prod_{a=1}^{j} mu(m,a).
        muprod.assign(static_cast<std::size_t>(Nt) + 1, one);
        for (long j = 1; j <= Nt; ++j)
            muprod[static_cast<std::size_t>(j)] =
                muprod[static_cast<std::size_t>(j - 1)] * mu(m, static_cast<std::size_t>(j - 1));

        Y.assign(total, zero);
        std::fill(k.begin(), k.end(), 0);
        bool more = true;
        while (more) {
            long j = 0;
            for (int v : k) j += v;
            const std::size_t ik = pop_index(k, prods);

            bool vanishes = false;
            T num = num_factorial<T>(static_cast<unsigned>(j));
            for (std::size_t r = 0; r < R && !vanishes; ++r) {
                if (k[r] > 0) {
                    if (Ls(m, r) == zero)
                        vanishes = true;
                    else
                        num *= num_pow_int(Ls(m, r), static_cast<unsigned>(k[r]));
                }
            }
            if (!vanishes) {
                T den = one;
                for (std::size_t r = 0; r < R; ++r) den *= num_factorial<T>(static_cast<unsigned>(k[r]));
                Y[ik] = num / den / muprod[static_cast<std::size_t>(j)];
            }
            more = next_pop(k, N);
        }

        Gs.assign(total, zero);
        std::fill(n.begin(), n.end(), 0);
        more = true;
        while (more) {
            const std::size_t idx = pop_index(n, prods);
            T acc = zero;
            std::fill(k.begin(), k.end(), 0);
            while (true) {
                std::size_t ik = 0, idiff = 0;
                for (std::size_t r = 0; r < R; ++r) {
                    ik += prods[r] * static_cast<std::size_t>(k[r]);
                    idiff += prods[r] * static_cast<std::size_t>(n[r] - k[r]);
                }
                if (!(Y[ik] == zero) && !(Gv[idiff] == zero)) acc += Y[ik] * Gv[idiff];
                bool carry = true;
                for (std::size_t r = 0; r < R; ++r) {
                    if (k[r] < n[r]) {
                        ++k[r];
                        carry = false;
                        break;
                    }
                    k[r] = 0;
                }
                if (carry) break;
            }
            Gs[idx] = acc;
            more = next_pop(n, N);
        }
        Gv.swap(Gs);
    }

    const T raw = Gv[total - 1];
    const double lG =
        num_traits<T>::log_as_double(raw) + static_cast<double>(Nt) * kscale * std::log(2.0);
    T Gn = raw;
    if constexpr (std::is_same<T, double>::value) {
        if (kscale != 0) Gn = std::ldexp(raw, static_cast<int>(Nt * static_cast<long>(kscale)));
    }
    return {Gn, lG};
}

/** Overload with all rates equal to one, i.e. every station a single server. */
template <class T>
NcResult<T> pfqn_gld(const Matrix<T>& L, const std::vector<int>& N) {
    long Nt = 0;
    for (int v : N)
        if (v > 0) Nt += v;
    Matrix<T> mu(L.rows(), static_cast<std::size_t>(Nt > 0 ? Nt : 1), num_traits<T>::from_int(1));
    return pfqn_gld(L, N, mu);
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_GLD_H
