/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_CA_H
#define LINE_API_PFQN_CA_H

/**
 * Convolution algorithm for the exact normalizing constant of a closed
 * product-form network (Buzen 1973, Reiser-Kobayashi 1975).
 *
 * Templated port of matlab/src/api/pfqn/pfqn_ca.m and
 * jar/src/main/java/jline/api/pfqn/nc/Pfqn_ca.java, cross-checked against
 * mp_pfqn's ca/convolution_multi_exact.c for the exact path.
 *
 * G_m(n) = G_{m-1}(n) + sum_r L(m,r) G_m(n - e_r), with G_0(n) the delay
 * balance function prod_r Z_r^{n_r}/n_r!.
 *
 * Scaling: in IEEE double the recursion overflows as soon as G(N) leaves the
 * double range, so the Lam (1982) dynamic scaling of the MATLAB implementation
 * is applied there. No other number type needs it: exact rationals have no
 * exponent range at all, and the high-precision binary floats have an exponent
 * range wide enough that G never leaves it in practice. The scaling is exact
 * either way, since dividing every demand by a power of two divides G(N) by
 * exactly that power raised to sum(N).
 */

#include <cmath>
#include <type_traits>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/util/population.h"

namespace line {
namespace pfqn {

/** Return value of the normalizing-constant family, mirroring Ret.pfqnNc. */
template <class T>
struct NcResult {
    T G;        ///< normalizing constant in the requested arithmetic
    double lG;  ///< log of the constant, always a double and always finite
};

namespace detail {

/** Delay balance function F_Z(n) = prod_r Z_r^{n_r} / n_r!. */
template <class T>
T pff_delay(const std::vector<T>& Z, const std::vector<int>& n) {
    int total = 0;
    for (int v : n) total += v;
    if (total == 0) return num_traits<T>::from_int(1);
    T f = num_traits<T>::from_int(1);
    for (std::size_t r = 0; r < n.size(); ++r) {
        if (n[r] == 0) continue;
        if (Z[r] == num_traits<T>::from_int(0)) return num_traits<T>::from_int(0);
        f *= num_pow_int(Z[r], static_cast<unsigned>(n[r])) /
             num_factorial<T>(static_cast<unsigned>(n[r]));
    }
    return f;
}

/**
 * Power-of-two scale factor centring log G near zero, from the largest state
 * term reachable by letting each class pick its own station or the delay.
 * Returns 0 for every arithmetic other than double, which needs no scaling.
 */
template <class T>
int scale_exponent(const Matrix<T>& L, const std::vector<int>& N, const std::vector<T>& Zsum) {
    if (!std::is_same<T, double>::value) return 0;
    const std::size_t M = L.rows(), R = L.cols();
    long Nt = 0;
    for (int v : N) Nt += v;
    // Each class independently takes whichever station -- or the delay -- gives it
    // its largest factor. The mixed state so named has term at least the product of
    // those factors, because a station holding several classes carries a multinomial
    // coefficient of at least one, so this is still a LOWER bound on log G. It
    // dominates the per-configuration maximum it replaces, which asked ONE station
    // (or the delay) to hold every class at once and so dropped the delay entirely
    // as soon as a single class had no think time. That collapse is what made the
    // scaling scale UP: on L=[1e-9,1], N=[99,1], Z=[1,0] the old estimate was the
    // all-at-the-queue -2051.6 against a true log G of -359.1, giving kscale=-30,
    // and Z/2^-30 = 1.07e9 overflowed the delay column Z^n/n! at n=[40,0].
    double lGest = 0.0;
    for (std::size_t r = 0; r < R; ++r) {
        if (N[r] <= 0) continue;
        double best = -std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < M; ++i) {
            double lir = num_traits<T>::to_double(L(i, r));
            if (lir > 0) best = std::max(best, N[r] * std::log(lir));
        }
        double zr = num_traits<T>::to_double(Zsum[r]);
        if (zr > 0) best = std::max(best, N[r] * std::log(zr) - std::lgamma(N[r] + 1.0));
        if (!std::isfinite(best)) {
            // no station and no delay can hold class r, so G(N) is exactly zero
            lGest = -std::numeric_limits<double>::infinity();
            break;
        }
        lGest += best;
    }
    if (!std::isfinite(lGest) || Nt == 0) return 0;
    return static_cast<int>(std::lround(lGest / (static_cast<double>(Nt) * std::log(2.0))));
}

}  // namespace detail

/**
 * @param L  (M x R) service demands, M queueing stations, R classes
 * @param N  (R) population per class
 * @param Z  (K x R) think times, summed over rows; may be empty
 */
template <class T>
NcResult<T> pfqn_ca(const Matrix<T>& L, const std::vector<int>& N, const Matrix<T>& Z) {
    const std::size_t M = L.rows();
    const std::size_t R = N.size();
    if (!L.empty() && L.cols() != R) throw InputError("pfqn_ca: L and N disagree on the class count");

    // Z summed over its rows, so a per-node think-time matrix is accepted.
    std::vector<T> Zsum(R, num_traits<T>::from_int(0));
    if (!Z.empty()) {
        if (Z.cols() != R) throw InputError("pfqn_ca: Z and N disagree on the class count");
        for (std::size_t k = 0; k < Z.rows(); ++k)
            for (std::size_t r = 0; r < R; ++r) Zsum[r] += Z(k, r);
    }

    long Nt = 0;
    bool negative = false;
    for (int v : N) {
        if (v < 0) negative = true;
        Nt += v;
    }
    if (negative) return {num_traits<T>::from_int(0), -std::numeric_limits<double>::infinity()};

    if (M == 0) {
        // Delay-only network: G = prod_r Z_r^{N_r}/N_r!.
        std::vector<int> n(N);
        T G = detail::pff_delay(Zsum, n);
        return {G, num_traits<T>::log_as_double(G)};
    }
    if (Nt == 0) return {num_traits<T>::from_int(1), 0.0};

    const int kscale = detail::scale_exponent(L, N, Zsum);
    Matrix<T> Ls = L;
    std::vector<T> Zs = Zsum;
    if constexpr (std::is_same<T, double>::value) {
        // exponent-only rescaling rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
        if (kscale != 0) {
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t r = 0; r < R; ++r) Ls(i, r) = std::ldexp(Ls(i, r), -kscale);
            for (std::size_t r = 0; r < R; ++r) Zs[r] = std::ldexp(Zs[r], -kscale);
        }
    }

    const std::vector<std::size_t> prods = plane_sizes(N);
    const std::size_t total = population_count(N);

    // G is (M+1) x total, laid out row-major with the station index outermost.
    std::vector<T> G(static_cast<std::size_t>(M + 1) * total, num_traits<T>::from_int(1));
    std::vector<int> n(R, 0);
    bool more = true;
    while (more) {
        const std::size_t idxn = pop_index(n, prods);
        G[idxn] = detail::pff_delay(Zs, n);
        for (std::size_t m = 1; m <= M; ++m) {
            T acc = G[(m - 1) * total + idxn];
            for (std::size_t r = 0; r < R; ++r)
                if (n[r] >= 1) acc += Ls(m - 1, r) * G[m * total + (idxn - prods[r])];
            G[m * total + idxn] = acc;
        }
        more = next_pop(n, N);
    }

    const T raw = G[M * total + (total - 1)];
    const double lG =
        num_traits<T>::log_as_double(raw) + static_cast<double>(Nt) * kscale * std::log(2.0);
    T Gn = raw;
    if constexpr (std::is_same<T, double>::value) {
        // exact exponent-adjustment recovery rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
        if (kscale != 0) Gn = std::ldexp(raw, static_cast<int>(static_cast<long>(Nt) * kscale));
    }
    return {Gn, lG};
}

/** Overload without think times. */
template <class T>
NcResult<T> pfqn_ca(const Matrix<T>& L, const std::vector<int>& N) {
    return pfqn_ca(L, N, Matrix<T>());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_CA_H
