/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_NCOI_H
#define LINE_API_PFQN_NCOI_H

/**
 * Normalizing constant of a closed network of ORDER-INDEPENDENT (OI) /
 * pass-and-swap stations with empty swap graph, plus one aggregated delay.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_ncoi.m.
 *
 * An OI station's total service rate mu_i(n) depends on the per-class
 * occupancy n only through which classes are present. Its balance function
 * satisfies the balanced-fairness recursion of Bonald and Proutiere (2003),
 *
 *   Phi_i(0) = 1,   Phi_i(n) = (1/mu_i(n)) sum_{r: n_r>0} Phi_i(n - e_r),
 *
 * and G(N) is the convolution of the per-station balance functions with the
 * multinomial delay factor prod_r Z_r^{n_r}/n_r!,
 *
 *   g_0(n) = F_Z(n),  g_i(n) = sum_{0<=x<=n} Phi_i(x) g_{i-1}(n-x),
 *   G(N) = g_K(N).
 *
 * This is a MACROSTATE routine: the balance functions and the convolution are
 * both tabulated over the count lattice 0 <= n <= N, never over orderings.
 * That is legitimate exactly because an OI rate is permutation-invariant, so
 * Phi(n) -- itself the sum of the ordered-prefix weights over all orderings of
 * the multiset n -- closes on the count vector. With a non-empty swap graph
 * the closure fails and pfqn_pas_nc (microstate) must be used instead.
 *
 * Cost: O(K R L) for the balance functions and O(K prod_r (N_r+1)(N_r+2)/2)
 * for the convolutions, with L = prod_r (N_r+1); that is the order of a
 * load-dependent Buzen convolution.
 *
 * Arithmetic: EXACT-CAPABLE. The routine performs only reciprocals, additions
 * and multiplications, plus the exact multinomial delay factor. Whether the
 * result is exact therefore depends only on the rate callables, which are
 * evaluated and never inspected. A nonpositive rate marks an occupancy the
 * station cannot serve: its balance value is zero, which prunes every ordering
 * through it, exactly as the reference does.
 */

#include <cstddef>
#include <functional>
#include <vector>

#include "line/api/pfqn/pfqn_ca.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** An OI station's total service rate as a function of the occupancy vector. */
template <class T>
using OiRate = std::function<T(const std::vector<int>&)>;

namespace detail {

/** Column-major lattice of the box 0 <= n <= N, with its strides. */
struct OiLattice {
    std::vector<int> dims;
    std::vector<int> strides;
    std::size_t ngrid;
    std::vector<std::vector<int>> counts;   ///< counts[k] is the point at index k
    std::vector<std::size_t> order;         ///< indices sorted by total population
};

inline OiLattice oi_lattice(const std::vector<int>& N) {
    const std::size_t R = N.size();
    OiLattice lat;
    lat.dims.resize(R);
    lat.strides.resize(R);
    std::size_t ngrid = 1;
    int total = 0;
    for (std::size_t r = 0; r < R; ++r) {
        lat.dims[r] = N[r] + 1;
        lat.strides[r] = static_cast<int>(ngrid);
        ngrid *= static_cast<std::size_t>(lat.dims[r]);
        total += N[r];
    }
    lat.ngrid = ngrid;
    lat.counts.assign(ngrid, std::vector<int>(R, 0));
    std::vector<std::vector<std::size_t>> byPop(static_cast<std::size_t>(total) + 1);
    for (std::size_t k = 0; k < ngrid; ++k) {
        std::size_t rest = k;
        int sum = 0;
        for (std::size_t r = 0; r < R; ++r) {
            lat.counts[k][r] = static_cast<int>(rest % static_cast<std::size_t>(lat.dims[r]));
            rest /= static_cast<std::size_t>(lat.dims[r]);
            sum += lat.counts[k][r];
        }
        byPop[static_cast<std::size_t>(sum)].push_back(k);
    }
    lat.order.reserve(ngrid);
    for (std::size_t s = 0; s < byPop.size(); ++s)
        for (std::size_t i = 0; i < byPop[s].size(); ++i) lat.order.push_back(byPop[s][i]);
    return lat;
}

/**
 * v-weighted balanced-fairness recursion over the count lattice, in increasing
 * population: Phi(n) = (1/mu(n)) sum_{r: n_r>0} v_r Phi(n - e_r). vis is the
 * station's class visit vector, whose geometric factor prod_r v_r^{n_r} carries
 * the OI-station visit ratio; empty means unit visits.
 */
template <class T>
std::vector<T> oi_nc_balance(const OiLattice& lat, const OiRate<T>& rate, std::size_t R,
                             const std::vector<T>& vis) {
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    std::vector<T> phi(lat.ngrid, zero);
    for (std::size_t t = 0; t < lat.order.size(); ++t) {
        const std::size_t k = lat.order[t];
        const std::vector<int>& n = lat.counts[k];
        int sum = 0;
        for (std::size_t r = 0; r < R; ++r) sum += n[r];
        if (sum == 0) {
            phi[k] = one;
            continue;
        }
        const T mun = rate(n);
        if (!(mun > zero)) continue;  // unreachable station state: zero balance
        T acc = zero;
        for (std::size_t r = 0; r < R; ++r)
            if (n[r] > 0)
                acc += (vis.empty() ? one : vis[r]) *
                       phi[k - static_cast<std::size_t>(lat.strides[r])];
        phi[k] = acc / mun;
    }
    return phi;
}

}  // namespace detail

/**
 * @param Z      (R) think-time demand of the aggregated delay node
 * @param N      (R) closed population, finite
 * @param mu     (K) OI rate callables, one per station; may be empty
 * @param visits (K x R) per-station class visit ratios weighting the balance
 *               recursion; empty for unit visits
 */
template <class T>
NcResult<T> pfqn_ncoi(const std::vector<T>& Z, const std::vector<int>& N,
                       const std::vector<OiRate<T>>& mu, const Matrix<T>& visits) {
    const std::size_t R = N.size();
    if (Z.size() != R) throw InputError("pfqn_ncoi: Z and N must have the same class count");
    if (!visits.empty() && (visits.rows() != mu.size() || visits.cols() != R))
        throw InputError("pfqn_ncoi: visits must be K x R");
    for (int v : N)
        if (v < 0) throw InputError("pfqn_ncoi: requires finite, nonnegative populations");
    for (std::size_t i = 0; i < mu.size(); ++i)
        if (!mu[i]) throw InputError("pfqn_ncoi: an OI rate callable is empty");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    if (R == 0) return {one, num_traits<T>::log_as_double(one)};

    const detail::OiLattice lat = detail::oi_lattice(N);

    // Delay balance function: the multinomial factor F_Z(n). A class with
    // population but no delay demand makes the state infeasible.
    std::vector<T> g(lat.ngrid, zero);
    for (std::size_t k = 0; k < lat.ngrid; ++k) {
        T f = one;
        bool feas = true;
        for (std::size_t r = 0; r < R; ++r) {
            const int nr = lat.counts[k][r];
            if (nr == 0) continue;
            if (!(Z[r] > zero)) {
                feas = false;
                break;
            }
            f *= num_pow_int(Z[r], static_cast<unsigned>(nr)) /
                 num_factorial<T>(static_cast<unsigned>(nr));
        }
        if (feas) g[k] = f;
    }

    // Convolve in one OI station at a time.
    std::vector<int> rem(R, 0), y(R, 0);
    for (std::size_t i = 0; i < mu.size(); ++i) {
        std::vector<T> vis;
        if (!visits.empty()) {
            vis.resize(R);
            for (std::size_t r = 0; r < R; ++r) vis[r] = visits(i, r);
        }
        const std::vector<T> phi = detail::oi_nc_balance<T>(lat, mu[i], R, vis);
        std::vector<T> gnext(lat.ngrid, zero);
        for (std::size_t kx = 0; kx < lat.ngrid; ++kx) {
            if (!(phi[kx] != zero)) continue;
            std::size_t base = 0;
            for (std::size_t r = 0; r < R; ++r) {
                rem[r] = N[r] - lat.counts[kx][r];
                y[r] = 0;
                base += static_cast<std::size_t>(lat.counts[kx][r]) *
                        static_cast<std::size_t>(lat.strides[r]);
            }
            // Odometer over the sub-box 0 <= y <= rem; lin() is linear, so the
            // target index of x + y is base + lin(y).
            for (;;) {
                std::size_t ylin = 0;
                for (std::size_t r = 0; r < R; ++r)
                    ylin += static_cast<std::size_t>(y[r]) * static_cast<std::size_t>(lat.strides[r]);
                gnext[base + ylin] += phi[kx] * g[ylin];
                std::size_t d = 0;
                while (d < R && y[d] == rem[d]) {
                    y[d] = 0;
                    ++d;
                }
                if (d == R) break;
                y[d] += 1;
            }
        }
        g.swap(gnext);
    }

    const T G = g[lat.ngrid - 1];
    return {G, num_traits<T>::log_as_double(G)};
}

/** Overload with unit visits. */
template <class T>
NcResult<T> pfqn_ncoi(const std::vector<T>& Z, const std::vector<int>& N,
                       const std::vector<OiRate<T>>& mu) {
    return pfqn_ncoi(Z, N, mu, Matrix<T>());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_NCOI_H
