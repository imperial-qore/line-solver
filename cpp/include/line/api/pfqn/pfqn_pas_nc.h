/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PAS_NC_H
#define LINE_API_PFQN_PAS_NC_H

/**
 * Normalizing constant G_C of one communicating class of a closed
 * PASS-AND-SWAP (P&S) network, plus one aggregated delay.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_pas_nc.m.
 *
 * With a non-empty swap graph the ordered-state chain is reducible (Comte and
 * Dorsman, 2021, arXiv:2009.12299): the recurrent communicating classes are the
 * placement-order-adhering sets and the product form pi(c) = prod_m Phi_m(c_m)/G_C
 * holds per class. This routine returns that per-class constant.
 *
 * Method. Build station M's chain head-first: appending class r at chain
 * position k = |occ|+1 is admissible iff no class already placed at that
 * station must come after r, and contributes the reciprocal OI prefix rate
 * 1/mu_M(occ+e_r); the chain may be finalized (recursing to station M-1) only
 * when occ is a placement-order ideal at full multiplicity. Once every P&S
 * station is peeled the residual population sits at the delay node with the
 * multinomial weight prod_r Z_r^{N_r}/N_r!.
 *
 * This is a MICROSTATE routine: it walks the ordered chains position by
 * position, because with a placement order the reachable set is a set of
 * ORDERINGS that does not collapse onto the count lattice. With an empty order
 * the node count is sum_{b<=N} C(|b|+M-1,M-1) |b|!/prod_r b_r!, factorial in
 * the total population -- use pfqn_ncoi for the plain OI case, which returns
 * the same constant on the count lattice.
 *
 * The occupancy shift is carried as an explicit vector per station rather than
 * by wrapping the callable in a new closure at every level as MATLAB does
 * (`shifted = @(state) active(state + e_r)`). The two are the same function;
 * the vector form avoids a closure chain whose depth is the total population.
 *
 * Arithmetic: EXACT-CAPABLE, on the same terms as pfqn_ncoi.
 */

#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_ca.h"
#include "line/api/pfqn/pfqn_ncoi.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace pfqn {

/**
 * Placement order of one station: prec[i][j] != 0 iff class i must be placed
 * before class j. This is the precedence closure of pas_placement, fed by the
 * global DAG of pas_swap2order.
 */
using PlacementOrder = std::vector<std::vector<int>>;

namespace detail {

/**
 * True iff occ is a placement-order ideal at full multiplicity: for every
 * i prec j, occ[j] > 0 requires occ[i] == Norig[i]. Reduces to support
 * downward-closure when Norig is all ones, and to "always true" for an empty
 * order.
 */
inline bool pas_nc_isideal(const std::vector<int>& occ, const PlacementOrder& prec,
                           const std::vector<int>& Norig) {
    if (prec.empty()) return true;
    const std::size_t R = occ.size();
    for (std::size_t i = 0; i < R; ++i)
        for (std::size_t j = 0; j < R; ++j)
            if (prec[i][j] != 0 && occ[j] > 0 && occ[i] < Norig[i]) return false;
    return true;
}

template <class T>
T pas_nc_rec(const std::vector<T>& Z, std::vector<int>& N, const std::vector<int>& Norig,
             const std::vector<OiRate<T>>& mu, const std::vector<PlacementOrder>& prec,
             std::size_t nsta, std::vector<int>& occ) {
    const std::size_t R = N.size();
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    if (nsta == 0) {
        // Only the aggregated delay is left: the multinomial weight.
        T f = one;
        for (std::size_t r = 0; r < R; ++r) {
            if (N[r] == 0) continue;
            if (!(Z[r] > zero)) return zero;  // population with no delay demand
            f *= num_pow_int(Z[r], static_cast<unsigned>(N[r])) /
                 num_factorial<T>(static_cast<unsigned>(N[r]));
        }
        return f;
    }

    const PlacementOrder& precm = prec[nsta - 1];

    // Step A: finalize this station's chain and peel to the next one, but only
    // if the accumulated occupancy is a placement-order ideal.
    T G = zero;
    if (pas_nc_isideal(occ, precm, Norig)) {
        std::vector<int> empty(R, 0);
        G = pas_nc_rec(Z, N, Norig, mu, prec, nsta - 1, empty);
    }

    // Step B: append one more class r at the next chain position.
    for (std::size_t r = 0; r < R; ++r) {
        if (N[r] == 0) continue;
        bool blocked = false;
        if (!precm.empty()) {
            for (std::size_t j = 0; j < R; ++j) {
                if (occ[j] > 0 && precm[r][j] != 0) {
                    blocked = true;
                    break;
                }
            }
        }
        if (blocked) continue;
        occ[r] += 1;
        const T rate = mu[nsta - 1](occ);
        if (rate > zero) {
            N[r] -= 1;
            G += pas_nc_rec(Z, N, Norig, mu, prec, nsta, occ) / rate;
            N[r] += 1;
        }
        occ[r] -= 1;
    }
    return G;
}

}  // namespace detail

/**
 * @param Z    (R) think-time demand of the aggregated delay node
 * @param N    (R) closed population, finite
 * @param mu   (M) P&S rate callables, one per station; may be empty
 * @param prec (M) placement orders, one per station, each R x R with
 *             prec[m][i][j] != 0 iff class i must be placed before class j at
 *             station m. An empty vector, or an empty matrix for a station,
 *             means no order there, so G is the plain OI constant. NOTE the
 *             orientation: around a cycle each downstream station traverses its
 *             chain in the opposite direction, so downstream stations take the
 *             TRANSPOSE of the upstream order. Passing the same matrix to both
 *             stations of a cycle silently returns a smaller, wrong G.
 */
template <class T>
NcResult<T> pfqn_pas_nc(const std::vector<T>& Z, const std::vector<int>& N,
                        const std::vector<OiRate<T>>& mu,
                        const std::vector<PlacementOrder>& prec) {
    const std::size_t R = N.size();
    if (Z.size() != R) throw InputError("pfqn_pas_nc: Z and N must have the same class count");
    for (int v : N)
        if (v < 0) throw InputError("pfqn_pas_nc: requires finite, nonnegative populations");
    for (std::size_t i = 0; i < mu.size(); ++i)
        if (!mu[i]) throw InputError("pfqn_pas_nc: a P&S rate callable is empty");

    std::vector<PlacementOrder> prel;
    if (prec.empty()) {
        prel.assign(mu.size(), PlacementOrder());
    } else if (prec.size() == 1 && mu.size() > 1) {
        prel.assign(mu.size(), prec[0]);
    } else {
        if (prec.size() != mu.size())
            throw InputError("pfqn_pas_nc: prec must supply one precedence matrix per station");
        prel = prec;
    }
    for (std::size_t m = 0; m < prel.size(); ++m) {
        if (prel[m].empty()) continue;
        if (prel[m].size() != R)
            throw InputError("pfqn_pas_nc: each precedence matrix must be R x R");
        for (std::size_t i = 0; i < R; ++i)
            if (prel[m][i].size() != R)
                throw InputError("pfqn_pas_nc: each precedence matrix must be R x R");
    }

    std::vector<int> Nw(N);
    std::vector<int> occ(R, 0);
    const T G = detail::pas_nc_rec(Z, Nw, N, mu, prel, mu.size(), occ);
    return {G, num_traits<T>::log_as_double(G)};
}

/** Plain OI case (no placement order); prefer pfqn_ncoi, which is cheaper. */
template <class T>
NcResult<T> pfqn_pas_nc(const std::vector<T>& Z, const std::vector<int>& N,
                        const std::vector<OiRate<T>>& mu) {
    return pfqn_pas_nc(Z, N, mu, std::vector<PlacementOrder>());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PAS_NC_H
