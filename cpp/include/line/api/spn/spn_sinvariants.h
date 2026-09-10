/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SPN_SPN_SINVARIANTS_H
#define LINE_API_SPN_SPN_SINVARIANTS_H

/**
 * Minimal-support S-invariants (P-invariants) of a stochastic Petri net, and
 * the load vector V = S m0.
 *
 * An S-invariant is a non-negative left null vector of the incidence matrix,
 * U' C = 0, so U' m is conserved by every firing. The minimal-support ones form
 * a basis of all of them (S. Balsamo, A. Marin, I. Stojic, FGCS 111 (2020),
 * Sec. 3.1) and are what the convolution algorithm `spn_conv` decomposes the
 * reachability set along; `spn_mdd` uses a single positive invariant for a much
 * weaker purpose, to bound each place a priori.
 *
 * FARKAS' ALGORITHM, on [C | I]: for each transition column in turn, keep the
 * rows that already annihilate it and add, for every pair of rows of opposite
 * sign in it, the positive combination that cancels it; then drop every row
 * whose support strictly contains another's, which is what leaves the minimal
 * supports. Rows are kept in integer arithmetic and divided by their gcd, so a
 * multiplicity is never lost to rounding and two invariants that differ only by
 * a positive scale are the same row.
 *
 * ARC MULTIPLICITIES MUST BE INTEGRAL. A fractional arc has no Petri-net
 * meaning and would make the gcd normalisation and the ILP-free convolution
 * both wrong, so it is refused rather than rounded.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <string>
#include <vector>

#include "line/lang/lang_types.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace spn {

/** The invariant basis of a net, in place-level coordinates. */
struct SpnInvariants {
    /** 1-based node indices of the places, in level order. */
    std::vector<std::size_t> places;
    /** S[i][p]: weight of place p in minimal-support invariant i. */
    std::vector<std::vector<long long>> S;
    /** V = S m0, the load vector. */
    std::vector<long long> V;
    /** The initial marking the load vector was taken against. */
    std::vector<long long> m0;
};

namespace detail {

inline long long spn_gcd(long long a, long long b) {
    a = a < 0 ? -a : a;
    b = b < 0 ? -b : b;
    while (b != 0) {
        const long long t = a % b;
        a = b;
        b = t;
    }
    return a;
}

/** An arc multiplicity, refused unless integral. */
inline long long spn_as_integer(double x, const char* what) {
    const double r = std::floor(x + 0.5);
    if (std::fabs(x - r) > 1e-9)
        throw InputError(std::string("spn_sinvariants: ") + what +
                         " is not integral; a fractional arc multiplicity has no Petri-net "
                         "meaning and no invariant basis over the integers");
    return static_cast<long long>(r);
}

/** True when the support of a is contained in the support of b. */
inline bool spn_support_subset(const std::vector<long long>& a, const std::vector<long long>& b,
                              std::size_t off, std::size_t n) {
    for (std::size_t k = 0; k < n; ++k)
        if (a[off + k] != 0 && b[off + k] == 0) return false;
    return true;
}

inline bool spn_support_equal(const std::vector<long long>& a, const std::vector<long long>& b,
                             std::size_t off, std::size_t n) {
    for (std::size_t k = 0; k < n; ++k)
        if ((a[off + k] != 0) != (b[off + k] != 0)) return false;
    return true;
}

}  // namespace detail

/**
 * Minimal-support S-invariants and the load vector of a net.
 *
 * @param sn a NetworkStruct holding Places and Transitions
 * @param init initial marking per place level; empty takes it from the
 *        reference station of each closed class, as `spn_mdd` does
 */
template <class T>
SpnInvariants spn_sinvariants(const qn::NetworkStruct<T>& sn,
                              const std::vector<double>& init = std::vector<double>()) {
    std::vector<std::size_t> places, transitions;
    for (std::size_t i = 1; i <= sn.nodes.size(); ++i) {
        if (sn.nodes[i - 1].nodetype == lang::NodeType::Place) places.push_back(i);
        else if (sn.nodes[i - 1].nodetype == lang::NodeType::Transition) transitions.push_back(i);
    }
    if (places.empty() || transitions.empty())
        throw InputError("spn_sinvariants: the model holds no Place or no Transition node");
    const std::size_t n = places.size();

    // ---- incidence matrix C[p][mode] = post - pre, one column per (transition, mode)
    std::vector<std::vector<long long>> C(n);
    std::size_t ncols = 0;
    for (std::size_t p = 0; p < n; ++p) C[p].clear();
    for (std::size_t t = 0; t < transitions.size(); ++t) {
        const typename std::map<std::size_t, qn::TransitionParam<T>>::const_iterator it =
            sn.transparam.find(transitions[t]);
        if (it == sn.transparam.end()) continue;
        const qn::TransitionParam<T>& tp = it->second;
        for (std::size_t m = 0; m < tp.nmodes; ++m) {
            // THE INCIDENCE IS OVER PLACES, so the arcs are summed over classes:
            // an S-invariant of the class-summed net is a genuine invariant of
            // the coloured one (every colour moves along the same arc), it is
            // just not the finest one -- the per-(place, class) invariants
            // refine it. Stated here rather than refused, because a weaker
            // invariant is still an invariant.
            const std::vector<T> en_t = qn::TransitionParam<T>::arc_total(tp.enabling, m);
            const std::vector<T> fi_t = qn::TransitionParam<T>::arc_total(tp.firing, m);
            for (std::size_t p = 0; p < n; ++p) {
                const std::size_t q = places[p] - 1;
                double pre = 0, post = 0;
                if (q < en_t.size()) pre = num_traits<T>::to_double(en_t[q]);
                if (q < fi_t.size()) post = num_traits<T>::to_double(fi_t[q]);
                C[p].push_back(detail::spn_as_integer(post, "a firing arc") -
                               detail::spn_as_integer(pre, "an enabling arc"));
            }
            ++ncols;
        }
    }

    // ---- Farkas on [C | I]: row p starts as (C[p], e_p)
    std::vector<std::vector<long long>> rows(n);
    for (std::size_t p = 0; p < n; ++p) {
        rows[p].assign(ncols + n, 0);
        for (std::size_t c = 0; c < ncols; ++c) rows[p][c] = C[p][c];
        rows[p][ncols + p] = 1;
    }
    for (std::size_t c = 0; c < ncols; ++c) {
        std::vector<std::vector<long long>> next;
        for (std::size_t r = 0; r < rows.size(); ++r)
            if (rows[r][c] == 0) next.push_back(rows[r]);
        for (std::size_t a = 0; a < rows.size(); ++a) {
            if (rows[a][c] <= 0) continue;
            for (std::size_t b = 0; b < rows.size(); ++b) {
                if (rows[b][c] >= 0) continue;
                const long long pa = rows[a][c], nb = -rows[b][c];
                const long long d = detail::spn_gcd(pa, nb);
                const long long fa = nb / d, fb = pa / d;
                std::vector<long long> combo(ncols + n, 0);
                long long g = 0;
                for (std::size_t k = 0; k < combo.size(); ++k) {
                    combo[k] = fa * rows[a][k] + fb * rows[b][k];
                    g = detail::spn_gcd(g, combo[k]);
                }
                if (g > 1)
                    for (std::size_t k = 0; k < combo.size(); ++k) combo[k] /= g;
                bool nonzero = false;
                for (std::size_t k = 0; k < n; ++k) nonzero = nonzero || combo[ncols + k] != 0;
                if (nonzero) next.push_back(combo);
            }
        }
        // support-minimality filter, applied at every step so the row set cannot
        // grow combinatorially on the way to the answer
        std::vector<std::vector<long long>> keep;
        for (std::size_t r = 0; r < next.size(); ++r) {
            bool dominated = false;
            for (std::size_t s = 0; s < next.size() && !dominated; ++s) {
                if (s == r) continue;
                if (!detail::spn_support_subset(next[s], next[r], ncols, n)) continue;
                const bool same = detail::spn_support_equal(next[s], next[r], ncols, n);
                if (!same || s < r) dominated = true;  // keep the first of equal supports
            }
            if (!dominated) keep.push_back(next[r]);
        }
        rows = keep;
    }

    SpnInvariants out;
    out.places = places;
    for (std::size_t r = 0; r < rows.size(); ++r) {
        bool nonneg = true;
        for (std::size_t k = 0; k < n; ++k) nonneg = nonneg && rows[r][ncols + k] >= 0;
        if (!nonneg) continue;  // an S-invariant is non-negative by definition
        std::vector<long long> y(n, 0);
        for (std::size_t k = 0; k < n; ++k) y[k] = rows[r][ncols + k];
        out.S.push_back(y);
    }

    // ---- initial marking and the load vector V = S m0
    out.m0.assign(n, 0);
    if (!init.empty()) {
        if (init.size() != n)
            throw InputError("spn_sinvariants: init must hold one token count per place");
        for (std::size_t p = 0; p < n; ++p)
            out.m0[p] = detail::spn_as_integer(init[p], "an initial marking");
    } else {
        for (std::size_t r = 0; r < sn.classes.size(); ++r) {
            const double njobs = sn.classes[r].population;
            if (!std::isfinite(njobs))
                throw UnsupportedError("spn_sinvariants: class " + std::to_string(r + 1) +
                                       " is open, so the net has no finite load vector");
            const std::size_t ref_node = sn.station_to_node[sn.classes[r].refstat - 1];
            for (std::size_t p = 0; p < n; ++p)
                if (places[p] == ref_node)
                    out.m0[p] += detail::spn_as_integer(njobs, "a class population");
        }
    }
    out.V.assign(out.S.size(), 0);
    for (std::size_t i = 0; i < out.S.size(); ++i)
        for (std::size_t p = 0; p < n; ++p) out.V[i] += out.S[i][p] * out.m0[p];
    return out;
}

}  // namespace spn
}  // namespace line

#endif  // LINE_API_SPN_SPN_SINVARIANTS_H
