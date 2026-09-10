/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_LQN_LSN_MAX_MULTIPLICITY_H
#define LINE_API_LQN_LSN_MAX_MULTIPLICITY_H

/**
 * Maximum sustainable multiplicity (concurrency level) of every element of a
 * layered software network.
 *
 * Templated port of matlab/src/api/lsn/lsn_max_multiplicity.m, cross-checked
 * against jar/src/main/java/jline/api/lsn/LsnMaxMultiplicity.java. It lives
 * under api/lqn/ because the port adds no api/lsn/ directory of its own; the
 * namespace is line::lsn, matching the MATLAB domain.
 *
 * Concurrency is propagated along the call graph in topological order (Kahn):
 * a reference task seeds its own multiplicity, an entry with open arrivals
 * seeds one thread, and every element passes on min(what reaches it, what it
 * can hold). A setup task is exempt from the caller bound: its instances
 * are provisioned by the platform rather than spawned by its callers, so it
 * passes on its declared multiplicity. A non-reference task with infinite
 * multiplicity ends up unbounded.
 *
 * SCOPE: the port takes the six plain-data fields the algorithm actually
 * reads -- the call graph, the multiplicities, the element types, the
 * reference and setup-task flags, and the per-entry open-arrival flag --
 * rather than a LayeredNetworkStruct, so no model layer is needed.
 *
 * DIVERGENCE, MATLAB vs JAR: MATLAB also seeds inflow(i) = 1 for an ENTRY that
 * has an open arrival (lsn.arrival{i} non-empty); the JAR omits that branch
 * entirely, so an open-arrival entry reachable from no reference task gets
 * outflow 0 there and 1 in MATLAB. This port follows MATLAB, the reference
 * implementation, and exposes the flag as entry_has_arrival.
 *
 * ARITHMETIC: comparisons and additions only, so a finite field computation,
 * exact in the exact instantiation. Infinite multiplicity is carried by an
 * explicit flag instead of a floating infinity, both because Rational has no
 * infinity and because Inf + Inf and min(Inf, Inf) are then decided by the
 * algorithm rather than by the number type.
 */

#include <cstddef>
#include <deque>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace lsn {

/** Element kinds, with the values of MATLAB's LayeredNetworkElement. */
enum class LsnElementType { HOST = 0, TASK = 1, ENTRY = 2, ACTIVITY = 3, CALL = 4 };

/** A multiplicity, possibly infinite; MATLAB's mult(i) = Inf. */
template <class T>
struct Multiplicity {
    bool infinite = false;
    T value = num_traits<T>::from_int(0);

    static Multiplicity finite(const T& v) {
        Multiplicity r;
        r.infinite = false;
        r.value = v;
        return r;
    }
    static Multiplicity of(long v) { return finite(num_traits<T>::from_int(v)); }
    static Multiplicity inf() {
        Multiplicity r;
        r.infinite = true;
        return r;
    }

    bool positive() const { return infinite || value > num_traits<T>::from_int(0); }
};

/** a + b, with infinity absorbing. */
template <class T>
Multiplicity<T> mult_add(const Multiplicity<T>& a, const Multiplicity<T>& b) {
    if (a.infinite || b.infinite) return Multiplicity<T>::inf();
    return Multiplicity<T>::finite(T(a.value + b.value));
}

/** min(a,b), with infinity as the top element. */
template <class T>
Multiplicity<T> mult_min(const Multiplicity<T>& a, const Multiplicity<T>& b) {
    if (a.infinite) return b;
    if (b.infinite) return a;
    return a.value < b.value ? a : b;
}

/** The plain-data fields of a layered software network read by the algorithm. */
template <class T>
struct LsnInput {
    Matrix<T> dag;                       ///< (n x n) call graph; an edge is a strictly positive entry
    std::vector<Multiplicity<T>> mult;   ///< (n) declared multiplicity; short vectors are padded with Inf
    std::vector<LsnElementType> type;    ///< (n) element kind
    std::vector<bool> isref;             ///< (n) reference task flag
    std::vector<bool> hassetup;        ///< (n) setup task flag; may be empty
    std::vector<bool> entry_has_arrival; ///< (n) entry with an open arrival; may be empty
};

namespace detail {

/**
 * Kahn topological sort of an adjacency matrix, port of matlab/util/kahn.m.
 * The queue is FIFO, as in both MATLAB and the JAR, so the order is identical.
 * A cycle leaves the order short, which the caller must reject: MATLAB then
 * indexes with a zero and the JAR runs off the end of its list, so neither
 * reference implementation defines an answer there.
 */
template <class T>
std::vector<std::size_t> kahn(const Matrix<T>& adj) {
    const std::size_t n = adj.rows();
    const T zero = num_traits<T>::from_int(0);
    std::vector<std::size_t> indeg(n, 0);
    for (std::size_t c = 0; c < n; ++c)
        for (std::size_t r = 0; r < n; ++r)
            if (adj(r, c) > zero) ++indeg[c];

    std::deque<std::size_t> q;
    for (std::size_t v = 0; v < n; ++v)
        if (indeg[v] == 0) q.push_back(v);

    std::vector<std::size_t> order;
    order.reserve(n);
    while (!q.empty()) {
        const std::size_t i = q.front();
        q.pop_front();
        order.push_back(i);
        for (std::size_t j = 0; j < n; ++j)
            if (adj(i, j) > zero && --indeg[j] == 0) q.push_back(j);
    }
    return order;
}

}  // namespace detail

/**
 * @param lsn the call graph and the per-element attributes
 * @return    (n) maximum multiplicity sustainable by each element
 */
template <class T>
std::vector<Multiplicity<T>> lsn_max_multiplicity(const LsnInput<T>& lsn) {
    const std::size_t n = lsn.dag.rows();
    if (lsn.dag.cols() != n) throw InputError("lsn_max_multiplicity: the call graph is not square");
    if (lsn.type.size() != n || lsn.isref.size() != n)
        throw InputError("lsn_max_multiplicity: type/isref disagree with the graph size");
    if (!lsn.hassetup.empty() && lsn.hassetup.size() != n)
        throw InputError("lsn_max_multiplicity: hassetup disagrees with the graph size");
    if (!lsn.entry_has_arrival.empty() && lsn.entry_has_arrival.size() != n)
        throw InputError("lsn_max_multiplicity: entry_has_arrival disagrees with the graph size");

    // an edge is any strictly positive weight, MATLAB's ag = lsn.dag > 0
    const T zero = num_traits<T>::from_int(0);
    Matrix<T> ag(n, n, zero);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            if (lsn.dag(i, j) > zero) ag(i, j) = num_traits<T>::from_int(1);

    const std::vector<std::size_t> order = detail::kahn(ag);
    if (order.size() != n)
        throw InputError("lsn_max_multiplicity: the call graph has a cycle, no topological order");

    // multiplicities shorter than the graph are unbounded, as MATLAB pads with Inf
    std::vector<Multiplicity<T>> mult = lsn.mult;
    if (mult.size() > n) throw InputError("lsn_max_multiplicity: more multiplicities than elements");
    while (mult.size() < n) mult.push_back(Multiplicity<T>::inf());

    std::vector<Multiplicity<T>> inflow(n, Multiplicity<T>::of(0));
    for (std::size_t i = 0; i < n; ++i) {
        if (lsn.type[i] == LsnElementType::TASK && lsn.isref[i]) {
            inflow[i] = mult[i];
        } else if (lsn.type[i] == LsnElementType::ENTRY && !lsn.entry_has_arrival.empty() &&
                   lsn.entry_has_arrival[i]) {
            // an entry fed by an open arrival needs at least one thread of its
            // parent task (MATLAB only; the JAR has no such branch)
            inflow[i] = Multiplicity<T>::of(1);
        }
    }

    std::vector<Multiplicity<T>> outflow(n, Multiplicity<T>::of(0));
    for (std::size_t k = 0; k < n; ++k) {
        const std::size_t i = order[k];
        const bool has_setup = !lsn.hassetup.empty() && lsn.hassetup[i];
        if (has_setup && inflow[i].positive()) {
            // setup-task multiplicity rationale: see _kb/03-api-layer.md (cpp port notes: lsn)
            outflow[i] = mult[i];
        } else {
            outflow[i] = mult_min(inflow[i], mult[i]);
        }
        for (std::size_t j = 0; j < n; ++j)
            if (j != i && ag(i, j) > zero) inflow[j] = mult_add(inflow[j], outflow[i]);
    }

    for (std::size_t i = 0; i < n; ++i)
        if (lsn.type[i] == LsnElementType::TASK && mult[i].infinite && !lsn.isref[i])
            outflow[i] = Multiplicity<T>::inf();

    return outflow;
}

}  // namespace lsn
}  // namespace line

#endif  // LINE_API_LQN_LSN_MAX_MULTIPLICITY_H
