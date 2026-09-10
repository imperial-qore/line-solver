/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_DAG_MAKESPAN_H
#define LINE_API_FJ_DAG_MAKESPAN_H

/**
 * Makespan of a task system with precedence constraints.
 *
 * Templated port of matlab/src/api/fj/fj_dag_makespan.m.
 *
 * Because the precedence relation is acyclic, so is the chain whose state is
 * the SET of completed tasks, and the makespan is swept level by level instead
 * of solved as a linear system. In a state with completed set S the eligible
 * tasks are those all of whose predecessors lie in S; task i among the k of
 * them completes at rate rate(i,k), so the state is held for M(S) = 1/T(S) and
 * moves to S + {i} with probability b(S,i) = rate(i,k)/T(S). Making the rate
 * depend on the concurrency is what couples the task system to the queueing
 * network underneath it.
 *
 *   p(R) = sum_{S -> R} p(S) b(S,R),
 *   D(R) = M(R) p(R) + sum_{S -> R} b(S,R) D(S),
 *
 * started at p(empty) = 1; the makespan is D at the fully completed state, and
 * the per-task initiation and completion times accumulate over the transitions
 * that start and that finish each task.
 */

#include <cstddef>
#include <vector>

#include "line/api/fj/fj_types.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace fj {

/** [C, I, Cend, E] of fj_dag_makespan. */
template <class T>
struct FJDagMakespanResult {
    T C;
    std::vector<T> I;
    std::vector<T> Cend;
    std::vector<T> E;
};

/**
 * @param pred n by n precedence relation, pred(i,j) nonzero when i precedes j
 * @param rate n by n table whose entry (i,k) is the rate of task i at concurrency k
 * @return     the makespan with the per-task initiation, completion and execution times
 */
template <class T>
FJDagMakespanResult<T> fj_dag_makespan(const Matrix<T>& pred, const Matrix<T>& rate) {
    const std::size_t n = pred.rows();
    if (pred.cols() != n) throw InputError("fj_dag_makespan: pred must be square");
    if (n < 1) throw InputError("fj_dag_makespan: at least one task is required");
    if (n > 20) throw InputError("fj_dag_makespan: the completed-set sweep enumerates 2^n states");
    if (rate.rows() != n || rate.cols() != n)
        throw InputError("fj_dag_makespan: the rate table must be n by n");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t k = 0; k < n; ++k)
            if (!(rate(i, k) > zero))
                throw InputError("fj_dag_makespan: all completion rates must be positive");

    // Predecessor masks, and a topological pass that rejects a cycle
    std::vector<std::size_t> predmask(n, 0);
    std::vector<long> indeg(n, 0);
    for (std::size_t j = 0; j < n; ++j)
        for (std::size_t i = 0; i < n; ++i)
            if (!(pred(i, j) == zero)) {
                predmask[j] |= (static_cast<std::size_t>(1) << i);
                ++indeg[j];
            }
    std::vector<bool> seen(n, false);
    std::size_t remaining = n;
    for (std::size_t pass = 0; pass < n; ++pass) {
        std::size_t pick = n;
        for (std::size_t i = 0; i < n; ++i)
            if (!seen[i] && indeg[i] == 0) { pick = i; break; }
        if (pick == n) break;
        seen[pick] = true;
        indeg[pick] = -1;
        for (std::size_t j = 0; j < n; ++j)
            if (!(pred(pick, j) == zero) && indeg[j] > 0) --indeg[j];
        --remaining;
    }
    if (remaining > 0)
        throw InputError("fj_dag_makespan: the precedence relation contains a cycle");

    const std::size_t nmask = static_cast<std::size_t>(1) << n;
    std::vector<std::size_t> eligmask(nmask, 0);
    std::vector<bool> closed(nmask, false);
    for (std::size_t mask = 0; mask < nmask; ++mask) {
        bool ok = true;
        std::size_t em = 0;
        for (std::size_t i = 0; i < n; ++i) {
            const std::size_t bit = static_cast<std::size_t>(1) << i;
            if (mask & bit) {
                // A completed task must have all of its predecessors completed
                if ((predmask[i] & mask) != predmask[i]) { ok = false; break; }
            } else if ((predmask[i] & mask) == predmask[i]) {
                em |= bit;
            }
        }
        closed[mask] = ok;
        if (ok) eligmask[mask] = em;
    }

    std::vector<T> p(nmask, zero), D(nmask, zero);
    p[0] = one;
    FJDagMakespanResult<T> out;
    out.I.assign(n, zero);
    out.Cend.assign(n, zero);
    out.E.assign(n, zero);

    for (std::size_t mask = 0; mask < nmask; ++mask) {
        if (!closed[mask]) continue;
        std::vector<std::size_t> elig;
        for (std::size_t i = 0; i < n; ++i)
            if (eligmask[mask] & (static_cast<std::size_t>(1) << i)) elig.push_back(i);
        const std::size_t k = elig.size();
        if (k == 0) continue;
        T Ttot = zero;
        for (std::size_t idx = 0; idx < k; ++idx) Ttot += rate(elig[idx], k - 1);
        // Holding time of this state, weighted by the probability of reaching it
        D[mask] += p[mask] / Ttot;
        for (std::size_t idx = 0; idx < k; ++idx) {
            const std::size_t i = elig[idx];
            const T b = rate(i, k - 1) / Ttot;
            const std::size_t nxt = mask | (static_cast<std::size_t>(1) << i);
            const T contrib = b * D[mask];
            p[nxt] += p[mask] * b;
            D[nxt] += contrib;
            // Task i completes on this transition
            out.Cend[i] += contrib;
            // Tasks that first become eligible on this transition start on it
            const std::size_t fresh = eligmask[nxt] & ~eligmask[mask];
            for (std::size_t j = 0; j < n; ++j)
                if (fresh & (static_cast<std::size_t>(1) << j)) out.I[j] += contrib;
        }
    }

    out.C = D[nmask - 1];
    for (std::size_t i = 0; i < n; ++i) out.E[i] = out.Cend[i] - out.I[i];
    return out;
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_DAG_MAKESPAN_H
