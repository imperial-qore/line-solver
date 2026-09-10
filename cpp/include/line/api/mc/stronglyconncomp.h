/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MC_STRONGLYCONNCOMP_H
#define LINE_API_MC_STRONGLYCONNCOMP_H

/**
 * Strongly connected components of a directed graph, and which of them are
 * recurrent (closed under the successor relation).
 *
 * Templated port of matlab/util/stronglyconncomp.m, the decomposition on which
 * dtmc_solve_reducible, ctmc_solve_reducible and
 * ctmc_solve_reducible_blkdecomp all rest. Tarjan's algorithm, with the
 * components renumbered by decreasing size exactly as MATLAB does (a stable
 * sort, so components of equal size keep their completion order), and a
 * component declared recurrent when no state in it has a successor outside it.
 *
 * Two details of the MATLAB version are reproduced deliberately. The depth
 * first search follows the COLUMNS of the adjacency matrix, i.e. the reversed
 * graph, while the recurrence test follows the ROWS; the component partition
 * is the same for a graph and its reverse, so this only affects the discovery
 * order, but reproducing it keeps the component numbering identical when sizes
 * tie. The recursion of the MATLAB original is replaced by an explicit stack,
 * which visits vertices in the same order and does not overflow on the tens of
 * thousands of states a lumped generator can carry.
 *
 * The computation is combinatorial: an entry only ever has its non-zero-ness
 * tested, no arithmetic is performed, so it is exact at every number type.
 */

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mc {

struct SccResult {
    /** Component index of each state, 1-based as in MATLAB (0 is never used). */
    std::vector<std::size_t> scc;
    /** recurrent[c-1] is true when component c has no edge leaving it. */
    std::vector<bool> recurrent;
    /** Member states of each component, ascending. */
    std::vector<std::vector<std::size_t>> members;
    std::size_t numSCC() const { return members.size(); }
};

/**
 * @param A adjacency matrix; an edge i -> j exists iff A(i,j) is non-zero
 */
template <class T>
SccResult stronglyconncomp(const Matrix<T>& A) {
    const std::size_t n = A.rows();
    if (A.cols() != n) throw InputError("stronglyconncomp: adjacency matrix is not square");
    const T zero = num_traits<T>::from_int(0);

    // Successor and predecessor lists. The search runs on predecessors (the
    // MATLAB find(e(:,i))), the recurrence test on successors.
    std::vector<std::vector<std::size_t>> pred(n), succ(n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            if (A(i, j) != zero) {
                succ[i].push_back(j);
                pred[j].push_back(i);
            }

    std::vector<std::size_t> index(n, 0), low(n, 0);
    std::vector<char> onstack(n, 0);
    std::vector<std::size_t> stack;
    std::vector<std::vector<std::size_t>> comps;
    std::size_t counter = 0;

    // Explicit emulation of the recursive Tarjan: each frame is a vertex and
    // the position reached in its predecessor list.
    std::vector<std::pair<std::size_t, std::size_t>> frames;
    for (std::size_t s = 0; s < n; ++s) {
        if (index[s] != 0) continue;
        frames.push_back(std::make_pair(s, static_cast<std::size_t>(0)));
        ++counter;
        index[s] = counter;
        low[s] = counter;
        stack.push_back(s);
        onstack[s] = 1;
        while (!frames.empty()) {
            const std::size_t u = frames.back().first;
            std::size_t& p = frames.back().second;
            if (p < pred[u].size()) {
                const std::size_t v = pred[u][p];
                ++p;
                if (index[v] == 0) {
                    ++counter;
                    index[v] = counter;
                    low[v] = counter;
                    stack.push_back(v);
                    onstack[v] = 1;
                    frames.push_back(std::make_pair(v, static_cast<std::size_t>(0)));
                } else if (onstack[v]) {
                    if (index[v] < low[u]) low[u] = index[v];
                }
                continue;
            }
            // u is finished: close a component or propagate its low link.
            if (low[u] == index[u]) {
                std::vector<std::size_t> comp;
                for (;;) {
                    const std::size_t w = stack.back();
                    stack.pop_back();
                    onstack[w] = 0;
                    comp.push_back(w);
                    if (w == u) break;
                }
                std::sort(comp.begin(), comp.end());
                comps.push_back(comp);
            }
            frames.pop_back();
            if (!frames.empty()) {
                const std::size_t parent = frames.back().first;
                if (low[u] < low[parent]) low[parent] = low[u];
            }
        }
    }

    // Renumber by decreasing size, stably, as MATLAB's sort(...,'descend') does.
    std::vector<std::size_t> order(comps.size());
    std::iota(order.begin(), order.end(), static_cast<std::size_t>(0));
    std::stable_sort(order.begin(), order.end(),
                     [&comps](std::size_t a, std::size_t b) {
                         return comps[a].size() > comps[b].size();
                     });

    SccResult r;
    r.members.resize(comps.size());
    r.scc.assign(n, 0);
    for (std::size_t k = 0; k < order.size(); ++k) {
        r.members[k] = comps[order[k]];
        for (std::size_t v : r.members[k]) r.scc[v] = k + 1;
    }

    r.recurrent.assign(comps.size(), true);
    for (std::size_t k = 0; k < r.members.size(); ++k) {
        for (std::size_t v : r.members[k]) {
            for (std::size_t w : succ[v])
                if (r.scc[w] != k + 1) {
                    r.recurrent[k] = false;
                    break;
                }
            if (!r.recurrent[k]) break;
        }
    }
    return r;
}

}  // namespace mc
}  // namespace line

#endif  // LINE_API_MC_STRONGLYCONNCOMP_H
