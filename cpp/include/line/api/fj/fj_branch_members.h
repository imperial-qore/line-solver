/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_FJ_BRANCH_MEMBERS_H
#define LINE_API_FJ_FJ_BRANCH_MEMBERS_H

/**
 * Activities belonging to each branch of an AND-join.
 *
 * Templated port of matlab/src/api/fj/fj_branch_members.m. Each immediate
 * predecessor of the join activity is the tail of one branch; the branch is
 * recovered by walking backwards along the activity graph until an activity
 * marked POST_AND is reached, that one being the head the AND-fork spawned.
 * Branches between a fork and its join are disjoint paths, so the walk is
 * unambiguous, and the guard on the number of steps bounds it by the activity
 * count.
 *
 * SCOPE. The reference takes a LayeredNetworkStruct, but it reads only four
 * plain numeric fields of it -- graph, ashift, nacts, actposttype -- and no
 * LayeredNetwork object, no cell array of processes and no derived index map.
 * Those four are taken here as explicit arguments, in a small view struct, so
 * the algorithm is ported in full without pulling the LQN object layer into
 * this tree. Nothing is stubbed: the walk, the activity-range test, the
 * merge/start stop condition and the guard are all as written.
 *
 * INDEXING. The reference is 1-based and its indices are absolute LQN element
 * indices, activities occupying (ashift, ashift + nacts]. The port keeps that
 * convention exactly, so an index that appears in the result can be compared
 * against a MATLAB one without a shift: joinIdx, ashift and every returned
 * index are 1-based absolute indices.
 *
 * ARITHMETIC. Only the test graph(i, j) > 0 touches the number type, so this
 * is a structural algorithm and is instantiated at T = Rational as well.
 */

#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace fj {

/** MATLAB's ActivityPrecedenceType codes, as stored in lqn.actposttype. */
enum ActivityPrecedenceCode {
    APC_PRE_SEQ = 1,
    APC_PRE_AND = 2,
    APC_PRE_OR = 3,
    APC_POST_SEQ = 11,
    APC_POST_AND = 12,
    APC_POST_OR = 13,
    APC_POST_LOOP = 14,
    APC_POST_CACHE = 15
};

/**
 * The four LayeredNetworkStruct fields fj_branch_members reads.
 *
 * graph is the full (nidx x nidx) call/precedence graph in absolute 1-based
 * indices, actposttype is indexed the same way, and activities are the indices
 * in (ashift, ashift + nacts].
 */
template <class T>
struct LqnBranchView {
    Matrix<T> graph;
    std::size_t ashift = 0;
    std::size_t nacts = 0;
    std::vector<int> actposttype;  ///< absolute-index vector, 1-based reading
};

/**
 * Branch membership of an AND-join.
 *
 * @param lqn      the four fields listed above
 * @param joinaidx absolute 1-based index of the AND-join activity
 * @return one vector per branch, head last: entry 0 is the tail (the immediate
 *         predecessor of the join) and the last entry is the branch head, which
 *         is the order the reference builds the chain in
 */
template <class T>
std::vector<std::vector<std::size_t>> fj_branch_members(const LqnBranchView<T>& lqn,
                                                        std::size_t joinaidx) {
    const std::size_t nidx = lqn.graph.rows();
    if (lqn.graph.cols() != nidx) throw InputError("fj_branch_members: graph is not square");
    if (joinaidx < 1 || joinaidx > nidx)
        throw InputError("fj_branch_members: the join index is out of the graph");
    if (lqn.actposttype.size() < lqn.ashift + lqn.nacts)
        throw InputError("fj_branch_members: actposttype is shorter than the activity range");
    const T zero = num_traits<T>::from_int(0);

    std::vector<std::vector<std::size_t>> members;
    for (std::size_t tail = 1; tail <= nidx; ++tail) {
        if (!(lqn.graph(tail - 1, joinaidx - 1) > zero)) continue;
        if (tail <= lqn.ashift || tail > lqn.ashift + lqn.nacts) continue;  // not an activity

        std::vector<std::size_t> chain;
        chain.push_back(tail);
        std::size_t cur = tail;
        std::size_t guard = 0;
        while (guard < lqn.nacts) {
            ++guard;
            if (lqn.actposttype[cur - 1] == APC_POST_AND) break;  // the branch head
            std::size_t prev = 0;
            std::size_t nprev = 0;
            for (std::size_t p = 1; p <= nidx; ++p) {
                if (!(lqn.graph(p - 1, cur - 1) > zero)) continue;
                if (p <= lqn.ashift || p > lqn.ashift + lqn.nacts) continue;
                ++nprev;
                prev = p;
            }
            if (nprev != 1) break;  // a merge, or the start of the graph
            cur = prev;
            chain.push_back(cur);
        }
        members.push_back(chain);
    }
    return members;
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_FJ_BRANCH_MEMBERS_H
