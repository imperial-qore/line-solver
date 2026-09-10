/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MDD_MDD_REACHSET_H
#define LINE_API_MDD_MDD_REACHSET_H

/**
 * Reachability set generation into a decision diagram.
 *
 * Port of matlab/src/api/mdd/mdd_reachset.m, jline.api.mdd.Mdd_reachset and
 * python/line_solver/api/mdd/reachset.py, after A.S. Miner, G. Ciardo,
 * "Efficient Reachability Set Generation and Storage Using Decision Diagrams",
 * ICATPN 1999, LNCS 1639, pp.6-25.
 *
 * This is the basic (explicit-frontier) realisation: a breadth-first search
 * enumerates successors while the MDD provides the O(K) membership test that
 * replaces the usual explicit visited hash. The stored set lives entirely in
 * the MDD (O(#nodes) memory); only the transient BFS frontier is held
 * explicitly. Symbolic image computation / saturation, which removes the
 * explicit frontier too, is the natural next step but is out of scope here and
 * is what caps the wall-clock saving -- the storage saving is real regardless.
 */

#include <cstddef>
#include <vector>

#include "line/api/mdd/mdd.h"
#include "line/api/mdd/mdd_types.h"

namespace line {
namespace mdd {

/**
 * Generate and store the reachability set into a quasi-reduced ordered MDD.
 *
 * @param domain per-level local-state counts, values 0..domain[k]-1
 * @param init the initial global state, 0-based local values
 * @param nextfun the next-state function
 * @return an MDD holding every state reachable from init
 */
inline MDD mdd_reachset(const std::vector<int>& domain, const std::vector<int>& init,
                        const MddNextState& nextfun) {
    MDD diagram(domain);
    diagram.insert(init);

    std::vector<std::vector<int>> frontier;
    frontier.push_back(init);
    std::size_t head = 0;
    while (head < frontier.size()) {
        const std::vector<int> s = frontier[head];
        ++head;
        if (nextfun) {
            const std::vector<std::vector<int>> successors = nextfun(s);
            for (std::size_t r = 0; r < successors.size(); ++r) {
                if (!diagram.member(successors[r])) {
                    diagram.insert(successors[r]);
                    frontier.push_back(successors[r]);
                }
            }
        }
        // drop already-expanded rows periodically to bound frontier memory
        if (head > 1024 && 2 * head > frontier.size()) {
            frontier.erase(frontier.begin(), frontier.begin() + static_cast<long>(head));
            head = 0;
        }
    }
    diagram.compact();  // reclaim the dead nodes left by the append-only inserts
    return diagram;
}

}  // namespace mdd
}  // namespace line

#endif  // LINE_API_MDD_MDD_REACHSET_H
