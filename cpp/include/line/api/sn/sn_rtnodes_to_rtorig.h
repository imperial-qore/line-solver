/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SN_SN_RTNODES_TO_RTORIG_H
#define LINE_API_SN_SN_RTNODES_TO_RTORIG_H

/**
 * Port of matlab/src/api/sn/sn_rtnodes_to_rtorig.m.
 *
 * `sn.rtnodes` is over every node the refresh built, INCLUDING the artificial
 * `CS_*` class-switch nodes that the routing expansion inserts. This recovers
 * the routing over the ORIGINAL nodes by taking the stochastic complement over
 * their rows, i.e. by eliminating the artificial ones rather than dropping
 * them -- dropping would lose the class switch the node performs.
 *
 * The original nodes are those BEFORE the first `CS_` node, because the
 * expansion appends. A model with no class switching therefore keeps every
 * node and the complement is the identity operation on `rtnodes`.
 *
 * ARITHMETIC: field. dtmc_stochcomp is a linear solve.
 */

#include <cstddef>
#include <string>
#include <vector>

#include "line/api/mc/dtmc_stochcomp.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"

namespace line {
namespace api {

/** What sn_rtnodes_to_rtorig returns: the flat matrix and its per-class-pair blocks. */
template <class T>
struct SnRtOrig {
    Matrix<T> rtorig;                        ///< (norig*nclasses) square
    std::vector<std::vector<Matrix<T>>> cells;  ///< cells[r][s] is (norig x norig)
    std::size_t norig = 0;                   ///< nodes kept, the `csshift` of the reference
};

template <class T>
SnRtOrig<T> sn_rtnodes_to_rtorig(const qn::NetworkStruct<T>& sn) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t K = sn.nclasses, I = sn.nodes.size();
    SnRtOrig<T> out;
    std::size_t csshift = I;
    for (std::size_t a = 0; a < I; ++a)
        if (sn.nodes[a].name.compare(0, 3, "CS_") == 0) {
            csshift = a;
            break;
        }
    out.norig = csshift;
    std::vector<std::size_t> keep;
    keep.reserve(csshift * K);
    for (std::size_t a = 0; a < csshift; ++a)
        for (std::size_t r = 0; r < K; ++r) keep.push_back(a * K + r);
    if (sn.rtnodes.rows() == I * K && !keep.empty()) {
        out.rtorig = mc::dtmc_stochcomp(sn.rtnodes, keep);
    } else {
        out.rtorig = Matrix<T>(csshift * K, csshift * K, zero);
    }
    out.cells.assign(K, std::vector<Matrix<T>>(K, Matrix<T>(csshift, csshift, zero)));
    for (std::size_t a = 0; a < csshift; ++a) {
        if (sn.nodes[a].nodetype == qn::NodeType::Sink) continue;
        for (std::size_t b = 0; b < csshift; ++b)
            for (std::size_t r = 0; r < K; ++r)
                for (std::size_t s = 0; s < K; ++s) {
                    const double v = num_traits<T>::to_double(out.rtorig(a * K + r, b * K + s));
                    if (v != v) continue;  // the reference zeroes NaN before the copy
                    out.cells[r][s](a, b) = out.rtorig(a * K + r, b * K + s);
                }
    }
    return out;
}

}  // namespace api
}  // namespace line

#endif  // LINE_API_SN_SN_RTNODES_TO_RTORIG_H
