/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SN_SN_PRINT_ROUTING_MATRIX_H
#define LINE_API_SN_SN_PRINT_ROUTING_MATRIX_H

/**
 * Port of matlab/src/api/sn/sn_print_routing_matrix.m.
 *
 * The human-readable form of `sn.rtnodes`: one line per positive (node, class)
 * to (node, class) edge. A Cache's outgoing probability is reported as
 * "state-dependent" rather than as a number, because the number the refresh
 * left there is a placeholder the cache fixed point later replaces; a Sink has
 * no outgoing edge to report; and a class whose routing at that node is
 * DISABLED is skipped even though the expanded matrix may carry a value for it.
 *
 * Returned as a string rather than printed, so a caller can route it to a log,
 * a CLI or a test. `sn_print_routing_matrix` in the reference ends with a
 * newline and so does this.
 *
 * ARITHMETIC: field, but the rendering is in double.
 */

#include <cstddef>
#include <cstdio>
#include <string>

#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"

namespace line {
namespace api {

/**
 * @param onlyclass 1-based class index to restrict to, 0 for every class.
 *        The reference matches on the class NAME and keeps an edge when
 *        EITHER end names it, which is what the index test below reproduces.
 */
template <class T>
std::string sn_print_routing_matrix(const qn::NetworkStruct<T>& sn, std::size_t onlyclass = 0) {
    const std::size_t I = sn.nodes.size(), K = sn.nclasses;
    std::string out;
    if (sn.rtnodes.rows() != I * K) return out + "\n";
    for (std::size_t i = 0; i < I; ++i)
        for (std::size_t r = 0; r < K; ++r)
            for (std::size_t j = 0; j < I; ++j)
                for (std::size_t s = 0; s < K; ++s) {
                    if (!(num_traits<T>::to_double(sn.rtnodes(i * K + r, j * K + s)) > 0.0))
                        continue;
                    std::string pr;
                    if (sn.nodes[i].nodetype == qn::NodeType::Cache) {
                        pr = "state-dependent";
                    } else if (sn.nodes[i].nodetype == qn::NodeType::Sink) {
                        continue;
                    } else {
                        if (r < sn.nodes[i].routing.size() &&
                            sn.nodes[i].routing[r] == qn::RoutingStrategy::DISABLED)
                            continue;
                        char buf[64];
                        std::snprintf(buf, sizeof(buf), "%f",
                                      num_traits<T>::to_double(sn.rtnodes(i * K + r, j * K + s)));
                        pr = buf;
                    }
                    if (onlyclass != 0 && r + 1 != onlyclass && s + 1 != onlyclass) continue;
                    out += "\n" + sn.nodes[i].name + " [" + sn.classes[r].name + "] => " +
                           sn.nodes[j].name + " [" + sn.classes[s].name + "] : Pr=" + pr;
                }
    out += "\n";
    return out;
}

}  // namespace api
}  // namespace line

#endif  // LINE_API_SN_SN_PRINT_ROUTING_MATRIX_H
