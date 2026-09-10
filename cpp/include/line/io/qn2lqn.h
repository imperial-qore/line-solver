/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_IO_QN2LQN_H
#define LINE_IO_QN2LQN_H

/**
 * Convert a closed queueing network into the layered model used by SolverLQNS.
 *
 * Port of matlab/src/io/QN2LQN.m.  A chain becomes a reference task, every
 * Queue or Delay becomes an infinite-multiplicity server task, and the node
 * routing becomes activity precedences on the reference task.  Forks use an
 * AND fork and Joins an AND join; the other routing nodes use probabilistic OR
 * forks.
 *
 * C++ keeps visits on two explicit index spaces.  The MATLAB routine indexes
 * `sn.visits` with a node index, while this port uses `nodevisits`, the
 * node-indexed twin, so a Router/ClassSwitch before a station cannot shift the
 * service activity onto the wrong row.
 */

#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "line/lang/lqn/lqn_builder.h"
#include "line/lang/qn/network_struct.h"
#include "line/util/error.h"

namespace line {
namespace qn {
template <class T>
class Network;
}

namespace io {

namespace detail {

inline bool qn2lqn_routing_node(qn::NodeType t) {
    return t == qn::NodeType::Queue || t == qn::NodeType::Delay ||
           t == qn::NodeType::ClassSwitch || t == qn::NodeType::Router ||
           t == qn::NodeType::Logger || t == qn::NodeType::Fork ||
           t == qn::NodeType::Join;
}

inline std::string qn2lqn_indexed(const char* prefix, std::size_t i, std::size_t r) {
    return std::string(prefix) + std::to_string(i + 1) + "_" + std::to_string(r + 1);
}

inline std::string qn2lqn_pseudo(const char* prefix, std::size_t c, std::size_t i,
                                 std::size_t r) {
    return std::string(prefix) + "_" + std::to_string(c + 1) + "_" +
           std::to_string(i + 1) + "_" + std::to_string(r + 1);
}

}  // namespace detail

/** Port of MATLAB `QN2LQN(model)`, over the refreshed C++ NetworkStruct. */
template <class T>
lqn::LqnModel<T> qn2lqn(const qn::NetworkStruct<T>& sn) {
    typedef num_traits<T> nt;
    const T one = nt::from_int(1);
    const std::size_t I = sn.nodes.size();
    const std::size_t K = sn.nclasses;
    const std::size_t C = sn.nchains;

    if (I == 0 || K == 0 || C == 0)
        throw InputError("QN2LQN: the queueing network has no closed chain to convert");
    if (sn.nodevisits.size() != C || sn.inchain.size() != C)
        throw InputError("QN2LQN: the NetworkStruct visits or chains are not refreshed");
    for (std::size_t r = 0; r < K; ++r)
        if (!std::isfinite(sn.classes[r].population))
            throw UnsupportedError(
                "QN2LQN: open classes are not supported; SolverQNS sends open models to "
                "qnsolver directly");

    std::vector<std::size_t> chain_of(K, C);
    for (std::size_t c = 0; c < C; ++c)
        for (std::size_t x = 0; x < sn.inchain[c].size(); ++x) {
            const std::size_t r1 = sn.inchain[c][x];
            if (r1 == 0 || r1 > K)
                throw InputError("QN2LQN: chain " + std::to_string(c + 1) +
                                 " contains an out-of-range class index");
            if (chain_of[r1 - 1] != C)
                throw InputError("QN2LQN: class " + std::to_string(r1) +
                                 " belongs to more than one chain");
            chain_of[r1 - 1] = c;
        }
    for (std::size_t r = 0; r < K; ++r)
        if (chain_of[r] == C)
            throw InputError("QN2LQN: class " + std::to_string(r + 1) +
                             " belongs to no chain");

    lqn::LqnBuilder<T> b;
    const std::string pseudo_host = sn.name.empty() ? std::string("QN") : sn.name;
    b.processor(pseudo_host, std::numeric_limits<double>::infinity(), lang::SchedStrategy::INF);

    std::vector<std::string> ref_task(C), ref_entry(C);
    for (std::size_t c = 0; c < C; ++c) {
        double jobs = 0.0;
        for (std::size_t x = 0; x < sn.inchain[c].size(); ++x)
            jobs += sn.classes[sn.inchain[c][x] - 1].population;
        ref_task[c] = "RefTask_" + std::to_string(c + 1);
        ref_entry[c] = "Chain_" + std::to_string(c + 1);
        b.task(ref_task[c], jobs, lang::SchedStrategy::REF, pseudo_host);
        b.entry(ref_entry[c], ref_task[c]);
    }

    // MATLAB's cell arrays E/A/PA are represented by empty names for absent
    // elements.  A name is sufficient because LqnBuilder resolves handles by
    // name at finalization.
    std::vector<std::string> entries(I * K), service_acts(I * K);
    const auto ir = [K](std::size_t i, std::size_t r) { return i * K + r; };
    const auto cir = [I, K](std::size_t c, std::size_t i, std::size_t r) {
        return (c * I + i) * K + r;
    };

    for (std::size_t i = 0; i < I; ++i) {
        const qn::NodeDef& nd = sn.nodes[i];
        if (nd.nodetype != qn::NodeType::Queue && nd.nodetype != qn::NodeType::Delay)
            continue;
        if (nd.station == 0 || nd.station > sn.stations.size())
            throw InputError("QN2LQN: service node '" + nd.name + "' has no station row");
        const std::size_t ist = nd.station - 1;
        b.processor(nd.name, sn.stations[ist].nservers, sn.stations[ist].sched);
        const std::string task = "T_" + nd.name;
        b.task(task, std::numeric_limits<double>::infinity(), lang::SchedStrategy::INF, nd.name);
        for (std::size_t r = 0; r < K; ++r) {
            const std::size_t c = chain_of[r];
            if (sn.nodevisits[c].rows() != I || sn.nodevisits[c].cols() != K)
                throw InputError("QN2LQN: nodevisits has the wrong dimensions for chain " +
                                 std::to_string(c + 1));
            if (!(nt::to_double(sn.nodevisits[c](i, r)) > 0.0)) continue;
            if (ist >= sn.service.size() || r >= sn.service[ist].size() ||
                sn.service[ist][r].disabled)
                throw InputError("QN2LQN: visited station '" + nd.name + "', class '" +
                                 sn.classes[r].name + "' has no service distribution");
            entries[ir(i, r)] = detail::qn2lqn_indexed("E", i, r);
            service_acts[ir(i, r)] = detail::qn2lqn_indexed("Q", i, r);
            b.entry(entries[ir(i, r)], task);
            b.activity(service_acts[ir(i, r)], sn.service[ist][r], task);
            b.bound_to(service_acts[ir(i, r)], entries[ir(i, r)]);
            b.replies_to(service_acts[ir(i, r)], entries[ir(i, r)]);
        }
    }

    std::vector<std::string> pseudo(C * I * K);
    std::vector<bool> bound(C, false);
    std::vector<std::size_t> bound_i(C, 0), bound_r(C, 0);
    const auto incoming = [&](std::size_t i, std::size_t r) {
        const std::size_t col = i * K + r;
        if (sn.rtnodes.cols() <= col) return false;
        for (std::size_t row = 0; row < sn.rtnodes.rows(); ++row)
            if (nt::to_double(sn.rtnodes(row, col)) > 0.0) return true;
        return false;
    };

    for (std::size_t i = 0; i < I; ++i) {
        const qn::NodeType type = sn.nodes[i].nodetype;
        if (type == qn::NodeType::ClassSwitch || type == qn::NodeType::Router ||
            type == qn::NodeType::Logger || type == qn::NodeType::Fork ||
            type == qn::NodeType::Join) {
            const char* prefix = (type == qn::NodeType::Fork || type == qn::NodeType::Join)
                                     ? "FJ"
                                     : "CS";
            for (std::size_t r = 0; r < K; ++r) {
                if (!incoming(i, r)) continue;
                const std::size_t c = chain_of[r];
                pseudo[cir(c, i, r)] = detail::qn2lqn_pseudo(prefix, c, i, r);
                b.activity(pseudo[cir(c, i, r)], lang::Distrib<T>::immediate(), ref_task[c]);
            }
        } else if (type == qn::NodeType::Queue || type == qn::NodeType::Delay) {
            for (std::size_t r = 0; r < K; ++r) {
                const std::size_t c = chain_of[r];
                if (service_acts[ir(i, r)].empty()) continue;
                pseudo[cir(c, i, r)] = detail::qn2lqn_indexed("A", i, r);
                b.activity(pseudo[cir(c, i, r)], lang::Distrib<T>::immediate(), ref_task[c]);
                const std::size_t first_r = sn.inchain[c].front() - 1;
                const std::size_t refstat = sn.classes[first_r].refstat;
                if (refstat == 0 || refstat > sn.station_to_node.size())
                    throw InputError("QN2LQN: chain " + std::to_string(c + 1) +
                                     " has an invalid reference station");
                const std::size_t refnode = sn.station_to_node[refstat - 1] - 1;
                if (i == refnode && r == first_r) {
                    b.bound_to(pseudo[cir(c, i, r)], ref_entry[c]);
                    bound[c] = true;
                    bound_i[c] = i;
                    bound_r[c] = r;
                }
                b.sync_call(pseudo[cir(c, i, r)], entries[ir(i, r)], one);
            }
        } else if (type != qn::NodeType::Source && type != qn::NodeType::Sink) {
            throw UnsupportedError("QN2LQN: node '" + sn.nodes[i].name + "' has type " +
                                   std::string(lang::node_type_to_text(type)) +
                                   ", which the MATLAB conversion does not map");
        }
    }
    for (std::size_t c = 0; c < C; ++c)
        if (!bound[c])
            throw InputError("QN2LQN: chain " + std::to_string(c + 1) +
                             " has no visited reference station to bind to its entry");

    for (std::size_t c = 0; c < C; ++c) {
        for (std::size_t i = 0; i < I; ++i) {
            if (!detail::qn2lqn_routing_node(sn.nodes[i].nodetype)) continue;
            for (std::size_t rx = 0; rx < sn.inchain[c].size(); ++rx) {
                const std::size_t r = sn.inchain[c][rx] - 1;
                const std::string& pre = pseudo[cir(c, i, r)];
                if (pre.empty() || !incoming(i, r)) continue;
                std::vector<std::string> posts;
                std::vector<T> probs;
                for (std::size_t j = 0; j < I; ++j) {
                    if (!detail::qn2lqn_routing_node(sn.nodes[j].nodetype) ||
                        sn.nodes[j].nodetype == qn::NodeType::Join)
                        continue;
                    for (std::size_t sx = 0; sx < sn.inchain[c].size(); ++sx) {
                        const std::size_t s = sn.inchain[c][sx] - 1;
                        const T pr = sn.rtnodes(i * K + r, j * K + s);
                        if (!(nt::to_double(pr) > 0.0)) continue;
                        if (j == bound_i[c] && s == bound_r[c]) {
                            const std::string end = detail::qn2lqn_pseudo("End", c, i, r);
                            b.activity(end, lang::Distrib<T>::immediate(), ref_task[c]);
                            posts.push_back(end);
                        } else {
                            const std::string& post = pseudo[cir(c, j, s)];
                            if (post.empty())
                                throw InputError("QN2LQN: route from node '" + sn.nodes[i].name +
                                                 "' reaches a node/class with no pseudo-activity");
                            posts.push_back(post);
                        }
                        probs.push_back(pr);
                    }
                }
                if (posts.empty()) continue;
                if (sn.nodes[i].nodetype == qn::NodeType::Fork)
                    b.and_fork(pre, posts);
                else
                    b.or_fork(pre, posts, probs);
            }
        }
    }

    // Join destinations were deliberately omitted above.  They are the post
    // side of an AND join (or a serial edge for a degenerate one-input join).
    for (std::size_t c = 0; c < C; ++c)
        for (std::size_t j = 0; j < I; ++j) {
            if (sn.nodes[j].nodetype != qn::NodeType::Join) continue;
            for (std::size_t sx = 0; sx < sn.inchain[c].size(); ++sx) {
                const std::size_t s = sn.inchain[c][sx] - 1;
                const std::string& post = pseudo[cir(c, j, s)];
                if (post.empty()) continue;
                std::vector<std::string> pres;
                for (std::size_t i = 0; i < I; ++i) {
                    if (!detail::qn2lqn_routing_node(sn.nodes[i].nodetype)) continue;
                    for (std::size_t rx = 0; rx < sn.inchain[c].size(); ++rx) {
                        const std::size_t r = sn.inchain[c][rx] - 1;
                        if (nt::to_double(sn.rtnodes(i * K + r, j * K + s)) <= 0.0) continue;
                        const std::string& pre = pseudo[cir(c, i, r)];
                        if (!pre.empty()) pres.push_back(pre);
                    }
                }
                if (pres.size() > 1)
                    b.and_join(pres, post);
                else if (pres.size() == 1)
                    b.serial(pres.front(), post);
            }
        }

    // Finalize once here as a structural validation: every entry must have a
    // bound activity and every precedence/call target must exist.
    (void)b.build();
    return b.model();
}

/** MATLAB-compatible model-level entry point; refreshes the Network first. */
template <class T>
lqn::LqnModel<T> qn2lqn(qn::Network<T>& model) {
    return qn2lqn(model.get_struct());
}

}  // namespace io
}  // namespace line

#endif  // LINE_IO_QN2LQN_H
