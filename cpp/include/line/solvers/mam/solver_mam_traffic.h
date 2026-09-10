/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_MAM_SOLVER_MAM_TRAFFIC_H
#define LINE_SOLVERS_MAM_SOLVER_MAM_TRAFFIC_H

/**
 * Port of `solver_mam_traffic.m` and `solver_mam_traffic_mmap.m`: the traffic
 * step of the `dec.mmap` decomposition.
 *
 * WHAT IT COMPUTES. Given the DEPARTURE process of every class at every station
 * as a (D0,D1) pair, it produces the ARRIVAL process seen by every node. Each
 * station's per-class departures are superposed into one marked MAP, that MAP
 * is split along the routing (`npfqn_traffic_split_cs`, which also carries the
 * class switching), and the per-link flows arriving at a node are superposed
 * back (`npfqn_traffic_merge`). The result is the per-link traffic descriptor
 * table the outer fixed point of `solver_mam.m` iterates on.
 *
 * IT IS AN APPROXIMATION, and the approximation is not incidental. The exact
 * superposition of n marked MAPs is the Kronecker sum, whose order is the
 * PRODUCT of the operand orders, so a network of any size exhausts memory after
 * a handful of merges. Both the per-station class superposition here and
 * `npfqn_traffic_merge` therefore compress back to a bounded representation
 * (`mmap_compress`, an APH(2) mixture fit) whenever the order passes
 * `config.space_max`. Compression preserves the class probabilities and the
 * first three backward moments per class; it does NOT preserve the correlation
 * structure, and it does not preserve the per-class rates exactly. Nothing here
 * may be compared against an exact solver for equality: the descriptors are
 * moment-matched surrogates of the true superposed processes. What survives the
 * compression exactly is the SPLIT: `npfqn_traffic_split_cs` is linear in the
 * routing probabilities, so first-moment conservation across a split, and the
 * whole traffic table of a network that never trips `space_max` and never
 * merges more than one flow into a node, are identities.
 *
 * INDEXING. The reference works in an indexing over the NON-CLASS-SWITCH nodes
 * ("NCS"): the class-switch nodes are eliminated by taking the stochastic
 * complement of `sn.rtnodes` over the rows of the surviving nodes, so a switch
 * shows up only as a class change on the edges around it. `ARV` is returned
 * NODE-indexed and of length `nnodes` -- the reference preallocates it as
 * `cell(Inc,1)` and then assigns `ARV{ind}` at node indices, so MATLAB grows it
 * to `I`; the preallocation size is dead. An entry of order 0 is MATLAB's `[]`:
 * a class-switch node, or a Source, which has no arrivals to describe.
 *
 * REFERENCE DEFECT in `solver_mam_traffic.m`, line 104. When a node has no
 * incoming flow above `FineTol` the fallback reads `ARV{ind} = LINKS{jnd,1}`.
 * `jnd` is the loop variable left behind by the `for jnd=1:Inc` above it, so it
 * is `Inc`, and the column is the literal 1 rather than `inc`: the node is
 * handed the link from the LAST non-class-switch node to the FIRST one, which
 * is a flow between two other nodes entirely, or `[]` when that link was never
 * built. The FJ variant of the same file (`solver_mam_traffic_mmap.m`, lines
 * 179-188) writes the intended form -- any non-empty link INTO this node, and a
 * zero-rate MMAP when there is none -- and that is what this port does for both
 * entry points. Propagating the defect would attribute one node's traffic to
 * another, which no downstream consumer could detect.
 *
 * WHAT IS REFUSED BY NAME. The reference's node switch has a branch only for
 * Source, Delay and Queue (plus Fork and Join in the FJ variant). Any other
 * node that a flow passes THROUGH -- a Router, a Logger, a Cache, a Place, a
 * Transition, a Region -- silently emits no outgoing link there, so everything
 * downstream of it is described as receiving no traffic at all. That is not a
 * conservative approximation, it is a wrong answer with no symptom, so such a
 * model is refused by name rather than reproduced. A Sink is not refused: it
 * genuinely has no departures.
 *
 * WHAT IS NOT PORTED. `sn_build_fj_sync_map` and `mmap_max` have no home in
 * this tree yet -- the first belongs in `line/api/fj/`, the second in
 * `line/api/mam/` next to the rest of the M3A MMAP algebra -- and the FJ
 * variant is unusable without both. They are transcribed here under their
 * MATLAB names, exactly as `npfqn_traffic_split_cs.h` inlines `mmap_normalize`
 * for the same reason, and should be lifted into their own headers when one is
 * added. Neither is re-derived: both are line-by-line transcriptions.
 *
 * ARITHMETIC. The merge, the split and the synchronization are Kronecker
 * algebra and are exact in any field. Compression is not: `mmap_compress` fits
 * an APH(2) and needs square roots, so every call to it is behind the same
 * compile-time gate `npfqn_traffic_merge` uses, and the exact instantiation
 * refuses at run time only when compression is actually reached.
 */

#include <cstddef>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/api/mam/mmap_compress.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/api/mc/dtmc_stochcomp.h"
#include "line/api/npfqn/npfqn_traffic_merge.h"
#include "line/api/npfqn/npfqn_traffic_split_cs.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/mam/mam_types.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/** `DEP{i,r}`, the departure process of class r from i in (D0,D1) form. */
template <class T>
using DepTable = std::vector<std::vector<Map<T>>>;

/**
 * The fields of `options.config` the traffic step reads.
 *
 * `merge` carries MATLAB's `config.merge` and `config.compress` together, since
 * `npfqn_traffic_merge` already owns both.
 */
struct TrafficConfig {
    std::size_t space_max = 128;
    /** `config.fj_sync_q_len`, the join's synchronization queue; MATLAB's default. */
    std::size_t fj_sync_q_len = 2;
    npfqn::MergeConfig merge;
};

/** The traffic step's view of `SolverOptions('MAM')`. */
inline TrafficConfig traffic_config(const MamOptions& opt) {
    TrafficConfig c;
    c.space_max = opt.space_max;
    return c;
}

/**
 * `sn_build_fj_sync_map`: which incoming flows at a Join must be synchronized.
 *
 * `node_sync[j][i]` is the sync group of the flow from node i into node j, with
 * 0 meaning an independent flow. Both indices are 0-based over nodes; group ids
 * are 1-based, as in the reference.
 */
struct FjSyncMap {
    std::vector<std::vector<std::size_t>> node_sync;
    std::vector<std::size_t> fork_of_group;  ///< 1-based node index, by group
    std::vector<std::size_t> join_of_group;  ///< 1-based node index, by group
    std::size_t ngroups = 0;
};

namespace traffic_detail {

/**
 * `npfqn_traffic_split_cs` takes the MMAP as MATLAB's flat cell and the rest of
 * the port takes `mam::Mmap`, so the two representations have to meet here.
 */
template <class T>
npfqn::Mmap<T> to_flat(const Mmap<T>& m) {
    npfqn::Mmap<T> v;
    v.reserve(2 + m.classes());
    v.push_back(m.D0);
    v.push_back(m.D1);
    for (const Matrix<T>& Dc : m.Dc) v.push_back(Dc);
    return v;
}

template <class T>
Mmap<T> from_flat(const npfqn::Mmap<T>& v) {
    if (v.size() < 3) throw InputError("solver_mam_traffic: a split flow carries no marking");
    Mmap<T> m;
    m.D0 = v[0];
    m.D1 = v[1];
    m.Dc.assign(v.begin() + 2, v.end());
    return m;
}

/** See `npfqn_traffic_merge`'s compress_or_refuse: the same gate, named here. */
template <class T>
Mmap<T> compress_or_refuse(const Mmap<T>& m) {
    if constexpr (num_traits<T>::has_transcendental) {
        return mmap_compress(m, MmapCompressMethod::MixtureOrder1);
    } else {
        throw UnsupportedError(
            "solver_mam_traffic: the superposition passed config.space_max and must be compressed, "
            "which needs transcendental arithmetic (aph2_fit); raise space_max or run at double");
    }
}

/** An MMAP of `nclasses` classes that never fires, MATLAB's `{[0],[0],[0]}`. */
template <class T>
Mmap<T> silent_mmap(std::size_t nclasses) {
    const T zero = num_traits<T>::from_int(0);
    Mmap<T> m;
    m.D0 = Matrix<T>(1, 1, zero);
    m.D1 = Matrix<T>(1, 1, zero);
    m.Dc.assign(nclasses, Matrix<T>(1, 1, zero));
    return m;
}

/**
 * One cell of DEP as a single-class MMAP.
 *
 * MATLAB replaces an empty or NaN-carrying entry by the silent MMAP and
 * otherwise sets `D1^(1) = D1`, i.e. reads the renewal departure process as an
 * MMAP with one class. The NaN test is the reference's `any(any(isnan(D0)))`.
 */
template <class T>
Mmap<T> dep_cell(const DepTable<T>& DEP, std::size_t i, std::size_t r) {
    if (i >= DEP.size() || r >= DEP[i].size()) return silent_mmap<T>(1);
    const Map<T>& d = DEP[i][r];
    if (d.order() == 0) return silent_mmap<T>(1);
    for (std::size_t a = 0; a < d.D0.rows(); ++a)
        for (std::size_t b = 0; b < d.D0.cols(); ++b)
            if (npfqn::detail::num_isnan(d.D0(a, b))) return silent_mmap<T>(1);
    Mmap<T> m;
    m.D0 = d.D0;
    m.D1 = d.D1;
    m.Dc.assign(1, d.D1);
    return m;
}

/**
 * `sum(mmap_lambda(.))`, MATLAB's test for a link that carries traffic.
 *
 * A link whose D1 is identically zero is short-circuited to zero. Its rate is
 * zero by inspection -- lambda_c = theta D1^(c) e -- and computing theta would
 * mean solving a singular 0 x = 0, which is exactly the shape the silent MMAP
 * has. This is not a deviation: it is the same number, obtained without the
 * degenerate solve.
 */
template <class T>
T total_rate(const Mmap<T>& m) {
    const T zero = num_traits<T>::from_int(0);
    bool silent = true;
    for (std::size_t i = 0; i < m.D1.rows() && silent; ++i)
        for (std::size_t j = 0; j < m.D1.cols(); ++j)
            if (!(m.D1(i, j) == zero)) {
                silent = false;
                break;
            }
    if (silent) return zero;
    T s = zero;
    for (const T& v : mmap_lambda(m)) s += v;
    return s;
}

/** The node types the reference's switch has no branch for; see the header note. */
template <class T>
void reject_opaque_node(const qn::NetworkStruct<T>& sn, std::size_t a, bool fj_aware) {
    const qn::NodeType ty = sn.nodes[a].nodetype;
    const char* what = nullptr;
    switch (ty) {
        case qn::NodeType::Router: what = "Router"; break;
        case qn::NodeType::Logger: what = "Logger"; break;
        case qn::NodeType::Cache: what = "Cache"; break;
        case qn::NodeType::Place: what = "Place"; break;
        case qn::NodeType::Transition: what = "Transition"; break;
        case qn::NodeType::Region: what = "Region"; break;
        case qn::NodeType::Fork:
        case qn::NodeType::Join:
            if (fj_aware) return;
            throw UnsupportedError(
                std::string("solver_mam_traffic: node '") + sn.nodes[a].name +
                "' is a Fork or a Join, which the plain traffic step has no branch for; the "
                "fork-join variant is solver_mam_traffic_mmap, which synchronizes the join with "
                "mmap_max");
        default: return;
    }
    throw UnsupportedError(std::string("solver_mam_traffic: node '") + sn.nodes[a].name + "' is a " +
                           what +
                           " node. The reference's node switch has no branch for it, so it emits "
                           "no outgoing link and every node downstream of it is described as "
                           "receiving no traffic at all; that is a wrong traffic table with no "
                           "symptom, not a conservative approximation");
}

}  // namespace traffic_detail

/**
 * `mmap_max(MMAPa, MMAPb, k)`: the synchronization of two flows through a join
 * with a queue of length k on each side.
 *
 * The phase space is (phase of a) x (phase of b) x (lead), with `lead` running
 * over 2k+1 blocks: block 0 is "both streams matched", the ODD blocks 1..2k-1
 * are "a is ahead by 1..k" and the EVEN blocks 2..2k are "b is ahead by 1..k".
 * A stream that is k ahead is BLOCKED, which is why the two extreme diagonal
 * blocks carry only the other stream's hidden generator. An arrival is emitted
 * exactly when the lagging stream catches up, i.e. on every transition that
 * moves the lead towards 0, so the join's throughput is bounded by the slower
 * of the two inputs and falls short of it by the blocking probability.
 *
 * Transcription of matlab/lib/m3a/m3a/mmap/mmap_max.m. It belongs in
 * `line/api/mam/`; see the header note on why it is here.
 */
template <class T>
Mmap<T> mmap_max(const Mmap<T>& a, const Mmap<T>& b, std::size_t k) {
    if (k == 0) throw InputError("mmap_max: the synchronization queue must hold at least one job");
    if (a.classes() != b.classes())
        throw InputError("mmap_max: the two flows carry different numbers of classes");
    const T zero = num_traits<T>::from_int(0);
    const std::size_t na = a.order(), nbp = b.order(), n = na * nbp, nblk = 1 + 2 * k;
    const Matrix<T> Ia = eye<T>(na), Ib = eye<T>(nbp);
    const Matrix<T> A0B0 = krons(a.D0, b.D0);
    const Matrix<T> A1IB = kron(a.D1, Ib);
    const Matrix<T> IAB1 = kron(Ia, b.D1);
    const Matrix<T> IAB0 = kron(Ia, b.D0);
    const Matrix<T> A0IB = kron(a.D0, Ib);

    auto put = [n](Matrix<T>& M, std::size_t rb, std::size_t cb, const Matrix<T>& B) {
        for (std::size_t i = 0; i < B.rows(); ++i)
            for (std::size_t j = 0; j < B.cols(); ++j) M(rb * n + i, cb * n + j) = B(i, j);
    };

    Matrix<T> M0(n * nblk, n * nblk, zero);
    put(M0, 0, 0, A0B0);
    put(M0, 0, 1, A1IB);
    put(M0, 0, 2, IAB1);
    for (std::size_t blk = 1; blk + 2 <= 2 * k; ++blk) put(M0, blk, blk, A0B0);
    put(M0, 2 * k - 1, 2 * k - 1, IAB0);
    put(M0, 2 * k, 2 * k, A0IB);
    for (std::size_t i = 2; i <= k; ++i) {
        const std::size_t rb = 1 + 2 * (i - 2), cb = 3 + 2 * (i - 2);
        put(M0, rb, cb, A1IB);
        put(M0, rb + 1, cb + 1, IAB1);
    }

    // The emission matrices all share one pattern; only the two operand blocks
    // change, which is what makes sum_c Dc = D1 hold by linearity of kron.
    auto emissions = [&](const Matrix<T>& bArr, const Matrix<T>& aArr) {
        Matrix<T> M(n * nblk, n * nblk, zero);
        put(M, 1, 0, bArr);
        put(M, 2, 0, aArr);
        for (std::size_t i = 2; i <= k; ++i) {
            const std::size_t rb = 1 + 2 * (i - 1), cb = 1 + 2 * (i - 2);
            put(M, rb, cb, bArr);
            put(M, rb + 1, cb + 1, aArr);
        }
        return M;
    };

    Mmap<T> out;
    out.D0 = M0;
    out.D1 = emissions(IAB1, A1IB);
    out.Dc.reserve(a.classes());
    for (std::size_t c = 0; c < a.classes(); ++c)
        out.Dc.push_back(emissions(kron(Ia, b.Dc[c]), kron(a.Dc[c], Ib)));
    return out;
}

namespace traffic_detail {

/** The NCS bookkeeping both entry points open with. */
struct NcsIndex {
    std::vector<bool> is_ncs;             ///< by 0-based node
    std::vector<std::size_t> to_ncs;      ///< 1-based NCS index, 0 for a class switch
    std::vector<std::size_t> keep;        ///< 0-based rtnodes rows of the surviving nodes
    std::size_t inc_count = 0;
};

template <class T>
NcsIndex build_ncs_index(const qn::NetworkStruct<T>& sn) {
    const std::size_t I = sn.nof_nodes(), R = sn.nclasses;
    NcsIndex x;
    x.is_ncs.assign(I, false);
    x.to_ncs.assign(I, 0);
    for (std::size_t a = 0; a < I; ++a) {
        if (sn.nodes[a].nodetype == qn::NodeType::ClassSwitch) continue;
        x.is_ncs[a] = true;
        x.to_ncs[a] = ++x.inc_count;
        for (std::size_t r = 0; r < R; ++r) x.keep.push_back(a * R + r);
    }
    return x;
}

/** The per-station (or per-node) class superposition, bounded by space_max. */
template <class T>
Mmap<T> superpose_classes(const DepTable<T>& DEP, std::size_t row, std::size_t R,
                          const TrafficConfig& config) {
    Mmap<T> d = dep_cell(DEP, row, 0);
    for (std::size_t r = 1; r < R; ++r) {
        d = mmap_super(d, dep_cell(DEP, row, r));
        if (d.order() > config.space_max) d = compress_or_refuse(d);
    }
    return d;
}

/** Row `inc` of the complemented routing, reshaped into split_cs's (R x Inc R). */
template <class T>
Matrix<T> split_probabilities(const Matrix<T>& rtncs, std::size_t inc, std::size_t Inc,
                              std::size_t R) {
    Matrix<T> P(R, Inc * R, num_traits<T>::from_int(0));
    for (std::size_t r = 0; r < R; ++r)
        for (std::size_t jnc = 0; jnc < Inc; ++jnc)
            for (std::size_t s = 0; s < R; ++s)
                P(r, jnc * R + s) = rtncs(inc * R + r, jnc * R + s);
    return P;
}

/**
 * The reference's fallback for a node with no incoming flow, in the corrected
 * form the FJ variant writes; see the header note on the line-104 defect.
 */
template <class T>
Mmap<T> no_flow_arrival(const std::vector<std::vector<Mmap<T>>>& LINKS, std::size_t inc,
                        std::size_t R) {
    for (std::size_t jnc = 0; jnc < LINKS.size(); ++jnc)
        if (LINKS[jnc][inc].order() != 0) return LINKS[jnc][inc];
    return silent_mmap<T>(R);
}

}  // namespace traffic_detail

/**
 * `sn_build_fj_sync_map`.
 *
 * A node sits on the parallel path of a (Fork, Join) pair when the fork routes
 * to it and it routes to the join, both for at least one class; every such node
 * feeds one sync group at that join. Transcription of
 * matlab/src/api/fj/sn_build_fj_sync_map.m; it belongs in `line/api/fj/`.
 */
template <class T>
FjSyncMap sn_build_fj_sync_map(const qn::NetworkStruct<T>& sn) {
    const std::size_t I = sn.nof_nodes(), K = sn.nclasses;
    const T zero = num_traits<T>::from_int(0);
    FjSyncMap m;
    m.node_sync.assign(I, std::vector<std::size_t>(I, 0));
    for (const std::pair<std::size_t, std::size_t>& fjp : sn.fj) {
        const std::size_t f = fjp.first - 1, j = fjp.second - 1;
        const std::size_t gid = ++m.ngroups;
        m.fork_of_group.push_back(f + 1);
        m.join_of_group.push_back(j + 1);
        for (std::size_t a = 0; a < I; ++a) {
            if (a == f || a == j) continue;
            bool from_fork = false, to_join = false;
            for (std::size_t k = 0; k < K; ++k) {
                if (sn.rtnodes(f * K + k, a * K + k) > zero) from_fork = true;
                if (sn.rtnodes(a * K + k, j * K + k) > zero) to_join = true;
            }
            if (from_fork && to_join) m.node_sync[j][a] = gid;
        }
    }
    return m;
}

/**
 * Port of `solver_mam_traffic.m`.
 *
 * @param DEP    station-indexed: `DEP[ist][r]` is the class-r departure process
 *               of station `ist` (0-based), an empty Map standing for MATLAB's
 *               `[]`, "this class does not depart from here"
 * @param sn the refreshed network struct
 * @param config the `options.config` fields the traffic analyzer reads
 * @return       node-indexed arrival descriptors of length `nnodes`; an entry
 *               of order 0 is MATLAB's `[]`
 */
template <class T>
std::vector<Mmap<T>> solver_mam_traffic(const qn::NetworkStruct<T>& sn, const DepTable<T>& DEP,
                                        const TrafficConfig& config) {
    namespace td = traffic_detail;
    const std::size_t I = sn.nof_nodes(), R = sn.nclasses;
    if (R == 0) throw InputError("solver_mam_traffic: the model has no classes");

    const td::NcsIndex x = td::build_ncs_index(sn);
    const std::size_t Inc = x.inc_count;
    const Matrix<T> rtncs = mc::dtmc_stochcomp(sn.rtnodes, x.keep);

    // Outgoing flows: superpose each station's classes, then split along the
    // routing. Only the branches of the reference's switch produce links.
    std::vector<std::vector<Mmap<T>>> LINKS(Inc, std::vector<Mmap<T>>(Inc));
    for (std::size_t a = 0; a < I; ++a) {
        if (!x.is_ncs[a]) continue;
        const qn::NodeType ty = sn.nodes[a].nodetype;
        if (ty != qn::NodeType::Source && ty != qn::NodeType::Delay && ty != qn::NodeType::Queue) {
            td::reject_opaque_node(sn, a, false);
            continue;
        }
        const std::size_t inc = x.to_ncs[a] - 1;
        const std::size_t ist = sn.nodes[a].station;
        if (ist == 0)
            throw InputError("solver_mam_traffic: node '" + sn.nodes[a].name +
                             "' serves jobs but carries no station index");
        const Mmap<T> dep = td::superpose_classes(DEP, ist - 1, R, config);
        const Matrix<T> Psplit = td::split_probabilities(rtncs, inc, Inc, R);
        const std::vector<npfqn::Mmap<T>> F =
            npfqn::npfqn_traffic_split_cs(td::to_flat(dep), Psplit);
        for (std::size_t jnc = 0; jnc < Inc; ++jnc)
            LINKS[inc][jnc] = mmap_normalize(td::from_flat(F[jnc]));
    }

    // Incoming flows: superpose the links that carry traffic into each node.
    std::vector<Mmap<T>> ARV(I);
    for (std::size_t a = 0; a < I; ++a) {
        if (!x.is_ncs[a] || sn.nodes[a].nodetype == qn::NodeType::Source) continue;
        const std::size_t inc = x.to_ncs[a] - 1;
        std::vector<Mmap<T>> flows;
        for (std::size_t jnc = 0; jnc < Inc; ++jnc) {
            const Mmap<T>& lk = LINKS[jnc][inc];
            if (lk.order() == 0) continue;
            if (!(num_traits<T>::to_double(td::total_rate(lk)) > lang::GlobalConstants::FineTol))
                continue;
            flows.push_back(lk);
        }
        if (flows.size() > 1) ARV[a] = npfqn::npfqn_traffic_merge(flows, config.merge);
        else if (flows.size() == 1) ARV[a] = flows[0];
        else ARV[a] = td::no_flow_arrival(LINKS, inc, R);
    }
    return ARV;
}

/**
 * Port of `solver_mam_traffic_mmap.m`, the fork-join aware traffic step.
 *
 * It differs from `solver_mam_traffic` in three places: DEP is NODE-indexed,
 * Fork and Join nodes produce outgoing links, and the flows arriving at a join
 * along one sync group are combined with `mmap_max` rather than superposed --
 * a join fires when its slowest branch delivers, which superposition, being the
 * union of the two point processes, would get badly wrong.
 *
 * @param DEP        node-indexed: `DEP[ind][r]`, 0-based over nodes
 * @param fjSyncMap  from `sn_build_fj_sync_map`
 * @param sn the refreshed network struct
 * @param config the `options.config` fields the traffic analyzer reads
 */
template <class T>
std::vector<Mmap<T>> solver_mam_traffic_mmap(const qn::NetworkStruct<T>& sn, const DepTable<T>& DEP,
                                             const TrafficConfig& config,
                                             const FjSyncMap& fjSyncMap) {
    namespace td = traffic_detail;
    const std::size_t I = sn.nof_nodes(), R = sn.nclasses;
    if (R == 0) throw InputError("solver_mam_traffic_mmap: the model has no classes");
    if (fjSyncMap.node_sync.size() != I)
        throw InputError("solver_mam_traffic_mmap: the sync map is not indexed over the nodes");

    const td::NcsIndex x = td::build_ncs_index(sn);
    const std::size_t Inc = x.inc_count;
    const Matrix<T> rtncs = mc::dtmc_stochcomp(sn.rtnodes, x.keep);

    std::vector<std::vector<std::size_t>> syncNCS(Inc, std::vector<std::size_t>(Inc, 0));
    for (std::size_t a = 0; a < I; ++a) {
        if (!x.is_ncs[a]) continue;
        for (std::size_t b = 0; b < I; ++b)
            if (x.is_ncs[b] && fjSyncMap.node_sync[a][b] > 0)
                syncNCS[x.to_ncs[a] - 1][x.to_ncs[b] - 1] = fjSyncMap.node_sync[a][b];
    }

    std::vector<std::vector<Mmap<T>>> LINKS(Inc, std::vector<Mmap<T>>(Inc));
    for (std::size_t a = 0; a < I; ++a) {
        if (!x.is_ncs[a]) continue;
        const qn::NodeType ty = sn.nodes[a].nodetype;
        if (ty != qn::NodeType::Source && ty != qn::NodeType::Delay &&
            ty != qn::NodeType::Queue && ty != qn::NodeType::Fork && ty != qn::NodeType::Join) {
            td::reject_opaque_node(sn, a, true);
            continue;
        }
        const std::size_t inc = x.to_ncs[a] - 1;
        const Mmap<T> dep = td::superpose_classes(DEP, a, R, config);
        const Matrix<T> Psplit = td::split_probabilities(rtncs, inc, Inc, R);
        const std::vector<npfqn::Mmap<T>> F =
            npfqn::npfqn_traffic_split_cs(td::to_flat(dep), Psplit);
        for (std::size_t jnc = 0; jnc < Inc; ++jnc)
            LINKS[inc][jnc] = mmap_normalize(td::from_flat(F[jnc]));
    }

    std::vector<Mmap<T>> ARV(I);
    for (std::size_t a = 0; a < I; ++a) {
        if (!x.is_ncs[a] || sn.nodes[a].nodetype == qn::NodeType::Source) continue;
        const std::size_t inc = x.to_ncs[a] - 1;

        // Ordered by group id, which is what MATLAB's unique() delivers.
        std::map<std::size_t, std::vector<Mmap<T>>> groups;
        std::vector<Mmap<T>> independent;
        for (std::size_t jnc = 0; jnc < Inc; ++jnc) {
            const Mmap<T>& lk = LINKS[jnc][inc];
            if (lk.order() == 0) continue;
            if (!(num_traits<T>::to_double(td::total_rate(lk)) > lang::GlobalConstants::FineTol))
                continue;
            const std::size_t gid = syncNCS[inc][jnc];
            if (gid == 0) independent.push_back(lk);
            else groups[gid].push_back(lk);
        }

        std::vector<Mmap<T>> flows;
        for (const auto& g : groups) {
            Mmap<T> synced = g.second.front();
            for (std::size_t f = 1; f < g.second.size(); ++f) {
                synced = mmap_normalize(mmap_max(synced, g.second[f], config.fj_sync_q_len));
                if (synced.order() > config.space_max) synced = td::compress_or_refuse(synced);
            }
            flows.push_back(synced);
        }
        flows.insert(flows.end(), independent.begin(), independent.end());

        if (flows.size() > 1) ARV[a] = npfqn::npfqn_traffic_merge(flows, config.merge);
        else if (flows.size() == 1) ARV[a] = flows[0];
        else ARV[a] = td::no_flow_arrival(LINKS, inc, R);
    }
    return ARV;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_SOLVERS_MAM_SOLVER_MAM_TRAFFIC_H
