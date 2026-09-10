/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SN_SN_FJ_NODEVISITS_MMT_H
#define LINE_API_SN_SN_FJ_NODEVISITS_MMT_H

/**
 * Post-MMT node visits of a fork-join model.
 *
 * Port of the fork block of `matlab/src/lang/@MNetwork/refreshStruct.m`
 * (twins: the `sn.fj.any() && !isFJAugmented` tail of `Network.java:refreshStruct`,
 * and `network.py:_refresh_fork_join_nodevisits`).
 *
 * WHY THE ROUTING ANSWER IS NOT THE ANSWER. `sn_refresh_visits` solves one
 * traffic equation per chain, so a job that a Fork splits into siblings is
 * counted once: its blunt fork correction leaves every visited node at 1 and a
 * Join at its in-degree. That is the visit vector of ONE token, not of the
 * work the fork actually releases. The reference recovers the rest from the MMT
 * transformation, whose auxiliary open classes ARE the siblings: for every
 * auxiliary chain the transformed layer grows, the original class's node visits
 * become
 *
 *     V_orig(:,r) <- tasksPerLink * ( V_orig(:,r) + V_aux(:,r) )
 *
 * with `V_aux` read off the auxiliary chain and the transformed layer's own
 * Source, Sink and Fork rows zeroed (they carry the auxiliary arrival, not a
 * visit of the original class). On the two-branch closed fork-join this turns
 * (Think,F,Q1,Q2,J) = (1,1,1,1,2) into (1,2,1.5,1.5,3), which is what MATLAB,
 * the JAR and native Python all report.
 *
 * `V_orig` ON THE RIGHT IS THE PRE-CORRECTION SNAPSHOT, and the reference is a
 * value-semantics MATLAB struct where that happens for free: it reads `sn` and
 * writes `self.sn`. Two forks feeding the same class therefore do not compound
 * -- the last auxiliary chain wins, over the same untouched snapshot -- so this
 * port keeps the snapshot explicitly rather than accumulating in place.
 *
 * NOT APPLIED to a struct that came out of `fj_tag` (`isfjaugmented`): there the
 * Fork nodes are already resolved into tagged classes, and the reference skips
 * the correction for exactly that reason.
 *
 * ARITHMETIC: field. One multiplication and one addition per entry.
 */

#include <cstddef>
#include <string>
#include <vector>

#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/solvers/mva/fj_mmt.h"
#include "line/util/matrix.h"

namespace line {
namespace api {

/**
 * Rewrite `sn.nodevisits` with the MMT correction. A no-op on a model with no
 * fork-join pair, on a tag-augmented struct, and whenever the transformation
 * grows no auxiliary class.
 */
template <class T>
void sn_fj_nodevisits_mmt(qn::NetworkStruct<T>& sn) {
    if (sn.fj.empty() || sn.isfjaugmented) return;
    const mva::FjMmt<T> tr = mva::fj_mmt(sn);
    if (!tr.active()) return;
    const qn::NetworkStruct<T>& V = tr.V;
    if (V.nchains <= sn.nchains) return;

    // Every write reads the PRE-correction visits; see the header.
    const std::vector<Matrix<T>> X = sn.nodevisits;
    const std::size_t I = sn.nodes.size();

    // The transformed layer keeps the base nodes and appends its own Source and
    // Sink, so a row of V maps to a base row by NAME, and a row with no base
    // counterpart contributes nothing.
    std::vector<std::size_t> vrow_to_base(V.nodes.size(), 0);  // 1-based, 0 = absent
    for (std::size_t a = 0; a < V.nodes.size(); ++a) {
        const qn::NodeType nt = V.nodes[a].nodetype;
        if (nt == qn::NodeType::Source || nt == qn::NodeType::Sink ||
            nt == qn::NodeType::Fork)
            continue;  // the reference zeroes these rows before reading V_aux
        for (std::size_t i = 0; i < I; ++i)
            if (sn.nodes[i].name == V.nodes[a].name) {
                vrow_to_base[a] = i + 1;
                break;
            }
    }

    for (std::size_t nc = sn.nchains; nc < V.nchains; ++nc) {
        const std::vector<std::size_t>& aux = V.inchain[nc];
        if (aux.empty()) continue;
        const std::size_t a0 = aux[0];
        if (a0 >= tr.fjclassmap.size() || tr.fjclassmap[a0] == 0) continue;
        const std::size_t forkIdx = tr.fjforkmap[a0];
        if (forkIdx >= tr.forks.size()) continue;
        const T fanOut = num_traits<T>::from_double(tr.forks[forkIdx].fanOut);

        // the base chain of the original class this auxiliary chain mirrors
        std::size_t origChain = V.nchains;
        for (std::size_t c = 0; c < sn.nchains; ++c)
            for (std::size_t r : sn.inchain[c])
                if (r == tr.fjclassmap[a0]) origChain = c;
        if (origChain >= sn.nchains) continue;

        for (std::size_t jaux = 0; jaux < aux.size(); ++jaux) {
            const std::size_t a = aux[jaux];
            if (a >= tr.fjclassmap.size() || tr.fjclassmap[a] == 0) continue;
            const std::size_t r = tr.fjclassmap[a];  // original class mirrored by a

            std::vector<T> vaux(I, num_traits<T>::from_int(0));
            for (std::size_t b = 0; b < V.nodes.size(); ++b)
                if (vrow_to_base[b] != 0) vaux[vrow_to_base[b] - 1] = V.nodevisits[nc](b, a - 1);

            for (std::size_t i = 0; i < I; ++i)
                sn.nodevisits[origChain](i, r - 1) = T(fanOut * (X[origChain](i, r - 1) + vaux[i]));
        }
    }
}

}  // namespace api
}  // namespace line

#endif  // LINE_API_SN_SN_FJ_NODEVISITS_MMT_H
