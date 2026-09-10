/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FES_VALIDATE_H
#define LINE_API_FES_VALIDATE_H

/**
 * Input validation for Flow-Equivalent Server (FES) aggregation.
 *
 * Templated port of matlab/src/api/fes/fes_validate.m.
 *
 * Checks that the network structure and the station subset admit an FES
 * aggregation:
 *   - the model is closed (Norton's theorem is a closed-network result);
 *   - the subset is non-empty and a PROPER subset of the stations;
 *   - every index in the subset is a valid station;
 *   - every station in the subset is a Queue or a Delay, so that the isolated
 *     subnetwork has a load-dependent throughput table at all;
 *   - the subset holds no duplicates.
 *
 * The reference returns (isValid, errorMsg) rather than raising, because the
 * caller (ModelAdapter.aggregateFES) reports the message to the user; the port
 * keeps that contract. THE ORDER OF THE CHECKS IS PART OF THE CONTRACT: a
 * subset that is both too large and holds a Source would report the first
 * failure, and callers key on the message.
 *
 * NOTE the closed-model test is `sn_is_open_model`, which is true only when
 * EVERY class is open. A MIXED model therefore passes this validator; that is
 * the reference's behaviour and is reproduced rather than tightened, since
 * changing it would silently refuse models the MATLAB path accepts.
 *
 * Arithmetic: EXACT-CAPABLE. Structural checks only; no arithmetic on T.
 */

#include <cstddef>
#include <set>
#include <string>
#include <vector>

#include "line/api/sn/sn_predicates.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace fes {

/** The reference's (isValid, errorMsg) pair. */
struct FesValidateResult {
    bool isValid = false;
    std::string errorMsg;
};

/** NodeType.toText, so that the refusal messages read as the reference's do. */
inline std::string fes_node_type_text(qn::NodeType t) {
    switch (t) {
        case qn::NodeType::Region: return "Region";
        case qn::NodeType::Transition: return "Transition";
        case qn::NodeType::Place: return "Place";
        case qn::NodeType::Fork: return "Fork";
        case qn::NodeType::Router: return "Router";
        case qn::NodeType::Cache: return "Cache";
        case qn::NodeType::Logger: return "Logger";
        case qn::NodeType::ClassSwitch: return "ClassSwitch";
        case qn::NodeType::Delay: return "Delay";
        case qn::NodeType::Source: return "Source";
        case qn::NodeType::Queue: return "Queue";
        case qn::NodeType::Join: return "Join";
        case qn::NodeType::Sink: return "Sink";
    }
    throw InputError("fes_validate: unrecognized node type");
}

/**
 * @param sn            the network structure
 * @param subsetIndices 1-based station indices to aggregate, matching the
 *                      reference's index base so that the messages agree
 */
template <class T>
FesValidateResult fes_validate(const qn::NetworkStruct<T>& sn,
                               const std::vector<std::size_t>& subsetIndices) {
    FesValidateResult res;

    if (subsetIndices.empty()) {
        res.errorMsg = "Station subset cannot be empty.";
        return res;
    }
    if (api::sn_is_open_model(sn)) {
        res.errorMsg =
            "FES aggregation only applies to closed queueing networks. Model contains open "
            "classes.";
        return res;
    }

    const std::size_t M = sn.nstations;
    if (subsetIndices.size() >= M) {
        res.errorMsg =
            "Cannot aggregate all stations. The subset must be a proper subset of the network.";
        return res;
    }

    for (std::size_t i = 0; i < subsetIndices.size(); ++i) {
        const std::size_t stationIdx = subsetIndices[i];
        if (stationIdx < 1 || stationIdx > M) {
            res.errorMsg = "Station index " + std::to_string(stationIdx) +
                           " is invalid. Must be an integer between 1 and " + std::to_string(M) +
                           ".";
            return res;
        }
        const std::size_t nodeIdx = sn.station_to_node.at(stationIdx - 1);
        const qn::NodeType nodeType = sn.nodes.at(nodeIdx - 1).nodetype;
        if (nodeType != qn::NodeType::Queue && nodeType != qn::NodeType::Delay) {
            res.errorMsg = "Station " + std::to_string(stationIdx) + " has type " +
                           fes_node_type_text(nodeType) +
                           ". FES aggregation only supports Queue and Delay stations.";
            return res;
        }
    }

    const std::set<std::size_t> uniq(subsetIndices.begin(), subsetIndices.end());
    if (uniq.size() != subsetIndices.size()) {
        res.errorMsg = "Station subset contains duplicate stations.";
        return res;
    }

    res.isValid = true;
    return res;
}

}  // namespace fes
}  // namespace line

#endif  // LINE_API_FES_VALIDATE_H
