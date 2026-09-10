/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_LANG_QN_QN_LAYER_H
#define LINE_LANG_QN_QN_LAYER_H

/**
 * One SolverLN layer: a NetworkStruct plus the LQN annotations that say which
 * element of the layered model each station and class stands for.
 *
 * SCOPE. The queueing network itself, its refresh and its predicates live in
 * `network_struct.h` and are shared with every other solver; nothing here is
 * read by an algorithm. What a layer adds is the back-mapping SolverLN needs to
 * write a layer's metrics onto the tasks, entries, activities and calls of the
 * LayeredNetworkStruct it came from, which is the whole content of MATLAB's
 * `model.attribute` fields on a layer model.
 *
 * The shape SolverLN builds is narrow -- a client Delay, one or more replicas
 * of a server station, optionally a Source/Sink pair for open arrivals, and a
 * class per LQN element, wired by a class-switching routing matrix P{r,s}(i,j)
 * -- and `buildLayers` refuses by name when a layered model needs a construct
 * it cannot express, rather than building a network that silently omits it.
 */

#include <array>
#include <cstddef>
#include <utility>
#include <vector>

#include "line/lang/qn/network_struct.h"

namespace line {
namespace qn {

/**
 * A layer network: everything a NetworkStruct holds, plus the LQN back-mapping.
 *
 * `clientIdx` / `serverIdx` are the two stations every layer has by
 * construction (the callers' Delay and the served task's station); the `attr_*`
 * vectors pair a 1-based class index with the index of the LQN element it
 * stands for, in the LayeredNetworkStruct's own numbering.
 *
 * Under the SQUASHED layering (`flat`) a layer serves MANY elements at once, so
 * `serverIdx` is no longer the answer to "which station stands for element i":
 * `server_idx_of` is, and it is populated under both layerings so a consumer
 * can read it without knowing which one built the layer.
 */
template <class T>
class Layer : public NetworkStruct<T> {
public:
    std::size_t clientIdx = 0;  ///< 1-based station index of the client Delay, 0 = none
    std::size_t serverIdx = 0;  ///< 1-based station index of the server
    /** LQN element -> 1-based station index of its server here, 0 = not served here. */
    std::vector<std::size_t> server_idx_of;
    std::vector<std::size_t> host_stations;  ///< station indices of the processor servers
    std::vector<std::size_t> task_stations;  ///< station indices of the task servers
    bool flat = false;                       ///< true when built by the squashed layering
    std::vector<std::pair<std::size_t, std::size_t>> attr_tasks;       ///< (classIdx, tidx)
    std::vector<std::pair<std::size_t, std::size_t>> attr_entries;     ///< (classIdx, eidx)
    std::vector<std::pair<std::size_t, std::size_t>> attr_activities;  ///< (classIdx, aidx)
    /** (classIdx, cidx, callerActivity, calledEntry) */
    std::vector<std::array<std::size_t, 4>> attr_calls;
};

}  // namespace qn
}  // namespace line

#endif  // LINE_LANG_QN_QN_LAYER_H
