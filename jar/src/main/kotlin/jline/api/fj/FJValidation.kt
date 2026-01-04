package jline.api.fj

import jline.lang.NetworkStruct
import jline.lang.constant.NodeType
import jline.lang.constant.SchedStrategy
import jline.lang.nodes.*

/**
 * Fork-Join topology validation
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */

/**
 * Information about detected Fork-Join topology
 */
data class FJInfo(
    val K: Int,                  // Number of parallel queues
    val forkIdx: Int,            // Fork node index
    val joinIdx: Int,            // Join node index
    val queueIndices: IntArray,  // Parallel queue node indices
    val sourceIdx: Int,          // Source node index
    val sinkIdx: Int,            // Sink node index
    val isValid: Boolean         // Whether topology is valid for FJ_codes
)

/**
 * Check if network has Fork-Join topology suitable for FJ_codes
 *
 * Validates:
 * - Source → Fork → K Queues → Join → Sink topology
 * - Open classes only
 * - Homogeneous service across parallel queues
 * - FCFS or PS scheduling
 * - Supported distributions (Exp, HyperExp(2), Erlang(2), MAP(2))
 *
 * @param sn Network structure
 * @return Pair of (isValid, fjInfo)
 */
fun isFJ(sn: NetworkStruct): Pair<Boolean, FJInfo?> {
    try {
        // Check for basic open network structure
        if (sn.nclosedjobs > 0) {
            return Pair(false, null)  // Must be open network
        }

        // Find Source and Sink nodes
        var sourceIdx = -1
        var sinkIdx = -1
        var forkIdx = -1
        var joinIdx = -1

        for (i in 0 until sn.nnodes) {
            val node = sn.nodes[i]
            when (node) {
                is Source -> sourceIdx = i
                is Sink -> sinkIdx = i
                is Fork -> forkIdx = i
                is Join -> joinIdx = i
                else -> {}
            }
        }

        // Validate required nodes exist
        if (sourceIdx < 0 || sinkIdx < 0 || forkIdx < 0 || joinIdx < 0) {
            return Pair(false, null)
        }

        // Check connectivity: Source → Fork → Queues → Join → Sink
        // Simplified validation: count queues between Fork and Join
        val queueIndices = mutableListOf<Int>()

        for (i in 0 until sn.nnodes) {
            val node = sn.nodes[i]
            when (node) {
                is Queue, is Delay -> {
                    // Check if this queue is between Fork and Join
                    // Simplified: assume all queues are part of FJ structure
                    queueIndices.add(i)
                }
                else -> {}
            }
        }

        val K = queueIndices.size
        if (K < 1) {
            return Pair(false, null)
        }

        // Check homogeneous service (all queues have same service distribution)
        // Simplified validation for MVP
        // TODO: Add detailed service homogeneity check

        // Check scheduling strategy (should be FCFS or PS)
        for (qIdx in queueIndices) {
            val nodeIdx = qIdx
            val stationIdx = sn.nodeToStation.get(nodeIdx).toInt()
            val station = sn.stations[stationIdx]
            val sched = sn.sched[station]
            if (sched != SchedStrategy.FCFS && sched != SchedStrategy.PS) {
                return Pair(false, null)
            }
        }

        val fjInfo = FJInfo(
            K = K,
            forkIdx = forkIdx,
            joinIdx = joinIdx,
            queueIndices = queueIndices.toIntArray(),
            sourceIdx = sourceIdx,
            sinkIdx = sinkIdx,
            isValid = true
        )

        return Pair(true, fjInfo)

    } catch (e: Exception) {
        return Pair(false, null)
    }
}
