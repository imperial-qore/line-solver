package jline.api.fes

import jline.lang.Network
import jline.lang.nodes.Queue
import jline.lang.nodes.Station
import jline.util.matrix.Matrix

/**
 * Result of Flow-Equivalent Server (FES) aggregation
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
data class FESResult(
    /** New Network with FES replacing the subset */
    val fesModel: Network,

    /** Reference to the FES Queue station */
    val fesStation: Queue,

    /** Deaggregation information */
    val deaggInfo: FESDeaggInfo
)

/**
 * Information needed to deaggregate FES results back to original model
 */
data class FESDeaggInfo(
    /** Original model reference */
    val originalModel: Network,

    /** Original subset stations */
    val stationSubset: List<Station>,

    /** Original station indices in subset */
    val subsetIndices: IntArray,

    /** Original station indices in complement */
    val complementIndices: IntArray,

    /** Computed throughputs for all states (list per class) */
    val throughputTable: List<Matrix>,

    /** Per-class cutoffs used */
    val cutoffs: Matrix,

    /** Stochastic complement for subset routing */
    val stochCompSubset: Matrix,

    /** Stochastic complement for complement routing */
    val stochCompComplement: Matrix,

    /** Isolated subnetwork model */
    val isolatedModel: Network,

    /** FES node index in new model */
    val fesNodeIdx: Int
) {
    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (javaClass != other?.javaClass) return false

        other as FESDeaggInfo

        if (originalModel != other.originalModel) return false
        if (stationSubset != other.stationSubset) return false
        if (!subsetIndices.contentEquals(other.subsetIndices)) return false
        if (!complementIndices.contentEquals(other.complementIndices)) return false
        if (throughputTable != other.throughputTable) return false
        if (cutoffs != other.cutoffs) return false
        if (stochCompSubset != other.stochCompSubset) return false
        if (stochCompComplement != other.stochCompComplement) return false
        if (isolatedModel != other.isolatedModel) return false
        if (fesNodeIdx != other.fesNodeIdx) return false

        return true
    }

    override fun hashCode(): Int {
        var result = originalModel.hashCode()
        result = 31 * result + stationSubset.hashCode()
        result = 31 * result + subsetIndices.contentHashCode()
        result = 31 * result + complementIndices.contentHashCode()
        result = 31 * result + throughputTable.hashCode()
        result = 31 * result + cutoffs.hashCode()
        result = 31 * result + stochCompSubset.hashCode()
        result = 31 * result + stochCompComplement.hashCode()
        result = 31 * result + isolatedModel.hashCode()
        result = 31 * result + fesNodeIdx
        return result
    }
}
