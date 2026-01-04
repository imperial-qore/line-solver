package jline.lib.mom.solver

import org.apache.commons.math3.linear.RealMatrix

/**
 * Result container for MOM solver computations
 * 
 * @property X Throughput matrix (M×R) - throughput of each class at each station
 * @property Q Queue length matrix (M×R) - average queue length of each class at each station  
 * @property G Normalizing constant vector - normalizing constants for the queueing network
 */
data class MomSolverResult(
    val X: RealMatrix,
    val Q: RealMatrix,
    val G: DoubleArray
) {
    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (javaClass != other?.javaClass) return false

        other as MomSolverResult

        if (X != other.X) return false
        if (Q != other.Q) return false
        if (!G.contentEquals(other.G)) return false

        return true
    }

    override fun hashCode(): Int {
        var result = X.hashCode()
        result = 31 * result + Q.hashCode()
        result = 31 * result + G.contentHashCode()
        return result
    }
}