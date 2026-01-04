package jline.api.sn

import jline.lang.NetworkStruct
import jline.GlobalConstants
import jline.solvers.AvgHandle
import jline.util.matrix.Matrix

/**
 * Calculates the residence times at each station from the response times.
 *
 * @param sn      the NetworkStruct object for the queueing network model
 * @param RNclass the matrix of resp times
 * @param WH      residt handles for the solver
 * @return a matrix of residence times
 */

fun snGetResidTFromRespT(sn: NetworkStruct, RNclass: Matrix, WH: AvgHandle?): Matrix {
    val Vrows = sn.visits[0]!!.numRows
    val Vcols = sn.visits[0]!!.numCols
    val Vcells = sn.visits.size
    val V = Matrix(Vrows, Vcols)
    for (i in 0..<Vrows) {
        for (j in 0..<Vcols) {
            var tmpSum = 0.0
            for (k in 0..<Vcells) {
                tmpSum += sn.visits[k]!![i, j]
            }
            V[i, j] = tmpSum
        }
    }
    val M = sn.nstations
    val K = sn.nclasses
    val WNclass = RNclass.copy()
    WNclass.zero()
    
    for (ist in 0..<M) {
        for (k in 0..<K) {
            // Check if handle is disabled (following MATLAB logic)
            val isHandleDisabled = WH?.isEmpty ?: true || 
                WH?.get(sn.stations[ist])?.get(sn.jobclasses[k])?.isDisabled ?: true
            
            if (isHandleDisabled) {
                WNclass[ist, k] = Double.NaN
            } else if (!RNclass.isEmpty && RNclass[ist, k] > 0) {
                var c = -1
                for (chain in 0..<sn.chains.numRows) {
                    if (sn.chains[chain, k] > 0) {
                        c = chain
                        break
                    }
                }
                if (RNclass[ist, k] < GlobalConstants.FineTol) {
                    WNclass[ist, k] = RNclass[ist, k]
                } else {
                    val refClass = sn.refclass[0, c].toInt()
                    if (refClass >= 0) {
                        // If there is a reference class, use this:
                        WNclass[ist, k] = RNclass[ist, k] * V[ist, k] / V[sn.refstat[k, 0].toInt(), refClass]
                    } else {
                        val Vrow = sn.refstat[k, 0].toInt()
                        var Vsum = 0.0
                        for (col in 0..<sn.inchain[c]!!.numCols) {
                            Vsum += V[Vrow, sn.inchain[c]!![0, col].toInt()]
                        }
                        WNclass[ist, k] = RNclass[ist, k] * V[ist, k] / Vsum
                    }
                }
            }
        }
    }
    for (ist in 0..<M) {
        for (k in 0..<K) {
            if (java.lang.Double.isNaN(WNclass[ist, k]) || WNclass[ist, k] < 10 * GlobalConstants.FineTol || WNclass[ist, k] < GlobalConstants.Zero) {
                WNclass[ist, k] = 0
            }
        }
    }
    return WNclass
}
/**
 * Stochastic network GetResidTFromRespT algorithms
 */
@Suppress("unused")
class SngetresidtfromresptAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}