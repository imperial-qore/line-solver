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
    val M = sn.nstations
    val K = sn.nclasses
    val nstateful = sn.nstateful
    val Vcells = sn.visits.size

    // Create V as station-indexed matrix (M x K) by converting from stateful-indexed visits
    // sn.visits[chain] is stateful-indexed (nstateful x nclasses)
    // We need to sum visits at each station across all chains
    val V = Matrix(M, K)
    V.zero()

    for (chainIdx in 0..<Vcells) {
        val visits = sn.visits[chainIdx] ?: continue
        for (sf in 0..<nstateful) {
            // Convert stateful index to station index
            // Note: statefulToStation is a row vector (1 x nstateful), so access as [0, sf]
            // Note: NaN.toInt() returns 0 in Kotlin, so we must check for NaN first
            val stationIdxRaw = if (sn.statefulToStation != null && sf < sn.statefulToStation.numCols) {
                sn.statefulToStation[0, sf]
            } else {
                Double.NaN
            }
            // Only add if it's a valid station (skip NaN for non-station nodes like Router)
            if (!java.lang.Double.isNaN(stationIdxRaw)) {
                val stationIdx = stationIdxRaw.toInt()
                if (stationIdx >= 0 && stationIdx < M) {
                    for (r in 0..<K) {
                        if (r < visits.numCols) {
                            V[stationIdx, r] = V[stationIdx, r] + visits[sf, r]
                        }
                    }
                }
            }
        }
    }

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
                    // Java uses refclass=-1 to indicate "no specific reference class"
                    // (Java indices are 0-based, so -1 means unset)
                    // In that case, sum visits over all classes in the chain
                    if (refClass >= 0) {
                        // If there is a specific reference class (0-based index >= 0), use it:
                        val denom = V[sn.refstat[k, 0].toInt(), refClass]
                        WNclass[ist, k] = RNclass[ist, k] * V[ist, k] / denom
                    } else {
                        // No specific reference class - sum visits over all classes in the chain
                        val Vrow = sn.refstat[k, 0].toInt()
                        var Vsum = 0.0
                        for (col in 0..<sn.inchain[c]!!.numCols) {
                            val classIdx = sn.inchain[c]!![0, col].toInt()
                            Vsum += V[Vrow, classIdx]
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
