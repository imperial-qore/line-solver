/**
 * @file General-purpose MVA for mixed networks with multi-server stations
 * 
 * Provides comprehensive MVA solution for mixed queueing networks with multi-server stations.
 * Automatically selects appropriate solution method based on network characteristics,
 * supporting both load-independent and load-dependent multi-server configurations.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.mva

import jline.api.pfqn.ld.pfqn_mvald
import jline.api.pfqn.ld.pfqn_mvaldms
import jline.io.Ret
import jline.util.Maths
import jline.util.Utils
import jline.util.matrix.Matrix

/**
 * General purpose script to handle mixed Query Networks with multiserver nodes.
 *
 * @param lambda - arrival rates
 * @param L      - service demand matrix
 * @param N      - population vector
 * @param Z      - think time matrix
 * @param mi     - multiplicity vector
 * @param S      - server count for each station
 * @return - performance measures of the mixed query network.
 */

fun pfqn_mvams(lambda: Matrix, L: Matrix, N: Matrix, Z: Matrix = Matrix(1, L.numCols), mi: Matrix?, S: Matrix?): Ret.pfqnMVA {
    var Z = if (Z.isEmpty) Matrix(1, L.numCols) else Z.copy()  // Create local copy to avoid modifying the original
    var mi = mi?.copy()  // Create local copy to avoid modifying the original
    var S = S?.copy()  // Create local copy to avoid modifying the original
    val M = L.numRows
    val R = L.numCols
    var Ntot = 0.0
    var NhasInf = false
    for (i in 0..<N.numRows) {
        for (j in 0..<N.numCols) {
            val num = N[i, j]
            if (java.lang.Double.isFinite(num)) {
                Ntot += num
            } else if (Utils.isInf(num)) {
                NhasInf = true
            }
        }
    }
    val mu = Matrix.ones(M, Ntot.toInt())
    if (S == null) {
        S = Matrix.ones(M, 1)
    }
    if (mi == null) {
        mi = Matrix.ones(M, 1)
    }
    // Z is already handled at the beginning of the function
    for (ist in 0..<M) {
        var j = 0
        while (j < Ntot) {
            mu[ist, j] = Maths.min((j + 1).toDouble(), S!![ist])
            j++
        }
    }
    var maxS = 0.0
    var initialisedMax = false
    for (i in 0..<S!!.numRows) {
        for (j in 0..<S.numCols) {
            val num = S[i, j]
            if (java.lang.Double.isFinite(num)) {
                if (!initialisedMax || num > maxS) {
                    initialisedMax = true
                    maxS = num
                }
            }
        }
    }
    val returnObject: Ret.pfqnMVA
    if (maxS == 1.0) {
        // No multi-server nodes
        if (NhasInf) {
            // open or mixed model
            val retMVAMX = pfqn_mvamx(lambda, L, N, Z, mi)
            returnObject = Ret.pfqnMVA(retMVAMX.X, retMVAMX.Q, retMVAMX.U, retMVAMX.R, retMVAMX.lGN)
        } else {
            // closed model
            val retMVA = pfqn_mva(L, N, Z, mi)
            returnObject = Ret.pfqnMVA(retMVA.X, retMVA.Q, retMVA.U, retMVA.R, retMVA.lGN)
        }
    } else {
        // Multi-server nodes
        if (NhasInf) {
            // open or mixed model
            if (mi!!.elementMax() == 1.0) {
                val lG = Double.NaN
                val retMVALDMS = pfqn_mvaldms(lambda, L, N, Z, S)
                returnObject = Ret.pfqnMVA(retMVALDMS.X, retMVALDMS.Q, retMVALDMS.U, retMVALDMS.R, lG)
            } else {
                throw RuntimeException("pfqn_mvams: Queue replicas not available in exact MVA for mixed models.")
            }
        } else {
            val retMVALD = pfqn_mvald(L, N, Z, mu)
            val lG = retMVALD.lG[retMVALD.lG.size - 1]
            returnObject = Ret.pfqnMVA(retMVALD.X, retMVALD.Q, retMVALD.U, retMVALD.R, lG)
        }
    }
    return returnObject
}
/**
 * PFQN mvams algorithms
 */
@Suppress("unused")
class PfqnMvamsAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}