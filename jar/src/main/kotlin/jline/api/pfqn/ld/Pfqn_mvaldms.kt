/**
 * @file Load-dependent MVA wrapper for multi-server utilization adjustment
 * 
 * Provides wrapper functionality for load-dependent Mean Value Analysis with automatic
 * multi-server utilization adjustments. Handles proper utilization scaling for multi-server
 * stations in load-dependent closed queueing networks.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.ld

import jline.io.Ret
import jline.util.Maths
import jline.util.Utils
import jline.util.matrix.Matrix
import java.lang.Double

/**
 * Wrapper for pfqn_mvaldmx that adjusts utilizations to account for multiservers
 *
 * @param lambda - arrival rates
 * @param D      - service demand matrix
 * @param N      - population vector
 * @param Z      - think times
 * @param S      - servers at each station
 * @return - the performance measures for the given network.
 */

fun pfqn_mvaldms(lambda: Matrix, D: Matrix, N: Matrix, Z: Matrix, S: Matrix): Ret.pfqnMVA {
    var Z = Z
    val M = D.numRows
    val R = D.numCols
    var Nct = 0.0
    for (i in 0..<N.numRows) {
        for (j in 0..<N.numCols) {
            val num = N[i, j]
            if (Double.isFinite(num)) {
                Nct += num
            }
        }
    }
    val mu = Matrix.ones(M, Nct.toInt())
    for (ist in 0..<M) {
        for (j in 0..<mu.numCols) {
            mu[ist, j] = Maths.min((j + 1).toDouble(), S[ist])
        }
    }
    if (Z.isEmpty) {
        Z = Matrix(1, R)
    }
    val ret1 = pfqn_mvaldmx(lambda, D, N, Z, mu, S)
    val openClasses = ArrayList<Int>()
    val closedClasses = ArrayList<Int>()
    for (i in 0..<N.length()) {
        if (Utils.isInf(N[i])) {
            openClasses.add(i)
        } else {
            closedClasses.add(i)
        }
    }
    val XN = ret1.X
    val QN = ret1.Q
    val CN = ret1.R
    val lGN = ret1.lG
    val UN = Matrix(M, R)
    for (r in closedClasses) {
        for (ist in 0..<M) {
            UN[ist, r] = XN[r] * D[ist, r] / S[ist]
        }
    }
    for (r in openClasses) {
        for (ist in 0..<M) {
            UN[ist, r] = lambda[r] * D[ist, r] / S[ist]
        }
    }
    return Ret.pfqnMVA(XN, QN, UN, CN, lGN)
}
/**
 * PFQN mvaldms algorithms
 */
@Suppress("unused")
class PfqnMvaldmsAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}