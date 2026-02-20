/**
 * @file Mean Value Analysis for mixed open-closed queueing networks
 * 
 * Implements MVA for mixed networks containing both open and closed classes without multi-server
 * stations. Combines open class analysis with exact closed-class MVA using service demand inflation
 * to account for open class interference.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.mva

import jline.io.Ret
import jline.util.Utils
import jline.util.matrix.Matrix
import java.lang.Double

/**
 * Mean Value Analysis (MVA) method for open and mixed queueing networks with no multi-server nodes.
 *
 * @param lambda - arrival rates
 * @param D      - service demand matrix
 * @param N      - population vector
 * @param Z      - think times
 * @param mi     - multiplicity vector
 * @return - the performance measures for the given network.
 */

fun pfqn_mvamx(lambda: Matrix, D: Matrix, N: Matrix, Z: Matrix = Matrix(1, D.numCols), mi: Matrix?): Ret.pfqnMVA {
    var Z = if (Z.isEmpty) Matrix(1, D.numCols) else Z.copy()  // Create local copy to avoid modifying the original
    var mi = mi?.copy()  // Create local copy to avoid modifying the original
    for (i in 0..<lambda.numCols) {
        if (lambda[i] > 0 && N[i] > 0 && Double.isFinite(N[i])) {
            throw RuntimeException("pfqn_mvamx: Arrival rate cannot be specified on closed classes.")
        }
    }
    val M = D.numRows
    val R = D.numCols
    if (mi == null) {
        mi = Matrix.ones(M, 1)
    }
    val openClasses = ArrayList<Int>()
    val closedClasses = ArrayList<Int>()
    for (i in 0..<N.length()) {
        if (Utils.isInf(N[i])) {
            openClasses.add(i)
        } else {
            closedClasses.add(i)
        }
    }
    val XN = Matrix(1, R) // Throughputs
    val UN = Matrix(M, R) // Utilizations
    val CN = Matrix(M, R) // Residence times
    val QN = Matrix(M, R) // Queue lengths

    // Compute utilizations and throughputs for the open classes
    for (r in openClasses) {
        for (ist in 0..<M) {
            UN[ist, r] = lambda[r] * D[ist, r]
        }
        XN[0, r] = lambda[r]
    }

    val UNt = UN.sumRows()
    // Z is already handled at the beginning of the function
    val Dc = Matrix(D.numRows, closedClasses.size)
    val rep = UNt.repmat(1, closedClasses.size)
    for (i in 0..<Dc.numRows) {
        var j = 0
        for (closedClass in closedClasses) {
            Dc[i, j] = D[i, closedClass] / (1 - rep[i, j])
            j++
        }
    }
    val Nclosed = Matrix(1, closedClasses.size)
    val Zclosed = Matrix(1, closedClasses.size)
    var idx = 0
    for (closedClass in closedClasses) {
        Nclosed[0, idx] = N[closedClass]
        Zclosed[0, idx] = Z[closedClass]
        idx++
    }
    // Solve the closed model consisting of M centers and just the closed classes with teh inflated service demands.
    val ret1 = pfqn_mva(Dc, Nclosed, Zclosed, mi)
    for (i in closedClasses.indices) {
        XN[closedClasses[i]] = ret1.X[i]
        for (j in 0..<QN.numRows) {
            QN[j, closedClasses[i]] = ret1.Q[j, i]
        }
        for (j in 0..<CN.numRows) {
            CN[j, closedClasses[i]] = ret1.R[j, i]
        }
    }
    for (ist in 0..<M) {
        for (r in closedClasses) {
            UN[ist, r] = XN[r] * D[ist, r]
        }
    }
    // Compute the residence times and the queue lengths for the open classes.
    for (ist in 0..<M) {
        for (r in openClasses) {
            if (ret1.Q.isEmpty) {
                CN[ist, r] = D[ist, r] / (1 - UNt[ist])
            } else {
                CN[ist, r] = D[ist, r] * (1 + ret1.Q.sumRows(ist)) / (1 - UNt[ist])
            }
            QN[ist, r] = CN[ist, r] * XN[r]
        }
    }
    return Ret.pfqnMVA(XN, QN, UN, CN, ret1.lGN)
}
/**
 * PFQN mvamx algorithms
 */
@Suppress("unused")
class PfqnMvamxAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}