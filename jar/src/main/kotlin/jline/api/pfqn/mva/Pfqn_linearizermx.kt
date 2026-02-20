/**
 * @file Linearizer approximate MVA for mixed networks with multi-server stations
 * 
 * Implements linearizer-based approximation methods for mixed queueing networks with
 * multi-server stations. Provides multiple solution techniques including standard linearizer,
 * general-form linearizer, and extended general-form linearizer with automatic method selection.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.mva

import jline.io.Ret
import jline.lang.constant.SchedStrategy
import jline.util.Utils
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath

/**
 * Linearizer method for mixed models with multi-server stations
 *
 * @param lambda   - arrival rate of open classes
 * @param L        - the service demand matrix
 * @param N        - the population vector
 * @param Z        - the think times vector
 * @param nservers - number of servers at the stations
 * @param type     - scheduling strategy type
 * @param tol      - max tolerance admitted between successive iterations
 * @param maxiter  - maximum number of iterations
 * @param method   - solution method
 * @return - the performance metrics for this network
 */

fun pfqn_linearizermx(lambda: Matrix,
                      L: Matrix,
                      N: Matrix,
                      Z: Matrix,
                      nservers: Matrix,
                      type: Array<SchedStrategy>,
                      tol: Double,
                      maxiter: Int,
                      method: String): Ret.pfqnAMVA {
    var Z = Z
    val lambdaLocal = lambda.copy()  // Create local copy to avoid modifying the original
    val LLocal = L.copy()  // Create local copy to avoid modifying the original
    
    for (i in 0..<lambdaLocal.numCols) {
        if (lambdaLocal[i] > 0 && N[i] > 0 && java.lang.Double.isFinite(N[i])) {
            throw RuntimeException("pfqn_mvamx: Arrival rate cannot be specified on closed classes.")
        }
    }
    val M = LLocal.numRows
    val R = LLocal.numCols

    for (i in 0..<lambdaLocal.numRows) {
        for (j in 0..<lambdaLocal.numCols) {
            if (java.lang.Double.isNaN(lambdaLocal[i, j])) {
                lambdaLocal[i, j] = 0
            }
        }
    }
    for (i in 0..<LLocal.numRows) {
        for (j in 0..<LLocal.numCols) {
            if (java.lang.Double.isNaN(LLocal[i, j])) {
                LLocal[i, j] = 0
            }
        }
    }
    for (i in 0..<Z.numRows) {
        for (j in 0..<Z.numCols) {
            if (java.lang.Double.isNaN(Z[i, j])) {
                Z[i, j] = 0
            }
        }
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

    val XN = Matrix(1, R)
    val UN = Matrix(M, R)
    val WN = Matrix(M, R)
    val QN = Matrix(M, R)
    val CN = Matrix(1, R)

    for (r in openClasses) {
        for (i in 0..<M) {
            UN[i, r] = lambdaLocal[r] * LLocal[i, r]
        }
        XN[0, r] = lambdaLocal[r]
    }

    val UNt = UN.sumRows()

    Z = if (Z.isEmpty) {
        Matrix(1, R)
    } else {
        Z.sumCols()
    }
    val Dc = Matrix(LLocal.numRows, closedClasses.size)
    val rep = UNt.repmat(1, closedClasses.size)
    for (i in 0..<Dc.numRows) {
        var j = 0
        for (closedClass in closedClasses) {
            Dc[i, j] = LLocal[i, closedClass] / (1 - rep[i, j])
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

    val QNc: Matrix
    val UNc: Matrix
    val WNc: Matrix
    val CNc: Matrix
    val XNc: Matrix
    val totiter: Int

    if (nservers.elementMax() == 1.0) {
        var res: Ret.pfqnAMVA? = null
        if (method == "lin") {
            res = pfqn_linearizer(Dc, Nclosed, Zclosed, type, tol, maxiter)
        } else if (method == "gflin") {
            val linAlpha = 2.0
            res = pfqn_gflinearizer(Dc, Nclosed, Zclosed, type, tol, maxiter, linAlpha)
        } else {
            // Extended General Form Linearizer
            val alphaM = Matrix(1, Nclosed.numCols)
            /**
             * Gompertz function to approximate alpha as a function of the network population
             */
            for (i in 0..<Nclosed.numCols) {
                alphaM[i] = 0.6 + 1.4 * FastMath.exp(-8 * FastMath.exp(-0.8 * Nclosed[i]))
            }
            res = pfqn_egflinearizer(Dc, Nclosed, Zclosed, type, tol, maxiter, alphaM)
        }
        QNc = res.Q
        UNc = res.U
        WNc = res.R
        CNc = res.C
        XNc = res.X
        totiter = res.totiter
    } else {
        val typeMatrix: MutableList<SchedStrategy> = ArrayList()
        for (i in type.indices) {
            typeMatrix.add(type[i])
        }
        val res = pfqn_linearizerms(Dc, Nclosed, Zclosed, nservers, typeMatrix, tol, maxiter)
        QNc = res.Q
        UNc = res.U
        WNc = res.R
        CNc = res.C
        XNc = res.X
        totiter = res.totiter
    }

    for (i in closedClasses.indices) {
        XN[closedClasses[i]] = XNc[i]
        for (j in 0..<QN.numRows) {
            QN[j, closedClasses[i]] = QNc[j, i]
        }
        for (j in 0..<WN.numRows) {
            WN[j, closedClasses[i]] = WNc[j, i]
        }
        for (j in 0..<UN.numRows) {
            UN[j, closedClasses[i]] = UNc[j, i]
        }
        CN[closedClasses[i]] = CNc[i]
    }

    for (i in 0..<M) {
        for (r in closedClasses) {
            UN[i, r] = XN[r] * LLocal[i, r]
        }
    }

    for (i in 0..<M) {
        for (r in openClasses) {
            if (QNc.isEmpty) {
                WN[i, r] = LLocal[i, r] / (1 - UNt[i])
            } else {
                WN[i, r] = LLocal[i, r] * (1 + QNc.sumRows(i)) / (1 - UNt[i])
            }
            QN[i, r] = WN[i, r] * XN[r]
        }
    }
    for (r in openClasses) {
        CN[r] = WN.sumCols(r)
    }
    return Ret.pfqnAMVA(QN, UN, WN, null, CN, XN, totiter)
}
/**
 * PFQN linearizermx algorithms
 */
@Suppress("unused")
class PfqnLinearizermxAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}