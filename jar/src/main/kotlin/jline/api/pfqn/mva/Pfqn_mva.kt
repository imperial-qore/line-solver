/**
 * Mean Value Analysis for Product-Form Queueing Networks
 *
 * Implements the exact MVA algorithm for closed product-form networks with load-independent
 * stations. Computes exact performance measures including throughputs, queue lengths,
 * residence times, and utilizations using the MVA recursion.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva

import jline.api.pfqn.pfqn_unique
import jline.api.pfqn.pfqn_expand
import jline.api.pfqn.pfqn_combine_mi
import jline.io.Ret
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath

/**
 * Mean Value Analysis (MVA) Algorithm for closed Product-Form Queueing Networks. Exact solution is computed for
 * several performance measures.
 *
 * @param L  - service demand matrix
 * @param N  - population vector
 * @param Z  - think times
 * @param mi - multiplicity vector
 * @return - mean value of the computed performance measures, including residence times, throughputs, number
 * of jobs at a specific node, and utilizations.
 */
fun pfqn_mva(L: Matrix, N: Matrix, Z: Matrix = Matrix(1, N.numCols), mi: Matrix? = null): Ret.pfqnMVA {
    var N = N.copy()  // Create local copy to avoid modifying the original
    var Z = if (Z.isEmpty) Matrix(1, L.numCols) else Z
    var mi = mi
    val XN: Matrix // throughputs
    var QN: Matrix // queue lengths
    var UN: Matrix // utilizations
    var CN: Matrix // residence times
    var lGN = 0.0 // log of the normalizing constant
    var InfServ = 1.0
    if (Z.isEmpty && mi == null) {
        InfServ = 0.0
    }
    N = N.ceil()
    val M_original = L.numRows // Original number of stations
    val R = L.numCols // R Classes
    N = N.columnMajorOrder().transpose()

    // Normalize mi to row vector
    if (mi == null) {
        mi = Matrix.ones(1, M_original)
    } else {
        if (mi.numRows > 1) {
            mi = mi.transpose()
        }
    }

    // Detect and consolidate replicated stations
    val uniqueResult = pfqn_unique(L)
    val L_reduced = uniqueResult.L_unique
    val mapping = uniqueResult.mapping
    val M = L_reduced.numRows

    // Combine user-provided mi with detected multiplicity
    mi = pfqn_combine_mi(mi, mapping, M)
    if (!N.any()) {
        return Ret.pfqnMVA(Matrix(0, 0), Matrix(0, 0), Matrix(0, 0), Matrix(0, 0), 0.0)
    }
    val NR = N.length()
    if (R != NR) {
        throw RuntimeException("pfqn_mva: Demand matrix and population vector have different number of classes")
    }

    XN = Matrix(1, R)
    QN = Matrix(M, R)
    CN = Matrix(M, R)
    Z = if (InfServ == 1.0) {
        Z.columnMajorOrder().transpose()
    } else {
        Matrix(1, R)
    }
    val prods = Matrix(1, R - 1)
    for (w in 0..<R - 1) {
        // ones(1, R-(w+1)+1)
        // (w + 2) instead of (w + 1) because Java indexing starts at 0
        val o = Matrix.ones(1, R - (w + 2) + 1)
        // Addition: ones(1,R-(w+1)+1) + N(1, w+1:R)
        for (i in 0..<R - (w + 2) + 1) {
            o[0, i] = o[0, i] + N[0, w + 1 + i]
        }
        // Now take prod(o)
        prods[0, w] = o.elementMult()
    }

    var firstNonEmpty = R - 1
    while (firstNonEmpty >= 0 && N[0, firstNonEmpty] == 0.0) {
        firstNonEmpty--
    }
    val totpop = Matrix.ones(1, N.numCols).add(1.0, N).elementMult()
    var ctr = totpop
    val Q = Matrix(totpop.toInt(), M)
    var currentpop = 1

    val n = Matrix(1, R)
    n[0, firstNonEmpty] = 1
    while (ctr > 0) {
        var s = 0
        while (s < R) {
            var pos_n_1s = 0
            if (n[0, s] > 0) {
                n[0, s] = n[0, s] - 1
                pos_n_1s = n[0, R - 1].toInt()
                // for w=1:R-1
                var w = 0
                while (w < R - 1) {
                    pos_n_1s = (pos_n_1s + n[0, w] * prods[0, w]).toInt()
                    w++
                }
                n[0, s] = n[0, s] + 1
            }
            var CNtot = 0.0
            var i = 0
            // Compute the residence times. Compute the total residence time as well to avoid another iteration
            // through all the stations.
            while (i < M) {
                val Lis = L_reduced[i, s]
                CN[i, s] = Lis * (mi!![0, i] + Q[pos_n_1s, i])
                CNtot += CN[i, s]
                i++
            }
            // Compute the throughput for class s
            XN[0, s] = n[0, s] / (Z[0, s] + CNtot)
            i = 0
            // Compute the queue lengths
            while (i < M) {
                QN[i, s] = XN[0, s] * CN[i, s]
                Q[currentpop, i] = Q[currentpop, i] + QN[i, s]
                i++
            }
            s++
        }
        s = R - 1
        while (s >= 0 && (n[0, s] == N[0, s]) || s > firstNonEmpty) {
            s--
        }
        val nonZero = n.find()
        if (!nonZero.isEmpty) {
            var nonZeroIdx = nonZero.numRows - 1
            while (nonZeroIdx >= 0 && n[nonZero[nonZeroIdx, 0].toInt()] <= 0) {
                nonZeroIdx--
            }
            val last_nnz = nonZero[nonZeroIdx, 0].toInt()
            var sumn = 0.0
            var sumN = 0.0
            var sumnprime = 0.0
            for (i in 0..<last_nnz) {
                sumn += n[0, i]
                sumN += N[0, i]
            }
            for (i in last_nnz + 1..<R) {
                sumnprime += n[0, i]
            }
            if (sumn == sumN && sumnprime == 0.0) {
                val logX = FastMath.log(XN[0, last_nnz])
                lGN -= logX
            }
        }
        if (s == -1) {
            break
        }
        n[0, s] = n[0, s] + 1
        s++
        while (s < R) {
            n[0, s] = 0
            s++
        }
        ctr--
        currentpop++
    }

    UN = Matrix(M, R) // Utilizations
    for (m in 0..<M) {
        for (r in 0..<R) {
            UN[m, r] = XN[0, r] * L_reduced[m, r]
        }
    }

    var WN = Matrix(M, R) // Residence times
    for (m in 0..<M) {
        for (r in 0..<R) {
            WN[m, r] = QN[m, r] / XN[0, r]
        }
    }

    // Expand results back to original dimensions if stations were consolidated
    if (M < M_original) {
        val expanded = pfqn_expand(QN, UN, WN, mapping)
        QN = expanded.first
        UN = expanded.second
        WN = expanded.third
        // Also expand CN
        val expandedCN = pfqn_expand(CN, CN, CN, mapping)
        CN = expandedCN.first
    }

    return Ret.pfqnMVA(XN, QN, UN, WN, lGN)
}
/**
 * PFQN mva algorithms
 */
@Suppress("unused")
class PfqnMvaAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}