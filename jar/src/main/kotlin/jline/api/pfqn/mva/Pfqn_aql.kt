/**
 * @file Aggregate Queue Length approximate MVA for closed queueing networks
 * 
 * Implements the Aggregate Queue Length (AQL) approximation method for analyzing closed
 * product-form queueing networks. Provides efficient approximate solution for multi-class
 * systems using aggregate queue length representations and gamma correction factors.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.mva

import jline.io.Ret
import jline.util.matrix.Matrix

/**
 * AQL approximate MVA for a closed queueing network model.
 *
 * @param L       - demands at all stations
 * @param N       - number of jobs for each class
 * @param Z       - think time for each class
 * @param maxiter - maximum number of iterations
 * @return approximate performance metrics
 */

fun pfqn_aql(L: Matrix, N: Matrix, Z: Matrix = Matrix(1, L.numCols), maxiter: Int): Ret.pfqnAMVA {
    return pfqn_aql(L, N, Z, 1e-7, maxiter)
}

/**
 * AQL approximate MVA for a closed queueing network model.
 *
 * @param L - demands at all stations
 * @param N - number of jobs for each class
 * @param Z - think time for each class
 * @return approximate performance metrics
 */
fun pfqn_aql(L: Matrix, N: Matrix, Z: Matrix = Matrix(1, L.numCols), tol: Double = 1e-7, maxiter: Int = 1000): Ret.pfqnAMVA {
    var Z = Z
    val M = L.numRows
    val K = L.numCols

    Z = if (Z.isEmpty) Matrix(1, K) else Z

    val QN0 = Matrix(M, K)
    val value = N.sumRows(0) / M
    for (i in 0..<M) {
        for (j in 0..<K) {
            QN0[i, j] = value
        }
    }


    val Q: MutableList<Matrix> = ArrayList(K + 1)
    val R: MutableList<Matrix> = ArrayList(K + 1)
    val X: MutableList<Matrix> = ArrayList(K + 1)
    for (i in 0..K) {
        Q.add(Matrix(M, 1))
        R.add(Matrix(M, K))
        X.add(Matrix(1, K))
    }
    val gamma = Matrix(M, K)

    for (t in 0..K) {
        //Matrix n = oner(N, t);
        for (k in 0..<M) {
            Q[t][k, 0] = QN0[k, 0]
        }
    }

    var it = 0
    while (true) {
        val Q_olditer: MutableList<Matrix> = ArrayList()
        for (q in Q) {
            Q_olditer.add(q.copy())
        }
        it++

        for (t in 0..K) {
            val n = if (t > 0) {
                Matrix.oner(N, t - 1)
            } else {
                N
            }
            for (k in 0..<M) {
                for (s in 0..<K) {
                    val RValue = L[k, s] * (1 + (n.elementSum() - 1) * (Q[t][k, 0] / n.elementSum() - gamma[k, s]))
                    R[t][k, s] = RValue
                }
            }

            for (s in 0..<K) {
                val sumR = R[t].sumCols(s)
                val XValue = n[0, s] / (Z[0, s] + sumR)
                X[t][0, s] = XValue
            }

            for (k in 0..<M) {
                var QValue = 0.0
                for (s in 0..<K) {
                    QValue += X[t][0, s] * R[t][k, s]
                }
                Q[t][k, 0] = QValue
            }
        }

        for (k in 0..<M) {
            for (s in 0..<K) {
                val gammaValue = (Q[0][k, 0] / N.elementSum()) - (Q[s + 1][k, 0] / (N.elementSum() - 1))
                gamma[k, s] = gammaValue
            }
        }

        if (Matrix.maxAbsDiff(Q_olditer[0], Q[0]) < tol || it == maxiter) {
            break
        }
    }

    val XN = X[0]
    val RN = R[0]
    val UN = Matrix(M, K)
    val QN = Matrix(M, K)
    val AN = Matrix(M, K)
    val TN = Matrix(M, K)

    for (k in 0..<M) {
        for (s in 0..<K) {
            val UNValue = XN[0, s] * L[k, s]
            UN[k, s] = UNValue
            val QNValue = UNValue * (1 + Q[s + 1][k, 0])
            QN[k, s] = QNValue
            AN[k, s] = Q[s + 1][k, 0]
            TN[k, s] = AN[k, s]
        }
    }

    val CN = Matrix(1, K)
    for (r in 0..<K) {
        CN[0, r] = N[0, r] / XN[0, r]
    }

    return Ret.pfqnAMVA(QN, UN, RN, TN, CN, XN, it)
}
/**
 * PFQN aql algorithms
 */
@Suppress("unused")
class PfqnAqlAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}