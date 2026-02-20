/**
 * @file Bard-Schweitzer approximate Mean Value Analysis with priority support
 * 
 * Implements the classic Bard-Schweitzer approximate MVA algorithm for closed queueing networks
 * with optional weighted priority extensions. Provides efficient approximation for large
 * multi-class networks with support for both FCFS and PS scheduling disciplines.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.mva

import jline.io.Ret
import jline.lang.constant.SchedStrategy
import jline.util.Maths
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath
import java.util.*

/**
 * Bard-Schweitzer approximate mean value analysis algorithm
 */

fun pfqn_bs(L: Matrix, N: Matrix): Ret.pfqnAMVA {
    val Z = Matrix(N.numRows, N.numCols)
    return pfqn_bs(L, N, Z)
}

/**
 * Bard-Schweitzer approximate mean value analysis algorithm
 */
@JvmOverloads

fun pfqn_bs(L: Matrix,
            N: Matrix,
            Z: Matrix,
            tol: Double = 1.0e-6,
            maxiter: Int = 1000,
            QN0: Matrix? = null): Ret.pfqnAMVA {
    var QN0 = QN0
    val M = L.numRows
    if (QN0 == null || QN0.isEmpty) {
        QN0 = N.repmat(M, 1)
        for (i in 0..<QN0.numRows) {
            for (j in 0..<QN0.numCols) {
                QN0[i, j] = QN0[i, j] / M
            }
        }
    }
    val type = arrayOfNulls<SchedStrategy>(M)
    Arrays.fill(type, SchedStrategy.PS)
    return pfqn_bs(L, N, Z, tol, maxiter, QN0, type)
}

/**
 * Bard-Schweitzer approximate mean value analysis algorithm with weighted priorities
 *
 * @param L       - the service demand matrix
 * @param N       - the population vector
 * @param Z       - the think times vector
 * @param tol     - max tolerance admitted between successive iterations
 * @param maxiter - maximum number of iterations
 * @param QN0     - original queue lengths
 * @param weight  - weight matrix for priorities (MxR)
 * @return - the performance metrics for this network.
 */

fun pfqn_bs(L: Matrix,
            N: Matrix,
            Z: Matrix,
            tol: Double,
            maxiter: Int,
            QN0: Matrix?,
            weight: Matrix?): Ret.pfqnAMVA {
    val M = L.numRows
    val type = arrayOfNulls<SchedStrategy>(M)
    Arrays.fill(type, SchedStrategy.FCFS)
    return pfqn_bs(L, N, Z, tol, maxiter, QN0, type, weight)
}

/**
 * Bard-Schweitzer approximate mean value analysis algorithm
 *
 * @param L       - the service demand matrix
 * @param N       - the population vector
 * @param Z       - the think times vector
 * @param tol     - max tolerance admitted between successive iterations
 * @param maxiter - maximum number of iterations
 * @param QN0     - original queue lengths
 * @param type    - scheduling disciplines at each station
 * @return - the performance metrics for this network.
 */

fun pfqn_bs(L: Matrix,
            N: Matrix,
            Z: Matrix,
            tol: Double,
            maxiter: Int,
            QN0: Matrix?,
            type: Array<SchedStrategy?>): Ret.pfqnAMVA {
    return pfqn_bs(L, N, Z, tol, maxiter, QN0, type, null)
}

/**
 * Bard-Schweitzer approximate mean value analysis algorithm with optional weighted priorities
 *
 * @param L       - the service demand matrix
 * @param N       - the population vector
 * @param Z       - the think times vector
 * @param tol     - max tolerance admitted between successive iterations
 * @param maxiter - maximum number of iterations
 * @param QN0     - original queue lengths
 * @param type    - scheduling disciplines at each station
 * @param weight  - optional weight matrix for priorities (MxR)
 * @return - the performance metrics for this network.
 */

fun pfqn_bs(L: Matrix,
            N: Matrix,
            Z: Matrix,
            tol: Double,
            maxiter: Int,
            QN0: Matrix?,
            type: Array<SchedStrategy?>,
            weight: Matrix?): Ret.pfqnAMVA {
    var QN0 = QN0
    val M = L.numRows
    val R = L.numCols
    if (QN0 == null || QN0.isEmpty) {
        QN0 = N.repmat(M, 1)
        for (i in 0..<QN0.numRows) {
            for (j in 0..<QN0.numCols) {
                QN0[i, j] = QN0[i, j] / M
            }
        }
    } else {
        // Add small epsilon to avoid zero problems as in MATLAB
        for (i in 0..<QN0.numRows) {
            for (j in 0..<QN0.numCols) {
                QN0[i, j] = QN0[i, j] + FastMath.ulp(1.0)
            }
        }
    }
    
    // Initialize weight matrix if not provided
    val weightMatrix = weight ?: Matrix.ones(M, R)
    
    val CN = Matrix(M, R)
    val QN = QN0!!
    val XN = Matrix(1, R)
    val UN = Matrix(M, R)
    val relprio = Matrix(M, R)
    
    var it = 1
    while (it <= maxiter) {
        val QN_1 = Matrix(QN)
        
        // Calculate relative priorities if using weighted FCFS
        val useWeightedFCFS = weight != null
        if (useWeightedFCFS) {
            for (ist in 0..<M) {
                for (r in 0..<R) {
                    relprio[ist, r] = QN[ist, r] * weightMatrix[ist, r]
                }
            }
        }
        
        for (r in 0..<R) {
            for (ist in 0..<M) {
                CN[ist, r] = L[ist, r]
                if (L[ist, r] == 0.0) {
                    // 0 service demand at this station => this class does not visit the current node
                    continue
                }
                for (s in 0..<R) {
                    if (s != r) {
                        if (type[ist] == SchedStrategy.FCFS && useWeightedFCFS) {
                            // Weighted FCFS approximation
                            CN[ist, r] = CN[ist, r] + L[ist, s] * QN[ist, s] * relprio[ist, s] / relprio[ist, r]
                        } else if (type[ist] == SchedStrategy.FCFS) {
                            // Standard FCFS approximation
                            CN[ist, r] = CN[ist, r] + L[ist, s] * QN[ist, s]
                        } else {
                            // PS approximation
                            CN[ist, r] = CN[ist, r] + L[ist, r] * QN[ist, s]
                        }
                    } else {
                        if (type[ist] == SchedStrategy.FCFS && useWeightedFCFS) {
                            // Note: in MATLAB there's relprio(ist,s)/relprio(ist,r) but s==r here, so it's 1
                            CN[ist, r] = CN[ist, r] + L[ist, r] * QN[ist, r] * (N[r] - 1) / N[r]
                        } else {
                            CN[ist, r] = CN[ist, r] + L[ist, r] * QN[ist, r] * (N[r] - 1) / N[r]
                        }
                    }
                }
            }
            XN[r] = N[r] / (Z[r] + Matrix.extractColumn(CN, r, null).elementSum())
        }
        for (r in 0..<R) {
            for (ist in 0..<M) {
                QN[ist, r] = XN[r] * CN[ist, r]
            }
        }
        for (r in 0..<R) {
            for (ist in 0..<M) {
                UN[ist, r] = XN[r] * L[ist, r]
            }
        }
        var maxabs = Double.MIN_VALUE
        for (i in 0..<QN.numRows) {
            for (j in 0..<QN.numCols) {
                val absValue = FastMath.abs(1 - QN[i, j] / QN_1[i, j])
                maxabs = Maths.max(maxabs, absValue)
            }
        }
        if (maxabs < tol) {
            break
        }
        it++
    }
    val RN = XN.repmat(M, 1)
    for (i in 0..<RN.numRows) {
        for (j in 0..<RN.numCols) {
            RN[i, j] = QN[i, j] / RN[i, j]
        }
    }
    return Ret.pfqnAMVA(QN, UN, RN, null, CN, XN, it)
}
/**
 * PFQN bs algorithms
 */
@Suppress("unused")
class PfqnBsAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}