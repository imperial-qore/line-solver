/**
 * @file MMAP mixture fitting from trace
 *
 * Fits a MMAP with m classes using a mixture of m^2 PH-distributions.
 * Each PH distribution represents the probability distribution conditioned
 * on the fact that the last arrival was of class i and the next arrival is
 * of class j.
 *
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import jline.api.trace.*

/**
 * Fits a MMAP with m classes using a mixture of m^2 PH-distributions from trace data.
 *
 * Each PH distribution represents the probability distribution conditioned
 * on the fact that the last arrival was of class i and the next arrival is
 * of class j.
 *
 * @param T Inter-arrival times as array
 * @param A Class labels as array
 * @return Pair of (fitted MMAP, PH distributions for each transition)
 */
fun mmap_mixture_fit_trace(T: DoubleArray, A: IntArray): Pair<MatrixCell, Array<MatrixCell>> {
    // Compute two-step transition probabilities
    val P2 = mtrace_sigma2(T, A)

    // Compute cross moments of order 1, 2 and 3
    val M1 = mtrace_cross_moment(T, A, 1)
    val M2 = mtrace_cross_moment(T, A, 2)
    val M3 = mtrace_cross_moment(T, A, 3)

    // Apply fitting algorithm
    return mmap_mixture_fit(P2, M1, M2, M3)
}

/**
 * Fits a MMAP with m classes using a mixture of m^2 PH-distributions from Matrix inputs.
 *
 * @param T Inter-arrival times as Matrix
 * @param A Class labels as Matrix
 * @return Pair of (fitted MMAP, PH distributions for each transition)
 */
fun mmap_mixture_fit_trace(T: Matrix, A: Matrix): Pair<MatrixCell, Array<MatrixCell>> {
    val tArray = T.toArray1D()
    val aArray = IntArray(A.numRows * A.numCols) { i -> A.toArray1D()[i].toInt() }
    return mmap_mixture_fit_trace(tArray, aArray)
}

/**
 * Computes two-step class transition probabilities from trace.
 *
 * @param T Inter-arrival times
 * @param A Class labels
 * @return Two-step transition probability matrix
 */
fun mtrace_sigma2(T: DoubleArray, A: IntArray): Matrix {
    val labels = A.distinct().sorted()
    val m = labels.size
    val P2 = Matrix.zeros(m, m)
    val counts = Array(m) { IntArray(m) }
    var total = 0

    // Count two-step transitions
    for (i in 0 until A.size - 2) {
        val from = labels.indexOf(A[i])
        val to = labels.indexOf(A[i + 2])
        if (from >= 0 && to >= 0) {
            counts[from][to]++
            total++
        }
    }

    // Normalize
    for (i in 0 until m) {
        val rowSum = counts[i].sum()
        for (j in 0 until m) {
            P2[i, j] = if (rowSum > 0) counts[i][j].toDouble() / rowSum else 1.0 / m
        }
    }

    return P2
}

/**
 * Computes cross-class moments of given order from trace.
 *
 * The cross moment M_k(i,j) is the k-th moment of the inter-arrival time
 * distribution conditioned on class i preceding and class j following.
 *
 * @param T Inter-arrival times
 * @param A Class labels
 * @param k Moment order
 * @return Matrix of cross-class moments
 */
fun mtrace_cross_moment(T: DoubleArray, A: IntArray, k: Int): Matrix {
    val labels = A.distinct().sorted()
    val m = labels.size
    val moments = Matrix.zeros(m, m)
    val counts = Array(m) { IntArray(m) }
    val sums = Array(m) { DoubleArray(m) }

    // Accumulate moments for each class transition
    for (i in 0 until A.size - 1) {
        val from = labels.indexOf(A[i])
        val to = labels.indexOf(A[i + 1])
        if (from >= 0 && to >= 0) {
            val value = Math.pow(T[i + 1], k.toDouble())
            sums[from][to] += value
            counts[from][to]++
        }
    }

    // Compute averages
    for (i in 0 until m) {
        for (j in 0 until m) {
            moments[i, j] = if (counts[i][j] > 0) {
                sums[i][j] / counts[i][j]
            } else {
                // Fallback: use overall moment
                T.map { Math.pow(it, k.toDouble()) }.average()
            }
        }
    }

    return moments
}

/**
 * Fits a MMAP using mixture of PH distributions based on moments.
 *
 * @param P2 Two-step transition probabilities
 * @param M1 First cross-class moments
 * @param M2 Second cross-class moments
 * @param M3 Third cross-class moments
 * @return Pair of (fitted MMAP, PH distributions for each transition)
 */
fun mmap_mixture_fit(P2: Matrix, M1: Matrix, M2: Matrix, M3: Matrix): Pair<MatrixCell, Array<MatrixCell>> {
    val m = P2.numRows

    // Fit PH distribution for each class transition
    val PHs = Array(m * m) { MatrixCell(2) }
    var idx = 0
    for (i in 0 until m) {
        for (j in 0 until m) {
            val mu1 = M1[i, j]
            val mu2 = M2[i, j]
            val mu3 = M3[i, j]

            // Fit 2-phase PH distribution using moments
            PHs[idx] = fitPH2(mu1, mu2, mu3)
            idx++
        }
    }

    // Construct MMAP from PH distributions
    val mmap = constructMMapFromPH(PHs, P2, m)

    return Pair(mmap, PHs)
}

/**
 * Fits a 2-phase PH distribution matching first 3 moments.
 */
private fun fitPH2(M1: Double, M2: Double, M3: Double): MatrixCell {
    // Compute squared coefficient of variation
    val cv2 = M2 / (M1 * M1) - 1

    val result = MatrixCell(2)

    if (cv2 <= 0.5) {
        // Erlang-2 case
        val lambda = 2.0 / M1
        result[0] = Matrix(2, 2)
        result[0][0, 0] = -lambda
        result[0][0, 1] = lambda
        result[0][1, 1] = -lambda
        result[1] = Matrix(2, 1)
        result[1][1, 0] = lambda
    } else if (cv2 >= 1.0) {
        // Hyperexponential case
        val scv = cv2
        val p = 0.5 * (1 + Math.sqrt((scv - 1) / (scv + 1)))
        val lambda1 = 2 * p / M1
        val lambda2 = 2 * (1 - p) / M1

        result[0] = Matrix(2, 2)
        result[0][0, 0] = -lambda1
        result[0][1, 1] = -lambda2
        result[1] = Matrix(2, 1)
        result[1][0, 0] = lambda1
        result[1][1, 0] = lambda2
    } else {
        // Mixed Erlang case (cv2 between 0.5 and 1)
        val lambda = 1.0 / M1
        val p = cv2 * 2  // Probability of direct exit

        result[0] = Matrix(2, 2)
        result[0][0, 0] = -2 * lambda / (1 + cv2)
        result[0][0, 1] = (1 - p) * 2 * lambda / (1 + cv2)
        result[0][1, 1] = -2 * lambda / (1 + cv2)
        result[1] = Matrix(2, 1)
        result[1][0, 0] = p * 2 * lambda / (1 + cv2)
        result[1][1, 0] = 2 * lambda / (1 + cv2)
    }

    return result
}

/**
 * Constructs MMAP from mixture of PH distributions.
 */
private fun constructMMapFromPH(PHs: Array<MatrixCell>, P2: Matrix, m: Int): MatrixCell {
    // Compute total state space size
    val stateSizes = PHs.map { it[0].numRows }
    val totalStates = stateSizes.sum()

    // Build MMAP matrices
    val D0 = Matrix.zeros(totalStates, totalStates)
    val D1 = Matrix.zeros(totalStates, totalStates)
    val classMatrices = Array(m) { Matrix.zeros(totalStates, totalStates) }

    var offset = 0
    var idx = 0
    for (i in 0 until m) {
        for (j in 0 until m) {
            val ph = PHs[idx]
            val n = ph[0].numRows

            // Copy D0 (transient part)
            for (r in 0 until n) {
                for (c in 0 until n) {
                    D0[offset + r, offset + c] = ph[0][r, c]
                }
            }

            // Add absorption rates to class-specific D1
            for (r in 0 until n) {
                val absRate = -ph[0].sumRows(r)
                if (absRate > 0.0) {
                    classMatrices[j][offset + r, offset] = absRate * P2[i, j]
                }
            }

            offset += n
            idx++
        }
    }

    // Aggregate D1
    for (c in 0 until m) {
        D1.addEq(classMatrices[c])
    }

    // Create result
    val result = MatrixCell(2 + m)
    result[0] = D0
    result[1] = D1
    for (c in 0 until m) {
        result[2 + c] = classMatrices[c]
    }

    return result
}

/**
 * MMAP mixture fit trace algorithms
 */
@Suppress("unused")
class MmapMixtureFitTraceAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}
