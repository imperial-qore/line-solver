/**
 * @file M3PP trace-based counting process fitting
 *
 * Fits a M3PP(2,m) from trace data using counting process characteristics.
 * Supports multiple fitting methods: exact_delta, approx_delta, approx_cov, approx_ag.
 *
 * @since LINE 3.0
 */
package jline.api.mam.m3pp

import jline.GlobalConstants
import jline.VerboseLevel
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import jline.api.trace.*

/** Helper to check if verbose output is enabled */
private fun isVerbose(): Boolean = GlobalConstants.getVerbose() != VerboseLevel.SILENT

/**
 * Fits a M3PP(2,m) from trace data using counting process characteristics.
 *
 * @param T Inter-arrival times as array
 * @param A Class labels as array
 * @param method Fitting method: "exact_delta", "approx_delta", "approx_cov", or "approx_ag"
 * @param t1 Optional finite time scale
 * @param tinf Optional near-infinite time scale
 * @return Fitted M3PP(2,m)
 */
@JvmOverloads
fun m3pp2m_fitc_trace(
    T: DoubleArray,
    A: IntArray,
    method: String = "approx_delta",
    t1: Double? = null,
    tinf: Double? = null
): Array<Matrix> {
    val labels = A.distinct().sorted()
    val m = labels.size

    // Validate method for 2-class restriction
    if (method == "approx_cov" && m > 2) {
        throw IllegalArgumentException("Approximate covariance fitting only supported for two classes.")
    }

    val TC = cumulativeSum(T)

    // Default time scales
    val t1Used = t1 ?: (10 * T.average())
    val tinfUsed = tinf ?: maxOf(10 * t1Used, (TC.last() - TC.first()) / 100)
    val t2 = t1Used + T.average()
    val t3 = tinfUsed // Controls approximation accuracy

    if (isVerbose()) println("Computing counting process at resolution $t1Used")
    val mNt1 = mtrace_iat2counts(T, A, t1Used)
    val mNt2 = mNt1 // Same resolution
    if (isVerbose()) println("Computing counting process at resolution $tinfUsed")
    val mNtinf = mtrace_iat2counts(T, A, tinfUsed)
    val mNt3 = mNt1

    val Nt1 = sumRows(mNt1)
    val Nt2 = sumRows(mNt2)
    val Ntinf = sumRows(mNtinf)

    // Total rate
    val a = 1.0 / T.average()

    // Per-class rates
    val ai = DoubleArray(m) { i ->
        val classLabel = labels[i]
        a * A.count { it == classLabel } / A.size.toDouble()
    }

    if (isVerbose()) println("Rate: $a")

    // Joint-process characteristics
    val bt1 = variance(Nt1) / (a * t1Used)
    val bt2 = bt1
    val binf = variance(Ntinf) / (a * tinfUsed)

    // Third centered moment
    val meanNt2 = Nt2.average()
    val meanNt2sq = Nt2.map { it * it }.average()
    val meanNt2cub = Nt2.map { it * it * it }.average()
    val m3t2 = meanNt2cub - 3 * meanNt2sq * meanNt2 + 2 * meanNt2 * meanNt2 * meanNt2

    return when (method.lowercase()) {
        "exact_delta" -> {
            val dvt3 = computeTraceDelta(mNt3, m, labels)
            M3pp2m_fitc.m3pp2m_fitc(a, bt1, bt2, binf, m3t2, t1Used, t2, ai, dvt3, t3)
        }
        "approx_delta" -> {
            val dvt3 = computeTraceDelta(mNt3, m, labels)
            M3pp2m_fitc_approx.m3pp2m_fitc_approx(a, bt1, bt2, binf, m3t2, t1Used, t2, ai, dvt3, t3)
        }
        "approx_cov" -> {
            val totalVar = variance(sumRows(mNt3))
            val vi = DoubleArray(m) { i -> variance(getColumn(mNt3, i)) }
            val s = 0.5 * (totalVar - vi.sum())
            M3pp22_fitc_approx_cov.m3pp22_fitc_approx_cov(a, bt1, bt2, binf, m3t2, t1Used, t2, ai, s, t3)
        }
        "approx_ag" -> {
            val vt3 = DoubleArray(m) { i -> variance(getColumn(mNt3, i)) }
            val st3 = computeTraceCovariance(mNt3, m)
            val gt3 = DoubleArray(m) { i -> vt3[i] + st3[i] }
            M3pp2m_fitc_approx_ag.m3pp2m_fitc_approx_ag(a, bt1, bt2, binf, m3t2, t1Used, t2, ai, gt3, t3)
        }
        else -> throw IllegalArgumentException("Invalid method '$method'")
    }
}

/**
 * Fits a M3PP(2,m) from trace data using Matrix inputs.
 *
 * @param T Inter-arrival times as Matrix
 * @param A Class labels as Matrix
 * @param method Fitting method
 * @param t1 Optional finite time scale
 * @param tinf Optional near-infinite time scale
 * @return Fitted M3PP(2,m)
 */
@JvmOverloads
fun m3pp2m_fitc_trace(
    T: Matrix,
    A: Matrix,
    method: String = "approx_delta",
    t1: Double? = null,
    tinf: Double? = null
): Array<Matrix> {
    val tArray = T.toArray1D()
    val aArray = IntArray(A.numRows * A.numCols) { i -> A.toArray1D()[i].toInt() }
    return m3pp2m_fitc_trace(tArray, aArray, method, t1, tinf)
}

/**
 * Computes per-class variance difference for delta fitting from trace.
 */
private fun computeTraceDelta(mNt: Array<DoubleArray>, m: Int, labels: List<Int>): DoubleArray {
    val dvt = DoubleArray(m)
    for (i in 0 until m) {
        val classCol = getColumn(mNt, i)
        val otherSum = DoubleArray(mNt.size) { row ->
            var sum = 0.0
            for (j in 0 until m) {
                if (j != i) sum += mNt[row][j]
            }
            sum
        }
        dvt[i] = variance(classCol) - variance(otherSum)
    }
    return dvt
}

/**
 * Computes per-class covariance with complementary classes from trace.
 */
private fun computeTraceCovariance(mNt: Array<DoubleArray>, m: Int): DoubleArray {
    val st = DoubleArray(m)
    for (i in 0 until m) {
        val classCol = getColumn(mNt, i)
        val otherSum = DoubleArray(mNt.size) { row ->
            var sum = 0.0
            for (j in 0 until m) {
                if (j != i) sum += mNt[row][j]
            }
            sum
        }
        st[i] = covariance(classCol, otherSum)
    }
    return st
}

// Helper functions

private fun cumulativeSum(arr: DoubleArray): DoubleArray {
    val result = DoubleArray(arr.size)
    var sum = 0.0
    for (i in arr.indices) {
        sum += arr[i]
        result[i] = sum
    }
    return result
}

private fun sumRows(matrix: Array<DoubleArray>): DoubleArray {
    return DoubleArray(matrix.size) { row -> matrix[row].sum() }
}

private fun getColumn(matrix: Array<DoubleArray>, col: Int): DoubleArray {
    return DoubleArray(matrix.size) { row -> matrix[row][col] }
}

private fun variance(arr: DoubleArray): Double {
    if (arr.isEmpty()) return 0.0
    val mean = arr.average()
    return arr.map { (it - mean) * (it - mean) }.average()
}

private fun covariance(x: DoubleArray, y: DoubleArray): Double {
    if (x.size != y.size || x.isEmpty()) return 0.0
    val meanX = x.average()
    val meanY = y.average()
    var cov = 0.0
    for (i in x.indices) {
        cov += (x[i] - meanX) * (y[i] - meanY)
    }
    return cov / x.size
}

/**
 * Converts inter-arrival times to counts at given resolution.
 *
 * @param T Inter-arrival times
 * @param A Class labels
 * @param scale Time resolution
 * @return Matrix of counts per time bin per class
 */
fun mtrace_iat2counts(T: DoubleArray, A: IntArray, scale: Double): Array<DoubleArray> {
    val labels = A.distinct().sorted()
    val m = labels.size

    // Cumulative arrival times
    val cumT = cumulativeSum(T)
    val totalTime = cumT.last()
    val numBins = maxOf(1, (totalTime / scale).toInt())

    val counts = Array(numBins) { DoubleArray(m) }

    // Count arrivals in each bin for each class
    for (i in T.indices) {
        val arrivalTime = cumT[i]
        val binIndex = minOf((arrivalTime / scale).toInt(), numBins - 1)
        val classIndex = labels.indexOf(A[i])
        if (classIndex >= 0 && binIndex >= 0) {
            counts[binIndex][classIndex]++
        }
    }

    return counts
}

/**
 * M3PP 2m fitc trace algorithms
 */
@Suppress("unused")
class M3pp2mFitcTraceAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}
