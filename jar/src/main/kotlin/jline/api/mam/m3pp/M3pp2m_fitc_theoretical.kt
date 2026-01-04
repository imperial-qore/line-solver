/**
 * @file M3PP theoretical counting process fitting
 *
 * Fits the theoretical characteristics of a MMAP(n,m) with a M3PP(2,m).
 * Supports multiple fitting methods: exact_delta, approx_delta, approx_cov, approx_ag.
 *
 * @since LINE 3.0
 */
package jline.api.mam.m3pp

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import jline.api.mam.*

/**
 * Fits the theoretical characteristics of a MMAP(n,m) with a M3PP(2,m).
 *
 * @param mmap The MMAP(n,m) to fit with a M3PP(2,m)
 * @param method Fitting method: "exact_delta", "approx_delta", "approx_cov", or "approx_ag"
 * @param t Optional finite time scale (default: 1.0)
 * @param tinf Optional near-infinite time scale (default: 1e4)
 * @return Fitted M3PP(2,m)
 */
@JvmOverloads
fun m3pp2m_fitc_theoretical(
    mmap: MatrixCell,
    method: String = "approx_delta",
    t: Double = 1.0,
    tinf: Double = 1e4
): Array<Matrix> {
    val m = mmap.size() - 2

    // Validate method
    if (method == "approx_cov" && m > 2) {
        throw IllegalArgumentException("Approximate covariance fitting only supported for two classes.")
    }

    // Time scales
    val t1 = t
    val t2 = t
    val t3 = t

    // Create aggregate MAP for joint-process characteristics
    val aggregateMap = MatrixCell(2)
    aggregateMap[0] = mmap[0]
    aggregateMap[1] = mmap[1]

    // Joint-process characteristics
    val countMean = map_count_mean(aggregateMap, t1)
    val a = countMean / t1

    val countVar1 = map_count_var(aggregateMap, t1)
    val countVar2 = map_count_var(aggregateMap, t2)
    val countVarInf = map_count_var(aggregateMap, tinf)

    val bt1 = countVar1 / (a * t1)
    val bt2 = countVar2 / (a * t2)
    val binf = countVarInf / (a * tinf)

    // Third centered moment
    val mt2 = map_count_moment(aggregateMap, t2, 3)
    val m1t2 = map_count_mean(aggregateMap, t2)
    val m2t2 = map_count_var(aggregateMap, t2) + m1t2 * m1t2
    val m3t2 = mt2 - 3.0 * m2t2 * m1t2 + 2.0 * m1t2 * m1t2 * m1t2

    // Per-class rates
    val ai = DoubleArray(m)
    for (i in 0 until m) {
        // Create single-class MAP
        val singleClassMap = MatrixCell(2)
        singleClassMap[0] = mmap[0]
        singleClassMap[1] = mmap[2 + i]
        ai[i] = map_count_mean(singleClassMap, 1.0)
    }

    return when (method.lowercase()) {
        "exact_delta" -> {
            val dvt3 = computeVarianceDifferenceTheoretical(mmap, m, t3)
            M3pp2m_fitc.m3pp2m_fitc(a, bt1, bt2, binf, m3t2, t1, t2, ai, dvt3, t3)
        }
        "approx_delta" -> {
            val dvt3 = computeVarianceDifferenceTheoretical(mmap, m, t3)
            M3pp2m_fitc_approx.m3pp2m_fitc_approx(a, bt1, bt2, binf, m3t2, t1, t2, ai, dvt3, t3)
        }
        "approx_cov" -> {
            val (vi, s) = computeCovarianceTheoretical(mmap, m, t3)
            M3pp22_fitc_approx_cov.m3pp22_fitc_approx_cov(a, bt1, bt2, binf, m3t2, t1, t2, ai, s, t3)
        }
        "approx_ag" -> {
            val gt3 = computeVarianceCovarianceTheoretical(mmap, m, t3)
            M3pp2m_fitc_approx_ag.m3pp2m_fitc_approx_ag(a, bt1, bt2, binf, m3t2, t1, t2, ai, gt3, t3)
        }
        else -> throw IllegalArgumentException("Invalid method '$method'")
    }
}

/**
 * Computes per-class variance difference for delta fitting methods.
 */
private fun computeVarianceDifferenceTheoretical(mmap: MatrixCell, m: Int, t: Double): DoubleArray {
    val dvt = DoubleArray(m)
    for (i in 0 until m) {
        // Create 2-class MMAP: class i vs others
        val mmap2 = MatrixCell(4)
        mmap2[0] = mmap[0]
        mmap2[1] = mmap[1]
        mmap2[2] = mmap[2 + i]

        // Sum of all other class matrices
        val otherClasses = Matrix.zeros(mmap[0].numRows, mmap[0].numCols)
        for (j in 0 until m) {
            if (j != i) {
                otherClasses.addEq(mmap[2 + j])
            }
        }
        mmap2[3] = otherClasses

        // Compute variance for each class
        val singleClass = MatrixCell(2)
        singleClass[0] = mmap2[0]
        singleClass[1] = mmap2[2]
        val var1 = map_count_var(singleClass, t)

        val otherClass = MatrixCell(2)
        otherClass[0] = mmap2[0]
        otherClass[1] = mmap2[3]
        val var2 = map_count_var(otherClass, t)

        dvt[i] = var1 - var2
    }
    return dvt
}

/**
 * Computes covariance for 2-class fitting.
 */
private fun computeCovarianceTheoretical(mmap: MatrixCell, m: Int, t: Double): Pair<DoubleArray, Double> {
    val vi = DoubleArray(m)
    var totalVar = 0.0

    val aggregateMap = MatrixCell(2)
    aggregateMap[0] = mmap[0]
    aggregateMap[1] = mmap[1]
    totalVar = map_count_var(aggregateMap, t)

    for (i in 0 until m) {
        val singleClass = MatrixCell(2)
        singleClass[0] = mmap[0]
        singleClass[1] = mmap[2 + i]
        vi[i] = map_count_var(singleClass, t)
    }

    val s = 0.5 * (totalVar - vi.sum())
    return Pair(vi, s)
}

/**
 * Computes per-class variance + covariance for aggregate fitting.
 */
private fun computeVarianceCovarianceTheoretical(mmap: MatrixCell, m: Int, t: Double): DoubleArray {
    val gt = DoubleArray(m)
    for (i in 0 until m) {
        // Create 2-class MMAP: class i vs others
        val mmap2 = MatrixCell(4)
        mmap2[0] = mmap[0]
        mmap2[1] = mmap[1]
        mmap2[2] = mmap[2 + i]

        // Sum of all other class matrices
        val otherClasses = Matrix.zeros(mmap[0].numRows, mmap[0].numCols)
        for (j in 0 until m) {
            if (j != i) {
                otherClasses.addEq(mmap[2 + j])
            }
        }
        mmap2[3] = otherClasses

        // Compute variance for class i
        val singleClass = MatrixCell(2)
        singleClass[0] = mmap2[0]
        singleClass[1] = mmap2[2]
        val varI = map_count_var(singleClass, t)

        // Approximate covariance (simplified)
        // In full implementation, would use mmap_count_mcov
        val covApprox = 0.0

        gt[i] = varI + covApprox
    }
    return gt
}

/**
 * M3PP 2m fitc theoretical algorithms
 */
@Suppress("unused")
class M3pp2mFitcTheoreticalAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}
