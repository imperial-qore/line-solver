package jline.lib.fjcodes

import jline.util.Maths
import jline.util.matrix.Matrix
import kotlin.math.abs
import kotlin.math.exp

/**
 * Extract percentiles from Phase-Type distribution
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */

/**
 * Extract percentiles from a Phase-Type (PH) distribution
 *
 * Uses uniformization-based CDF computation to find response time percentiles.
 *
 * @param vector Initial probability vector (row vector)
 * @param matrix PH sub-generator matrix
 * @param pers Array of percentile levels (e.g., [0.50, 0.90, 0.95, 0.99])
 * @return Matrix with 2 columns: [percentile_level, percentile_value]
 */
fun returnPer(vector: Matrix, matrix: Matrix, pers: DoubleArray): Matrix {
    val m = matrix.getNumCols()

    // Compute mean response time: mean = -vector / Matrix
    // For PH distribution: mean = -alpha * S^(-1) * e
    val vectorDivMatrix = vector.rightMatrixDivide(matrix)
    val meanRT = -vectorDivMatrix.elementSum()

    // Uniformization parameter: c = max(-diag(Matrix))
    var c = 0.0
    for (i in 0 until m) {
        c = kotlin.math.max(c, -matrix.get(i, i))
    }

    // P_res = Matrix/c + I
    val P_res = matrix.scale(1.0 / c).add(1.0, Matrix.eye(m))

    // M = sum(vector / (I - P_res), 2)
    val ImMinusPres = Matrix.eye(m).add(-1.0, P_res)
    val M = vector.rightMatrixDivide(ImMinusPres).elementSum()

    // Compute coefficients ak using uniformization
    val a0 = vector.elementSum()
    var sum_a = a0
    var k = 0
    var vP = P_res.sumRows()  // Column vector of row sums
    val akList = mutableListOf<Double>()

    while (abs(sum_a - M) >= 1e-10) {
        k++
        val ak = vector.mult(vP).get(0, 0)
        akList.add(ak)
        sum_a += ak
        vP = P_res.mult(vP)
    }
    val K1 = k
    val ak = akList.toDoubleArray()

    // Compute percentiles
    val percentileRTs = Matrix(pers.size, 2)

    for (p in pers.indices) {
        percentileRTs.set(p, 0, pers[p])

        if (pers[p] < 1.0 - vector.elementSum()) {
            percentileRTs.set(p, 1, 0.0)
        } else {
            var maxTime = 3.0 * meanRT
            var flag = false

            while (!flag) {
                // Compute CDF at maxTime
                var pM = exp(-c * maxTime)
                var F = pM * a0

                for (ki in 0 until K1) {
                    pM = c * maxTime * pM / (ki + 1)
                    F += pM * ak[ki]
                }

                val perAchieved = 1.0 - F

                if (perAchieved < pers[p]) {
                    // Need larger time
                    maxTime += 0.5 * meanRT
                } else {
                    // Search backwards from maxTime with step 0.001
                    var tempPercentileRT = 0.0
                    var t = maxTime

                    while (t >= 0) {
                        pM = exp(-c * t)
                        F = pM * a0

                        for (ki in 0 until K1) {
                            pM = c * t * pM / (ki + 1)
                            F += pM * ak[ki]
                        }

                        val cdfAnalysis = 1.0 - F

                        if (cdfAnalysis < pers[p]) {
                            tempPercentileRT = t + 0.001
                            break
                        }

                        t -= 0.001
                    }

                    percentileRTs.set(p, 1, tempPercentileRT)
                    flag = true
                }
            }
        }
    }

    return percentileRTs
}

/**
 * Compute element-wise sum of all matrix elements
 */
private fun Matrix.elementSum(): Double {
    var sum = 0.0
    for (i in 0 until getNumRows()) {
        for (j in 0 until getNumCols()) {
            sum += get(i, j)
        }
    }
    return sum
}

/**
 * Compute row sums as column vector
 */
private fun Matrix.sumRows(): Matrix {
    val result = Matrix(getNumRows(), 1)
    for (i in 0 until getNumRows()) {
        var sum = 0.0
        for (j in 0 until getNumCols()) {
            sum += get(i, j)
        }
        result.set(i, 0, sum)
    }
    return result
}
