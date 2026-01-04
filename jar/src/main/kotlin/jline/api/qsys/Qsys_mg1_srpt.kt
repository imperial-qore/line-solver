/**
 * @file M/G/1 queueing system analysis with SRPT scheduling
 *
 * Implements analysis of M/G/1 queues with Shortest Remaining Processing Time (SRPT)
 * scheduling using the Schrage-Miller formula.
 *
 * For SRPT, the mean response time for a job of size x is:
 *   E[T(x)] = (lambda * integral_0^x t*F_bar(t)dt) / (1-rho(x))^2
 *             + integral_0^x dt/(1-rho(t))
 *
 * SRPT is "Sometimes Unfair" - fair for rho <= 1/2, unfair for higher loads.
 *
 * References:
 *   - L. E. Schrage and L. W. Miller, "The queue M/G/1 with the shortest
 *     remaining processing time discipline", Operations Research, 14:670-684, 1966.
 *   - A. Wierman and M. Harchol-Balter, "Classifying scheduling policies with
 *     respect to unfairness in an M/GI/1", SIGMETRICS 2003.
 *
 * @since LINE 3.0
 */
package jline.api.qsys

import jline.io.Ret.qsys_prio
import jline.util.matrix.Matrix
import kotlin.math.exp

/**
 * Analyzes an M/G/1 queueing system with SRPT (Shortest Remaining Processing Time) scheduling.
 *
 * For multiple classes with different service times, SRPT effectively gives preemptive
 * priority to smaller jobs. Classes are internally sorted by mean service time.
 *
 * @param lambda Matrix (column vector) of arrival rates per class.
 * @param mu     Matrix (column vector) of service rates per class.
 * @param cs     Matrix (column vector) of coefficients of variation per class.
 * @return qsys_prio containing:
 *   - W (Matrix): Vector of mean response times per class
 *   - rho (Double): Overall system utilization
 * @throws IllegalArgumentException if inputs have different lengths or non-positive values
 * @throws IllegalStateException if system is unstable (rho >= 1)
 */
fun qsys_mg1_srpt(lambda: Matrix, mu: Matrix, cs: Matrix): qsys_prio {
    val lambdaArr = lambda.toArray1D()
    val muArr = mu.toArray1D()
    val csArr = cs.toArray1D()

    require(lambdaArr.size == muArr.size && lambdaArr.size == csArr.size) {
        "lambda, mu, and cs must have the same length"
    }

    val K = lambdaArr.size

    for (i in 0 until K) {
        require(lambdaArr[i] > 0.0) { "lambda[$i] must be positive" }
        require(muArr[i] > 0.0) { "mu[$i] must be positive" }
        require(csArr[i] >= 0.0) { "cs[$i] must be non-negative" }
    }

    // Compute mean service times and sort by ascending order
    val meanService = DoubleArray(K) { i -> 1.0 / muArr[i] }
    val sortIdx = meanService.indices.sortedBy { meanService[it] }

    // Reorder parameters
    val lambdaSorted = sortIdx.map { lambdaArr[it] }.toDoubleArray()
    val muSorted = sortIdx.map { muArr[it] }.toDoubleArray()
    val csSorted = sortIdx.map { csArr[it] }.toDoubleArray()

    // Compute utilizations
    val rhoI = DoubleArray(K) { i -> lambdaSorted[i] / muSorted[i] }
    val rhoTotal = rhoI.sum()

    check(rhoTotal < 1.0) { "System is unstable: utilization rho = $rhoTotal >= 1" }

    // Check if all exponential (cs = 1)
    val isExponential = csArr.all { kotlin.math.abs(it - 1.0) < 1e-10 }

    val WSorted = if (isExponential || K == 1) {
        qsysMg1SrptExp(lambdaSorted, muSorted)
    } else {
        qsysMg1SrptGeneral(lambdaSorted, muSorted, csSorted)
    }

    // Restore original ordering
    val unsortIdx = IntArray(K)
    for (i in sortIdx.indices) {
        unsortIdx[sortIdx[i]] = i
    }
    val WArr = DoubleArray(K) { i -> WSorted[unsortIdx[i]] }

    // Compute rhohat
    var Q = 0.0
    for (i in 0 until K) {
        Q += lambdaArr[i] * WArr[i]
    }
    val rhohat = Q / (1.0 + Q)

    return qsys_prio(Matrix(WArr), rhohat)
}

/**
 * SRPT for exponential service - uses preemptive priority formula.
 */
private fun qsysMg1SrptExp(lambda: DoubleArray, mu: DoubleArray): DoubleArray {
    val K = lambda.size
    val rhoI = DoubleArray(K) { i -> lambda[i] / mu[i] }

    // Mean residual service time: E[R] = sum_i lambda_i / mu_i^2
    var E_R = 0.0
    for (i in 0 until K) {
        E_R += lambda[i] / (mu[i] * mu[i])
    }

    val W = DoubleArray(K)
    for (k in 0 until K) {
        var rhoPrev = 0.0
        for (i in 0 until k) {
            rhoPrev += rhoI[i]
        }
        val rhoCurr = rhoPrev + rhoI[k]

        val Wq = E_R / ((1.0 - rhoPrev) * (1.0 - rhoCurr))
        W[k] = Wq + 1.0 / mu[k]
    }

    return W
}

/**
 * SRPT for general service distributions using numerical integration.
 */
private fun qsysMg1SrptGeneral(lambda: DoubleArray, mu: DoubleArray, cs: DoubleArray): DoubleArray {
    val K = lambda.size
    val W = DoubleArray(K)

    for (k in 0 until K) {
        val x = 1.0 / mu[k]

        // Truncated load from smaller classes
        var rhoX = 0.0
        for (i in 0..k) {
            rhoX += lambda[i] / mu[i]
        }

        if (rhoX >= 1.0) {
            W[k] = Double.POSITIVE_INFINITY
            continue
        }

        // Numerator: sum of lambda_i / mu_i^2 for smaller classes
        var numerator = 0.0
        for (i in 0..k) {
            numerator += lambda[i] / (mu[i] * mu[i])
        }

        // Second integral approximation using quadrature
        var secondTerm = 0.0
        val nSteps = 100
        val dt = x / nSteps
        for (step in 1..nSteps) {
            val t = (step - 0.5) * dt
            var rhoT = 0.0
            for (i in 0 until K) {
                rhoT += lambda[i] / mu[i] * (1.0 - exp(-mu[i] * t) * (1.0 + mu[i] * t))
            }
            if (rhoT < 1.0) {
                secondTerm += dt / (1.0 - rhoT)
            }
        }

        W[k] = numerator / ((1.0 - rhoX) * (1.0 - rhoX)) + secondTerm
    }

    return W
}

/**
 * Queueing system M/G/1 with SRPT scheduling algorithms.
 */
@Suppress("unused")
class QsysMg1SrptAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}
