/**
 * @file M/G/1 queueing system analysis with PSJF scheduling
 *
 * Implements analysis of M/G/1 queues with Preemptive Shortest Job First (PSJF)
 * scheduling. Under PSJF, priority is based on original job size (not remaining).
 *
 * For PSJF, the mean response time for a job of size x is:
 *   E[T(x)]^PSJF = (lambda * integral_0^x t^2*f(t)dt) / (2*(1-rho(x))^2)
 *                  + x / (1 - rho(x))
 *
 * PSJF is "Always Unfair" - some job sizes are treated unfairly under all loads.
 *
 * References:
 *   - A. Wierman and M. Harchol-Balter, "Classifying scheduling policies with
 *     respect to unfairness in an M/GI/1", SIGMETRICS 2003, Section 3.2.
 *   - L. Kleinrock, "Queueing Systems, Volume II: Computer Applications", 1976.
 *
 * @since LINE 3.0
 */
package jline.api.qsys

import jline.io.Ret.qsys_prio
import jline.util.matrix.Matrix

/**
 * Analyzes an M/G/1 queueing system with PSJF (Preemptive Shortest Job First) scheduling.
 *
 * Under PSJF, jobs with smaller original sizes always preempt jobs with larger sizes.
 * This differs from SRPT where priority is based on remaining processing time.
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
fun qsys_mg1_psjf(lambda: Matrix, mu: Matrix, cs: Matrix): qsys_prio {
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

    // Sort classes by mean service time (ascending)
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

    // Compute response times using PSJF formula
    val WSorted = DoubleArray(K)

    for (k in 0 until K) {
        val x = 1.0 / muSorted[k]  // Mean service time

        // Truncated load
        var rhoX = 0.0
        for (i in 0..k) {
            rhoX += rhoI[i]
        }

        // Truncated second moment
        var m2X = 0.0
        for (i in 0..k) {
            // E[S_i^2] = (1 + cs_i^2) / mu_i^2
            val ES2i = (1.0 + csSorted[i] * csSorted[i]) / (muSorted[i] * muSorted[i])
            m2X += lambdaSorted[i] * ES2i
        }

        if (rhoX >= 1.0) {
            WSorted[k] = Double.POSITIVE_INFINITY
        } else {
            val waitingTerm = m2X / (2.0 * (1.0 - rhoX) * (1.0 - rhoX))
            val serviceTerm = x / (1.0 - rhoX)
            WSorted[k] = waitingTerm + serviceTerm
        }
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
 * Queueing system M/G/1 with PSJF scheduling algorithms.
 */
@Suppress("unused")
class QsysMg1PsjfAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}
