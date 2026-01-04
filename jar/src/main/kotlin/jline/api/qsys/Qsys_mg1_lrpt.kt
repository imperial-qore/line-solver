/**
 * @file M/G/1 queueing system analysis with LRPT scheduling
 *
 * Implements analysis of M/G/1 queues with Longest Remaining Processing Time (LRPT)
 * scheduling. Under LRPT, jobs with the longest remaining time share the processor.
 *
 * For LRPT, the expected slowdown for a job of size x is:
 *   E[S(x)]^LRPT = 1/(1-rho) + lambda*E[X^2] / (2*x*(1-rho)^2)
 *
 * Therefore the expected response time is:
 *   E[T(x)]^LRPT = x/(1-rho) + lambda*E[X^2] / (2*(1-rho)^2)
 *
 * LRPT is "Always Unfair" - E[S(y)]^LRPT > 1/(1-rho) for all finite job sizes y.
 * The slowdown is monotonically decreasing in x, converging to 1/(1-rho).
 *
 * References:
 *   - A. Wierman and M. Harchol-Balter, "Classifying scheduling policies with
 *     respect to unfairness in an M/GI/1", SIGMETRICS 2003, Section 3.2.
 *   - M. Harchol-Balter, K. Sigman, and A. Wierman, "Asymptotic convergence
 *     of scheduling policies with respect to slowdown", Performance Evaluation, 2002.
 *
 * @since LINE 3.0
 */
package jline.api.qsys

import jline.io.Ret.qsys_prio
import jline.util.matrix.Matrix

/**
 * Analyzes an M/G/1 queueing system with LRPT (Longest Remaining Processing Time) scheduling.
 *
 * Under LRPT, jobs with the longest remaining processing time share the processor
 * evenly. This prioritizes large jobs over small jobs.
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
fun qsys_mg1_lrpt(lambda: Matrix, mu: Matrix, cs: Matrix): qsys_prio {
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

    // Compute utilizations
    val rhoI = DoubleArray(K) { i -> lambdaArr[i] / muArr[i] }
    val rhoTotal = rhoI.sum()

    check(rhoTotal < 1.0) { "System is unstable: utilization rho = $rhoTotal >= 1" }

    // Compute overall second moment of service time
    // E[X^2] = sum_i (lambda_i / lambda_total) * E[S_i^2]
    val lambdaTotal = lambdaArr.sum()
    val p = DoubleArray(K) { i -> lambdaArr[i] / lambdaTotal }

    var EX2 = 0.0
    for (i in 0 until K) {
        // E[S_i^2] = (1 + cs_i^2) / mu_i^2
        val ES2i = (1.0 + csArr[i] * csArr[i]) / (muArr[i] * muArr[i])
        EX2 += p[i] * ES2i
    }

    // Compute response times using LRPT formula
    // E[T(x)] = x/(1-rho) + lambda*E[X^2] / (2*(1-rho)^2)
    val W = DoubleArray(K)

    val term2 = lambdaTotal * EX2 / (2.0 * (1.0 - rhoTotal) * (1.0 - rhoTotal))

    for (k in 0 until K) {
        val x = 1.0 / muArr[k]  // Mean service time
        val term1 = x / (1.0 - rhoTotal)
        W[k] = term1 + term2
    }

    // Compute rhohat
    var Q = 0.0
    for (i in 0 until K) {
        Q += lambdaArr[i] * W[i]
    }
    val rhohat = Q / (1.0 + Q)

    return qsys_prio(Matrix(W), rhohat)
}

/**
 * Queueing system M/G/1 with LRPT scheduling algorithms.
 */
@Suppress("unused")
class QsysMg1LrptAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}
