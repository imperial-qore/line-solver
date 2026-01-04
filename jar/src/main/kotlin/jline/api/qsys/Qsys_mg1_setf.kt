/**
 * @file M/G/1 queueing system analysis with SETF scheduling
 *
 * Implements analysis of M/G/1 queues with SETF (Shortest Elapsed Time First)
 * scheduling - the non-preemptive version of FB/LAS.
 *
 * Under SETF:
 * - Jobs are ordered by their attained service (elapsed processing time)
 * - The job with the least attained service has highest priority
 * - However, once a job begins service, it runs to completion (non-preemptive)
 *
 * The difference from FB/LAS is that arriving jobs with less attained service
 * must wait until the currently serving job completes, rather than preempting it.
 *
 * For SETF, the mean response time follows a modified FB formula:
 *   E[T(x)]^SETF = E[T(x)]^FB + E[R] / (1 - rho_x)
 *
 * where E[R] is the mean residual service time.
 *
 * References:
 *   - M. Nuyens and A. Wierman, "The Foreground-Background queue: A survey",
 *     Performance Evaluation, 2008.
 *   - A. Wierman and M. Harchol-Balter, "Classifying scheduling policies with
 *     respect to unfairness in an M/GI/1", SIGMETRICS 2003.
 *
 * @since LINE 3.0
 */
package jline.api.qsys

import jline.io.Ret.qsys_prio
import jline.util.matrix.Matrix
import kotlin.math.exp
import kotlin.math.min

/**
 * Analyzes an M/G/1 queueing system with SETF (non-preemptive FB) scheduling.
 *
 * SETF is the non-preemptive variant of FB/LAS. Priority is based on
 * attained service, but once a job starts it runs to completion.
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
fun qsys_mg1_setf(lambda: Matrix, mu: Matrix, cs: Matrix): qsys_prio {
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

    // Compute mean residual service time for the mixture distribution
    // E[R] = sum_i (lambda_i / lambda_total) * E[S_i^2] / (2 * E[S_i])
    val lambdaTotal = lambdaArr.sum()
    var meanResidual = 0.0
    for (i in 0 until K) {
        val pI = lambdaArr[i] / lambdaTotal
        val meanS = 1.0 / muArr[i]
        val meanS2 = (1.0 + csArr[i] * csArr[i]) / (muArr[i] * muArr[i])
        meanResidual += pI * meanS2 / (2.0 * meanS)
    }

    // Compute response times using SETF formula
    val W = DoubleArray(K)

    for (k in 0 until K) {
        val x = 1.0 / muArr[k]  // Mean service time for this class

        // rho_x = lambda * integral_0^x F_bar(t)dt
        var rhoX = 0.0
        for (i in 0 until K) {
            val integralFbar = if (kotlin.math.abs(csArr[i] - 1.0) < 1e-10) {
                // Exponential: integral_0^x exp(-mu*t)dt = (1 - exp(-mu*x)) / mu
                (1.0 - exp(-muArr[i] * x)) / muArr[i]
            } else {
                // Approximation for non-exponential
                min(x, 1.0 / muArr[i])
            }
            rhoX += lambdaArr[i] * integralFbar
        }

        // Numerator: lambda * integral_0^x t*F_bar(t)dt
        var numerator = 0.0
        for (i in 0 until K) {
            val integralTFbar = if (kotlin.math.abs(csArr[i] - 1.0) < 1e-10) {
                // Exponential: integral_0^x t*exp(-mu*t)dt = (1 - exp(-mu*x)*(1+mu*x)) / mu^2
                val muI = muArr[i]
                (1.0 - exp(-muI * x) * (1.0 + muI * x)) / (muI * muI)
            } else {
                // Approximation for non-exponential
                min(x * x / 2.0, 1.0 / (muArr[i] * muArr[i]))
            }
            numerator += lambdaArr[i] * integralTFbar
        }

        if (rhoX >= 1.0) {
            W[k] = Double.POSITIVE_INFINITY
        } else {
            // FB waiting term
            val fbWaitingTerm = numerator / ((1.0 - rhoX) * (1.0 - rhoX))
            // FB service term (slowdown)
            val fbServiceTerm = x / (1.0 - rhoX)
            // Non-preemptive penalty: residual service time adjusted
            val npPenalty = meanResidual / (1.0 - rhoX)

            W[k] = fbWaitingTerm + fbServiceTerm + npPenalty
        }
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
 * Queueing system M/G/1 with SETF (non-preemptive FB) scheduling algorithms.
 */
@Suppress("unused")
class QsysMg1SetfAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}
