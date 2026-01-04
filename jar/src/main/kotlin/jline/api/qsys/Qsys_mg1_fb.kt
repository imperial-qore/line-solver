/**
 * @file M/G/1 queueing system analysis with FB/LAS scheduling
 *
 * Implements analysis of M/G/1 queues with Feedback (FB) scheduling, also known as
 * Least Attained Service (LAS) or Shortest Elapsed Time (SET).
 *
 * Under FB/LAS, the job with the least attained service (smallest age) receives
 * priority. This is an age-based policy where priority depends on how much service
 * a job has received, not its original or remaining size.
 *
 * For FB, the mean response time for a job of size x is:
 *   E[T(x)]^FB = (lambda * integral_0^x t*F_bar(t)dt) / (1-rho_x)^2 + x / (1-rho_x)
 *
 * where rho_x = lambda * integral_0^x F_bar(t)dt
 *
 * FB is "Always Unfair" but approximates SRPT for heavy-tailed distributions
 * and is practical since job sizes need not be known in advance.
 *
 * References:
 *   - A. Wierman and M. Harchol-Balter, "Classifying scheduling policies with
 *     respect to unfairness in an M/GI/1", SIGMETRICS 2003, Section 3.3.
 *   - L. Kleinrock, "Queueing Systems, Volume II: Computer Applications", 1976.
 *
 * @since LINE 3.0
 */
package jline.api.qsys

import jline.io.Ret.qsys_prio
import jline.util.matrix.Matrix
import kotlin.math.exp
import kotlin.math.min

/**
 * Analyzes an M/G/1 queueing system with FB (Feedback/LAS) scheduling.
 *
 * Under FB/LAS, the job with the least attained service gets priority.
 * This is practical when job sizes are unknown since it only requires
 * tracking how much service each job has received.
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
fun qsys_mg1_fb(lambda: Matrix, mu: Matrix, cs: Matrix): qsys_prio {
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

    // Compute response times using FB formula
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
            val waitingTerm = numerator / ((1.0 - rhoX) * (1.0 - rhoX))
            val serviceTerm = x / (1.0 - rhoX)
            W[k] = waitingTerm + serviceTerm
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
 * Queueing system M/G/1 with FB/LAS scheduling algorithms.
 */
@Suppress("unused")
class QsysMg1FbAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}
