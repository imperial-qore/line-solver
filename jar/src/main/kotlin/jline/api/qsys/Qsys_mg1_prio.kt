/**
 * @file M/G/1 queueing system analysis with non-preemptive priorities
 *
 * Implements analysis of M/G/1 queues with Head-of-Line (non-preemptive) priority
 * scheduling. Extends the Pollaczek-Khinchine formula to handle multiple priority
 * classes using conservation laws for non-preemptive priorities.
 *
 * For K priority classes (class 1 = highest priority), the mean response time for
 * class k is:
 *   W_k = B_0 / ((1 - sum_{i=1}^{k-1} rho_i) * (1 - sum_{i=1}^{k} rho_i)) + 1/mu_k
 *
 * Where:
 *   - rho_i = lambda_i / mu_i (utilization of class i)
 *   - B_0 = sum_{i=1}^{K} lambda_i * E[S_i^2] / 2
 *   - E[S_i^2] = (1 + cs_i^2) / mu_i^2
 *   - cs_i = coefficient of variation for class i service time
 *
 * @since LINE 3.0
 */
package jline.api.qsys

import jline.io.Ret.qsys_prio
import jline.util.matrix.Matrix

/**
 * Analyzes an M/G/1 queueing system with non-preemptive (Head-of-Line) priorities.
 *
 * Computes per-class mean response times and overall utilization for a system with
 * multiple priority classes using the Pollaczek-Khinchine formula extended for
 * non-preemptive priority scheduling.
 *
 * @param lambda Matrix (column vector) of arrival rates per priority class. Class 1 = highest priority.
 * @param mu     Matrix (column vector) of service rates per priority class.
 * @param cs     Matrix (column vector) of coefficients of variation per priority class.
 * @return qsys_prio containing:
 *   - W (Matrix): Vector of mean response times per priority class
 *   - rho (Double): Overall system utilization
 * @throws IllegalArgumentException if inputs have different lengths or non-positive values
 * @throws IllegalStateException if system is unstable (rho >= 1)
 */
fun qsys_mg1_prio(lambda: Matrix, mu: Matrix, cs: Matrix): qsys_prio {
    // Convert Matrix to arrays for computation
    val lambdaArr = lambda.toArray1D()
    val muArr = mu.toArray1D()
    val csArr = cs.toArray1D()

    // Validate input lengths
    require(lambdaArr.size == muArr.size && lambdaArr.size == csArr.size) {
        "lambda, mu, and cs must have the same length: " +
                "got ${lambdaArr.size}, ${muArr.size}, ${csArr.size}"
    }

    val K = lambdaArr.size

    // Validate positive values
    for (i in 0 until K) {
        require(lambdaArr[i] > 0.0) { "lambda[$i] must be positive, got ${lambdaArr[i]}" }
        require(muArr[i] > 0.0) { "mu[$i] must be positive, got ${muArr[i]}" }
        require(csArr[i] > 0.0) { "cs[$i] must be positive, got ${csArr[i]}" }
    }

    // Compute per-class utilizations
    val rho_i = DoubleArray(K) { i -> lambdaArr[i] / muArr[i] }

    // Overall utilization
    val rho = rho_i.sum()

    // Stability check
    check(rho < 1.0) {
        "System is unstable: utilization rho = $rho >= 1"
    }

    // Compute mean second moment of service time (used in waiting time formula)
    // B_0 = sum_i lambda_i * E[S_i^2] / 2
    // where E[S_i^2] = (1 + cs_i^2) / mu_i^2
    var B_0 = 0.0
    for (i in 0 until K) {
        B_0 += lambdaArr[i] * (1.0 + csArr[i] * csArr[i]) / (muArr[i] * muArr[i])
    }
    B_0 /= 2.0

    // Compute per-class waiting times (in queue, not including service)
    val W_q_arr = DoubleArray(K)

    for (k in 0 until K) {
        // Cumulative sum of utilizations for higher priority classes (0 to k-1)
        var rho_prev = 0.0
        for (i in 0 until k) {
            rho_prev += rho_i[i]
        }

        // Cumulative sum including current class (0 to k)
        var rho_curr = rho_prev + rho_i[k]

        // Waiting time (in queue) for class k
        // W_q_k = B_0 / ((1 - rho_prev) * (1 - rho_curr))
        W_q_arr[k] = B_0 / ((1.0 - rho_prev) * (1.0 - rho_curr))
    }

    // Response time = waiting time + service time
    // W_k = W_q_k + 1/mu_k
    val W_arr = DoubleArray(K)
    for (i in 0 until K) {
        W_arr[i] = W_q_arr[i] + 1.0 / muArr[i]
    }

    // Compute rhohat = Q/(1+Q) to match qsys convention
    // Q is the mean number of customers in the system (Little's law)
    // Q = sum_k lambda_k * W_k
    var Q = 0.0
    for (i in 0 until K) {
        Q += lambdaArr[i] * W_arr[i]
    }
    val rhohat = Q / (1.0 + Q)

    // Convert result array to Matrix and return
    val W = Matrix(W_arr)
    return qsys_prio(W, rhohat)
}

/**
 * Queueing system M/G/1 with priorities algorithms.
 */
@Suppress("unused")
class QsysMg1PrioAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}
