/**
 * @file MX/M/1 queueing system analysis (batch arrivals)
 *
 * Implements analytical solutions for the MX/M/1 queue with batch Poisson arrivals
 * and exponential service times. Uses the Pollaczek-Khinchine formula extended
 * to batch arrivals.
 *
 * In an MX/M/1 queue:
 * - Batches arrive according to a Poisson process with rate λ_b
 * - Each batch has random size X with distribution specified by PMF
 * - Service times are exponentially distributed with rate μ
 * - Single server with FCFS discipline
 *
 * @since LINE 3.0
 */
package jline.api.qsys

import jline.io.Ret.qsys

/**
 * Analyzes an MX/M/1 queueing system using batch arrival moments.
 *
 * Uses the correct MX/M/1 formula that accounts for both queueing delay
 * and internal batch delay (jobs within the same batch wait for each other):
 *
 * - Effective arrival rate: λ = λ_b * E[X]
 * - Utilization: ρ = λ / μ
 * - Mean waiting time in queue:
 *   W_q = ρ/(μ(1-ρ)) + (E[X²] - E[X])/(2μE[X](1-ρ))
 *        \_________/   \__________________________/
 *         M/M/1 term      Internal batch delay
 * - Mean time in system: W = W_q + 1/μ
 *
 * @param lambdaBatch Batch arrival rate (λ_b)
 * @param mu Service rate
 * @param meanBatchSize Mean batch size E[X]
 * @param secondMomentBatchSize Second moment of batch size E[X²]
 * @return qsys containing average time in system (W) and utilization (rho)
 */
fun qsys_mxm1(lambdaBatch: Double, mu: Double, meanBatchSize: Double, secondMomentBatchSize: Double): qsys {
    // Effective arrival rate
    val lambda = lambdaBatch * meanBatchSize
    // Utilization
    val rho = lambda / mu

    if (rho >= 1.0) {
        return qsys(Double.POSITIVE_INFINITY, rho)
    }

    // Mean waiting time in queue for MX/M/1
    // Two components:
    // 1. Standard M/M/1 queueing delay: ρ / (μ * (1 - ρ))
    // 2. Internal batch delay (waiting for other jobs in same batch):
    //    (E[X²] - E[X]) / (2 * μ * E[X] * (1 - ρ))
    val mm1Term = rho / (mu * (1.0 - rho))
    val batchTerm = (secondMomentBatchSize - meanBatchSize) / (2.0 * mu * meanBatchSize * (1.0 - rho))
    val Wq = mm1Term + batchTerm

    // Mean time in system
    val W = Wq + 1.0 / mu

    return qsys(W, rho)
}

/**
 * Analyzes an MX/M/1 queueing system using batch sizes and PMF.
 *
 * @param lambdaBatch Batch arrival rate (λ_b)
 * @param mu Service rate
 * @param batchSizes Array of possible batch sizes
 * @param pmf Probability mass function for batch sizes (will be normalized)
 * @return qsys containing average time in system (W) and utilization (rho)
 */
fun qsys_mxm1(lambdaBatch: Double, mu: Double, batchSizes: IntArray, pmf: DoubleArray): qsys {
    require(batchSizes.size == pmf.size) { "Batch sizes and PMF must have same length" }

    // Normalize PMF
    val sum = pmf.sum()
    val normPmf = pmf.map { it / sum }

    // Compute E[X] and E[X²]
    var meanBatchSize = 0.0
    var secondMomentBatchSize = 0.0
    for (i in batchSizes.indices) {
        meanBatchSize += batchSizes[i] * normPmf[i]
        secondMomentBatchSize += batchSizes[i] * batchSizes[i] * normPmf[i]
    }

    return qsys_mxm1(lambdaBatch, mu, meanBatchSize, secondMomentBatchSize)
}

/**
 * Analyzes an MX/M/1 queueing system using mean and variance of batch size.
 *
 * @param lambdaBatch Batch arrival rate (λ_b)
 * @param mu Service rate
 * @param meanBatchSize Mean batch size E[X]
 * @param varianceBatchSize Variance of batch size Var(X)
 * @return qsys containing average time in system (W) and utilization (rho)
 */
fun qsys_mxm1_var(lambdaBatch: Double, mu: Double, meanBatchSize: Double, varianceBatchSize: Double): qsys {
    // E[X²] = Var(X) + E[X]²
    val secondMomentBatchSize = varianceBatchSize + meanBatchSize * meanBatchSize
    return qsys_mxm1(lambdaBatch, mu, meanBatchSize, secondMomentBatchSize)
}

/**
 * Queueing system MX/M/1 algorithms (batch arrivals)
 */
@Suppress("unused")
class QsysMxm1Algo {
    companion object {
        // Class documentation marker for Dokka
    }
}
