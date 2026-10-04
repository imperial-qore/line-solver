/**
 * @file MX/M/1 queueing system analysis (batch arrivals)
 *
 * Implements analytical solutions for the MX/M/1 queue with batch Poisson arrivals
 * and exponential service times.
 *
 * @since LINE 3.0
 */
package jline.api.qsys;

import jline.io.Ret;

public final class Qsys_mxm1 {
    private Qsys_mxm1() {}

    /**
     * Analyzes an MX/M/1 queueing system using batch arrival moments.
     */
    public static Ret.qsys qsys_mxm1(double lambdaBatch, double mu, double meanBatchSize, double secondMomentBatchSize) {
        // Effective arrival rate
        double lambda = lambdaBatch * meanBatchSize;
        // Utilization
        double rho = lambda / mu;

        if (rho >= 1.0) {
            return new Ret.qsys(Double.POSITIVE_INFINITY, rho);
        }

        // Mean waiting time in queue for MX/M/1
        double mm1Term = rho / (mu * (1.0 - rho));
        double batchTerm = (secondMomentBatchSize - meanBatchSize) / (2.0 * mu * meanBatchSize * (1.0 - rho));
        double Wq = mm1Term + batchTerm;

        // Mean time in system
        double W = Wq + 1.0 / mu;

        return new Ret.qsys(W, rho);
    }

    /**
     * Analyzes an MX/M/1 queueing system using batch sizes and PMF.
     */
    public static Ret.qsys qsys_mxm1(double lambdaBatch, double mu, int[] batchSizes, double[] pmf) {
        if (batchSizes.length != pmf.length) {
            throw new IllegalArgumentException("Batch sizes and PMF must have same length");
        }

        // Normalize PMF
        double sum = 0.0;
        for (double p : pmf) sum += p;
        double[] normPmf = new double[pmf.length];
        for (int i = 0; i < pmf.length; i++) {
            normPmf[i] = pmf[i] / sum;
        }

        // Compute E[X] and E[X^2]
        double meanBatchSize = 0.0;
        double secondMomentBatchSize = 0.0;
        for (int i = 0; i < batchSizes.length; i++) {
            meanBatchSize += batchSizes[i] * normPmf[i];
            secondMomentBatchSize += batchSizes[i] * batchSizes[i] * normPmf[i];
        }

        return qsys_mxm1(lambdaBatch, mu, meanBatchSize, secondMomentBatchSize);
    }

    /**
     * Analyzes an MX/M/1 queueing system using mean and variance of batch size.
     */
    public static Ret.qsys qsys_mxm1_var(double lambdaBatch, double mu, double meanBatchSize, double varianceBatchSize) {
        // E[X^2] = Var(X) + E[X]^2
        double secondMomentBatchSize = varianceBatchSize + meanBatchSize * meanBatchSize;
        return qsys_mxm1(lambdaBatch, mu, meanBatchSize, secondMomentBatchSize);
    }
}
