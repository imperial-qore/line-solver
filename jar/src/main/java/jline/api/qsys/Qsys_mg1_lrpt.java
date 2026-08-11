/**
 * @file M/G/1 queueing system analysis with LRPT scheduling
 *
 * Implements analysis of M/G/1 queues with Longest Remaining Processing Time (LRPT)
 * scheduling. Under LRPT, jobs with the longest remaining time share the processor.
 *
 * @since LINE 3.0
 */
package jline.api.qsys;

import jline.io.Ret;
import jline.util.matrix.Matrix;

public final class Qsys_mg1_lrpt {
    private Qsys_mg1_lrpt() {}

    /**
     * Analyzes an M/G/1 queueing system with LRPT (Longest Remaining Processing Time) scheduling.
     */
    public static Ret.qsys_prio qsys_mg1_lrpt(Matrix lambda, Matrix mu, Matrix cs) {
        double[] lambdaArr = lambda.toArray1D();
        double[] muArr = mu.toArray1D();
        double[] csArr = cs.toArray1D();

        if (lambdaArr.length != muArr.length || lambdaArr.length != csArr.length) {
            throw new IllegalArgumentException("lambda, mu, and cs must have the same length");
        }

        int K = lambdaArr.length;

        for (int i = 0; i < K; i++) {
            if (!(lambdaArr[i] > 0.0)) throw new IllegalArgumentException("lambda[" + i + "] must be positive");
            if (!(muArr[i] > 0.0)) throw new IllegalArgumentException("mu[" + i + "] must be positive");
            if (!(csArr[i] >= 0.0)) throw new IllegalArgumentException("cs[" + i + "] must be non-negative");
        }

        // Compute utilizations
        double[] rhoI = new double[K];
        double rhoTotal = 0.0;
        for (int i = 0; i < K; i++) {
            rhoI[i] = lambdaArr[i] / muArr[i];
            rhoTotal += rhoI[i];
        }

        if (!(rhoTotal < 1.0)) {
            throw new IllegalStateException("System is unstable: utilization rho = " + rhoTotal + " >= 1");
        }

        // Compute overall second moment of service time
        double lambdaTotal = 0.0;
        for (double l : lambdaArr) lambdaTotal += l;
        double[] p = new double[K];
        for (int i = 0; i < K; i++) p[i] = lambdaArr[i] / lambdaTotal;

        double EX2 = 0.0;
        for (int i = 0; i < K; i++) {
            double ES2i = (1.0 + csArr[i] * csArr[i]) / (muArr[i] * muArr[i]);
            EX2 += p[i] * ES2i;
        }

        // Compute response times using LRPT formula
        double[] W = new double[K];
        double term2 = lambdaTotal * EX2 / (2.0 * (1.0 - rhoTotal) * (1.0 - rhoTotal));

        for (int k = 0; k < K; k++) {
            double x = 1.0 / muArr[k];
            double term1 = x / (1.0 - rhoTotal);
            W[k] = term1 + term2;
        }

        // Compute rhohat
        double Q = 0.0;
        for (int i = 0; i < K; i++) {
            Q += lambdaArr[i] * W[i];
        }
        double rhohat = Q / (1.0 + Q);

        return new Ret.qsys_prio(new Matrix(W), rhohat);
    }
}
