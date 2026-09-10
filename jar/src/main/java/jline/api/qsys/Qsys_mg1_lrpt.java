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

        // The reference has TWO branches and only the exponential one was ported
        // here, so a class with cs != 1 was answered with the exponential
        // formula: on lambda=[0.3,0.2], mu=[1,0.5], cs=[1,sqrt(0.5)] that read
        // W = [13.33, 5.33] against MATLAB's [7.111, 3.333].
        double[] W;
        boolean allExp = true;
        for (int i = 0; i < K; i++) {
            if (Math.abs(csArr[i] - 1.0) >= 1e-6) {
                allExp = false;
                break;
            }
        }
        if (allExp) {
            W = new double[K];
            double term2 = lambdaTotal * EX2 / (2.0 * (1.0 - rhoTotal) * (1.0 - rhoTotal));
            for (int k = 0; k < K; k++) {
                W[k] = (1.0 / muArr[k]) / (1.0 - rhoTotal) + term2;
            }
        } else {
            W = lrptGeneral(lambdaArr, muArr);
        }

        // Compute rhohat
        double Q = 0.0;
        for (int i = 0; i < K; i++) {
            Q += lambdaArr[i] * W[i];
        }
        double rhohat = Q / (1.0 + Q);

        return new Ret.qsys_prio(new Matrix(W), rhohat);
    }

    /**
     * The general (non-exponential) branch: MATLAB's class-based preemptive
     * priority surrogate, with the classes ordered by DECREASING mean service
     * time, which is the order LRPT serves them in.
     *
     * <pre>
     *   W_q(k) = (sum_{j&lt;=k} lambda_j/mu_j^2) / ((1 - rho_{&lt;k})(1 - rho_{&lt;=k}))
     *   W(k)   = W_q(k) + 1/mu_k
     * </pre>
     *
     * <p>The sort must be STABLE, as MATLAB's descending sort is, or two classes
     * of equal mean size swap places and the cumulative sums change.
     */
    private static double[] lrptGeneral(double[] lambdaArr, double[] muArr) {
        int K = lambdaArr.length;
        Integer[] idx = new Integer[K];
        for (int i = 0; i < K; i++) {
            idx[i] = Integer.valueOf(i);
        }
        final double[] meanService = new double[K];
        for (int i = 0; i < K; i++) {
            meanService[i] = 1.0 / muArr[i];
        }
        java.util.Arrays.sort(idx, new java.util.Comparator<Integer>() {
            @Override
            public int compare(Integer a, Integer b) {
                return Double.compare(meanService[b.intValue()], meanService[a.intValue()]);
            }
        });

        double[] W = new double[K];
        double rhoPrev = 0.0;
        double ERk = 0.0;
        for (int k = 0; k < K; k++) {
            int r = idx[k].intValue();
            double rhoCurr = rhoPrev + lambdaArr[r] / muArr[r];
            ERk += lambdaArr[r] / (muArr[r] * muArr[r]);
            W[r] = ERk / ((1.0 - rhoPrev) * (1.0 - rhoCurr)) + 1.0 / muArr[r];
            rhoPrev = rhoCurr;
        }
        return W;
    }
}
