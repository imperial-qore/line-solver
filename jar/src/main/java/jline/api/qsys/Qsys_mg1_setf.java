/**
 * @file M/G/1 queueing system analysis with SETF scheduling
 *
 * @since LINE 3.0
 */
package jline.api.qsys;

import jline.io.Ret;
import jline.util.matrix.Matrix;

public final class Qsys_mg1_setf {
    private Qsys_mg1_setf() {}

    /**
     * Analyzes an M/G/1 queueing system with SETF (non-preemptive FB) scheduling.
     */
    public static Ret.qsys_prio qsys_mg1_setf(Matrix lambda, Matrix mu, Matrix cs) {
        double[] lambdaArr = lambda.toArray1D();
        double[] muArr = mu.toArray1D();
        double[] csArr = cs.toArray1D();

        if (!(lambdaArr.length == muArr.length && lambdaArr.length == csArr.length)) {
            throw new IllegalArgumentException("lambda, mu, and cs must have the same length");
        }

        int K = lambdaArr.length;

        for (int i = 0; i < K; i++) {
            if (!(lambdaArr[i] > 0.0)) {
                throw new IllegalArgumentException("lambda[" + i + "] must be positive");
            }
            if (!(muArr[i] > 0.0)) {
                throw new IllegalArgumentException("mu[" + i + "] must be positive");
            }
            if (!(csArr[i] >= 0.0)) {
                throw new IllegalArgumentException("cs[" + i + "] must be non-negative");
            }
        }

        double[] rhoI = new double[K];
        double rhoTotal = 0.0;
        for (int i = 0; i < K; i++) {
            rhoI[i] = lambdaArr[i] / muArr[i];
            rhoTotal += rhoI[i];
        }

        if (!(rhoTotal < 1.0)) {
            throw new IllegalStateException("System is unstable: utilization rho = " + rhoTotal + " >= 1");
        }

        double lambdaTotal = 0.0;
        for (int i = 0; i < K; i++) {
            lambdaTotal += lambdaArr[i];
        }
        double meanResidual = 0.0;
        for (int i = 0; i < K; i++) {
            double pI = lambdaArr[i] / lambdaTotal;
            double meanS = 1.0 / muArr[i];
            double meanS2 = (1.0 + csArr[i] * csArr[i]) / (muArr[i] * muArr[i]);
            meanResidual += pI * meanS2 / (2.0 * meanS);
        }

        double[] W = new double[K];

        for (int k = 0; k < K; k++) {
            double x = 1.0 / muArr[k];

            double rhoX = 0.0;
            for (int i = 0; i < K; i++) {
                double integralFbar;
                if (Math.abs(csArr[i] - 1.0) < 1e-10) {
                    integralFbar = (1.0 - Math.exp(-muArr[i] * x)) / muArr[i];
                } else {
                    integralFbar = Math.min(x, 1.0 / muArr[i]);
                }
                rhoX += lambdaArr[i] * integralFbar;
            }

            double numerator = 0.0;
            for (int i = 0; i < K; i++) {
                double integralTFbar;
                if (Math.abs(csArr[i] - 1.0) < 1e-10) {
                    double muI = muArr[i];
                    integralTFbar = (1.0 - Math.exp(-muI * x) * (1.0 + muI * x)) / (muI * muI);
                } else {
                    integralTFbar = Math.min(x * x / 2.0, 1.0 / (muArr[i] * muArr[i]));
                }
                numerator += lambdaArr[i] * integralTFbar;
            }

            if (rhoX >= 1.0) {
                W[k] = Double.POSITIVE_INFINITY;
            } else {
                double fbWaitingTerm = numerator / ((1.0 - rhoX) * (1.0 - rhoX));
                double fbServiceTerm = x / (1.0 - rhoX);
                double npPenalty = meanResidual / (1.0 - rhoX);

                W[k] = fbWaitingTerm + fbServiceTerm + npPenalty;
            }
        }

        double Q = 0.0;
        for (int i = 0; i < K; i++) {
            Q += lambdaArr[i] * W[i];
        }
        double rhohat = Q / (1.0 + Q);

        return new Ret.qsys_prio(new Matrix(W), rhohat);
    }
}
