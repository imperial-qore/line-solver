/**
 * @file M/G/1 queueing system analysis with non-preemptive priorities
 *
 * @since LINE 3.0
 */
package jline.api.qsys;

import jline.io.Ret;
import jline.util.matrix.Matrix;

public final class Qsys_mg1_prio {
    private Qsys_mg1_prio() {}

    /**
     * Analyzes an M/G/1 queueing system with non-preemptive (Head-of-Line) priorities.
     */
    public static Ret.qsys_prio qsys_mg1_prio(Matrix lambda, Matrix mu, Matrix cs) {
        double[] lambdaArr = lambda.toArray1D();
        double[] muArr = mu.toArray1D();
        double[] csArr = cs.toArray1D();

        if (lambdaArr.length != muArr.length || lambdaArr.length != csArr.length) {
            throw new IllegalArgumentException(
                "lambda, mu, and cs must have the same length: got "
                + lambdaArr.length + ", " + muArr.length + ", " + csArr.length);
        }

        int K = lambdaArr.length;

        for (int i = 0; i < K; i++) {
            if (lambdaArr[i] <= 0.0) throw new IllegalArgumentException("lambda[" + i + "] must be positive, got " + lambdaArr[i]);
            if (muArr[i] <= 0.0) throw new IllegalArgumentException("mu[" + i + "] must be positive, got " + muArr[i]);
            if (csArr[i] <= 0.0) throw new IllegalArgumentException("cs[" + i + "] must be positive, got " + csArr[i]);
        }

        double[] rho_i = new double[K];
        double rho = 0.0;
        for (int i = 0; i < K; i++) {
            rho_i[i] = lambdaArr[i] / muArr[i];
            rho += rho_i[i];
        }

        if (rho >= 1.0) {
            throw new IllegalStateException("System is unstable: utilization rho = " + rho + " >= 1");
        }

        double B_0 = 0.0;
        for (int i = 0; i < K; i++) {
            B_0 += lambdaArr[i] * (1.0 + csArr[i] * csArr[i]) / (muArr[i] * muArr[i]);
        }
        B_0 /= 2.0;

        double[] W_q_arr = new double[K];
        for (int k = 0; k < K; k++) {
            double rho_prev = 0.0;
            for (int i = 0; i < k; i++) {
                rho_prev += rho_i[i];
            }
            double rho_curr = rho_prev + rho_i[k];
            W_q_arr[k] = B_0 / ((1.0 - rho_prev) * (1.0 - rho_curr));
        }

        double[] W_arr = new double[K];
        for (int i = 0; i < K; i++) {
            W_arr[i] = W_q_arr[i] + 1.0 / muArr[i];
        }

        double Q = 0.0;
        for (int i = 0; i < K; i++) {
            Q += lambdaArr[i] * W_arr[i];
        }
        double rhohat = Q / (1.0 + Q);

        return new Ret.qsys_prio(new Matrix(W_arr), rhohat);
    }
}
