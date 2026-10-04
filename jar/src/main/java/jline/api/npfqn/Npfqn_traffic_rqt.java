/**
 * @file Effective arrival processes of a network under the Robust Queueing calculus
 *
 * Network characterization of Robust Queueing Theory, per C. Bandi,
 * D. Bertsimas, N. Youssef (2015), "Robust Queueing Theory", Operations
 * Research 63(3), 676-700, Theorems 4-7 and 10.
 *
 * @since LINE 3.0
 */
package jline.api.npfqn;

import jline.util.matrix.Matrix;

import org.apache.commons.math3.util.FastMath;

public final class Npfqn_traffic_rqt {
    private Npfqn_traffic_rqt() {}

    /**
     * Effective arrival process perceived at each node of a single-class open
     * queueing network under the Robust Queueing Theory calculus.
     *
     * The characterization composes three operators: passage through a queue
     * with adversarial servers leaves the uncertainty set unchanged (robust
     * Burke, Theorem 4), superposition merges sets by Theorem 5, and thinning by
     * a fraction f scales the rate by f and the variability by f^(-1/alpha)
     * (Theorem 6). The resulting equations are
     * lambda_j = lambda0_j + sum_i lambda_i f_ij and
     * Gamma_j = (1/lambda_j) [ 1{a0_j=ab_j} (lambda0_j Gamma0_j)^(p_j)
     * + sum_i 1{ab_i=ab_j} (lambda_i Gamma_i)^(p_i) f_ij ]^(1/p_j), with
     * p_j = ab_j/(ab_j-1) and ab_j the minimum tail coefficient among the
     * streams feeding j: the heaviest tail upstream dominates.
     *
     * <p>Both are solved exactly rather than iteratively. The rate equations are
     * the usual traffic equations, and in the variables z_j = (lambda_j
     * Gamma_j)^(p_j) the variability equations are linear as well, so each is one
     * linear system; ab is obtained by propagating the minimum to a fixed point.
     *
     * @param lambda0 external arrival rate at each node, 0 where there is none
     * @param Gamma0  variability parameter of each external arrival process,
     *                which for a renewal stream is the interarrival standard
     *                deviation
     * @param alpha0  tail coefficient in (1,2] of each external arrival process
     * @param F       routing probability matrix, F(i,j) = fraction of the jobs
     *                leaving node i that go to node j (row sums &lt;= 1)
     * @return array {lambda, Gamma, alpha} of length-J vectors: the effective
     *         arrival rate, variability parameter and tail coefficient at each
     *         node
     */
    public static double[][] npfqn_traffic_rqt(double[] lambda0, double[] Gamma0, double[] alpha0, Matrix F) {
        int J = lambda0.length;

        // traffic equations
        Matrix eyeJ = Matrix.eye(J);
        Matrix Xi = eyeJ.sub(F.transpose()).inv();   // (I - F')^{-1}
        double[] lambda = new double[J];
        for (int j = 0; j < J; j++) {
            double s = 0;
            for (int i = 0; i < J; i++) {
                s += Xi.get(j, i) * lambda0[i];
            }
            lambda[j] = FastMath.abs(s) < 1e-14 ? 0.0 : s;
        }

        // effective tail coefficient: the minimum propagated along the routing graph
        double[] alpha = new double[J];
        for (int j = 0; j < J; j++) {
            alpha[j] = lambda0[j] > 0 ? alpha0[j] : Double.POSITIVE_INFINITY;
        }
        for (int it = 0; it < J; it++) {
            boolean changed = false;
            for (int j = 0; j < J; j++) {
                for (int i = 0; i < J; i++) {
                    if (F.get(i, j) > 0 && lambda[i] > 0 && alpha[i] < alpha[j]) {
                        alpha[j] = alpha[i];
                        changed = true;
                    }
                }
            }
            if (!changed) {
                break;
            }
        }
        for (int j = 0; j < J; j++) {
            if (!Double.isFinite(alpha[j])) {
                alpha[j] = 2.0;   // an unreachable node keeps the light-tailed default
            }
        }

        // variability equations, linear in z_j = (lambda_j Gamma_j)^(p_j)
        double[] p = new double[J];
        for (int j = 0; j < J; j++) {
            p[j] = alpha[j] / (alpha[j] - 1);
        }
        Matrix z0 = new Matrix(J, 1);
        for (int j = 0; j < J; j++) {
            if (lambda0[j] > 0 && FastMath.abs(alpha0[j] - alpha[j]) < 1e-12) {
                z0.set(j, 0, FastMath.pow(lambda0[j] * Gamma0[j], p[j]));
            }
        }
        Matrix A = new Matrix(J, J);
        for (int i = 0; i < J; i++) {
            for (int j = 0; j < J; j++) {
                if (F.get(i, j) > 0 && FastMath.abs(alpha[i] - alpha[j]) < 1e-12) {
                    A.set(i, j, F.get(i, j));
                }
            }
        }
        Matrix z = eyeJ.sub(A.transpose()).inv().mult(z0);

        double[] Gamma = new double[J];
        for (int j = 0; j < J; j++) {
            double zj = FastMath.max(z.get(j, 0), 0.0);
            if (lambda[j] > 0) {
                Gamma[j] = FastMath.pow(zj, 1.0 / p[j]) / lambda[j];
            }
        }
        return new double[][]{lambda, Gamma, alpha};
    }
}
