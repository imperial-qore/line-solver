/**
 * @file 1-Limited Polling System Analysis
 *
 * Implements analysis algorithms for 1-limited polling systems where the
 * server serves at most one customer per visit to each queue. Provides
 * performance evaluation for multi-queue polling systems with limited service.
 *
 * @since LINE 3.0
 */
package jline.api.polling;

import jline.api.mam.Map_lambda;
import jline.api.mam.Map_mean;
import jline.api.mam.Map_moment;
import jline.api.mam.Map_var;
import jline.util.matrix.MatrixCell;

public final class Polling_qsys_1limited {
    private Polling_qsys_1limited() {}

    /**
     * Computes the exact mean waiting time solution for a polling system with open arrivals.
     * The system assumes that all queues use 1-limited service discipline.
     *
     * Reference: O. J. Boxma and B. Meister, "Waiting-time approximations for cyclic-service
     * systems with switch-over times", SIGMETRICS/PERFORMANCE '86, pp. 254-262.
     *
     * @param arvMAPs    arrival process MAPs
     * @param svcMAPs    service process MAPs
     * @param switchMAPs switching times MAPs
     * @return mean waiting times for each queue in the system
     */
    public static double[] polling_qsys_1limited(MatrixCell[] arvMAPs, MatrixCell[] svcMAPs, MatrixCell[] switchMAPs) {
        int n = arvMAPs.length;

        double[] lambda = new double[n];
        double[] b = new double[n];
        double[] b2 = new double[n];
        double[] rho1 = new double[n];
        double[] r1 = new double[n];
        double[] delta2 = new double[n];
        double rho = 0.0;
        double r = 0.0;
        double d = 0.0;

        for (int i = 0; i < n; i++) {
            lambda[i] = Map_lambda.map_lambda(arvMAPs[i]);
            b[i] = Map_mean.map_mean(svcMAPs[i]);
            b2[i] = Map_moment.map_moment(svcMAPs[i], 2);
            rho1[i] = lambda[i] * b[i];
            rho += rho1[i];
            r1[i] = Map_mean.map_mean(switchMAPs[i]);
            r += r1[i];
            delta2[i] = Map_var.map_var(switchMAPs[i]);
            d += delta2[i];
        }

        double[] W = new double[n];
        double sumOfSquares = 0.0;
        for (int i = 0; i < rho1.length; i++) {
            sumOfSquares += rho1[i] * rho1[i];
        }
        double sumProduct = 0.0;
        for (int i = 0; i < lambda.length; i++) {
            sumProduct += lambda[i] * b2[i];
        }

        for (int i = 0; i < n; i++) {
            W[i] = (1 - rho + rho1[i]) / (1 - rho - lambda[i] * r);

            W[i] = W[i] * (1 - rho) / ((1 - rho) * rho + sumOfSquares);
            W[i] = W[i] * (rho / (2 * (1 - rho)) * sumProduct + rho * d / 2 / r
                    + r / (2 * (1 - rho)) * rho1[i] * (1 + rho1[i]));
        }
        return W;
    }
}
