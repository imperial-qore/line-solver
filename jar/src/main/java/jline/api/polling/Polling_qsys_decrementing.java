/**
 * @file Decrementing (Semiexhaustive) Polling System Analysis
 *
 * @since LINE 3.0
 */
package jline.api.polling;

import jline.api.mam.*;
import jline.util.matrix.MatrixCell;

/**
 * Exact mean waiting-time analysis for a symmetric decrementing (semiexhaustive)
 * polling system with Poisson arrivals and general service and switchover times.
 * <p>
 * The decrementing discipline serves a queue until the number of jobs present drops
 * to one less than the number found at the polling instant. For a symmetric system of
 * N queues the mean message waiting time is given in closed form by Pittel (1973) and
 * Takagi (1984); see also Takagi, "Queuing Analysis of Polling Models", ACM Computing
 * Surveys 20(1), 1988, eq. (28):
 * <pre>
 *   E[W] = delta2/(2 r)
 *        + ( N lambda b2 (1 - lambda r) + (r + lambda delta2)(N - rho) )
 *          / ( 2 ( 1 - rho - lambda r (N - rho) ) )
 * </pre>
 * where lambda is the per-queue arrival rate, b and b2 the first and second moments of
 * the service time, r and delta2 the mean and variance of the per-queue switchover time,
 * and rho = N lambda b the total offered load. No exact closed-form expression for the
 * individual E[W_i] is known for asymmetric decrementing systems, so this analysis is
 * restricted to the symmetric case.
 */
public final class Polling_qsys_decrementing {
    private Polling_qsys_decrementing() {}

    /**
     * Computes mean waiting times for a symmetric decrementing polling system.
     *
     * @param arvMAPs    per-class arrival MAPs (assumed Poisson; only the rate is used)
     * @param svcMAPs    per-class service MAPs (first and second moments are used)
     * @param switchMAPs per-class switchover MAPs (mean and variance are used)
     * @return array of mean waiting times, one per class (all equal in the symmetric case)
     * @throws RuntimeException if the per-class parameters are not symmetric
     */
    public static double[] polling_qsys_decrementing(MatrixCell[] arvMAPs,
                                                     MatrixCell[] svcMAPs,
                                                     MatrixCell[] switchMAPs) {
        int n = arvMAPs.length;

        double[] lambda = new double[n];
        double[] b = new double[n];
        double[] b2 = new double[n];
        double[] r1 = new double[n];
        double[] delta2 = new double[n];

        for (int i = 0; i < n; i++) {
            lambda[i] = Map_lambda.map_lambda(arvMAPs[i]);
            b[i] = Map_mean.map_mean(svcMAPs[i]);
            b2[i] = Map_moment.map_moment(svcMAPs[i], 2);
            r1[i] = Map_mean.map_mean(switchMAPs[i]);
            delta2[i] = Map_var.map_var(switchMAPs[i]);
        }

        // The exact result of Takagi (1984) requires a symmetric system: identical
        // arrival, service and switchover parameters across all queues.
        double tol = 1e-6;
        for (int i = 1; i < n; i++) {
            if (relDiff(lambda[i], lambda[0]) > tol
                    || relDiff(b[i], b[0]) > tol
                    || relDiff(b2[i], b2[0]) > tol
                    || relDiff(r1[i], r1[0]) > tol
                    || relDiff(delta2[i], delta2[0]) > tol) {
                throw new RuntimeException("MVA analysis for decrementing polling is only "
                        + "available for symmetric systems (identical arrival, service and "
                        + "switchover parameters across all queues).");
            }
        }

        double N = n;
        double lam = lambda[0];
        double b2s = b2[0];
        double r = r1[0];
        double d2 = delta2[0];
        double rho = N * lam * b[0];

        double denom = 2.0 * (1.0 - rho - lam * r * (N - rho));
        if (denom <= 0.0) {
            throw new RuntimeException("Decrementing polling system is unstable: "
                    + "rho + lambda*r*(N-rho) >= 1.");
        }

        // Mean residual switchover time; zero when there is no switchover (r = 0).
        double residualSwitchover = (r > 0.0) ? d2 / (2.0 * r) : 0.0;

        double W = residualSwitchover
                + (N * lam * b2s * (1.0 - lam * r) + (r + lam * d2) * (N - rho)) / denom;

        double[] result = new double[n];
        for (int i = 0; i < n; i++) {
            result[i] = W;
        }
        return result;
    }

    private static double relDiff(double a, double bb) {
        double scale = Math.max(Math.abs(a), Math.abs(bb));
        if (scale < 1e-30) {
            return 0.0;
        }
        return Math.abs(a - bb) / scale;
    }
}
