/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.npfqn;

import java.util.ArrayList;
import java.util.List;
import java.util.function.DoubleUnaryOperator;

import jline.api.qsys.Qsys_gtmtst_fluid;
import jline.api.qsys.QsysTvFluidResult;

/**
 * Time-varying open network of many-server fluid queues with abandonment.
 *
 * <p>Each queue is the Gt/Mt/st+GI fluid queue of {@link Qsys_gtmtst_fluid}; the
 * departure flow of queue i is routed to queue j with the proportion P[i][j],
 * and whatever is left leaves the network.
 *
 * <p>THE NETWORK IS A FIXED POINT. The total arrival rate of queue j is
 * lambda_j(t) = lambda_j^0(t) + sum_i sigma_i(t)P_ij(t) with sigma_i = mu_i B_i
 * the service completion rate (eqs. 23-24), and sigma_i itself depends on
 * lambda_i. The iteration starts from the external rates alone and adds one more
 * traversal of the network per round, so the nth iterate is the fluid that has
 * made n transitions; the map is a monotone contraction, so the rates increase
 * to the fixed point rather than oscillating.
 *
 * <p>Only the SERVICE COMPLETION flow is routed. Abandoning fluid leaves the
 * network, which is what makes the traffic equations linear in sigma.
 *
 * <p>Port of MATLAB npfqn_gtmtst_fluid.m.
 *
 * <p>Reference: Y. Liu, W. Whitt (2014). Algorithms for time-varying networks of
 * many-server fluid queues. INFORMS Journal on Computing 26(1), 59-73.
 *
 * @since LINE 3.1.0
 */
public final class Npfqn_gtmtst_fluid {

    private Npfqn_gtmtst_fluid() {
    }

    /** Result of the network solve. */
    public static class Result {
        /** The time grid. */
        public final double[] times;
        /** The per-queue trajectories. */
        public final List<QsysTvFluidResult> queues;
        /** The converged total arrival rates, one row per queue. */
        public final double[][] arrivalRates;
        /** Iterations of the traffic-rate fixed point. */
        public final int iterations;
        /** Sup-norm change at the last iteration. */
        public final double residual;

        public Result(double[] times, List<QsysTvFluidResult> queues, double[][] arrivalRates,
                      int iterations, double residual) {
            this.times = times;
            this.queues = queues;
            this.arrivalRates = arrivalRates;
            this.iterations = iterations;
            this.residual = residual;
        }
    }

    /**
     * @param lambdaFuns    external arrival rate of each queue
     * @param sFuns         staffing of each queue
     * @param muFuns        service rate of each queue
     * @param patienceCcdfs patience ccdf of each queue
     * @param P             routing proportions, substochastic
     * @param T             horizon
     * @param dt            grid step; non-positive takes T/2000
     * @param B0            initial fluid in service, or null for empty
     * @param w0            initial boundary waiting times, or null for empty
     * @param tol           sup-norm tolerance on the arrival-rate iteration
     * @param maxIter       cap on the iterations
     * @return the network solution
     */
    public static Result npfqn_gtmtst_fluid(DoubleUnaryOperator[] lambdaFuns,
                                            DoubleUnaryOperator[] sFuns,
                                            DoubleUnaryOperator[] muFuns,
                                            DoubleUnaryOperator[] patienceCcdfs, double[][] P,
                                            double T, double dt, double[] B0, double[] w0,
                                            double tol, int maxIter) {
        int m = lambdaFuns.length;
        if (sFuns.length != m || muFuns.length != m || patienceCcdfs.length != m) {
            throw new RuntimeException("npfqn_gtmtst_fluid: every queue needs an arrival rate, a "
                    + "staffing, a service rate and a patience law");
        }
        if (T <= 0) {
            throw new RuntimeException("npfqn_gtmtst_fluid: the horizon T must be positive");
        }
        if (!(dt > 0)) {
            dt = T / 2000.0;
        }
        int n = (int) Math.round(T / dt) + 1;
        double[] t = new double[n];
        for (int i = 0; i < n; i++) {
            t[i] = T * i / (n - 1.0);
        }
        if (P.length != m) {
            throw new RuntimeException("npfqn_gtmtst_fluid: the routing matrix must be m x m");
        }
        for (int i = 0; i < m; i++) {
            if (P[i].length != m) {
                throw new RuntimeException("npfqn_gtmtst_fluid: the routing matrix must be m x m");
            }
            double row = 0.0;
            for (int j = 0; j < m; j++) {
                if (P[i][j] < -1e-12) {
                    throw new RuntimeException(
                            "npfqn_gtmtst_fluid: the routing matrix must be non-negative");
                }
                row += P[i][j];
            }
            if (row > 1 + 1e-9) {
                throw new RuntimeException(
                        "npfqn_gtmtst_fluid: the routing matrix must be substochastic");
            }
        }
        double[] b0 = B0 != null ? B0 : new double[m];
        double[] ww0 = w0 != null ? w0 : new double[m];

        double[][] ext = new double[m][n];
        for (int i = 0; i < m; i++) {
            for (int k = 0; k < n; k++) {
                ext[i][k] = lambdaFuns[i].applyAsDouble(t[k]);
            }
        }
        double[][] lam = new double[m][n];
        for (int i = 0; i < m; i++) {
            System.arraycopy(ext[i], 0, lam[i], 0, n);
        }

        List<QsysTvFluidResult> queues = new ArrayList<QsysTvFluidResult>();
        double residual = Double.POSITIVE_INFINITY;
        int iter = 0;
        for (iter = 1; iter <= maxIter; iter++) {
            queues = new ArrayList<QsysTvFluidResult>();
            double[][] sigma = new double[m][n];
            for (int i = 0; i < m; i++) {
                final double[] row = lam[i];
                final double[] grid = t;
                DoubleUnaryOperator fun = new DoubleUnaryOperator() {
                    @Override
                    public double applyAsDouble(double u) {
                        return Qsys_gtmtst_fluid.interp(grid, row, u);
                    }
                };
                QsysTvFluidResult r = Qsys_gtmtst_fluid.qsys_gtmtst_fluid(fun, sFuns[i], muFuns[i],
                        patienceCcdfs[i], T, dt, b0[i], ww0[i], null, null, null);
                queues.add(r);
                System.arraycopy(r.sigma, 0, sigma[i], 0, n);
            }
            double[][] newlam = new double[m][n];
            double diff = 0.0;
            for (int j = 0; j < m; j++) {
                for (int k = 0; k < n; k++) {
                    double acc = ext[j][k];
                    for (int i = 0; i < m; i++) {
                        acc += sigma[i][k] * P[i][j];
                    }
                    newlam[j][k] = acc;
                    diff = Math.max(diff, Math.abs(acc - lam[j][k]));
                }
            }
            lam = newlam;
            residual = diff;
            if (residual < tol) {
                break;
            }
        }
        return new Result(t, queues, lam, Math.min(iter, maxIter), residual);
    }
}
