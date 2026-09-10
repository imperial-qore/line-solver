/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.examples.java.advanced;

import jline.api.pfqn.Pfqn_busyp;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.processes.Exp;
import jline.solvers.ldes.LDESOptions;
import jline.solvers.ldes.SolverLDES;
import jline.util.matrix.Matrix;

/**
 * Busy period of a subnetwork, exactly and by simulation.
 *
 * The busy period of order n for a set of stations is the time from the instant a
 * job entering the set finds n-1 jobs in it up to the next instant when fewer than
 * n remain. For large n it is a heavy-traffic period of the subnetwork, so the
 * whole family describes how long congestion of a given depth persists.
 *
 * {@link Pfqn_busyp} evaluates the mean value analysis of H. Daduna, "Busy Periods
 * for Subnetworks in Stochastic Networks: Mean Value Analysis", J. ACM 35(3), 1988,
 * from the normalizing constants of the subnetwork and of its complement.
 * {@link SolverLDES#getAvgBusyPeriod} measures the same quantity along the simulated
 * sample path, for a set of stations and optionally for a single class.
 */
public class BusyPeriodExample {

    public static void main(String[] args) {
        int N = 5;
        double[] rate = {1.5, 0.9, 2.0};
        double[][] Pc = {{0, 0.6, 0.4}, {0.7, 0, 0.3}, {0.5, 0.5, 0}};

        Network model = new Network("busyPeriodModel");
        Queue[] q = new Queue[3];
        for (int i = 0; i < 3; i++) {
            q[i] = new Queue(model, "Queue" + (i + 1), SchedStrategy.FCFS);
        }
        ClosedClass jobclass = new ClosedClass(model, "Class1", N, q[0], 0);
        for (int i = 0; i < 3; i++) {
            q[i].setService(jobclass, new Exp(rate[i]));
        }
        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                if (Pc[i][j] > 0) {
                    routingMatrix.set(jobclass, jobclass, q[i], q[j], Pc[i][j]);
                }
            }
        }
        model.link(routingMatrix);

        // visit ratios: the stochastic solution of x*P = x
        Matrix alpha = visits(Pc);
        Matrix mu = new Matrix(3, N);
        for (int i = 0; i < 3; i++) {
            for (int k = 0; k < N; k++) {
                mu.set(i, k, rate[i]);
            }
        }
        Matrix P = new Matrix(3, 3);
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                P.set(i, j, Pc[i][j]);
            }
        }

        LDESOptions options = new LDESOptions();
        options.samples = 2000000;
        options.seed(23000);
        options.busyPeriodOrders = 5;
        options.busyPeriodSubnet(0, 1);
        SolverLDES solver = new SolverLDES(model, options);

        int[][] subnets = {{0}, {1}, {2}, {0, 1}};
        System.out.printf("%-14s %10s %10s %10s %10s%n",
                "subnetwork", "n", "exact", "simulated", "rel.err");
        for (int[] subnet : subnets) {
            for (int n : new int[]{1, 3, 5}) {
                double exact = Pfqn_busyp.pfqn_busyp(alpha, mu, P, N, subnet, n, null);
                double sim = solver.getAvgBusyPeriod(subnet, -1, n);
                System.out.printf("%-14s %10d %10.4f %10.4f %10.2e%n",
                        java.util.Arrays.toString(subnet), n, exact, sim,
                        Math.abs(sim - exact) / exact);
            }
        }
    }

    /** Stochastic solution of x*P = x by power iteration. */
    private static Matrix visits(double[][] P) {
        int J = P.length;
        double[] x = new double[J];
        java.util.Arrays.fill(x, 1.0 / J);
        for (int it = 0; it < 20000; it++) {
            double[] y = new double[J];
            for (int i = 0; i < J; i++) {
                for (int j = 0; j < J; j++) {
                    y[j] += x[i] * P[i][j];
                }
            }
            double s = 0;
            for (int j = 0; j < J; j++) {
                s += y[j];
            }
            for (int j = 0; j < J; j++) {
                y[j] /= s;
            }
            x = y;
        }
        Matrix alpha = new Matrix(1, J);
        for (int j = 0; j < J; j++) {
            alpha.set(0, j, x[j]);
        }
        return alpha;
    }
}
