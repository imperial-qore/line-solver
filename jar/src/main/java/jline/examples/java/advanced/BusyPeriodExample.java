/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.examples.java.advanced;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.processes.Exp;
import jline.solvers.ldes.LDESOptions;
import jline.solvers.ldes.SolverLDES;
import jline.solvers.nc.SolverNC;

/**
 * Busy period of a subnetwork, exactly and by simulation.
 *
 * The busy period of order n for a set of stations is the time from the instant a
 * job entering the set finds n-1 jobs in it up to the next instant when fewer than
 * n remain. For large n it is a heavy-traffic period of the subnetwork, so the
 * whole family describes how long congestion of a given depth persists.
 *
 * {@link SolverNC#getAvgBusyPeriod} evaluates the mean value analysis of H. Daduna,
 * "Busy Periods for Subnetworks in Stochastic Networks: Mean Value Analysis",
 * J. ACM 35(3), 1988, from the normalizing constants of the subnetwork and of its
 * complement. {@link SolverLDES#getAvgBusyPeriod} measures the same quantity along
 * the simulated sample path, for a set of stations and optionally for a single class.
 */
public class BusyPeriodExample {

    /**
     * Mean busy period of order n for a subnetwork, exactly and by simulation
     * (busyp_subnetwork.m).
     *
     * <p>Station indexes are ZERO-based here and one-based in the reference, so
     * the subnetwork {@code [0, 1]} below is the reference's {@code [1 2]}.</p>
     */
    public static void busyp_subnetwork() {
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

        int[][] subnets = {{0}, {1}, {2}, {0, 1}};
        int[] orders = {1, 3, 5};

        SolverNC ncSolver = new SolverNC(model);
        System.out.printf("%nMean busy period of order n (exact, SolverNC)%n");
        System.out.printf("%-18s %10s %10s %10s%n", "subnetwork", "n=1", "n=3", "n=5");
        for (int[] subnet : subnets) {
            double[] b = ncSolver.getAvgBusyPeriod(subnet, orders);
            System.out.printf("%-18s %10.4f %10.4f %10.4f%n",
                    java.util.Arrays.toString(subnet), b[0], b[1], b[2]);
        }

        LDESOptions options = new LDESOptions();
        options.samples = 2000000;
        options.seed(23000);
        options.busyPeriodOrders = 5;
        options.busyPeriodSubnet(0, 1);
        SolverLDES ldesSolver = new SolverLDES(model, options);
        System.out.printf("%nMean busy period of order n (measured, SolverLDES)%n");
        System.out.printf("%-18s %10s %10s %10s%n", "subnetwork", "n=1", "n=3", "n=5");
        for (int[] subnet : subnets) {
            double[] b = new double[orders.length];
            for (int k = 0; k < orders.length; k++) {
                b[k] = ldesSolver.getAvgBusyPeriod(subnet, -1, orders[k]);
            }
            System.out.printf("%-18s %10.4f %10.4f %10.4f%n",
                    java.util.Arrays.toString(subnet), b[0], b[1], b[2]);
        }
    }

    public static void main(String[] args) {
        busyp_subnetwork();
    }

}
