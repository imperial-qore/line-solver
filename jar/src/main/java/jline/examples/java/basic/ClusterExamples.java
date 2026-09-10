/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java.basic;

import jline.gen.Cluster;
import jline.lang.Network;
import jline.lang.constant.RoutingStrategy;
import jline.lang.constant.SchedStrategy;
import jline.solvers.NetworkAvgTable;
import jline.solvers.ldes.SolverLDES;
import jline.solvers.mva.MVA;
import jline.solvers.mva.SolverMVA;
import jline.solvers.ssa.SolverSSA;

import java.util.Map;
import java.util.Scanner;

/**
 * Demonstrates cluster modelling APIs:
 * <ul>
 *   <li>The static factories {@code Network.cluster*}</li>
 *   <li>The {@link Cluster} builder with comparison and sweep helpers</li>
 * </ul>
 */
public class ClusterExamples {

    private static final Scanner scanner = new Scanner(System.in);

    private static void pauseForUser() {
        if (System.console() == null) {
            System.out.println("\n[Running in non-interactive mode, continuing...]");
            return;
        }
        System.out.println("\nPress Enter to continue to next example...");
        try { scanner.nextLine(); } catch (Exception e) { /* ignore */ }
    }

    /** Solves a four-PS-server open cluster with random dispatching. */
    public static void ex1_basic() {
        System.out.println("\n--- ex1_basic: 4 PS servers, random dispatching ---");
        Network model = ClusterModel.cl_basic();
        new MVA(model).getAvgTable().print();
    }

    /** Compares random vs round-robin vs JSQ dispatching with the builder API. */
    public static void ex2_compare_dispatching() {
        System.out.println("\n--- ex2_compare_dispatching: RAND vs RROBIN vs JSQ ---");
        Cluster cluster = new Cluster().setNumStations(4).setArrivalRate(1.0).setServiceRate(0.4)
                .setScheduling(SchedStrategy.PS);
        // RROBIN and JSQ routing are outside the MVA feature set, so the
        // dispatching comparison uses the LDES simulator, which supports all three
        Map<RoutingStrategy, NetworkAvgTable> results = cluster.compareDispatching(
                SolverLDES.class,
                RoutingStrategy.RAND,
                RoutingStrategy.RROBIN,
                RoutingStrategy.JSQ);
        for (Map.Entry<RoutingStrategy, NetworkAvgTable> e : results.entrySet()) {
            System.out.println("\nDispatching = " + e.getKey());
            e.getValue().print();
        }
    }

    /** Three FCFS servers with non-uniform service rates dispatched JSQ. */
    public static void ex3_heterogeneous() {
        System.out.println("\n--- ex3_heterogeneous: fast/medium/slow servers ---");
        Network model = ClusterModel.cl_heterogeneous();
        // JSQ routing is state-dependent and outside the MVA feature set;
        // the LDES simulator supports it
        new SolverLDES(model, "seed", 23000).getAvgTable().print();
    }

    /** Closed-network variant: 8 jobs cycling between Think and 3 PS servers. */
    public static void ex4_closed() {
        System.out.println("\n--- ex4_closed: 8-job closed cluster ---");
        Network model = ClusterModel.cl_closed();
        new MVA(model).getAvgTable().print();
    }

    /** Two-class open cluster (interactive + batch) on two PS servers. */
    public static void ex5_multiclass() {
        System.out.println("\n--- ex5_multiclass: two classes ---");
        Network model = ClusterModel.cl_multiclass();
        new MVA(model).getAvgTable().print();
    }

    /** Mixed variant: an open class and a closed class share the same two servers. */
    public static void ex7_mixed() {
        System.out.println("\n--- ex7_mixed: open + closed classes on shared servers ---");
        Network model = ClusterModel.cl_mixed();
        new MVA(model).getAvgTable().print();
    }

    /** Sweeps the arrival rate to see how response time grows toward saturation. */
    public static void ex6_sweep() {
        System.out.println("\n--- ex6_sweep: arrival-rate sweep ---");
        Cluster cluster = new Cluster().setNumStations(2).setArrivalRate(0.1).setServiceRate(1.0)
                .setScheduling(SchedStrategy.PS)
                .setDispatching(RoutingStrategy.RAND);
        Map<Double, NetworkAvgTable> sweep = cluster.sweepArrivalRate(
                new double[]{0.2, 0.5, 0.9, 1.5}, SolverMVA.class);
        for (Map.Entry<Double, NetworkAvgTable> e : sweep.entrySet()) {
            System.out.println("\nlambda = " + e.getKey());
            e.getValue().print();
        }
    }

    /**
     * Sweeps the number of parallel servers, the reference's {@code cl_stations}.
     *
     * The cluster is declared at one station and 1.6 arrivals against a service
     * rate of 1.0 -- unstable as declared, which is the point: the sweep is over
     * the server count that makes it stable, and utilization falls as 1/m while
     * the response time collapses once m passes the offered load.
     */
    public static void cl_stations() {
        System.out.println("\n--- cl_stations: sweep the number of servers ---");
        Cluster cluster = new Cluster().setNumStations(1).setArrivalRate(1.6).setServiceRate(1.0)
                .setScheduling(SchedStrategy.PS)
                .setDispatching(RoutingStrategy.RAND);
        Map<Integer, NetworkAvgTable> sweep = cluster.sweepNumStations(
                new int[]{2, 3, 4, 6}, SolverMVA.class);
        for (Map.Entry<Integer, NetworkAvgTable> e : sweep.entrySet()) {
            System.out.println("\nstations = " + e.getKey());
            e.getValue().print();
        }
    }

    /**
     * Compares scheduling disciplines on one cluster, the reference's
     * {@code cl_scheduling}.
     *
     * Service is hyperexponential (SCV 4), which is where the disciplines part
     * company: PS is insensitive to the service distribution beyond its mean,
     * FCFS is not, so the queue lengths differ even though the load does not.
     * The simulator is used rather than MVA because FCFS with non-exponential
     * service is outside the product-form assumptions.
     */
    public static void cl_scheduling() {
        System.out.println("\n--- cl_scheduling: FCFS against PS at SCV 4 ---");
        Cluster cluster = new Cluster().setNumStations(3).setArrivalRate(0.9).setServiceRate(0.5)
                .setDispatching(RoutingStrategy.RAND);
        cluster.setServiceSCV(4.0);
        Map<SchedStrategy, NetworkAvgTable> results = cluster.compareScheduling(
                model -> new SolverSSA(model, "seed", 23000, "samples", 20000).getAvgTable(),
                SchedStrategy.FCFS, SchedStrategy.PS);
        for (Map.Entry<SchedStrategy, NetworkAvgTable> e : results.entrySet()) {
            System.out.println("\nScheduling = " + e.getKey());
            e.getValue().print();
        }
    }

    public static void main(String[] args) {
        ex1_basic();
        pauseForUser();
        ex2_compare_dispatching();
        pauseForUser();
        ex3_heterogeneous();
        pauseForUser();
        ex4_closed();
        pauseForUser();
        ex5_multiclass();
        pauseForUser();
        ex7_mixed();
        pauseForUser();
        ex6_sweep();
        pauseForUser();
        cl_stations();
        pauseForUser();
        cl_scheduling();
    }
}
