/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java.solvers;

import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.processes.Exp;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.solvers.NetworkAvgTable;
import jline.solvers.ctmc.SolverCTMC;

import java.util.List;

/**
 * The solver-comparison examples of {@code matlab/examples/solvers}.
 *
 * <p>These exist to show a solver against a closed form rather than against
 * another solver: each builds a model whose exact answer is known on paper, so
 * the printed table can be read as a correctness statement and not only as a
 * number. The CTMC is the reference solver throughout, since it enumerates the
 * chain rather than approximating it.
 */
public class SolverExamples {

    /** The index of {@code station} in an average table, or -1 if it is absent. */
    private static int stationIndex(NetworkAvgTable table, String station) {
        List<String> names = table.getStationNames();
        for (int i = 0; i < names.size(); i++) {
            if (names.get(i).equals(station)) {
                return i;
            }
        }
        return -1;
    }

    private static void report(NetworkAvgTable table, String station, String label) {
        int i = stationIndex(table, station);
        if (i < 0) {
            System.out.println("\n" + label + ": not present in the table");
            return;
        }
        System.out.println("\n" + label + ":");
        System.out.printf("  Queue Length:   %.6f%n", table.getQLen().get(i));
        System.out.printf("  Utilization:    %.6f%n", table.getUtil().get(i));
        System.out.printf("  Response Time:  %.6f%n", table.getRespT().get(i));
        System.out.printf("  Throughput:     %.6f%n", table.getTput().get(i));
    }

    /**
     * Two M/M/1 queues in series, solved exactly.
     *
     * <p>Burke's theorem is what makes the closed form available: the departure
     * stream of a stationary M/M/1 queue is Poisson at the arrival rate, so the
     * second station sees the same 0.6 the first one does and the two stations
     * are independent M/M/1 queues. Each therefore has queue length
     * {@code rho/(1-rho)} -- 1.5 at station one, 1.0 at station two -- and the
     * CTMC below reproduces both without being told any of it.
     */
    public static void ctmc_tandem_mm1() {
        System.out.println("========================================================================");
        System.out.println("Tandem M/M/1 Queues (Series of Two Stations)");
        System.out.println("========================================================================");

        Network model = new Network("tandem_mm1");
        Source source = new Source(model, "source");
        Queue queue1 = new Queue(model, "queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "queue2", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "sink");
        OpenClass jobclass = new OpenClass(model, "jobs");

        source.setArrival(jobclass, Exp.fitMean(1 / 0.6));
        queue1.setService(jobclass, Exp.fitMean(1.0));
        queue2.setService(jobclass, Exp.fitMean(1 / 1.2));

        RoutingMatrix R = model.initRoutingMatrix();
        R.set(jobclass, jobclass, source, queue1, 1.0);
        R.set(jobclass, jobclass, queue1, queue2, 1.0);
        R.set(jobclass, jobclass, queue2, sink, 1.0);
        model.link(R);

        System.out.println("\nSolving tandem M/M/1 system with CTMC...");
        System.out.println("Parameters:");
        System.out.println("  Station 1: lambda=0.6, mu=1.0, rho=0.6");
        System.out.println("  Station 2: lambda=0.6, mu=1.2, rho=0.5\n");

        NetworkAvgTable avgTable = new SolverCTMC(model, "exact").getAvgTable();
        System.out.println("CTMC Results for Tandem M/M/1:");
        avgTable.print();

        System.out.println("\n========================================================================");
        System.out.println("CTMC Results Summary:");
        System.out.println("========================================================================");
        report(avgTable, "queue1", "Queue 1");
        report(avgTable, "queue2", "Queue 2");
    }

    public static void main(String[] args) {
        ctmc_tandem_mm1();
    }
}
