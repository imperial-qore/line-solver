/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java.basic;

import java.util.ArrayList;
import java.util.List;

import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.lang.processes.MAP;
import jline.lang.processes.MAPt;
import jline.solvers.ldes.SolverLDES;
import jline.util.matrix.Matrix;

/**
 * Open queue with a MAP_t SERVICE process, simulated by LDES.
 *
 * <p>A MAP_t service is the only process LDES simulates that is both non-renewal and
 * non-stationary. Segment k covers [breakpoints[k], breakpoints[k+1]) and carries the pair
 * (D0[k], D1[k]); D1 fires a service COMPLETION and D0 only moves the modulating phase. The
 * engine walks that phase process forward from the instant the job ENTERS SERVICE -- not from
 * the instant it arrived, which may fall under a different segment -- and resumes from the phase
 * the previous completion left behind, so successive service times are correlated exactly as for
 * an ordinary MAP.
 *
 * <p>Because the walk is a wall-clock one, it is exact only where service runs continuously at
 * unit rate once started. LDES therefore refuses a MAP_t (or PH_t) service under processor
 * sharing, under a preemptive discipline, with load dependence and with heterogeneous server
 * pools, rather than returning a number from a sample path that models none of those. Use INF or
 * a non-preemptive FCFS/LCFS-family discipline, as here.
 */
public class OpenMAPtServiceModel {

    private static Matrix matrix(double[][] entries) {
        Matrix out = new Matrix(entries.length, entries[0].length);
        for (int i = 0; i < entries.length; i++) {
            for (int j = 0; j < entries[0].length; j++) {
                out.set(i, j, entries[i][j]);
            }
        }
        return out;
    }

    /**
     * Open M/MAP_t/1/5 whose Queue serves under a three-segment 2-phase MAP_t.
     *
     * <p>The schedule is an MMPP whose two environment states hold the same sojourn rates
     * throughout while their completion intensities change from segment to segment, so the
     * server alternates between a fast and a slow mode AND the pair of modes itself changes with
     * the clock. The segments are held for 1, 1 and 2 time units, repeating with period 4.
     *
     * @return configured open network with a MAP_t service process
     */
    public static Network oqn_mapt_service() {
        Network model = new Network("model");

        // Block 1: nodes
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        // Block 2: classes
        OpenClass jobclass = new OpenClass(model, "OpenClass", 0);

        double[] breakpoints = new double[]{0.0, 1.0, 2.0, 4.0};
        double[][] mu = new double[][]{{2.0, 8.0, 4.0}, {4.0, 2.0, 8.0}};
        double[] gamma = new double[]{1.0, 2.0};
        List<Matrix> d0 = new ArrayList<Matrix>();
        List<Matrix> d1 = new ArrayList<Matrix>();
        for (int k = 0; k < breakpoints.length - 1; k++) {
            d0.add(matrix(new double[][]{
                {-gamma[0] - mu[0][k], 1.0},
                {2.0, -gamma[1] - mu[1][k]}}));
            d1.add(matrix(new double[][]{{mu[0][k], 0.0}, {0.0, mu[1][k]}}));
        }

        source.setArrival(jobclass, new Exp(10));
        queue.setService(jobclass, new MAPt(breakpoints, d0, d1, true));
        queue.setCapacity(5);
        queue.setNumberOfServers(1);

        // Block 3: topology
        model.link(Network.serialRouting(source, queue, sink));

        return model;
    }

    /**
     * The same queue served by a MAP_t whose segments are IDENTICAL, which carries no time
     * dependence and must therefore reproduce {@link #homogeneousReference()}.
     *
     * @return the degenerate constant-schedule model
     */
    public static Network constantSchedule() {
        Network model = new Network("flat");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass jobclass = new OpenClass(model, "OpenClass", 0);
        List<Matrix> d0 = new ArrayList<Matrix>();
        List<Matrix> d1 = new ArrayList<Matrix>();
        for (int k = 0; k < 2; k++) {
            d0.add(matrix(new double[][]{{-6.0, 1.0}, {2.0, -12.0}}));
            d1.add(matrix(new double[][]{{5.0, 0.0}, {0.0, 10.0}}));
        }
        source.setArrival(jobclass, new Exp(1));
        queue.setService(jobclass, new MAPt(new double[]{0.0, 0.25, 0.5}, d0, d1, true));
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    /**
     * The ordinary MAP the constant schedule above must agree with.
     *
     * @return the time-homogeneous reference model
     */
    public static Network homogeneousReference() {
        Network model = new Network("homog");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass jobclass = new OpenClass(model, "OpenClass", 0);
        source.setArrival(jobclass, new Exp(1));
        queue.setService(jobclass, new MAP(
                matrix(new double[][]{{-6.0, 1.0}, {2.0, -12.0}}),
                matrix(new double[][]{{5.0, 0.0}, {0.0, 10.0}})));
        model.link(Network.serialRouting(source, queue, sink));
        return model;
    }

    /**
     * Solves the MAP_t service model with the LDES simulation engine and prints the average
     * performance table, then the constant-schedule degeneracy beside its MAP reference.
     *
     * @param args command line arguments (not used)
     * @throws Exception if the solver encounters an error
     */
    public static void main(String[] args) throws Exception {
        // The offered load exceeds what the server can clear, so the finite buffer drops: Tput
        // falls short of the arrival rate and QLen sits just under the capacity. Both are
        // signatures that the service process is being simulated -- a service that sampled as
        // zero would report QLen = Util = RespT = 0 with Tput equal to the arrival rate.
        new SolverLDES(oqn_mapt_service(), "samples", 1000000, "seed", 1234)
                .getAvgTable().print();

        System.out.println("\nconstant schedule vs ordinary MAP (must agree):");
        new SolverLDES(constantSchedule(), "samples", 400000, "seed", 1234)
                .getAvgTable().print();
        new SolverLDES(homogeneousReference(), "samples", 400000, "seed", 1234)
                .getAvgTable().print();
    }
}
