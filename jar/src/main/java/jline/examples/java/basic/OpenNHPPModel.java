/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java.basic;

import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.lang.processes.NHPP;
import jline.solvers.ldes.SolverLDES;

/**
 * Open queueing network with an NHPP (cyclic) arrival process.
 *
 * <p>NHPP is a non-homogeneous Poisson process with a piecewise-constant
 * intensity: segment i covers [breakpoints[i], breakpoints[i+1]) and carries
 * rate rates[i]. With cyclic=true the schedule repeats with period
 * breakpoints[n]-breakpoints[0], giving a cyclic Poisson process. Only the LDES
 * simulation engine honours the exact schedule; every analytical solver rejects
 * a model using it via the standard unsupported-feature check.
 */
public class OpenNHPPModel {

    /**
     * Open M/M/1 whose Source uses a cyclic NHPP arrival stream. Rates 2, 8, 4
     * are held for 3, 1, 2 time units respectively, repeating cyclically, so the
     * breakpoints are cumsum([0, 3, 1, 2]) = [0, 3, 4, 6].
     *
     * @return configured open network with an NHPP arrival process
     */
    public static Network oqn_nhpp() {
        Network model = new Network("model");

        // Block 1: nodes
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        // Block 2: classes
        OpenClass jobclass = new OpenClass(model, "OpenClass", 0);

        source.setArrival(jobclass, new NHPP(new double[]{0, 3, 4, 6},
                new double[]{2, 8, 4}, true));
        queue.setService(jobclass, new Exp(10));

        // Block 3: topology
        model.link(Network.serialRouting(source, queue, sink));

        return model;
    }

    /**
     * Solves the NHPP model with the LDES simulation engine and prints the
     * average performance table.
     *
     * @param args command line arguments (not used)
     * @throws Exception if the solver encounters an error
     */
    public static void main(String[] args) throws Exception {
        Network model = oqn_nhpp();
        SolverLDES solver = new SolverLDES(model, "samples", 100000, "seed", 1234);
        solver.getAvgTable().print();
    }
}
