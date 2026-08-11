/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.nc;

import jline.VerboseLevel;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * End-to-end tests of SolverNC method "mem" (Maximum Entropy Method,
 * Kouvatsos 1994) on Markovian networks with exact references: M/M/1,
 * M/M/3 with a downstream Delay (Erlang-C and IS exact) and an M/M/1
 * with Bernoulli feedback (Jackson exact).
 */
public class SolverNCMemTest {

    private static final double TOL = 1e-6;

    @Test
    public void testMemMM1() {
        Network model = new Network("mm1");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass oclass = new OpenClass(model, "Class1");
        source.setArrival(oclass, new Exp(0.6));
        queue.setService(oclass, new Exp(1.0));
        model.link(Network.serialRouting(source, queue, sink));

        SolverNC solver = new SolverNC(model, "method", "mem", "verbose", VerboseLevel.SILENT);
        solver.getAvgTable();
        NCResult res = (NCResult) solver.result;
        assertEquals(1.5, res.QN.get(1, 0), TOL);
        assertEquals(0.6, res.UN.get(1, 0), TOL);
        assertEquals(2.5, res.RN.get(1, 0), TOL);
        assertEquals(0.6, res.TN.get(0, 0), TOL); // source row reports arrivals
    }

    @Test
    public void testMemMMcAndDelay() {
        // M/M/3 (lambda=2, mu=1; Erlang-C exact L=26/9) followed by an IS
        // delay (exact L=1)
        Network model = new Network("mmc_is");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        queue.setNumberOfServers(3);
        Delay delay = new Delay(model, "Delay");
        Sink sink = new Sink(model, "Sink");
        OpenClass oclass = new OpenClass(model, "Class1");
        source.setArrival(oclass, new Exp(2.0));
        queue.setService(oclass, new Exp(1.0));
        delay.setService(oclass, new Exp(2.0));
        model.link(Network.serialRouting(source, queue, delay, sink));

        SolverNC solver = new SolverNC(model, "method", "mem", "verbose", VerboseLevel.SILENT);
        solver.getAvgTable();
        NCResult res = (NCResult) solver.result;
        assertEquals(26.0 / 9.0, res.QN.get(1, 0), TOL);
        assertEquals(2.0 / 3.0, res.UN.get(1, 0), TOL);
        assertEquals(1.0, res.QN.get(2, 0), TOL);
    }

    @Test
    public void testMemFeedback() {
        // M/M/1 with 50% Bernoulli feedback: Jackson-exact L=1 at rho=0.5;
        // throughput is visit-inclusive (2.0)
        Network model = new Network("fb");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass oclass = new OpenClass(model, "Class1");
        source.setArrival(oclass, new Exp(1.0));
        queue.setService(oclass, new Exp(4.0));
        RoutingMatrix routing = model.initRoutingMatrix();
        routing.set(oclass, source, queue, 1.0);
        routing.set(oclass, queue, queue, 0.5);
        routing.set(oclass, queue, sink, 0.5);
        model.link(routing);

        SolverNC solver = new SolverNC(model, "method", "mem", "verbose", VerboseLevel.SILENT);
        solver.getAvgTable();
        NCResult res = (NCResult) solver.result;
        assertEquals(1.0, res.QN.get(1, 0), TOL);
        assertEquals(0.5, res.UN.get(1, 0), TOL);
        assertEquals(2.0, res.TN.get(1, 0), TOL);
    }
}
