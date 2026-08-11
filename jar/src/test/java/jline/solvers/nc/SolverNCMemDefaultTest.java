/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.nc;

import jline.VerboseLevel;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.ClosedClass;
import jline.lang.processes.APH;
import jline.lang.processes.Exp;
import jline.lang.processes.HyperExp;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Tests that SolverNC method "default" keeps the exact normalizing-constant
 * path and never engages the Maximum Entropy Method (Kouvatsos 1994) on its
 * own. MEM is available only on explicit request (method="mem").
 */
public class SolverNCMemDefaultTest {

    private static final double TOL = 1e-6;

    @Test
    public void testDefaultKeepsPathForBurstyLowLoad() {
        // Bursty arrivals (scv=64) at rho=0.33: the GE bound is loose at low
        // load, so the selector keeps the native product-form path
        Network model = new Network("bursty");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass oclass = new OpenClass(model, "Class1");
        source.setArrival(oclass, APH.fitMeanAndSCV(3.0, 64.0));
        queue.setService(oclass, new Exp(1.0));
        model.link(Network.serialRouting(source, queue, sink));

        SolverNC solver = new SolverNC(model, "verbose", VerboseLevel.SILENT);
        solver.getAvg();
        NCResult res = (NCResult) solver.result;
        if ("mem".equals(res.method)) {
            throw new AssertionError("bursty low-load open model must not route to mem by default");
        }
    }

    @Test
    public void testDefaultKeepsPathForHyperexponentialClosed() {
        // Closed cyclic with H2 service (scv=4): the native path is closer
        // in this regime, so the selector must not engage MEM
        Network model = new Network("h2c");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.FCFS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.FCFS);
        ClosedClass cclass = new ClosedClass(model, "C", 3, q1);
        q1.setService(cclass, HyperExp.fitMeanAndSCV(1.0, 4.0));
        q2.setService(cclass, new Exp(2.0));
        model.link(Network.serialRouting(q1, q2));

        SolverNC solver = new SolverNC(model, "verbose", VerboseLevel.SILENT);
        solver.getAvg();
        NCResult res = (NCResult) solver.result;
        if ("mem".equals(res.method)) {
            throw new AssertionError("hyperexponential closed model must not route to mem by default");
        }
    }

    @Test
    public void testDefaultKeepsProductFormPath() {
        // M/M/1 is product form: the default must keep the exact
        // normalizing-constant path (and still return the exact answer)
        Network model = new Network("mm1");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass oclass = new OpenClass(model, "Class1");
        source.setArrival(oclass, new Exp(0.6));
        queue.setService(oclass, new Exp(1.0));
        model.link(Network.serialRouting(source, queue, sink));

        SolverNC solver = new SolverNC(model, "verbose", VerboseLevel.SILENT);
        solver.getAvg();
        NCResult res = (NCResult) solver.result;
        assertEquals(1.5, res.QN.get(1, 0), TOL);
        if ("mem".equals(res.method)) {
            throw new AssertionError("product-form open model must not route to mem by default");
        }
    }
}
