/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.nc;

import jline.VerboseLevel;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * End-to-end test of the mixed open/closed Maximum Entropy Method in
 * SolverNC: the composition of the open (Section 3.2) and closed
 * (Section 3.3) Kouvatsos algorithms is exact in the BCMP product-form
 * limit, where it matches exact mixed MVA (references cross-checked
 * against MATLAB/python and SolverCTMC).
 */
public class SolverNCMemMixedTest {

    private static final double TOL = 1e-3;

    @Test
    public void testMemMixedProductForm() {
        // Open class (Poisson 0.3) traverses Q1->Q2; closed class (N=2)
        // cycles Q1<->Q2; exponential services. Exact mixed MVA:
        // L=[1.0822 1.5252; 0.2603 0.4748], X=[0.3, 0.6249]
        Network model = new Network("mix");
        Source source = new Source(model, "Source");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.PS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.PS);
        Sink sink = new Sink(model, "Sink");
        OpenClass oclass = new OpenClass(model, "O");
        ClosedClass cclass = new ClosedClass(model, "C", 2, q1);
        source.setArrival(oclass, new Exp(0.3));
        q1.setService(oclass, new Exp(1.0));
        q1.setService(cclass, new Exp(1.0));
        q2.setService(oclass, new Exp(2.0));
        q2.setService(cclass, new Exp(2.0));
        RoutingMatrix routing = model.initRoutingMatrix();
        routing.set(oclass, oclass, source, q1, 1.0);
        routing.set(oclass, oclass, q1, q2, 1.0);
        routing.set(oclass, oclass, q2, sink, 1.0);
        routing.set(cclass, cclass, q1, q2, 1.0);
        routing.set(cclass, cclass, q2, q1, 1.0);
        model.link(routing);

        SolverNC solver = new SolverNC(model, "method", "mem", "verbose", VerboseLevel.SILENT);
        solver.getAvg();
        NCResult res = (NCResult) solver.result;
        assertEquals("mem", res.method);
        assertEquals(1.0822, res.QN.get(1, 0), TOL);
        assertEquals(1.5252, res.QN.get(1, 1), TOL);
        assertEquals(0.2603, res.QN.get(2, 0), TOL);
        assertEquals(0.4748, res.QN.get(2, 1), TOL);
        assertEquals(0.3, res.XN.get(0, 0), TOL);
        assertEquals(0.6249, res.XN.get(0, 1), TOL);
    }
}
