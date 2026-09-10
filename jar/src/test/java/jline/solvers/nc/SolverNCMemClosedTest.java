/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.nc;

import jline.VerboseLevel;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Exp;
import jline.lang.processes.HyperExp;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

/**
 * End-to-end tests of the closed-network Maximum Entropy Method
 * (Kouvatsos 1994, Section 3.3) in SolverNC: the closed ME product form
 * reduces to the exact BCMP solution on Markovian single-class networks
 * (cyclic and repairmen models); multiclass and GE-type cases are pinned
 * to the cross-codebase MEM reference values.
 */
public class SolverNCMemClosedTest {

    private static final double TOL = 1e-3;

    @Test
    public void testMemClosedCyclic() {
        // Cyclic 2-station M/M/1, N=3: exact L=[2.2667, 0.7333], X=0.9333
        Network model = new Network("cyc");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.FCFS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.FCFS);
        ClosedClass cclass = new ClosedClass(model, "C", 3, q1);
        q1.setService(cclass, new Exp(1.0));
        q2.setService(cclass, new Exp(2.0));
        model.link(Network.serialRouting(q1, q2));

        SolverNC solver = new SolverNC(model, "method", "mem", "verbose", VerboseLevel.SILENT);
        solver.getAvg();
        NCResult res = (NCResult) solver.result;
        assertEquals(2.2667, res.QN.get(0, 0), TOL);
        assertEquals(0.7333, res.QN.get(1, 0), TOL);
        assertEquals(0.9333, res.TN.get(0, 0), TOL);
    }

    @Test
    public void testMemClosedRepairmen() {
        // Delay (mu=1) + queue (mu=1.25), N=4: exact L=[1.2132, 2.7868]
        Network model = new Network("rep");
        Delay d1 = new Delay(model, "D1");
        Queue q2 = new Queue(model, "Q2", SchedStrategy.FCFS);
        ClosedClass cclass = new ClosedClass(model, "C", 4, d1);
        d1.setService(cclass, new Exp(1.0));
        q2.setService(cclass, new Exp(1.25));
        model.link(Network.serialRouting(d1, q2));

        SolverNC solver = new SolverNC(model, "method", "mem", "verbose", VerboseLevel.SILENT);
        solver.getAvg();
        NCResult res = (NCResult) solver.result;
        assertEquals(1.2132, res.QN.get(0, 0), TOL);
        assertEquals(2.7868, res.QN.get(1, 0), TOL);
    }

    @Test
    public void testMemClosedH2Reference() {
        // Cyclic with H2 service (scv=4), N=3: MEM cross-codebase reference
        // values (MATLAB/python parity: L1=2.018982, X=0.820026)
        Network model = new Network("h2c");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.FCFS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.FCFS);
        ClosedClass cclass = new ClosedClass(model, "C", 3, q1);
        q1.setService(cclass, HyperExp.fitMeanAndSCV(1.0, 4.0));
        q2.setService(cclass, new Exp(2.0));
        model.link(Network.serialRouting(q1, q2));

        SolverNC solver = new SolverNC(model, "method", "mem", "verbose", VerboseLevel.SILENT);
        solver.getAvg();
        NCResult res = (NCResult) solver.result;
        assertEquals(2.018982, res.QN.get(0, 0), 1e-4);
        assertEquals(0.820026, res.TN.get(0, 0), 1e-4);
    }

    @Test
    public void testMemClosedTwoClassReference() {
        // Two-class cyclic with class-dependent rates, N=[2,2]: MEM
        // cross-codebase reference values (0.1% off exact CTMC on QLen)
        Network model = new Network("mc");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.FCFS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.FCFS);
        ClosedClass c1 = new ClosedClass(model, "C1", 2, q1);
        ClosedClass c2 = new ClosedClass(model, "C2", 2, q1);
        q1.setService(c1, new Exp(1.0));
        q1.setService(c2, new Exp(0.8));
        q2.setService(c1, new Exp(2.0));
        q2.setService(c2, new Exp(1.5));
        RoutingMatrix routing = model.initRoutingMatrix();
        routing.set(c1, c1, q1, q2, 1.0);
        routing.set(c1, c1, q2, q1, 1.0);
        routing.set(c2, c2, q1, q2, 1.0);
        routing.set(c2, c2, q2, q1, 1.0);
        model.link(routing);

        SolverNC solver = new SolverNC(model, "method", "mem", "verbose", VerboseLevel.SILENT);
        solver.getAvg();
        NCResult res = (NCResult) solver.result;
        assertEquals(1.5690, res.QN.get(0, 0), TOL);
        assertEquals(1.5434, res.QN.get(0, 1), TOL);
        assertEquals(0.4446, res.XN.get(0, 0), TOL);
        assertEquals(0.4140, res.XN.get(0, 1), TOL);
    }

    @Test
    public void testMemClosedRejectsMultiserver() {
        // Closed MEM builds on G/G/1 and G/G/inf only: finite multiserver
        // stations must be rejected
        Network model = new Network("ms");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.FCFS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.FCFS);
        q2.setNumberOfServers(2);
        ClosedClass cclass = new ClosedClass(model, "C", 3, q1);
        q1.setService(cclass, new Exp(1.0));
        q2.setService(cclass, new Exp(2.0));
        model.link(Network.serialRouting(q1, q2));

        SolverNC solver = new SolverNC(model, "method", "mem", "verbose", VerboseLevel.SILENT);
        assertThrows(RuntimeException.class, new org.junit.jupiter.api.function.Executable() {
            public void execute() {
                solver.getAvg();
            }
        });
    }
}
