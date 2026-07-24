/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api;

import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.constant.SchedStrategy;
import jline.lang.processes.Det;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.solvers.mam.SolverMAM;
import jline.solvers.mva.SolverMVA;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Regression tests for the qsys / MAM / MVA dispatch paths added in the
 * Hillier-Yu benchmark cases: M/G/c (PK, Crommelin), G/M/c (Smith, matrix-geometric).
 *
 * Each test builds a single-class open Source-Queue-Sink network at rho=0.99
 * and asserts the expected Lq within 5e-4 absolute. Reference values agree
 * with the analytical / matrix-geometric formulas and across all three LINE
 * codebases (MATLAB, JAR, Python native).
 */
public class HillierYuQsysTest {

    private static final double TOL = 5e-4;

    private static Network buildSerial(String name, int c) {
        Network m = new Network(name);
        Source src = new Source(m, "Source");
        Queue q = new Queue(m, "Queue1", SchedStrategy.FCFS);
        q.setNumberOfServers(c);
        Sink sink = new Sink(m, "Sink");
        OpenClass cls = new OpenClass(m, "Class1", 0);
        m.link(Network.serialRouting(src, q, sink));
        return m;
    }

    private static double lqMam(Network m, int c) throws Exception {
        SolverMAM solver = new SolverMAM(m);
        Matrix Q = solver.getAvgQLen();
        Matrix U = solver.getAvgUtil();
        return Q.get(1, 0) - c * U.get(1, 0);
    }

    private static double lqMvaExact(Network m, int c) throws Exception {
        SolverMVA solver = new SolverMVA(m, "exact");
        Matrix Q = solver.getAvgQLen();
        Matrix U = solver.getAvgUtil();
        return Q.get(1, 0) - c * U.get(1, 0);
    }

    private static OpenClass openClassOf(Network m) {
        for (jline.lang.JobClass jc : m.getClasses()) {
            if (jc instanceof OpenClass) return (OpenClass) jc;
        }
        throw new IllegalStateException("No OpenClass in model");
    }

    private static void rebuildArrival(Network m, int c, jline.lang.processes.Distribution arrival,
                                       jline.lang.processes.Distribution service) {
        // Helper: rebuild the model from scratch to set arrival/service.
        // (Network does not allow mutating existing classes after link.)
        // Not used here — each test builds its own model directly.
    }

    // ---------------- M/M/k via MVA exact ---------------------------------
    @Test
    public void testMM4MvaExact() throws Exception {
        Network m = new Network("MM4");
        Source src = new Source(m, "Source");
        Queue q = new Queue(m, "Queue1", SchedStrategy.FCFS);
        q.setNumberOfServers(4);
        Sink sink = new Sink(m, "Sink");
        OpenClass cls = new OpenClass(m, "Class1", 0);
        src.setArrival(cls, new Exp(3.96));
        q.setService(cls, new Exp(1.0));
        m.link(Network.serialRouting(src, q, sink));
        assertEquals(96.8126, lqMvaExact(m, 4), TOL);
    }

    // ---------------- M/E3/1 via MVA exact (PK) ---------------------------
    @Test
    public void testME3_1MvaExact() throws Exception {
        Network m = new Network("ME3_1");
        Source src = new Source(m, "Source");
        Queue q = new Queue(m, "Queue1", SchedStrategy.FCFS);
        Sink sink = new Sink(m, "Sink");
        OpenClass cls = new OpenClass(m, "Class1", 0);
        src.setArrival(cls, new Exp(0.99));
        q.setService(cls, Erlang.fitMeanAndOrder(1.0, 3));
        m.link(Network.serialRouting(src, q, sink));
        assertEquals(65.3400, lqMvaExact(m, 1), TOL);
    }

    // ---------------- M/E3/c via MAM (rate-scaling + surrogate) -----------
    @Test
    public void testME3_3Mam() throws Exception {
        Network m = new Network("ME3_3");
        Source src = new Source(m, "Source");
        Queue q = new Queue(m, "Queue1", SchedStrategy.FCFS);
        q.setNumberOfServers(3);
        Sink sink = new Sink(m, "Sink");
        OpenClass cls = new OpenClass(m, "Class1", 0);
        src.setArrival(cls, new Exp(2.97));
        q.setService(cls, Erlang.fitMeanAndOrder(1.0, 3));
        m.link(Network.serialRouting(src, q, sink));
        assertEquals(65.3400, lqMam(m, 3), TOL);
    }

    // ---------------- M/D/c via MAM (Crommelin embedded DTMC) -------------
    @Test
    public void testMD3MamCrommelin() throws Exception {
        Network m = new Network("MD_3");
        Source src = new Source(m, "Source");
        Queue q = new Queue(m, "Queue1", SchedStrategy.FCFS);
        q.setNumberOfServers(3);
        Sink sink = new Sink(m, "Sink");
        OpenClass cls = new OpenClass(m, "Class1", 0);
        src.setArrival(cls, new Exp(2.97));
        q.setService(cls, new Det(1.0));
        m.link(Network.serialRouting(src, q, sink));
        assertEquals(48.6193, lqMam(m, 3), TOL);
    }

    // ---------------- D/M/c via MAM (Smith embedded DTMC) -----------------
    @Test
    public void testDM4MamSmith() throws Exception {
        Network m = new Network("DM_4");
        Source src = new Source(m, "Source");
        Queue q = new Queue(m, "Queue1", SchedStrategy.FCFS);
        q.setNumberOfServers(4);
        Sink sink = new Sink(m, "Sink");
        OpenClass cls = new OpenClass(m, "Class1", 0);
        src.setArrival(cls, new Det(1.0 / 3.96));
        q.setService(cls, new Exp(1.0));
        m.link(Network.serialRouting(src, q, sink));
        assertEquals(47.8380, lqMam(m, 4), TOL);
    }

    // ---------------- E3/M/c via MAM (matrix-geometric) -------------------
    private double erlangMcMam(int c) throws Exception {
        Network m = new Network("E3M" + c);
        Source src = new Source(m, "Source");
        Queue q = new Queue(m, "Queue1", SchedStrategy.FCFS);
        q.setNumberOfServers(c);
        Sink sink = new Sink(m, "Sink");
        OpenClass cls = new OpenClass(m, "Class1", 0);
        src.setArrival(cls, Erlang.fitMeanAndOrder(1.0 / (c * 0.99), 3));
        q.setService(cls, new Exp(1.0));
        m.link(Network.serialRouting(src, q, sink));
        return lqMam(m, c);
    }

    @Test public void testE3M1Mam()  throws Exception { assertEquals(65.1206, erlangMcMam(1),  TOL); }
    @Test public void testE3M2Mam()  throws Exception { assertEquals(64.7159, erlangMcMam(2),  TOL); }
    @Test public void testE3M3Mam()  throws Exception { assertEquals(64.4035, erlangMcMam(3),  TOL); }
    @Test public void testE3M4Mam()  throws Exception { assertEquals(64.1398, erlangMcMam(4),  TOL); }
    @Test public void testE3M5Mam()  throws Exception { assertEquals(63.9075, erlangMcMam(5),  TOL); }
    @Test public void testE3M8Mam()  throws Exception { assertEquals(63.3256, erlangMcMam(8),  TOL); }
    @Test public void testE3M10Mam() throws Exception { assertEquals(62.9986, erlangMcMam(10), TOL); }

    // ---------------- E3/M/1 via MVA exact (PH/M/1 path) ------------------
    @Test
    public void testE3M1MvaExact() throws Exception {
        Network m = new Network("E3M1_MVA");
        Source src = new Source(m, "Source");
        Queue q = new Queue(m, "Queue1", SchedStrategy.FCFS);
        Sink sink = new Sink(m, "Sink");
        OpenClass cls = new OpenClass(m, "Class1", 0);
        src.setArrival(cls, Erlang.fitMeanAndOrder(1.0 / 0.99, 3));
        q.setService(cls, new Exp(1.0));
        m.link(Network.serialRouting(src, q, sink));
        assertEquals(65.1206, lqMvaExact(m, 1), TOL);
    }
}
