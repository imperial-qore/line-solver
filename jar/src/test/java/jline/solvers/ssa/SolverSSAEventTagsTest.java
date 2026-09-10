/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.ssa;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.util.matrix.Matrix;

/**
 * The derived START and PREEMPT rates estimated over the SAMPLE PATH.
 *
 * <p>SolverCTMC reports the exact value of the same two quantities; this is the
 * simulated estimator of them, and the point of the tests here is that the two
 * SSA engines answer alike. They must, or the counters would depend on
 * {@code options.method}: the serial engine accumulates the tag RATE of the
 * enabled transitions while the NRM engine COUNTS the events it fires and
 * divides by the simulated time. The two estimators agree in the limit and
 * differ by simulation error at any finite budget.
 *
 * <p>The oracle is exact rather than recorded. At a lossless station with no
 * in-service abandonment
 *
 * <pre>
 *     startRate(i,r) == throughput(i,r) + preemptRate(i,r)
 * </pre>
 *
 * which collapses to startRate == throughput at a non-preemptive station. A
 * Source CREATES jobs rather than admitting them to service, so it seizes
 * nothing and its row of both derived rates is zero.
 *
 * <p>The MATLAB twin is line-test.git/test/testsCTMC/test_ctmc_event_filtration.m, the Python
 * twin is python/tests/test_ctmc_event_filtration.py and the C++ twin is
 * cpp/tests/test_ctmc_event_filtration.cpp.
 */
public class SolverSSAEventTagsTest {

    private static final int SAMPLES = 100000;
    private static final int SEED = 23000;

    /** Source -&gt; Queue(FCFS) -&gt; Sink, one open class. */
    private static Network openSingleClass() {
        Network model = new Network("mm1");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass cls = new OpenClass(model, "Class1");
        source.setArrival(cls, new Exp(0.5));
        queue.setService(cls, new Exp(1.0));
        model.link(model.serialRouting(source, queue, sink));
        return model;
    }

    /** Source -&gt; Queue(FCFSPRPRIO) -&gt; Sink, an urgent and a normal open class. */
    private static Network openTwoClassPrio() {
        return openTwoClassPrio(0.4);
    }

    /** As above at a chosen per-class arrival rate. */
    private static Network openTwoClassPrio(double lambda) {
        Network model = new Network("mm1prio");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFSPRPRIO);
        Sink sink = new Sink(model, "Sink");
        OpenClass urgent = new OpenClass(model, "Urgent", 0);
        OpenClass normal = new OpenClass(model, "Normal", 1);
        source.setArrival(urgent, new Exp(lambda));
        source.setArrival(normal, new Exp(lambda));
        queue.setService(urgent, new Exp(1.0));
        queue.setService(normal, new Exp(1.0));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(urgent, urgent, model.serialRouting(source, queue, sink));
        P.set(normal, normal, model.serialRouting(source, queue, sink));
        model.link(P);
        return model;
    }

    private static SolverSSA run(Network model, String method) {
        SolverSSA solver = new SolverSSA(model);
        solver.options.method = method;
        solver.options.samples = SAMPLES;
        solver.options.seed = SEED;
        solver.getAvgTable();
        return solver;
    }

    /** Both engines: at an FCFS station every start is a later departure and nothing is displaced. */
    @Test
    public void bothEnginesAgreeAtANonPreemptiveStation() throws Exception {
        String[] methods = {"serial", "nrm"};
        for (String method : methods) {
            SolverSSA solver = run(openSingleClass(), method);
            Matrix TN = solver.getAvgTput();
            Matrix startN = solver.getStartRate();
            Matrix preemptN = solver.getPreemptRate();

            assertTrue(startN.get(1, 0) > 0, method + ": services must actually start");
            assertEquals(TN.get(1, 0), startN.get(1, 0), 5e-3,
                    method + ": one service starts per departure");
            assertEquals(0.0, preemptN.get(1, 0), 1e-12, method + ": FCFS displaces nothing");
            assertEquals(0.0, startN.get(0, 0), 1e-12, method + ": the Source seizes nothing");
        }
    }

    /**
     * The sampled path of a preempt-resume priority station satisfies the
     * identity, displaces the low class and never the high one, and starts the
     * high class once per arrival -- the last an oracle that does not read the
     * throughput.
     */
    @Test
    public void aPreemptiveStationReportsPreemptions() throws Exception {
        SolverSSA solver = run(openTwoClassPrio(), "serial");
        Matrix TN = solver.getAvgTput();
        Matrix startN = solver.getStartRate();
        Matrix preemptN = solver.getPreemptRate();

        for (int r = 0; r < 2; r++) {
            assertEquals(TN.get(1, r) + preemptN.get(1, r), startN.get(1, r), 1e-2,
                    "startRate == TN + preemptRate, class " + r);
        }
        assertEquals(0.4, startN.get(1, 0), 2e-2, "the high class starts once per arrival");
        assertEquals(0.0, preemptN.get(1, 0), 1e-12, "the high class is never displaced");
        assertTrue(preemptN.get(1, 1) > 0, "the low class must be displaced sometimes");
    }

    /**
     * A REGRESSION PIN for the state-re-encoding defect fixed on 2026-08-14.
     *
     * <p>{@code update_paddings_dense} widens the recorded history whenever a
     * station's local row grows mid-run. It took {@code deltalen} as a BOOLEAN
     * and always inserted exactly one row, however far the block had actually
     * widened. At a gain of one the two agree, so every per-class-count buffer
     * was unaffected; the preemptive families store (class,phase) PAIRS and gain
     * two, so their recorded history was re-encoded and states differing only in
     * which class held the server were conflated.
     *
     * <p>Both classes are lossless here, so each must carry its own arrival rate
     * whatever the discipline does -- an exact statement, not a recorded one.
     * LIGHTLY LOADED ON PURPOSE: at a heavier load the state-space cutoff
     * truncates and the exact answer falls below lambda, so pinning to lambda
     * there would pass on a coincidence. At lambda = 0.08 the exact CTMC gives
     * [0.08 0.08] to six digits, and this is also the regime where the defect
     * was starkest -- MATLAB's twin reported [0.018 0.142].
     */
    @Test
    public void aPreemptiveStationConservesEachClassArrivalRate() throws Exception {
        final double lambda = 0.08;
        SolverSSA solver = run(openTwoClassPrio(lambda), "serial");
        Matrix TN = solver.getAvgTput();
        assertEquals(lambda, TN.get(1, 0), 5e-3, "the high class departs as fast as it arrives");
        assertEquals(lambda, TN.get(1, 1), 5e-3, "so does the low class");
    }
}
