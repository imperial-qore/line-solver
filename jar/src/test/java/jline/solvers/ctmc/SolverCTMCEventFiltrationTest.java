/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.ctmc;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.EventType;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.util.matrix.Matrix;

/**
 * The derived START and PREEMPT filtrations of SolverCTMC.
 *
 * <p>The identity that pins them is
 *
 * <pre>
 *     startRate(i,r) == throughput(i,r) + preemptRate(i,r)
 * </pre>
 *
 * at a lossless station with no in-service abandonment: every job starts
 * service once per entry into a server, and every preemption is followed by
 * exactly one later resume or restart. At a non-preemptive station it collapses
 * to startRate == throughput, which is an exact oracle rather than a
 * regression-recorded number.
 *
 * <p>The tags must also leave the chain alone: they are annotations on arcs the
 * generator already carries, so the generator and the mean measures must be
 * what they were before the annotation existed.
 */
public class SolverCTMCEventFiltrationTest {

    private static final double TOL = 1e-9;

    /** Source -&gt; Queue(sched) -&gt; Sink, one open class. */
    private static Network openSingleClass(SchedStrategy sched, int nservers) {
        Network model = new Network("mm1");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", sched);
        Sink sink = new Sink(model, "Sink");
        OpenClass cls = new OpenClass(model, "Class1");
        source.setArrival(cls, new Exp(0.5));
        queue.setService(cls, new Exp(1.0));
        if (nservers != 1) {
            queue.setNumberOfServers(nservers);
        }
        model.link(model.serialRouting(source, queue, sink));
        return model;
    }

    /** Source -&gt; Queue(sched) -&gt; Sink, an urgent and a normal open class. */
    private static Network openTwoClassPrio(SchedStrategy sched) {
        Network model = new Network("mm1prio");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", sched);
        Sink sink = new Sink(model, "Sink");
        OpenClass urgent = new OpenClass(model, "Urgent", 0);
        OpenClass normal = new OpenClass(model, "Normal", 1);
        source.setArrival(urgent, new Exp(0.4));
        source.setArrival(normal, new Exp(0.4));
        queue.setService(urgent, new Exp(1.0));
        queue.setService(normal, new Exp(1.0));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(urgent, urgent, model.serialRouting(source, queue, sink));
        P.set(normal, normal, model.serialRouting(source, queue, sink));
        model.link(P);
        return model;
    }

    /**
     * FCFS, PS and a two-server FCFS all promote without ever displacing anyone,
     * so the identity collapses to startRate == throughput at each of them.
     */
    @Test
    public void nonPreemptiveStartsOneServicePerDeparture() throws Exception {
        Network[] models = {
            openSingleClass(SchedStrategy.FCFS, 1),
            openSingleClass(SchedStrategy.PS, 1),
            openSingleClass(SchedStrategy.FCFS, 2),
        };
        for (int i = 0; i < models.length; i++) {
            SolverCTMC solver = new SolverCTMC(models[i]);
            solver.getOptions().cutoff = Matrix.singleton(4);
            Matrix TN = solver.getAvgTput();
            Matrix startN = solver.getStartRate();
            Matrix preemptN = solver.getPreemptRate();

            assertEquals(TN.get(1, 0), startN.get(1, 0), TOL,
                    "model " + i + ": one service starts per departure");
            assertEquals(0.0, preemptN.get(1, 0), TOL, "model " + i + ": nothing is displaced");
            // a Source CREATES jobs rather than admitting them to service, so it
            // seizes nothing and its row of both derived rates is zero
            assertEquals(0.0, startN.get(0, 0), TOL, "model " + i + ": the Source seizes nothing");
            assertEquals(0.0, preemptN.get(0, 0), TOL, "model " + i + ": the Source displaces nothing");
        }
    }

    /**
     * Preempt-resume and preempt-independent differ in the phase the displaced
     * job resumes in, not in how often it is displaced, so the preemption rate
     * is the SAME number.
     */
    @Test
    public void resumeAndIndependentReportTheSamePreemptionRate() throws Exception {
        SolverCTMC pr = new SolverCTMC(openTwoClassPrio(SchedStrategy.FCFSPRPRIO));
        pr.getOptions().cutoff = Matrix.singleton(3);
        SolverCTMC pi = new SolverCTMC(openTwoClassPrio(SchedStrategy.FCFSPIPRIO));
        pi.getOptions().cutoff = Matrix.singleton(3);
        Matrix prRate = pr.getPreemptRate();
        Matrix piRate = pi.getPreemptRate();
        for (int r = 0; r < 2; r++) {
            assertEquals(prRate.get(1, r), piRate.get(1, r), TOL,
                    "PR and PI displace equally often, class " + r);
        }
        assertTrue(prRate.get(1, 1) > 0, "the lower-priority class must be displaced sometimes");
    }

    /** Two priority classes at a preempt-resume station: the identity, and a positive preemption rate. */
    @Test
    public void fcfsPrPrioSatisfiesTheStartIdentity() throws Exception {
        Network model = new Network("mm1prio");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFSPRPRIO);
        Sink sink = new Sink(model, "Sink");
        OpenClass urgent = new OpenClass(model, "Urgent", 0);
        OpenClass normal = new OpenClass(model, "Normal", 1);
        source.setArrival(urgent, new Exp(0.4));
        source.setArrival(normal, new Exp(0.4));
        queue.setService(urgent, new Exp(1.0));
        queue.setService(normal, new Exp(1.0));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(urgent, urgent, model.serialRouting(source, queue, sink));
        P.set(normal, normal, model.serialRouting(source, queue, sink));
        model.link(P);

        SolverCTMC solver = new SolverCTMC(model);
        solver.getOptions().cutoff = Matrix.singleton(3);
        Matrix TN = solver.getAvgTput();
        Matrix startN = solver.getStartRate();
        Matrix preemptN = solver.getPreemptRate();

        for (int r = 0; r < 2; r++) {
            assertEquals(TN.get(1, r) + preemptN.get(1, r), startN.get(1, r), 1e-9,
                    "startRate == TN + preemptRate at station 1, class " + r);
        }
        assertEquals(0.0, preemptN.get(1, 0), TOL, "the urgent class is never displaced");
        assertTrue(preemptN.get(1, 1) > 0, "the lower-priority class must be displaced sometimes");
    }

    /** Closed cyclic model: the identity holds at both stations. */
    @Test
    public void closedCyclicSatisfiesTheStartIdentity() throws Exception {
        Network model = new Network("cyclic");
        Delay think = new Delay(model, "Think");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        ClosedClass jobs = new ClosedClass(model, "Jobs", 3, think);
        think.setService(jobs, new Exp(1.0));
        queue.setService(jobs, new Exp(2.0));
        model.link(model.serialRouting(think, queue));

        SolverCTMC solver = new SolverCTMC(model);
        Matrix TN = solver.getAvgTput();
        Matrix startN = solver.getStartRate();
        Matrix preemptN = solver.getPreemptRate();

        for (int i = 0; i < 2; i++) {
            assertEquals(TN.get(i, 0) + preemptN.get(i, 0), startN.get(i, 0), TOL,
                    "startRate == TN + preemptRate at station " + i);
        }
    }

    /** The filtration is reachable by event type and refuses a synchronization type. */
    @Test
    public void eventFiltrationServesTheDerivedEventsOnly() throws Exception {
        Network model = new Network("mm1");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass cls = new OpenClass(model, "Class1");
        source.setArrival(cls, new Exp(0.5));
        queue.setService(cls, new Exp(1.0));
        model.link(model.serialRouting(source, queue, sink));

        SolverCTMC solver = new SolverCTMC(model);
        solver.getOptions().cutoff = Matrix.singleton(3);
        Matrix[][] startFilt = solver.getEventFiltration(EventType.START);
        assertNotNull(startFilt, "the START filtration must be present");
        assertNotNull(startFilt[1][0], "station 1, class 0 must carry a filtration matrix");

        boolean refused = false;
        try {
            solver.getEventFiltration(EventType.DEP);
        } catch (IllegalArgumentException e) {
            refused = true;
        }
        assertTrue(refused, "a synchronization type must be refused, not silently served");
    }
}
