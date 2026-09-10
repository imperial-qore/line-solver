/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.mam;

import jline.VerboseLevel;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.solvers.SolverOptions;
import jline.solvers.mva.SolverMVA;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * The bgchain method on a PURELY CLOSED model, where it is the degenerate case
 * of the same construction rather than a separate algorithm.
 *
 * <p>With no open class the fixed point has nothing to iterate: cshare never
 * leaves its initial min(e,c), which is exactly the capacity the closed jobs
 * hold when no open work competes for it, so the background chain alone answers
 * and it answers with the EXACT closed CTMC at chain granularity. On a
 * product-form model that is exact MVA to machine precision, which is what the
 * first two tests assert -- not an approximation tolerance.</p>
 *
 * <p>The remaining tests pin the DEFAULT chooser: bgchain takes a closed model
 * whose chain fits in bgstates_max and whose service laws it represents exactly,
 * and stands aside otherwise. The background chain is built from the MEAN
 * service time alone, so a non-exponential law at a discipline that is not
 * insensitive to it belongs to mna, which carries the phase-type instead.</p>
 */
public class MamBgchainClosedTest {

    private static final double EXACT_TOL = 1e-9;

    /** Delay + mq PS queues in tandem, one closed class of n jobs. */
    private static Network tandem(int mq, int n, SchedStrategy sched, boolean erlang) {
        Network model = new Network("closedTandem");
        Delay delay = new Delay(model, "Think");
        Queue[] q = new Queue[mq];
        for (int i = 0; i < mq; i++) {
            q[i] = new Queue(model, "Q" + (i + 1), sched);
        }
        ClosedClass cclass = new ClosedClass(model, "C", n, delay);
        delay.setService(cclass, new Exp(1.0));
        for (int i = 0; i < mq; i++) {
            double rate = 1.0 + 0.3 * (i + 1);
            if (erlang) {
                q[i].setService(cclass, Erlang.fitMeanAndOrder(1.0 / rate, 3L));
            } else {
                q[i].setService(cclass, new Exp(rate));
            }
        }
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(cclass, cclass, delay, q[0], 1.0);
        for (int i = 0; i < mq - 1; i++) {
            P.set(cclass, cclass, q[i], q[i + 1], 1.0);
        }
        P.set(cclass, cclass, q[mq - 1], delay, 1.0);
        model.link(P);
        return model;
    }

    private static SolverOptions opts(String method) {
        SolverOptions options = SolverMAM.defaultOptions();
        options.method = method;
        options.verbose = VerboseLevel.SILENT;
        return options;
    }

    /** Explicitly requested bgchain must ANSWER a purely closed model, exactly. */
    @Test
    public void testClosedTandemMatchesExactMva() {
        int[] queues = {2, 4, 8};
        for (int mq : queues) {
            Network model = tandem(mq, 5, SchedStrategy.PS, false);
            Matrix QN = new SolverMAM(model, opts("bgchain")).getAvgQLen();
            Matrix QMVA = new SolverMVA(tandem(mq, 5, SchedStrategy.PS, false)).getAvgQLen();
            double held = 0.0;
            for (int i = 0; i < QN.getNumRows(); i++) {
                held += QN.get(i, 0);
                assertEquals(QMVA.get(i, 0), QN.get(i, 0), EXACT_TOL,
                        "queue length at station " + i + " with " + mq + " queues");
            }
            assertEquals(5.0, held, EXACT_TOL, "the closed population must be conserved");
        }
    }

    /** FCFS with exponential service is the same exact chain. */
    @Test
    public void testClosedFcfsTandemMatchesExactMva() {
        Network model = tandem(3, 4, SchedStrategy.FCFS, false);
        Matrix QN = new SolverMAM(model, opts("bgchain")).getAvgQLen();
        Matrix QMVA = new SolverMVA(tandem(3, 4, SchedStrategy.FCFS, false)).getAvgQLen();
        double held = 0.0;
        for (int i = 0; i < QN.getNumRows(); i++) {
            held += QN.get(i, 0);
            assertEquals(QMVA.get(i, 0), QN.get(i, 0), EXACT_TOL, "queue length at station " + i);
        }
        assertEquals(4.0, held, EXACT_TOL, "the closed population must be conserved");
    }

    /** The DEFAULT chooser takes bgchain on a closed model whose chain fits. */
    @Test
    public void testClosedDefaultIsBgchain() {
        Network model = tandem(4, 5, SchedStrategy.PS, false);
        SolverMAM solver = new SolverMAM(model, opts("default"));
        solver.getAvgQLen();
        // The JAR labels a resolved default as "default/<method>".
        assertTrue(solver.result.method.endsWith("bgchain"),
                "a closed model the background chain covers must default to it, got "
                        + solver.result.method);
    }

    /**
     * A non-exponential law under FCFS is NOT insensitive, and the background
     * chain carries only the mean, so the default must stand aside for mna.
     */
    @Test
    public void testClosedNonExponentialFcfsDoesNotDefaultToBgchain() {
        Network model = tandem(2, 4, SchedStrategy.FCFS, true);
        SolverMAM solver = new SolverMAM(model, opts("default"));
        solver.getAvgQLen();
        assertTrue(!solver.result.method.endsWith("bgchain"),
                "Erlang service under FCFS must not be answered by the mean-only chain, got "
                        + solver.result.method);
    }

    /** A population too large for bgstates_max must fall back, not refuse. */
    @Test
    public void testLargeClosedPopulationFallsBack() {
        Network model = tandem(6, 40, SchedStrategy.PS, false);
        SolverMAM solver = new SolverMAM(model, opts("default"));
        Matrix QN = solver.getAvgQLen();
        assertTrue(!solver.result.method.endsWith("bgchain"),
                "a background chain above bgstates_max must not be built as the default, got "
                        + solver.result.method);
        assertTrue(QN.getNumRows() > 0, "the fallback must still answer");
    }

    /** Delay -> Q -> Delay for two closed chains of 2 jobs each. */
    private static Network twoChainCycle(double rate1, double rate2, SchedStrategy sched) {
        Network model = new Network("twoChain");
        Delay delay = new Delay(model, "D");
        Queue queue = new Queue(model, "Q1", sched);
        ClosedClass c1 = new ClosedClass(model, "C1", 2, delay);
        ClosedClass c2 = new ClosedClass(model, "C2", 2, delay);
        delay.setService(c1, new Exp(1.0));
        delay.setService(c2, new Exp(1.0));
        queue.setService(c1, new Exp(rate1));
        queue.setService(c2, new Exp(rate2));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(c1, c1, delay, queue, 1.0);
        P.set(c1, c1, queue, delay, 1.0);
        P.set(c2, c2, delay, queue, 1.0);
        P.set(c2, c2, queue, delay, 1.0);
        model.link(P);
        return model;
    }

    /**
     * The chain splits a station's capacity in proportion to the job COUNTS,
     * i.e. service in random order. Under FCFS that is exact only at one shared
     * rate: with Exp(3) against Exp(0.8) it reads 25.2% off SolverCTMC, so the
     * default must stand aside.
     */
    @Test
    public void testClassDependentFcfsRatesDoNotDefaultToBgchain() {
        SolverMAM solver = new SolverMAM(twoChainCycle(3.0, 0.8, SchedStrategy.FCFS), opts("default"));
        solver.getAvgQLen();
        assertTrue(!solver.result.method.endsWith("bgchain"),
                "class-dependent FCFS rates must not be answered in random order, got "
                        + solver.result.method);
    }

    /** One shared rate under FCFS is the exact chain again. */
    @Test
    public void testClassIndependentFcfsRatesDefaultToBgchain() {
        SolverMAM solver = new SolverMAM(twoChainCycle(1.5, 1.5, SchedStrategy.FCFS), opts("default"));
        solver.getAvgQLen();
        assertTrue(solver.result.method.endsWith("bgchain"),
                "class-independent FCFS rates are exact here, got " + solver.result.method);
    }

    /** PS splits by DEMAND, not by count, so class-dependent rates cost nothing. */
    @Test
    public void testClassDependentRatesUnderPsDefaultToBgchain() {
        SolverMAM solver = new SolverMAM(twoChainCycle(3.0, 0.8, SchedStrategy.PS), opts("default"));
        solver.getAvgQLen();
        assertTrue(solver.result.method.endsWith("bgchain"),
                "PS is insensitive to the rate split, got " + solver.result.method);
    }

    /** A purely OPEN model has no closed population vector to build the chain from. */
    @Test
    public void testPurelyOpenModelIsRefused() {
        Network model = new Network("open");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Q1", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass open = new OpenClass(model, "Open");
        source.setArrival(open, new Exp(0.5));
        queue.setService(open, new Exp(1.0));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(open, open, source, queue, 1.0);
        P.set(open, open, queue, sink, 1.0);
        model.link(P);

        boolean refused = false;
        try {
            new SolverMAM(model, opts("bgchain")).getAvgQLen();
        } catch (RuntimeException e) {
            refused = true;
        }
        assertTrue(refused, "bgchain must refuse a model with no closed class");
    }
}
