/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mva;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.HashSet;
import java.util.Set;

import org.junit.jupiter.api.Test;

import jline.VerboseLevel;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Exp;
import jline.solvers.ctmc.SolverCTMC;
import jline.util.matrix.Matrix;

/**
 * The closed-population AMVA family on a CLASS-SWITCHING model.
 *
 * <p>THE DEFECT THIS PINS. These sixteen algorithms recur on a population
 * vector, and the recursion presumes each entry of it is CONSERVED: the
 * arrival-instant estimate E[Q(N - 1_r)] is only meaningful if removing a
 * customer of r leaves a network of the same shape. Under class switching a job
 * CHANGES CLASS as it moves, so no per-class population is conserved -- the
 * conserved quantity is the CHAIN population. Solver_amva builds its whole
 * product-form branch out of SnGetProductFormChainParams and deaggregates
 * through SnDeaggregateChainResults at the end, and every arm honoured that
 * except the two Schmidt ones, which reached past it for {@code sn.rates} and
 * {@code sn.njobs}. On the model below they therefore solved a different
 * network, and their per-class answer was then read column by column as if it
 * were per-chain and deaggregated a second time.</p>
 *
 * <p>WHAT IS ASSERTED, and why none of it is a number read back out of the
 * implementation: the exact answer comes from SolverCTMC, which solves the
 * generator as written and never goes near these kernels; every approximation
 * must land near it; and they must not all land on the SAME number, because
 * sixteen differently-derived estimators agreeing bit for bit is the signature
 * of the defect rather than of accuracy.</p>
 */
public class SolverMVAClassSwitchingTest {

    /** The closed-population algorithms, less "sqni"; see the note on it below. */
    private static final String[] FAMILY = {
            "bs", "aql", "qsa", "tay", "scat", "lcp", "chow", "pamb", "pami",
            "pamt", "clust", "dmlin", "ab", "schmidt", "schmidt-ext"};

    /**
     * Delay -> Queue -> Delay with the class relabelled on each hop.
     *
     * <p>One chain of 2 jobs spread over two classes: C1 is served at the delay
     * and switches to C2 on the way to the queue, C2 switches back on the way
     * out. C2 therefore has a population of 0 of its own while the chain holds
     * 2.</p>
     *
     * @param sched the discipline of the queueing station
     * @return the model
     */
    private static Network classSwitching(SchedStrategy sched) {
        Network model = new Network("cs");
        Delay d = new Delay(model, "D");
        Queue q = new Queue(model, "Q", sched);
        ClosedClass c1 = new ClosedClass(model, "C1", 2, d);
        ClosedClass c2 = new ClosedClass(model, "C2", 0, d);
        d.setService(c1, new Exp(1.0));
        d.setService(c2, new Exp(1.0));
        q.setService(c1, new Exp(2.0));
        q.setService(c2, new Exp(3.0));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(c1, c2, d, q, 1.0);
        P.set(c2, c1, q, d, 1.0);
        model.link(P);
        return model;
    }

    /** Total mean queue length at station {@code i}, summed over the classes. */
    private static double qlenAt(Matrix QN, int i) {
        double s = 0.0;
        for (int r = 0; r < QN.getNumCols(); r++) {
            s += QN.get(i, r);
        }
        return s;
    }

    private static Matrix mvaQlen(Network model, String method) {
        SolverMVA s = new SolverMVA(model, "method", method, "verbose", VerboseLevel.SILENT);
        s.getAvg();
        return s.result.QN;
    }

    private static Matrix exactQlen(Network model) {
        SolverCTMC s = new SolverCTMC(model, "verbose", VerboseLevel.SILENT);
        s.getAvg();
        return s.result.QN;
    }

    @Test
    public void theClosedPopulationFamilySolvesAClassSwitchingModel() {
        // The reference, from a solver that never touches these kernels.
        Matrix exact = exactQlen(classSwitching(SchedStrategy.PS));
        double eq0 = qlenAt(exact, 0);
        double eq1 = qlenAt(exact, 1);
        assertEquals(1.4117647058823530, eq0, 1e-9);
        assertEquals(0.5882352941176471, eq1, 1e-9);

        Set<Long> distinct = new HashSet<Long>();
        for (String name : FAMILY) {
            Matrix q = mvaQlen(classSwitching(SchedStrategy.PS), name);
            double q0 = qlenAt(q, 0);
            double q1 = qlenAt(q, 1);
            // N = 2 jobs are somewhere, whatever approximation is used.
            assertEquals(2.0, q0 + q1, 1e-6, "method '" + name + "' lost the population");
            // None may be off by the whole queue: [2, 0] is 0.59 out at both
            // stations, which is the entire content of the queue row.
            assertTrue(q1 > 0.4, "method '" + name + "' left the queue empty");
            assertEquals(eq0, q0, 0.25, "method '" + name + "' is not near the exact answer");
            assertEquals(eq1, q1, 0.25, "method '" + name + "' is not near the exact answer");
            distinct.add(Long.valueOf(Math.round(q1 * 1e9)));
        }
        // Sixteen differently-derived estimators cannot all land on one number.
        assertTrue(distinct.size() > 3,
                "the family collapsed onto " + distinct.size() + " distinct answers");
    }

    @Test
    public void theSameClassSwitchingModelServedFcfs() {
        // FCFS is where the extended Schmidt correction is formed, so the two
        // Schmidt arms take a different route through the kernel here.
        Matrix exact = exactQlen(classSwitching(SchedStrategy.FCFS));
        double eq0 = qlenAt(exact, 0);
        double eq1 = qlenAt(exact, 1);
        String[] names = {"schmidt", "schmidt-ext", "ab", "bs"};
        for (String name : names) {
            Matrix q = mvaQlen(classSwitching(SchedStrategy.FCFS), name);
            assertEquals(eq0, qlenAt(q, 0), 0.25, "method '" + name + "'");
            assertEquals(eq1, qlenAt(q, 1), 0.25, "method '" + name + "'");
        }
    }

    @Test
    public void sqniReportsItsOwnClosedFormOnTheClassSwitchingModel() {
        // "sqni" is held out of the accuracy case above on purpose. Its closed
        // form reports Q = N - X Z off a square-root estimate of X, and on this
        // model that estimate is the saturation value X = 2, which drives the
        // queue term to exactly zero. That is the algorithm and not the defect --
        // the reference formula in pfqn_sqni.m gives the same 2 and 0 -- so what
        // is pinned here is that it still conserves the population and still runs.
        Matrix q = mvaQlen(classSwitching(SchedStrategy.PS), "sqni");
        assertEquals(2.0, qlenAt(q, 0) + qlenAt(q, 1), 1e-6);
        assertEquals(2.0, qlenAt(q, 0), 1e-6);
    }
}
