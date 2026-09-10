/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.ctmc;

import jline.VerboseLevel;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Exp;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 * {@code options.config.chain_aggregation} on SolverCTMC, the first solver
 * consumer of the two class-level ModelAdapter transforms.
 *
 * <p>ModelAdapter.aggregateChains collapses every chain onto a single class and
 * SnDeaggregateChainResults maps chain-level metrics back through alpha. Both
 * existed in all four codebases with NO solver consumer at all: the transform
 * was exercised by examples and tests only, so nothing in the solver stack
 * depended on it and a defect in it could not surface as a wrong answer.
 *
 * <p>The oracle is an identity rather than a golden. On a PRODUCT-FORM model the
 * chain is the unit MVA and convolution already solve in, so the aggregation is
 * exact and the aggregated solve must reproduce the exact multiclass CTMC table
 * station by station and class by class. That is a statement about the transform
 * and the deaggregation together, and it is what fails if either mis-derives
 * alpha.
 */
public class ChainAggregationCTMCTest {

    private static SolverOptions silent() {
        SolverOptions options = new SolverOptions(SolverType.CTMC);
        options.verbose = VerboseLevel.SILENT;
        return options;
    }

    /** Delay -> Q1 -> Q2 -> Delay, the class switch on the Q1 -> Q2 link. */
    private static Network twoClassOneChain(int n) {
        Network model = new Network("agg");
        Delay d = new Delay(model, "Think");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.PS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.PS);
        ClosedClass c1 = new ClosedClass(model, "C1", n, d);
        ClosedClass c2 = new ClosedClass(model, "C2", 0, d);
        d.setService(c1, new Exp(1.0));
        d.setService(c2, new Exp(1.0));
        q1.setService(c1, new Exp(2.0));
        q1.setService(c2, new Exp(2.0));
        q2.setService(c1, new Exp(3.0));
        q2.setService(c2, new Exp(3.0));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(c1, c1, d, q1, 1.0);
        P.set(c1, c2, q1, q2, 1.0);
        P.set(c2, c1, q2, d, 1.0);
        model.link(P);
        return model;
    }

    private static Network oneClass() {
        Network model = new Network("plain");
        Delay d = new Delay(model, "Think");
        Queue q = new Queue(model, "Q1", SchedStrategy.PS);
        ClosedClass c = new ClosedClass(model, "C1", 3, d);
        d.setService(c, new Exp(1.0));
        q.setService(c, new Exp(2.0));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(c, c, model.serialRouting(d, q));
        model.link(P);
        return model;
    }

    private static SolverOptions aggregating() {
        SolverOptions options = silent();
        options.config.put("chain_aggregation", Boolean.TRUE);
        return options;
    }

    @Test
    public void theModelReallyHasTwoClassesInOneChain() {
        Network model = twoClassOneChain(3);
        assertEquals(2, model.getStruct(true).nclasses);
        assertEquals(1, model.getStruct(true).nchains);
    }

    @Test
    public void aggregatedSolveReproducesTheExactTable() {
        SolverCTMC exact = new SolverCTMC(twoClassOneChain(3), silent());
        SolverCTMC agg = new SolverCTMC(twoClassOneChain(3), aggregating());
        Matrix qe = exact.getAvgQLen();
        Matrix qa = agg.getAvgQLen();
        Matrix ue = exact.getAvgUtil();
        Matrix ua = agg.getAvgUtil();
        Matrix te = exact.getAvgTput();
        Matrix ta = agg.getAvgTput();
        assertEquals(qe.getNumRows(), qa.getNumRows());
        assertEquals(qe.getNumCols(), qa.getNumCols());
        for (int i = 0; i < qe.getNumRows(); i++) {
            for (int k = 0; k < qe.getNumCols(); k++) {
                assertEquals(qe.get(i, k), qa.get(i, k), 1e-9, "QLen at " + i + "," + k);
                assertEquals(ue.get(i, k), ua.get(i, k), 1e-9, "Util at " + i + "," + k);
                assertEquals(te.get(i, k), ta.get(i, k), 1e-9, "Tput at " + i + "," + k);
            }
        }
    }

    @Test
    public void flowIsConservedThroughTheDeaggregation() {
        SolverCTMC agg = new SolverCTMC(twoClassOneChain(3), aggregating());
        Matrix TN = agg.getAvgTput();
        double ref = 0;
        for (int k = 0; k < TN.getNumCols(); k++) ref += TN.get(0, k);
        assertTrue(ref > 0);
        for (int i = 1; i < TN.getNumRows(); i++) {
            double t = 0;
            for (int k = 0; k < TN.getNumCols(); k++) t += TN.get(i, k);
            assertEquals(ref, t, 1e-9, "throughput at station " + i);
        }
    }

    @Test
    public void populationIsConservedThroughTheDeaggregation() {
        SolverCTMC agg = new SolverCTMC(twoClassOneChain(5), aggregating());
        Matrix QN = agg.getAvgQLen();
        double total = 0;
        for (int i = 0; i < QN.getNumRows(); i++)
            for (int k = 0; k < QN.getNumCols(); k++) total += QN.get(i, k);
        assertEquals(5.0, total, 1e-9);
    }

    @Test
    public void theOptionIsDeclinedWhereTheTransformIsTheIdentity() {
        // nchains == nclasses, so the guard sends the model down the ordinary path.
        Matrix qe = new SolverCTMC(oneClass(), silent()).getAvgQLen();
        Matrix qa = new SolverCTMC(oneClass(), aggregating()).getAvgQLen();
        for (int i = 0; i < qe.getNumRows(); i++)
            assertEquals(qe.get(i, 0), qa.get(i, 0), 1e-12);
    }
}
