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
 * {@code options.config.fes_stations} on SolverCTMC, the first solver consumer
 * of {@code FESAggregator.aggregateFES}.
 *
 * <p>Flow-equivalent aggregation collapses a station subset into one
 * load-dependent station whose rate is the isolated subnetwork's throughput. The
 * transform has existed in all four codebases with NO solver consumer at all: it
 * was exercised by examples and tests only, so nothing in the solver stack
 * depended on it. This is that consumer, and the collapsed stations' own metrics
 * come back through the Chandy-Herzog-Woo conditional sum
 * E[Q_i] = sum_n P(N_fes = n) * Q_i(n).
 *
 * <p>The oracle is an identity, not a golden. On a product-form model the
 * decomposition is EXACT, so the reduced solve plus the conditioning must
 * reproduce the full chain's table station by station, INCLUDING the collapsed
 * stations. A wrong FES rate moves the surviving stations, a wrong conditional
 * sum moves only the collapsed ones, and a wrong visit ratio moves only their
 * throughput.
 */
public class FesAggregationCTMCTest {

    /** Think -> Q1 -> Q2 -> Q3 -> Think, a closed product-form cycle. */
    private static Network cycle(int n) {
        Network model = new Network("fes");
        Delay d = new Delay(model, "Think");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.PS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.PS);
        Queue q3 = new Queue(model, "Q3", SchedStrategy.PS);
        ClosedClass c = new ClosedClass(model, "C1", n, d);
        d.setService(c, new Exp(1.0));
        q1.setService(c, new Exp(2.0));
        q2.setService(c, new Exp(3.0));
        q3.setService(c, new Exp(4.0));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(c, c, model.serialRouting(d, q1, q2, q3));
        model.link(P);
        return model;
    }

    private static SolverOptions opts(int[] fes) {
        SolverOptions options = new SolverOptions(SolverType.CTMC);
        options.verbose = VerboseLevel.SILENT;
        if (fes != null) {
            options.config.put("fes_stations", fes);
        }
        return options;
    }

    @Test
    public void reducedSolveReproducesTheExactTable() {
        SolverCTMC exact = new SolverCTMC(cycle(4), opts(null));
        SolverCTMC fes = new SolverCTMC(cycle(4), opts(new int[]{2, 3}));
        Matrix qe = exact.getAvgQLen();
        Matrix ue = exact.getAvgUtil();
        Matrix te = exact.getAvgTput();
        Matrix qf = fes.getAvgQLen();
        Matrix uf = fes.getAvgUtil();
        Matrix tf = fes.getAvgTput();
        assertEquals(qe.getNumRows(), qf.getNumRows());
        for (int i = 0; i < qe.getNumRows(); i++) {
            assertEquals(qe.get(i, 0), qf.get(i, 0), 1e-9, "QLen at station " + i);
            assertEquals(ue.get(i, 0), uf.get(i, 0), 1e-9, "Util at station " + i);
            assertEquals(te.get(i, 0), tf.get(i, 0), 1e-9, "Tput at station " + i);
        }
    }

    @Test
    public void populationIsConservedByTheConditionalSplit() {
        SolverCTMC fes = new SolverCTMC(cycle(4), opts(new int[]{2, 3}));
        Matrix QN = fes.getAvgQLen();
        double total = 0;
        for (int i = 0; i < QN.getNumRows(); i++) total += QN.get(i, 0);
        assertEquals(4.0, total, 1e-9);
    }

    @Test
    public void aDifferentSubsetGivesTheSameAnswer() {
        // The choice of subset is the caller's and changes only which stations
        // are enumerated, so an exact decomposition must be invariant to it.
        Matrix qe = new SolverCTMC(cycle(4), opts(null)).getAvgQLen();
        Matrix a = new SolverCTMC(cycle(4), opts(new int[]{1, 2})).getAvgQLen();
        Matrix b = new SolverCTMC(cycle(4), opts(new int[]{2, 3})).getAvgQLen();
        for (int i = 0; i < qe.getNumRows(); i++) {
            assertEquals(qe.get(i, 0), a.get(i, 0), 1e-9, "subset {1,2} at station " + i);
            assertEquals(qe.get(i, 0), b.get(i, 0), 1e-9, "subset {2,3} at station " + i);
        }
    }

    @Test
    public void aSubsetThatSavesNothingIsRefusedByName() {
        assertThrows(RuntimeException.class,
                () -> new SolverCTMC(cycle(4), opts(new int[]{2})).getAvgQLen());
        assertThrows(RuntimeException.class,
                () -> new SolverCTMC(cycle(4), opts(new int[]{0, 1, 2, 3})).getAvgQLen());
        assertThrows(RuntimeException.class,
                () -> new SolverCTMC(cycle(4), opts(new int[]{2, 99})).getAvgQLen());
    }
}
