/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api;

import jline.api.npfqn.Npfqn_dps_morrison;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.solvers.nc.SolverNC;
import jline.solvers.nc.analyzers.Solver_nc_dps_analyzer;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Morrison's heavy-usage expansion for a closed think+DPS network, and its NC solver route.
 *
 * <p>The pinned values are the MATLAB reference (matlab/src/api/npfqn/npfqn_dps_morrison.m) and
 * agree with the C++ and python ports to the digits written here; the same models appear in the
 * test of each codebase, so a divergence in any of them fails somewhere.</p>
 *
 * <p>Reference: J.A. Morrison, Queueing Systems 9 (1991) 191-214.</p>
 */
public class DpsMorrisonTest {

    private static final double TOL = 1e-9;

    // M2: two classes, K = 100 each, rho = 0.9 -- the regime the expansion is derived for.
    private static final double[] M2_N = {100.0, 100.0};
    private static final double[] M2_Z = {1.0, 0.5};
    private static final double[] M2_S = {0.006, 0.0015};
    private static final double[] M2_W = {1.0, 4.0};

    // M3: three classes, unequal populations and weights, rho = 1.06 (appendix A's saturated side).
    private static final double[] M3_N = {20.0, 12.0, 28.0};
    private static final double[] M3_Z = {1.0, 0.5, 2.0};
    private static final double[] M3_S = {0.02, 0.01, 0.03};
    private static final double[] M3_W = {1.0, 4.0, 2.0};

    private static Matrix row(double[] v) {
        Matrix m = new Matrix(1, v.length);
        for (int i = 0; i < v.length; i++) {
            m.set(0, i, v[i]);
        }
        return m;
    }

    private static Npfqn_dps_morrison.Result run(double[] n, double[] z, double[] s, double[] w) {
        return Npfqn_dps_morrison.npfqn_dps_morrison(row(n), row(z), row(s), row(w));
    }

    @Test
    public void testTwoClassesMatchesReference() {
        Npfqn_dps_morrison.Result res = run(M2_N, M2_Z, M2_S, M2_W);
        assertEquals(0.9, res.rho, 1e-12);
        assertEquals(0.1, res.a, 1e-12);
        assertEquals(4.053964005705, res.Q.get(0, 0), TOL);
        assertEquals(0.930378969124, res.Q.get(0, 1), TOL);
        assertEquals(0.042452089190, res.R.get(0, 0), TOL);
        assertEquals(0.004666835854, res.R.get(0, 1), TOL);
        assertEquals(95.946035994295, res.X.get(0, 0), 1e-8);
        assertEquals(198.139242061752, res.X.get(0, 1), 1e-8);
        assertEquals(4.373155763276, res.Qlead.get(0, 0), TOL);
        assertEquals(0.546644470409, res.Qlead.get(0, 1), TOL);
        assertEquals(0.084297520661, res.sigma.get(0, 0), TOL);
        assertEquals(-0.337190082645, res.sigma.get(0, 1), TOL);
    }

    @Test
    public void testThreeClassesMatchesReference() {
        Npfqn_dps_morrison.Result res = run(M3_N, M3_Z, M3_S, M3_W);
        assertEquals(1.06, res.rho, 1e-12);
        assertEquals(3.213021210443, res.Q.get(0, 0), TOL);
        assertEquals(0.918773941158, res.Q.get(0, 1), TOL);
        assertEquals(2.645208279899, res.Q.get(0, 2), TOL);
        assertEquals(0.208463525546, res.R.get(0, 0), TOL);
        assertEquals(0.039776387080, res.R.get(0, 1), TOL);
        assertEquals(0.202390704352, res.R.get(0, 2), TOL);
    }

    /**
     * Throughput closes the think station EXACTLY, which is why the analyzers report R = Q/T.
     * Morrison's RESULT 2 (4.17) is the EXPANDED ratio (4.11)/(4.15), so it differs from Q/X by
     * higher-order terms -- 22% at K=6, 0.5% at K=100, 0.09% at K=400. The two agree only
     * asymptotically, and that convergence is what is asserted here.
     */
    @Test
    public void testLittlesLawAndTheTwoResponseTimeForms() {
        Npfqn_dps_morrison.Result res = run(M3_N, M3_Z, M3_S, M3_W);
        for (int j = 0; j < M3_N.length; j++) {
            assertEquals((M3_N[j] - res.Q.get(0, j)) / M3_Z[j], res.X.get(0, j), 1e-10);
        }
        double g6 = gap(run(new double[]{6.0, 6.0}, M2_Z, new double[]{0.075, 0.0375}, M2_W));
        double g100 = gap(run(M2_N, M2_Z, M2_S, M2_W));
        double g400 = gap(run(new double[]{400.0, 400.0}, M2_Z, new double[]{0.0015, 0.000375}, M2_W));
        assertTrue(g400 < g100, g400 + " !< " + g100);
        assertTrue(g100 < g6, g100 + " !< " + g6);
        assertTrue(g100 < 0.01, "gap at K=100 is " + g100);
    }

    private static double gap(Npfqn_dps_morrison.Result k) {
        double g = 0;
        for (int j = 0; j < k.R.getNumCols(); j++) {
            g = Math.max(g, Math.abs(k.R.get(0, j) - k.Q.get(0, j) / k.X.get(0, j)) / k.R.get(0, j));
        }
        return g;
    }

    /**
     * Equal weights make the network product-form, so Morrison's (4.18)-(4.21) collapse: D = B,
     * Q = C, M = H = I, L = J, hence R = S = V = U = 0, A = -K, delta = 0 and sigma = 0.
     */
    @Test
    public void testEqualWeightsCollapseTheCorrection() {
        Npfqn_dps_morrison.Result res = run(M2_N, M2_Z, M2_S, new double[]{1.0, 1.0});
        assertEquals(res.cB, res.cD, 1e-9 * res.cB);
        assertEquals(res.cC, res.cQ, 1e-9 * res.cC);
        assertEquals(res.cH, res.cI, 1e-9 * res.cH);
        assertEquals(res.cH, res.cM, 1e-9 * res.cH);
        assertEquals(res.cJ, res.cL, 1e-9 * res.cL);
        assertEquals(0.0, res.cR, 1e-9);
        assertEquals(0.0, res.cS, 1e-9);
        assertEquals(0.0, res.cV, 1e-9);
        assertEquals(0.0, res.cU, 1e-9);
        assertEquals(-res.cK, res.cA, 1e-9 * res.cK);
        assertEquals(0.0, res.delta, 1e-9);
        assertEquals(0.0, res.sigma.get(0, 0), 1e-9);
        assertEquals(0.0, res.sigma.get(0, 1), 1e-9);
        assertEquals(3.737598812660, res.Q.get(0, 0), TOL);
        assertEquals(1.937595955666, res.Q.get(0, 1), TOL);
    }

    /**
     * An asymptotic expansion has to get relatively better as the populations grow at fixed usage.
     * Both models run at rho = 0.9; the exact values are from the CTMC of eq. (2.1).
     */
    @Test
    public void testAccuracyImprovesWithPopulation() {
        double[] exactSmall = {1.178868777433, 0.746116560046};
        double[] exactBig = {4.082236085180, 0.862829749692};
        Npfqn_dps_morrison.Result small =
                run(new double[]{6.0, 6.0}, M2_Z, new double[]{0.075, 0.0375}, M2_W);
        Npfqn_dps_morrison.Result big = run(M2_N, M2_Z, M2_S, M2_W);
        for (int j = 0; j < 2; j++) {
            double eSmall = Math.abs(small.Q.get(0, j) - exactSmall[j]) / exactSmall[j];
            double eBig = Math.abs(big.Q.get(0, j) - exactBig[j]) / exactBig[j];
            assertTrue(eBig < eSmall, "class " + j + ": " + eBig + " !< " + eSmall);
            assertTrue(eBig < 0.1, "class " + j + ": " + eBig);
        }
    }

    private static Network thinkDps(int n1, int n2, double s1, double s2, double w2) {
        Network model = new Network("morrison");
        Delay delay = new Delay(model, "Think");
        Queue queue = new Queue(model, "DPS", SchedStrategy.DPS);
        ClosedClass c1 = new ClosedClass(model, "Class1", n1, delay);
        ClosedClass c2 = new ClosedClass(model, "Class2", n2, delay);
        delay.setService(c1, new Exp(1.0));
        delay.setService(c2, new Exp(2.0));
        queue.setService(c1, new Exp(1.0 / s1));
        queue.setService(c2, new Exp(1.0 / s2));
        queue.setSchedStrategyPar(c1, 1.0);
        queue.setSchedStrategyPar(c2, w2);
        model.link(model.serialRouting(delay, queue));
        return model;
    }

    /** SolverNC answers the think+DPS shape by Morrison, and reports no normalizing constant. */
    @Test
    public void testSolverNcRoutesToMorrison() throws Exception {
        Network model = thinkDps(100, 100, 0.006, 0.0015, 4.0);
        SolverNC solver = new SolverNC(model);
        Matrix q = solver.getAvgQLen();
        // station order is Delay, Queue: the DPS row carries Morrison's answer
        assertEquals(4.053964005705, q.get(1, 0), 1e-6);
        assertEquals(0.930378969124, q.get(1, 1), 1e-6);
        // population is conserved exactly by the analyzer
        assertEquals(100.0, q.get(0, 0) + q.get(1, 0), 1e-9);
        assertEquals(100.0, q.get(0, 1) + q.get(1, 1), 1e-9);
        // not product-form: no normalizing constant
        assertTrue(Double.isNaN(((jline.solvers.nc.NCResult) solver.result).logNormConstAggr()));
    }

    /** "morrison" is advertised, and every other NC method is refused on that shape. */
    @Test
    public void testMethodGate() {
        Network model = thinkDps(100, 100, 0.006, 0.0015, 4.0);
        boolean advertised = false;
        for (String m : new SolverNC(model).listValidMethods()) {
            if ("morrison".equals(m)) {
                advertised = true;
            }
        }
        assertTrue(advertised, "morrison must be advertised by listValidMethods");

        SolverNC exact = new SolverNC(thinkDps(100, 100, 0.006, 0.0015, 4.0), "exact");
        assertThrows(Exception.class, exact::getAvgQLen,
                "an NC method that cannot see the DPS weights must be refused, not answered");
    }

    /** A DPS station outside Morrison's shape is refused rather than approximated. */
    @Test
    public void testDpsOutsideShapeIsRefused() {
        // non-exponential service at the DPS station
        Network model = new Network("erlangDps");
        Delay delay = new Delay(model, "Think");
        Queue queue = new Queue(model, "DPS", SchedStrategy.DPS);
        ClosedClass c1 = new ClosedClass(model, "Class1", 10, delay);
        ClosedClass c2 = new ClosedClass(model, "Class2", 10, delay);
        delay.setService(c1, new Exp(1.0));
        delay.setService(c2, new Exp(2.0));
        queue.setService(c1, Erlang.fitMeanAndOrder(0.006, 3));
        queue.setService(c2, new Exp(1.0 / 0.0015));
        queue.setSchedStrategyPar(c1, 1.0);
        queue.setSchedStrategyPar(c2, 4.0);
        model.link(model.serialRouting(delay, queue));
        SolverNC solver = new SolverNC(model);
        assertThrows(Exception.class, solver::getAvgQLen,
                "a DPS model outside the think+DPS shape must be refused");
    }

    private static Network plainPs() {
        Network model = new Network("ps");
        Delay delay = new Delay(model, "Think");
        Queue queue = new Queue(model, "PS", SchedStrategy.PS);
        ClosedClass c1 = new ClosedClass(model, "Class1", 6, delay);
        ClosedClass c2 = new ClosedClass(model, "Class2", 6, delay);
        delay.setService(c1, new Exp(1.0));
        delay.setService(c2, new Exp(2.0));
        queue.setService(c1, new Exp(3.0));
        queue.setService(c2, new Exp(5.0));
        model.link(model.serialRouting(delay, queue));
        return model;
    }

    /**
     * "morrison" is in listValidMethods, so without an explicit arm it would pass the method gate,
     * match no route, and fall through to the ordinary normalizing-constant path -- answering a
     * product-form model UNDER THE CALLER'S LABEL. It must be refused by name instead.
     */
    @Test
    public void testMorrisonNamedOnANonDpsModelIsRefused() {
        SolverNC solver = new SolverNC(plainPs(), "morrison");
        assertThrows(Exception.class, solver::getAvgQLen);
    }

    /** The analyzer re-checks the shape rather than trusting its caller. */
    @Test
    public void testAnalyzerRefusesAModelOffTheShape() {
        Network model = plainPs();
        assertThrows(RuntimeException.class,
                () -> Solver_nc_dps_analyzer.solver_nc_dps_analyzer(model.getStruct(), null));
    }

    /** The gate must not disturb the ordinary product-form path. */
    @Test
    public void testDefaultStillSolvesThePsModel() throws Exception {
        Matrix q = new SolverNC(plainPs()).getAvgQLen();
        assertTrue(q.get(1, 0) > 0 && q.get(1, 1) > 0, "PS model must still be solved");
    }

    /** The kernel rejects inputs it has no derivation for rather than returning a number. */
    @Test
    public void testInvalidInputsRejected() {
        assertThrows(RuntimeException.class,
                () -> run(new double[]{100.0, 100.0}, M2_Z, M2_S, new double[]{1.0, -4.0}));
        assertThrows(RuntimeException.class,
                () -> run(new double[]{0.0, 100.0}, M2_Z, M2_S, M2_W));
        assertFalse(Double.isNaN(run(M2_N, M2_Z, M2_S, M2_W).Q.get(0, 0)));
    }
}
