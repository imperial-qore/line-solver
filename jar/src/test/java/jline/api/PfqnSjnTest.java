/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api;

import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import jline.GlobalConstants;
import jline.VerboseLevel;

import jline.api.pfqn.mva.Pfqn_amvasjn;
import jline.api.pfqn.mva.Pfqn_mvasjn;
import jline.api.pfqn.mva.Pfqn_bs;
import jline.api.pfqn.mva.Pfqn_mva;
import jline.api.pfqn.mva.SjnOptions;
import jline.api.pfqn.mva.SjnStarvationException;
import jline.io.Ret;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.JobClass;
import jline.lang.ClosedClass;
import jline.lang.processes.Exp;
import jline.lang.processes.HyperExp;
import jline.solvers.SolverOptions;
import jline.solvers.mva.SolverMVA;
import jline.util.matrix.Matrix;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Validates the shortest-job-next (SJN/SJF) station of Kant (1992):
 * {@link Pfqn_mvasjn} over the population lattice and {@link Pfqn_amvasjn} through its Schweitzer
 * fixed point.
 *
 * <p>The mathematics is anchored to results that do not come from the paper: a network with no SJN
 * station must reproduce exact MVA and Bard-Schweitzer, a single-job closed network cannot queue
 * at all, Little's law must hold at every population including where the utilization cap binds,
 * and SJN must beat any size-blind discipline on mean response time. The literal values are the
 * MATLAB outputs of pfqn_mvasjn.m and pfqn_amvasjn.m, so a divergence between the codebases shows
 * up here. Mirrors line-test.git/test/testsAPI/test_pfqn_sjn.m.</p>
 */
public class PfqnSjnTest {

    private static VerboseLevel savedVerbose;

    @BeforeAll
    public static void silenceStateWarning() {
        // The SJF stations warn on every initial state built (no discrete descriptor).
        savedVerbose = GlobalConstants.getVerbose();
        GlobalConstants.Verbose = VerboseLevel.SILENT;
    }

    @AfterAll
    public static void restoreVerbose() {
        GlobalConstants.Verbose = savedVerbose;
    }

    private static Matrix demands() {
        Matrix L = new Matrix(3, 1);
        L.set(0, 0, 0.125);
        L.set(1, 0, 0.100);
        L.set(2, 0, 0.050);
        return L;
    }

    private static Matrix scalar(double v) {
        Matrix m = new Matrix(1, 1);
        m.set(0, 0, v);
        return m;
    }

    private static Matrix scv(double cv2) {
        Matrix s = new Matrix(3, 1);
        s.set(0, 0, cv2);
        s.set(1, 0, 1.0);
        s.set(2, 0, 1.0);
        return s;
    }

    private static final int[] SJN0 = new int[] {0};

    @Test
    public void testNoSjnMatchesExactMva() {
        // With no SJN station the recursion is the Reiser-Lavenberg one.
        Matrix L = demands();
        int[] pops = new int[] {1, 4, 9};
        for (int p = 0; p < pops.length; p++) {
            Pfqn_mvasjn.Result r = Pfqn_mvasjn.pfqn_mvasjn(L, scalar(pops[p]), scalar(1.0), null,
                    null, null, null);
            Ret.pfqnMVA ref = Pfqn_mva.pfqn_mva(L, scalar(pops[p]), scalar(1.0));
            assertEquals(ref.X.get(0), r.X.get(0, 0), 1e-11);
            for (int m = 0; m < 3; m++) {
                assertEquals(ref.Q.get(m), r.Q.get(m, 0), 1e-11);
            }
        }
    }

    @Test
    public void testNoSjnMatchesBardSchweitzer() {
        // With no SJN station the fixed point is the Bard-Schweitzer one.
        Matrix L = demands();
        Pfqn_mvasjn.Result r = Pfqn_amvasjn.pfqn_amvasjn(L, scalar(8), scalar(1.0), null, null,
                null, null);
        Ret.pfqnAMVA ref = Pfqn_bs.pfqn_bs(L, scalar(8), scalar(1.0));
        assertEquals(ref.X.get(0), r.X.get(0, 0), 1e-5);
        for (int m = 0; m < 3; m++) {
            assertEquals(ref.Q.get(m), r.Q.get(m, 0), 1e-5);
        }
    }

    @Test
    public void testSingleJobCannotQueue() {
        // One job in a closed network never waits, whatever the size distribution.
        Matrix L = demands();
        double[] cv2s = new double[] {0.25, 1.0, 4.0};
        for (int c = 0; c < cv2s.length; c++) {
            Pfqn_mvasjn.Result a = Pfqn_mvasjn.pfqn_mvasjn(L, scalar(1), scalar(2.0), scv(cv2s[c]),
                    SJN0, null, null);
            Pfqn_mvasjn.Result b = Pfqn_amvasjn.pfqn_amvasjn(L, scalar(1), scalar(2.0),
                    scv(cv2s[c]), SJN0, null, null);
            double X = 1.0 / (2.0 + 0.275);
            assertEquals(X, a.X.get(0, 0), 1e-10);
            assertEquals(X, b.X.get(0, 0), 1e-10);
            for (int m = 0; m < 3; m++) {
                assertEquals(L.get(m), a.C.get(m, 0), 1e-10);
                assertEquals(L.get(m), b.C.get(m, 0), 1e-10);
            }
        }
    }

    @Test
    public void testPopulationConservation() {
        // Little's law over the whole network, including where the cap binds. The cap acts on the
        // waiting time and the throughput is derived from it, so no jobs are lost. Only the fixed
        // point may cap: the lattice refuses the starvation regime instead.
        Matrix L = demands();
        double[] cv2s = new double[] {0.5, 1.0, 2.0};
        int[] pops = new int[] {2, 7, 15, 40};
        for (int c = 0; c < cv2s.length; c++) {
            for (int p = 0; p < pops.length; p++) {
                final int pop = pops[p];
                final double cv2 = cv2s[c];
                if (pop == 40) {
                    assertThrows(SjnStarvationException.class, () -> Pfqn_mvasjn.pfqn_mvasjn(L,
                            scalar(pop), scalar(1.0), scv(cv2), SJN0, null, null));
                } else {
                    Pfqn_mvasjn.Result a = Pfqn_mvasjn.pfqn_mvasjn(L, scalar(pop), scalar(1.0),
                            scv(cv2), SJN0, null, null);
                    double asum = a.X.get(0, 0) * 1.0;
                    for (int m = 0; m < 3; m++) {
                        asum += a.Q.get(m, 0);
                    }
                    assertEquals(pop, asum, 1e-8);
                }
                Pfqn_mvasjn.Result b = Pfqn_amvasjn.pfqn_amvasjn(L, scalar(pops[p]), scalar(1.0),
                        scv(cv2s[c]), SJN0, null, null);
                double sum = b.X.get(0, 0) * 1.0;
                for (int m = 0; m < 3; m++) {
                    sum += b.Q.get(m, 0);
                }
                assertEquals(pops[p], sum, 1e-8);
            }
        }
    }

    @Test
    public void testSjnBeatsSizeBlindScheduling() {
        // SJN minimises the mean response time among non-preemptive disciplines. The comparison is
        // only meaningful at CV^2 = 1, where the product-form solution describes the SAME service
        // distribution: product form is insensitive.
        Matrix L = demands();
        int[] pops = new int[] {4, 6, 9};
        for (int p = 0; p < pops.length; p++) {
            Pfqn_mvasjn.Result a = Pfqn_mvasjn.pfqn_mvasjn(L, scalar(pops[p]), scalar(1.0), scv(1.0),
                    SJN0, null, null);
            Ret.pfqnMVA ref = Pfqn_mva.pfqn_mva(L, scalar(pops[p]), scalar(1.0));
            assertTrue(a.X.get(0, 0) > ref.X.get(0), "SJN must carry more work than product form");
            assertTrue(a.C.get(0, 0) < ref.R.get(0), "SJN must respond faster than product form");
        }
    }

    @Test
    public void testMatlabParity() {
        // Literal values of the MATLAB pfqn_mvasjn.m and pfqn_amvasjn.m at N = 6, Z = 1.
        Matrix L = demands();
        double[] cv2s = new double[] {0.5, 1.0, 2.0, 4.0};
        double[] Xref = new double[] {4.298002, 4.256655, 4.172968, 4.021977};
        double[] Rref = new double[] {0.184871, 0.199049, 0.228550, 0.284778};
        for (int c = 0; c < cv2s.length; c++) {
            Pfqn_mvasjn.Result a = Pfqn_mvasjn.pfqn_mvasjn(L, scalar(6), scalar(1.0), scv(cv2s[c]),
                    SJN0, null, null);
            assertEquals(Xref[c], a.X.get(0, 0), 1e-5);
            assertEquals(Rref[c], a.C.get(0, 0), 1e-5);
        }
        Pfqn_mvasjn.Result b = Pfqn_amvasjn.pfqn_amvasjn(L, scalar(6), scalar(1.0), scv(2.0), SJN0,
                null, null);
        assertEquals(4.162429, b.X.get(0, 0), 1e-5);
        assertEquals(0.227866, b.C.get(0, 0), 1e-5);
    }

    @Test
    public void testUtilizationCapIsEnforced() {
        // A saturated SJN station has no solution to the open-form response time equation. The
        // utilization law caps it in the FIXED POINT, which re-derives the profile from the capped
        // state and converges, so it returns instead of failing. The lattice cannot: the same
        // factor would rescale the conditional waiting time profile its next population step reads
        // back, so it refuses and the analyzer re-solves with the fixed point.
        Matrix L = demands();
        int[] pops = new int[] {24, 40, 80};
        for (int p = 0; p < pops.length; p++) {
            Pfqn_mvasjn.Result a = Pfqn_amvasjn.pfqn_amvasjn(L, scalar(pops[p]), scalar(1.0),
                    scv(2.0), SJN0, null, null);
            assertTrue(a.U.get(0, 0) <= 0.999 + 1e-9, "the cap must hold at N=" + pops[p]);
            final int pop = pops[p];
            assertThrows(SjnStarvationException.class, () -> Pfqn_mvasjn.pfqn_mvasjn(L,
                    scalar(pop), scalar(1.0), scv(2.0), SJN0, null, null));
        }
        Pfqn_mvasjn.Result a24 = Pfqn_amvasjn.pfqn_amvasjn(L, scalar(24), scalar(1.0), scv(2.0),
                SJN0, null, null);
        assertEquals(7.992000, a24.X.get(0, 0), 1e-5);
        assertEquals(0.999000, a24.U.get(0, 0), 1e-8);
        SjnOptions opt = new SjnOptions();
        opt.umax = 0.9;
        Pfqn_mvasjn.Result low = Pfqn_amvasjn.pfqn_amvasjn(L, scalar(40), scalar(1.0), scv(2.0),
                SJN0, null, opt);
        assertTrue(low.U.get(0, 0) <= 0.9 + 1e-9);
    }

    @Test
    public void testFixedPointTracksTheLattice() {
        // Away from saturation the Schweitzer closure must stay close to the lattice it replaces.
        Matrix L = demands();
        double[] cv2s = new double[] {0.5, 1.0, 2.0, 4.0};
        int[] pops = new int[] {3, 6};
        for (int c = 0; c < cv2s.length; c++) {
            for (int p = 0; p < pops.length; p++) {
                Pfqn_mvasjn.Result a = Pfqn_mvasjn.pfqn_mvasjn(L, scalar(pops[p]), scalar(1.0),
                        scv(cv2s[c]), SJN0, null, null);
                Pfqn_mvasjn.Result b = Pfqn_amvasjn.pfqn_amvasjn(L, scalar(pops[p]), scalar(1.0),
                        scv(cv2s[c]), SJN0, null, null);
                assertEquals(a.X.get(0, 0), b.X.get(0, 0), 0.01 * a.X.get(0, 0));
                assertEquals(a.C.get(0, 0), b.C.get(0, 0), 0.01 * a.C.get(0, 0));
            }
        }
    }

    @Test
    public void testGridRefinementConverges() {
        // The quadrature must be converged at the default grid.
        Matrix L = demands();
        SjnOptions o32 = new SjnOptions();
        o32.ns = 32;
        SjnOptions o64 = new SjnOptions();
        o64.ns = 64;
        Pfqn_mvasjn.Result a = Pfqn_mvasjn.pfqn_mvasjn(L, scalar(6), scalar(1.0), scv(2.0), SJN0,
                null, o32);
        Pfqn_mvasjn.Result b = Pfqn_mvasjn.pfqn_mvasjn(L, scalar(6), scalar(1.0), scv(2.0), SJN0,
                null, o64);
        assertEquals(a.X.get(0, 0), b.X.get(0, 0), 1e-4 * a.X.get(0, 0));
        assertEquals(a.C.get(0, 0), b.C.get(0, 0), 1e-4 * a.C.get(0, 0));
    }

    @Test
    public void testOddGridRejected() {
        // Composite Simpson integrates over panels of two subdivisions.
        SjnOptions odd = new SjnOptions();
        odd.ns = 31;
        assertThrows(IllegalArgumentException.class, new org.junit.jupiter.api.function.Executable() {
            @Override
            public void execute() {
                Pfqn_mvasjn.pfqn_mvasjn(demands(), scalar(4), scalar(1.0), scv(1.0), SJN0, null,
                        odd);
            }
        });
    }

    /** Closed model with one SJF queue, shared by the dispatch tests. */
    private static Network dispatchModel() {
        Network model = new Network("sjn_dispatch");
        Delay delay = new Delay(model, "Think");
        Queue q1 = new Queue(model, "SJN", SchedStrategy.SJF);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.PS);
        ClosedClass cl = new ClosedClass(model, "C1", 6, delay);
        delay.setService(cl, new Exp(1.0));
        q1.setService(cl, HyperExp.fitMeanAndSCVBalanced(0.125, 2.0));
        q2.setService(cl, new Exp(1 / 0.100));
        model.link(Network.serialRouting(cl, delay, q1, q2));
        return model;
    }

    private static Matrix dispatchDemands() {
        Matrix L = new Matrix(2, 1);
        L.set(0, 0, 0.125);
        L.set(1, 0, 0.100);
        return L;
    }

    private static Matrix dispatchScv() {
        Matrix s = new Matrix(2, 1);
        s.set(0, 0, 2.0);
        s.set(1, 0, 1.0);
        return s;
    }

    @Test
    public void testSolverMvaAmvaMethod() {
        // Asking for 'amva' must reach the fixed point, not the lattice.
        SolverMVA solver = new SolverMVA(dispatchModel(), "amva");
        Matrix Ts = solver.getAvgTput();
        Matrix Rs = solver.getAvgRespT();
        // the analyzer hands options.iter_tol to the fixed point, so the direct call must use it
        SjnOptions opt = new SjnOptions();
        opt.tol = SolverMVA.defaultOptions().iter_tol;
        opt.iterMax = SolverMVA.defaultOptions().iter_max;
        Pfqn_mvasjn.Result f = Pfqn_amvasjn.pfqn_amvasjn(dispatchDemands(), scalar(6), scalar(1.0),
                dispatchScv(), SJN0, null, opt);
        assertEquals(f.X.get(0, 0), Ts.get(1, 0), 1e-8);
        assertEquals(f.C.get(0, 0), Rs.get(1, 0), 1e-8);
        // the two routes must differ, or the test would pass on the lattice
        Pfqn_mvasjn.Result e = Pfqn_mvasjn.pfqn_mvasjn(dispatchDemands(), scalar(6), scalar(1.0),
                dispatchScv(), SJN0, null, null);
        assertNotEquals(e.X.get(0, 0), f.X.get(0, 0));
        assertNotEquals(e.C.get(0, 0), f.C.get(0, 0));
    }

    @Test
    public void testSolverMvaLatticeThreshold() {
        // Past config.sjn_lattice_max the default dispatch switches to the fixed point.
        SolverOptions big = SolverMVA.defaultOptions();
        big.config.sjn_lattice_max = 1e5;
        SolverMVA sbig = new SolverMVA(dispatchModel(), big);
        Pfqn_mvasjn.Result e = Pfqn_mvasjn.pfqn_mvasjn(dispatchDemands(), scalar(6), scalar(1.0),
                dispatchScv(), SJN0, null, null);
        assertEquals(e.X.get(0, 0), sbig.getAvgTput().get(1, 0), 1e-8);
        assertEquals(e.C.get(0, 0), sbig.getAvgRespT().get(1, 0), 1e-8);
        SolverOptions small = SolverMVA.defaultOptions();
        small.config.sjn_lattice_max = 2.0;   // the 7-state lattice no longer fits
        SolverMVA ssmall = new SolverMVA(dispatchModel(), small);
        SjnOptions opt = new SjnOptions();
        opt.tol = small.iter_tol;
        opt.iterMax = small.iter_max;
        Pfqn_mvasjn.Result f = Pfqn_amvasjn.pfqn_amvasjn(dispatchDemands(), scalar(6), scalar(1.0),
                dispatchScv(), SJN0, null, opt);
        assertEquals(f.X.get(0, 0), ssmall.getAvgTput().get(1, 0), 1e-8);
        assertEquals(f.C.get(0, 0), ssmall.getAvgRespT().get(1, 0), 1e-8);
    }

    @Test
    public void testMulticlassPopulationConservation() {
        // The same invariant per class in a two-class model, pooled and prioritised.
        Matrix L = new Matrix(2, 2);
        L.set(0, 0, 0.125);
        L.set(0, 1, 0.060);
        L.set(1, 0, 0.100);
        L.set(1, 1, 0.080);
        Matrix N = new Matrix(1, 2);
        N.set(0, 0, 4);
        N.set(0, 1, 3);
        Matrix Z = new Matrix(1, 2);
        Z.set(0, 0, 1);
        Z.set(0, 1, 1);
        Matrix s = new Matrix(2, 2);
        s.set(0, 0, 2);
        s.set(0, 1, 1);
        s.set(1, 0, 1);
        s.set(1, 1, 1);
        Pfqn_mvasjn.Result pooled = Pfqn_mvasjn.pfqn_mvasjn(L, N, Z, s, SJN0, null, null);
        SjnOptions prio = new SjnOptions();
        prio.prio = new int[] {1, 2};
        Pfqn_mvasjn.Result byprio = Pfqn_mvasjn.pfqn_mvasjn(L, N, Z, s, SJN0, null, prio);
        for (int r = 0; r < 2; r++) {
            double sp = pooled.X.get(0, r) * Z.get(0, r);
            double sq = byprio.X.get(0, r) * Z.get(0, r);
            for (int m = 0; m < 2; m++) {
                sp += pooled.Q.get(m, r);
                sq += byprio.Q.get(m, r);
            }
            assertEquals(N.get(0, r), sp, 1e-8);
            assertEquals(N.get(0, r), sq, 1e-8);
        }
    }

    /**
     * Two classes in one chain, both traversing every station and switching at Q2.
     *
     * <p>Disjoint class routing would make ST(i,k) equal STchain(i,c) at every station with
     * nonzero alpha, so the two deaggregation routes would coincide and the model could not
     * tell them apart.</p>
     */
    private static Network multiclassChainModel() {
        Network model = new Network("sjn_multiclass_chain");
        Delay think = new Delay(model, "Think");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.SJF);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.PS);
        q1.setNumberOfServers(1);
        q2.setNumberOfServers(1);
        ClosedClass ca = new ClosedClass(model, "A", 4, think, 0);
        ClosedClass cb = new ClosedClass(model, "B", 0, think, 0);
        think.setService(ca, Exp.fitMean(1.0));
        think.setService(cb, Exp.fitMean(2.0));
        q1.setService(ca, Exp.fitMean(0.5));
        q1.setService(cb, Exp.fitMean(0.9));
        q2.setService(ca, Exp.fitMean(0.7));
        q2.setService(cb, Exp.fitMean(1.25));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(ca, ca, think, q1, 1.0);
        P.set(ca, ca, q1, q2, 1.0);
        P.set(ca, cb, q2, think, 1.0);
        P.set(cb, cb, think, q1, 1.0);
        P.set(cb, cb, q1, q2, 1.0);
        P.set(cb, ca, q2, think, 1.0);
        model.link(P);
        return model;
    }

    @Test
    public void testMulticlassChainDeaggregation() {
        // On a multiclass chain Q and U must be rebuilt from Rchain, not scaled from Qchain.
        // solver_mva_sjn_analyzer.m:138 passes [] in the Qchain and Uchain slots, so
        // sn_deaggregate_chain_results rebuilds them per class. Passing the chain matrices
        // instead scales by alpha alone and drops the ST(i,k)/STchain(i,c) weighting, which
        // collapses both classes onto the chain average.
        //
        // The literals are the MATLAB ground truth for this model, which has
        // STchain = [1.5; 0.7; 0.975] against class demands [1.0 2.0], [0.5 0.9] and
        // [0.7 1.25], so no station has ST(i,k) equal to STchain(i,c).
        SolverMVA solver = new SolverMVA(multiclassChainModel(), "sjn.mva");
        Matrix Q = solver.getAvgQLen();
        Matrix R = solver.getAvgRespT();
        // a delay station holds each class for its own mean service time, never the chain
        // average, so this alone separates the two deaggregation routes
        assertEquals(1.0, R.get(0, 0), 1e-9);
        assertEquals(2.0, R.get(0, 1), 1e-9);
        assertEquals(0.4161596548, Q.get(0, 0), 1e-8);
        assertEquals(0.8323193095, Q.get(0, 1), 1e-8);
        assertEquals(0.3480101577, Q.get(1, 0), 1e-8);
        assertEquals(0.6264182838, Q.get(1, 1), 1e-8);
        assertEquals(0.6379306748, Q.get(2, 0), 1e-8);
        assertEquals(1.1391619194, Q.get(2, 1), 1e-8);
        double total = 0;
        for (int i = 0; i < Q.getNumRows(); i++) {
            for (int r = 0; r < Q.getNumCols(); r++) {
                total += Q.get(i, r);
            }
        }
        assertEquals(4.0, total, 1e-8);
        // the classes must not share a queue length, which is what the chain-matrix route gives
        assertNotEquals(Q.get(1, 0), Q.get(1, 1), 1e-6);
    }

    @Test
    public void testMultiServerQueueIsRefused() {
        // solver_mva_sjn_analyzer.m:44-48 refuses a PS station with nservers ~= 1, because the
        // remaining stations are solved with the single-server MVA equation.
        //
        // This does not separate the two station-classification routes in Java, only in MATLAB
        // and Python. The distinguishing case is nservers infinite with sched still PS, which
        // Java cannot build: Network.java:2931-2933 rewrites sched to INF whenever nservers is
        // Integer.MAX_VALUE, so the struct never presents that pair. MATLAB keeps PS with
        // nservers Inf and refuses at :46.
        final Network model = new Network("sjn_multiserver_ps");
        Delay think = new Delay(model, "Think");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.SJF);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.PS);
        q1.setNumberOfServers(1);
        q2.setNumberOfServers(2);
        ClosedClass cl = new ClosedClass(model, "A", 3, think, 0);
        think.setService(cl, Exp.fitMean(1.0));
        q1.setService(cl, Exp.fitMean(0.5));
        q2.setService(cl, Exp.fitMean(0.7));
        model.link(Network.serialRouting(cl, think, q1, q2));
        assertThrows(RuntimeException.class, new org.junit.jupiter.api.function.Executable() {
            public void execute() {
                new SolverMVA(model, "sjn.mva").getAvgQLen();
            }
        });
    }

    @Test
    public void testEmptyChainMetricsAreZero() {
        // A chain with no jobs carries no load, per solver_mva_sjn_analyzer.m:131-135.
        // This is a regression guard, not evidence for the chain-level zeroing: it passes
        // without it too, because an empty chain has Xchain == 0 and the final non-finite
        // sweep in snDeaggregateChainResults already flattens the NaN that Qchain / Tchain
        // leaves behind. The analyzer still zeroes explicitly, to match MATLAB and to keep a
        // NaN out of the deaggregation input.
        Network model = new Network("sjn_empty_chain");
        Delay think = new Delay(model, "Think");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.SJF);
        q1.setNumberOfServers(1);
        ClosedClass ca = new ClosedClass(model, "A", 3, think, 0);
        ClosedClass cb = new ClosedClass(model, "B", 0, think, 0);
        think.setService(ca, Exp.fitMean(1.0));
        think.setService(cb, Exp.fitMean(2.0));
        q1.setService(ca, Exp.fitMean(0.5));
        q1.setService(cb, Exp.fitMean(0.9));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(ca, ca, think, q1, 1.0);
        P.set(ca, ca, q1, think, 1.0);
        P.set(cb, cb, think, q1, 1.0);
        P.set(cb, cb, q1, think, 1.0);
        model.link(P);
        SolverMVA solver = new SolverMVA(model, "sjn.mva");
        Matrix Q = solver.getAvgQLen();
        Matrix U = solver.getAvgUtil();
        Matrix T = solver.getAvgTput();
        double totalA = 0;
        for (int i = 0; i < Q.getNumRows(); i++) {
            assertEquals(0.0, Q.get(i, 1), 0.0);
            assertEquals(0.0, U.get(i, 1), 0.0);
            assertEquals(0.0, T.get(i, 1), 0.0);
            totalA += Q.get(i, 0);
        }
        assertEquals(3.0, totalA, 1e-8);
    }
}
