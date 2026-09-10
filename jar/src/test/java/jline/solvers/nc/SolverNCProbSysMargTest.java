/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.nc;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.io.Ret.ProbabilityResult;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Exp;
import jline.lang.state.State;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.mva.SolverMVA;
import jline.api.pfqn.Pfqn_jointmarg;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Tests {@link SolverNC#getProbSysMarg(Matrix)}, the joint law of the per-station
 * TOTAL queue lengths evaluated through the Calame permanent identity.
 *
 * <p>The assertions run from weakest to strongest: the law must normalize, its
 * first moment must reproduce the exact CTMC queue lengths, and it must keep both
 * of those properties on the zero-bearing demand matrices (a class that skips a
 * station) and on several infinite servers, each of which contributes its own
 * 1/n_j!. The last two tests fix the refusals: an approximate permanent engine on
 * a structurally zero demand matrix, and the getter on a solver that does not
 * implement it.</p>
 */
public class SolverNCProbSysMargTest {

    private static final double TOL = 1e-9;

    @BeforeAll
    public static void setUp() {
        GlobalConstants.setVerbose(VerboseLevel.SILENT);
    }

    /** Delay + two PS queues, 2 classes, N = (2,1). Dense demand matrix. */
    private static Network denseModel() {
        Network model = new Network("marg_dense");
        Delay d = new Delay(model, "Delay");
        Queue q1 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue q2 = new Queue(model, "Queue2", SchedStrategy.PS);
        ClosedClass c1 = new ClosedClass(model, "Class1", 2, d, 0);
        ClosedClass c2 = new ClosedClass(model, "Class2", 1, d, 0);
        d.setService(c1, Exp.fitMean(0.7));
        d.setService(c2, Exp.fitMean(1.3));
        q1.setService(c1, Exp.fitMean(1.5));
        q1.setService(c2, Exp.fitMean(0.8));
        q2.setService(c1, Exp.fitMean(0.9));
        q2.setService(c2, Exp.fitMean(1.1));
        RoutingMatrix rm = model.initRoutingMatrix();
        rm.set(c1, c1, d, q1, 1.0);
        rm.set(c1, c1, q1, q2, 1.0);
        rm.set(c1, c1, q2, d, 1.0);
        rm.set(c2, c2, d, q1, 1.0);
        rm.set(c2, c2, q1, q2, 1.0);
        rm.set(c2, c2, q2, d, 1.0);
        model.link(rm);
        return model;
    }

    /**
     * Same three stations, but class 2 never visits Queue2, so its demand there is
     * a structural zero rather than a small number.
     */
    private static Network zeroDemandModel() {
        Network model = new Network("marg_zero");
        Delay d = new Delay(model, "Delay");
        Queue q1 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue q2 = new Queue(model, "Queue2", SchedStrategy.PS);
        ClosedClass c1 = new ClosedClass(model, "Class1", 2, d, 0);
        ClosedClass c2 = new ClosedClass(model, "Class2", 1, d, 0);
        d.setService(c1, Exp.fitMean(0.7));
        d.setService(c2, Exp.fitMean(1.3));
        q1.setService(c1, Exp.fitMean(1.5));
        q1.setService(c2, Exp.fitMean(0.8));
        q2.setService(c1, Exp.fitMean(0.9));
        q2.setService(c2, Exp.fitMean(1.1));
        RoutingMatrix rm = model.initRoutingMatrix();
        rm.set(c1, c1, d, q1, 1.0);
        rm.set(c1, c1, q1, q2, 1.0);
        rm.set(c1, c1, q2, d, 1.0);
        // Class 2 bypasses Queue2 entirely.
        rm.set(c2, c2, d, q1, 1.0);
        rm.set(c2, c2, q1, d, 1.0);
        model.link(rm);
        return model;
    }

    /** TWO infinite servers plus one queue: each delay carries its own 1/n_j!. */
    private static Network twoDelayModel() {
        Network model = new Network("marg_twodelay");
        Delay d1 = new Delay(model, "Delay1");
        Delay d2 = new Delay(model, "Delay2");
        Queue q = new Queue(model, "Queue1", SchedStrategy.PS);
        ClosedClass c1 = new ClosedClass(model, "Class1", 2, d1, 0);
        ClosedClass c2 = new ClosedClass(model, "Class2", 1, d1, 0);
        d1.setService(c1, Exp.fitMean(0.7));
        d1.setService(c2, Exp.fitMean(1.3));
        d2.setService(c1, Exp.fitMean(1.1));
        d2.setService(c2, Exp.fitMean(0.6));
        q.setService(c1, Exp.fitMean(1.5));
        q.setService(c2, Exp.fitMean(0.8));
        RoutingMatrix rm = model.initRoutingMatrix();
        rm.set(c1, c1, d1, d2, 1.0);
        rm.set(c1, c1, d2, q, 1.0);
        rm.set(c1, c1, q, d1, 1.0);
        rm.set(c2, c2, d1, d2, 1.0);
        rm.set(c2, c2, d2, q, 1.0);
        rm.set(c2, c2, q, d1, 1.0);
        model.link(rm);
        return model;
    }

    /** Every composition of {@code total} into {@code M} nonnegative parts. */
    private static List<int[]> compositions(int M, int total) {
        List<int[]> out = new ArrayList<int[]>();
        compositionsRec(M, total, new int[M], 0, out);
        return out;
    }

    private static void compositionsRec(int M, int left, int[] acc, int pos, List<int[]> out) {
        if (pos == M - 1) {
            acc[pos] = left;
            out.add(acc.clone());
            return;
        }
        for (int k = 0; k <= left; k++) {
            acc[pos] = k;
            compositionsRec(M, left - k, acc, pos + 1, out);
        }
    }

    private static Matrix rowVector(int[] v) {
        Matrix m = new Matrix(1, v.length);
        for (int i = 0; i < v.length; i++) {
            m.set(0, i, v[i]);
        }
        return m;
    }

    /** Sum over all total states, and the induced mean queue length per station. */
    private static double[] lawMoments(Network model, int M, int total, double[] mean) {
        SolverNC solver = new SolverNC(model);
        double mass = 0;
        for (int[] n : compositions(M, total)) {
            ProbabilityResult pr = solver.getProbSysMarg(rowVector(n));
            double p = pr.probability.get(0);
            assertTrue(p >= -TOL, "a probability must not be negative, was " + p + " at " + java.util.Arrays.toString(n));
            mass += p;
            for (int i = 0; i < M; i++) {
                mean[i] += n[i] * p;
            }
        }
        return new double[]{mass};
    }

    /** Station totals from the exact CTMC solution, classes summed out. */
    private static double[] ctmcTotals(Network model, int M) {
        Matrix qlen = new SolverCTMC(model).getAvgQLen();
        double[] out = new double[M];
        for (int i = 0; i < M; i++) {
            double s = 0;
            for (int r = 0; r < qlen.getNumCols(); r++) {
                s += qlen.get(i, r);
            }
            out[i] = s;
        }
        return out;
    }

    /**
     * The law must normalize and its first moment must be the exact one. Both are
     * checked on the same sweep: a normalizing bug and a per-state bug are
     * distinguishable only if the mean is tested too.
     */
    @Test
    public void testDenseModelNormalizesAndMatchesCTMC() {
        Network model = denseModel();
        double[] mean = new double[3];
        double mass = lawMoments(model, 3, 3, mean)[0];
        assertEquals(1.0, mass, 1e-12, "the total states must carry all the probability");

        double[] exact = ctmcTotals(denseModel(), 3);
        for (int i = 0; i < 3; i++) {
            assertEquals(exact[i], mean[i], 1e-9,
                    "E[n] at station " + (i + 1) + " must match the CTMC");
        }
    }

    /** A structural zero must not disturb either property under the exact engine. */
    @Test
    public void testZeroDemandNormalizesAndMatchesCTMC() {
        Network model = zeroDemandModel();
        double[] mean = new double[3];
        double mass = lawMoments(model, 3, 3, mean)[0];
        assertEquals(1.0, mass, 1e-12, "a zero demand must not leak probability");

        double[] exact = ctmcTotals(zeroDemandModel(), 3);
        for (int i = 0; i < 3; i++) {
            assertEquals(exact[i], mean[i], 1e-9,
                    "E[n] at station " + (i + 1) + " must match the CTMC with a zero demand");
        }
    }

    /**
     * With two infinite servers the identity needs one 1/n_j! per delay. Dividing
     * once, as a single aggregated delay row would, leaves the law unnormalized.
     */
    @Test
    public void testTwoDelaysNormalizeAndMatchCTMC() {
        Network model = twoDelayModel();
        double[] mean = new double[3];
        double mass = lawMoments(model, 3, 3, mean)[0];
        assertEquals(1.0, mass, 1e-12, "each infinite server must contribute its own 1/n_j!");

        double[] exact = ctmcTotals(twoDelayModel(), 3);
        for (int i = 0; i < 3; i++) {
            assertEquals(exact[i], mean[i], 1e-9,
                    "E[n] at station " + (i + 1) + " must match the CTMC with two delays");
        }
    }

    /** An empty network, n = 0 everywhere, is the only state when N = 0. */
    @Test
    public void testEmptyPopulationIsCertain() {
        Network model = new Network("marg_empty");
        Delay d = new Delay(model, "Delay");
        Queue q = new Queue(model, "Queue1", SchedStrategy.PS);
        ClosedClass c1 = new ClosedClass(model, "Class1", 2, d, 0);
        ClosedClass c2 = new ClosedClass(model, "Class2", 0, d, 0);
        d.setService(c1, Exp.fitMean(0.7));
        d.setService(c2, Exp.fitMean(1.3));
        q.setService(c1, Exp.fitMean(1.5));
        q.setService(c2, Exp.fitMean(0.8));
        RoutingMatrix rm = model.initRoutingMatrix();
        rm.set(c1, c1, d, q, 1.0);
        rm.set(c1, c1, q, d, 1.0);
        rm.set(c2, c2, d, q, 1.0);
        rm.set(c2, c2, q, d, 1.0);
        model.link(rm);

        SolverNC solver = new SolverNC(model);
        // A zero-population class leaves 2 jobs to place over 2 stations.
        double mass = 0;
        for (int[] n : compositions(2, 2)) {
            mass += solver.getProbSysMarg(rowVector(n)).probability.get(0);
        }
        assertEquals(1.0, mass, 1e-12, "a zero-population class must not change the law");

        // A state whose total does not match the population is impossible.
        assertEquals(0.0, solver.getProbSysMarg(rowVector(new int[]{1, 0})).probability.get(0), 0.0,
                "a state whose total differs from sum(N) has probability 0");
    }


    /** The 3-station 2-class demand matrix the saddle-point engine is calibrated on. */
    private static Matrix spmDemands() {
        double[][] v = {{0.286, 0.437}, {1.001, 0.782}, {0.294, 0.633}};
        Matrix L = new Matrix(3, 2);
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 2; j++) {
                L.set(i, j, v[i][j]);
            }
        }
        return L;
    }

    private static Matrix rowOf(int... v) {
        Matrix m = new Matrix(1, v.length);
        for (int i = 0; i < v.length; i++) {
            m.set(0, i, v[i]);
        }
        return m;
    }

    /**
     * The saddle point is asymptotic in the COLUMN multiplicities, which here are
     * the class populations, so its relative error must fall like 1/min(N) rather
     * than staying flat as the sampling and mean-field engines do.
     */
    @Test
    public void testSpmEngineErrorFallsAsThePopulationsGrow() {
        Matrix L = spmDemands();
        double previous = Double.POSITIVE_INFINITY;
        for (int k = 1; k <= 4; k++) {
            Matrix N = rowOf(k, k);
            double worst = 0.0, massExact = 0.0, massSpm = 0.0;
            for (int i = 0; i <= 2 * k; i++) {
                for (int j = 0; i + j <= 2 * k; j++) {
                    Matrix nvec = rowOf(i, j, 2 * k - i - j);
                    double ex = Pfqn_jointmarg.pfqn_jointmarg(nvec, L, N, null, null, "exact").pjoint;
                    double sp = Pfqn_jointmarg.pfqn_jointmarg(nvec, L, N, null, null, "spm").pjoint;
                    massExact += ex;
                    massSpm += sp;
                    if (ex > 1e-12) {
                        // the expansion overestimates the permanent, hence the probability
                        assertTrue(sp >= ex * (1.0 - 1e-9),
                                "the saddle point must not fall below the exact probability");
                        worst = Math.max(worst, Math.abs(sp - ex) / ex);
                    }
                }
            }
            assertEquals(1.0, massExact, 1e-9, "the exact law must normalize");
            assertTrue(worst < 0.2 / k, "the error must track 1/(8 min N), was " + worst);
            assertTrue(worst < previous, "the error must fall as the populations grow");
            assertTrue(massSpm > 1.0, "a uniform overestimate puts the mass above one");
            previous = worst;
        }
    }

    /**
     * The bias is nearly the same at every state, so a caller sweeping the lattice
     * and renormalizing to sum to one keeps an order of magnitude less of it.
     */
    @Test
    public void testSpmEngineBiasCancelsUnderRenormalization() {
        Matrix L = spmDemands();
        double previous = Double.POSITIVE_INFINITY;
        for (int k = 1; k <= 3; k++) {
            Matrix N = rowOf(k, k);
            List<Double> exact = new ArrayList<Double>();
            List<Double> spm = new ArrayList<Double>();
            for (int i = 0; i <= 2 * k; i++) {
                for (int j = 0; i + j <= 2 * k; j++) {
                    Matrix nvec = rowOf(i, j, 2 * k - i - j);
                    exact.add(Pfqn_jointmarg.pfqn_jointmarg(nvec, L, N, null, null, "exact").pjoint);
                    spm.add(Pfqn_jointmarg.pfqn_jointmarg(nvec, L, N, null, null, "spm").pjoint);
                }
            }
            double total = 0.0, raw = 0.0, tvd = 0.0;
            for (int t = 0; t < spm.size(); t++) {
                total += spm.get(t);
            }
            for (int t = 0; t < spm.size(); t++) {
                if (exact.get(t) > 1e-12) {
                    raw = Math.max(raw, Math.abs(spm.get(t) - exact.get(t)) / exact.get(t));
                }
                tvd += Math.abs(spm.get(t) / total - exact.get(t));
            }
            tvd *= 0.5;
            assertTrue(tvd < 0.1 * raw, "renormalization must remove most of the bias");
            assertTrue(tvd < previous, "the renormalized error must fall with the population");
            previous = tvd;
        }
    }

    /** The engine must be reachable by name from the solver getter. */
    @Test
    public void testSpmEngineIsReachableFromTheSolver() {
        SolverNC solver = new SolverNC(denseModel());
        Matrix nvec = rowVector(new int[]{1, 1, 1});
        double exact = solver.getProbSysMarg(nvec).probability.get(0);
        double approx = solver.getProbSysMarg(nvec, "spm").probability.get(0);
        assertTrue(approx > 0.0, "the saddle-point engine must answer");
        assertTrue(Math.abs(approx - exact) < 0.5 * exact,
                "the saddle-point engine must stay in the neighbourhood of the exact value");
    }

    /**
     * The approximate engines need full support: Sinkhorn stalls without it and
     * the Bethe gap is a state-dependent lower bound that does not cancel under
     * normalization. They must refuse by name rather than floor the zero.
     */
    @Test
    public void testApproximateEnginesRefuseAStructuralZero() {
        SolverNC solver = new SolverNC(zeroDemandModel());
        // Put at least one job at Queue2, so the zero-demand column reaches A.
        Matrix nvec = rowVector(new int[]{1, 1, 1});
        for (String engine : new String[]{"spm", "bethe", "heur", "huberlaw", "adapart"}) {
            RuntimeException e = assertThrows(RuntimeException.class,
                    () -> solver.getProbSysMarg(nvec, engine),
                    "engine '" + engine + "' must refuse a demand matrix with a structural zero");
            String msg = e.getMessage();
            assertNotNull(msg, "the refusal must carry a message");
            assertTrue(msg.contains("zero"),
                    "the refusal must name the zero entry, was: " + msg);
        }
        // The exact engine is unaffected on the same state.
        assertTrue(solver.getProbSysMarg(nvec).probability.get(0) > 0,
                "the exact engine must still answer on a zero-bearing matrix");
    }

    /** Delay + PS queue + FCFS queue, 2 classes: one shared server, one buffered. */
    private static Network startedModel() {
        Network model = new Network("started");
        Delay d = new Delay(model, "Delay");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.PS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.FCFS);
        ClosedClass c1 = new ClosedClass(model, "C1", 2, d, 0);
        ClosedClass c2 = new ClosedClass(model, "C2", 1, d, 0);
        d.setService(c1, Exp.fitMean(0.7));
        d.setService(c2, Exp.fitMean(1.3));
        q1.setService(c1, Exp.fitMean(1.5));
        q1.setService(c2, Exp.fitMean(0.8));
        q2.setService(c1, Exp.fitMean(0.9));
        q2.setService(c2, Exp.fitMean(1.1));
        RoutingMatrix rm = model.initRoutingMatrix();
        rm.set(c1, c1, d, q1, 1.0);
        rm.set(c1, c1, q1, q2, 1.0);
        rm.set(c1, c1, q2, d, 1.0);
        rm.set(c2, c2, d, q1, 1.0);
        rm.set(c2, c2, q1, q2, 1.0);
        rm.set(c2, c2, q2, d, 1.0);
        model.link(rm);
        return model;
    }

    /** A state space as a nested int list, so an expected set reads literally. */
    private static List<List<Integer>> rows(Matrix space) {
        List<List<Integer>> out = new ArrayList<List<Integer>>();
        for (int i = 0; i < space.getNumRows(); i++) {
            List<Integer> r = new ArrayList<Integer>();
            for (int j = 0; j < space.getNumCols(); j++) {
                r.add((int) Math.round(space.get(i, j)));
            }
            out.add(r);
        }
        return out;
    }

    private static List<Integer> row(int... v) {
        List<Integer> r = new ArrayList<Integer>();
        for (int x : v) {
            r.add(x);
        }
        return r;
    }

    /**
     * The union over every (n,s) pair consistent with the two totals, pinned
     * against the MATLAB reference. Class 2 has population 1, so the split
     * [0,2] must be dropped BEFORE the per-class builder is asked for it: an
     * empty local space is absorbed by the cartesian product rather than
     * annihilating it, so the job would silently disappear.
     */
    @Test
    public void testFromMargAndStartedMatchesTheReference() {
        Network model = startedModel();

        // A shared server holds every job present, so the state carries n, not s.
        assertEquals(Arrays.asList(row(2, 0), row(1, 1)),
                rows(State.fromMargAndStarted(model, 1, 2, 1)),
                "the PS station must carry the totals, and never the split [0,2]");

        // An ordered buffer carries the waiting class tags, the server block the
        // started counts in phase one: [buffer | srv(C1) | srv(C2)].
        assertEquals(Arrays.asList(row(2, 1, 0), row(1, 1, 0), row(1, 0, 1)),
                rows(State.fromMargAndStarted(model, 2, 2, 1)),
                "the FCFS station must emit the three (n,s) splits in reference order");

        // Both totals zero is the one empty state, at the discipline's own width.
        assertEquals(Arrays.asList(row(0, 0, 0)),
                rows(State.fromMargAndStarted(model, 2, 0, 0)),
                "an empty station keeps the width the decoder expects");
    }

    /** The class-summed marginal, with the same capped split enumeration. */
    @Test
    public void testFromMargMatchesTheReference() {
        Network model = startedModel();
        assertEquals(Arrays.asList(row(2, 0), row(1, 1)), rows(State.fromMarg(model, 1, 2)));
        assertEquals(Arrays.asList(row(0, 0)), rows(State.fromMarg(model, 1, 0)));
    }

    /** Every other solver refuses the getter with an actionable message. */
    @Test
    public void testOtherSolversRefuseTheGetter() {
        SolverMVA solver = new SolverMVA(denseModel());
        RuntimeException e = assertThrows(RuntimeException.class,
                () -> solver.getProbSysMarg(rowVector(new int[]{1, 1, 1})));
        assertTrue(e.getMessage().contains("getProbSysMarg is not supported"),
                "message must name the unsupported getter, was: " + e.getMessage());
    }

    /** A multiserver station has no n_i!, so the identity is refused by name. */
    @Test
    public void testMultiserverIsRefused() {
        Network model = new Network("marg_multiserver");
        Delay d = new Delay(model, "Delay");
        Queue q = new Queue(model, "Queue1", SchedStrategy.FCFS);
        q.setNumberOfServers(2);
        ClosedClass c1 = new ClosedClass(model, "Class1", 3, d, 0);
        d.setService(c1, Exp.fitMean(0.7));
        q.setService(c1, Exp.fitMean(1.5));
        RoutingMatrix rm = model.initRoutingMatrix();
        rm.set(c1, c1, d, q, 1.0);
        rm.set(c1, c1, q, d, 1.0);
        model.link(rm);

        SolverNC solver = new SolverNC(model);
        RuntimeException e = assertThrows(RuntimeException.class,
                () -> solver.getProbSysMarg(rowVector(new int[]{1, 2})));
        assertTrue(e.getMessage().contains("multiserver"),
                "message must name the multiserver station, was: " + e.getMessage());
    }
}
