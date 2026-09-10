/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.fj;

import jline.VerboseLevel;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Fork;
import jline.lang.nodes.Join;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.lang.constant.SolverType;
import jline.solvers.Solver;
import jline.solvers.SolverOptions;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.mva.SolverMVA;
import jline.solvers.nc.SolverNC;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Regression tests for the solver-agnostic fork-join driver
 * ({@link FJFixedPoint}) and for the SolverNC fork-join route built on it.
 *
 * <p>The MMT transformation replaces every fork by a router and every join by a
 * zero-service delay, carrying the parallelism on auxiliary open classes. The
 * fixed point around it consumes only the metric matrices of an inner solve, so
 * it was extracted out of MVARunner and is now driven by both MVA and NC. The
 * tests assert (i) that the extraction left SolverMVA unchanged on fork-join
 * and on ordinary models, (ii) that SolverNC now solves fork-join models and
 * agrees with SolverMVA where both are exact, and (iii) that the closed-model
 * route neither degenerates on the auxiliary classes nor drifts from the MVA
 * route, measured against the exact SolverCTMC solution of the same model.
 *
 * <p>Java mirror of line-test.git/test/testsFJ/test_fj_driver_nc.m.
 */
public class FJDriverNCTest {

    /**
     * MVA/NC defaults, not the generic ones. Solver.defaultOptions() leaves
     * iter_tol at the generic 1e-4, while the MVA branch of SolverOptions sets
     * 1e-6; the fork-join fixed point stops on iter_tol, so the generic value
     * left the MMT iterate ~1e-4 short and the JAR reported queue lengths that
     * MATLAB, native Python and this same JAR under its own MVA defaults do not
     * produce. NC inherits the generic numeric controls, so one helper serves
     * both and keeps the three codebases on a single golden.
     */
    private static SolverOptions silent() {
        SolverOptions options = new SolverOptions(SolverType.MVA);
        options.verbose = VerboseLevel.SILENT;
        return options;
    }

    /**
     * CTMC needs its own solver defaults, not the generic ones: only the CTMC
     * branch of SolverOptions sets config.hide_immediate, and without it the
     * vanishing states of the tag-augmented copy are kept in the chain with a
     * large-but-finite Immediate rate, biasing every marginal by O(mu/Immediate).
     */
    private static SolverOptions silentCtmc() {
        SolverOptions options = new SolverOptions(SolverType.CTMC);
        options.verbose = VerboseLevel.SILENT;
        return options;
    }

    private static Network buildOpenFJ() {
        Network model = new Network("fjopen");
        Source source = new Source(model, "Source");
        Fork fork = new Fork(model, "Fork");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.FCFS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.FCFS);
        Join join = new Join(model, "Join", fork);
        Sink sink = new Sink(model, "Sink");
        OpenClass jobclass = new OpenClass(model, "C1", 0);
        source.setArrival(jobclass, new Exp(0.5));
        q1.setService(jobclass, new Exp(2.0));
        q2.setService(jobclass, new Exp(3.0));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobclass, jobclass, source, fork, 1.0);
        P.set(jobclass, jobclass, fork, q1, 1.0);
        P.set(jobclass, jobclass, fork, q2, 1.0);
        P.set(jobclass, jobclass, q1, join, 1.0);
        P.set(jobclass, jobclass, q2, join, 1.0);
        P.set(jobclass, jobclass, join, sink, 1.0);
        model.link(P);
        return model;
    }

    private static Network buildClosedFJ(int N) {
        Network model = new Network("fjclosed");
        Delay delay = new Delay(model, "Think");
        Fork fork = new Fork(model, "Fork");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.PS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.PS);
        Join join = new Join(model, "Join", fork);
        ClosedClass jobclass = new ClosedClass(model, "C1", N, delay);
        delay.setService(jobclass, new Exp(1.0));
        q1.setService(jobclass, new Exp(2.0));
        q2.setService(jobclass, new Exp(2.5));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobclass, jobclass, delay, fork, 1.0);
        P.set(jobclass, jobclass, fork, q1, 1.0);
        P.set(jobclass, jobclass, fork, q2, 1.0);
        P.set(jobclass, jobclass, q1, join, 1.0);
        P.set(jobclass, jobclass, q2, join, 1.0);
        P.set(jobclass, jobclass, join, delay, 1.0);
        model.link(P);
        return model;
    }

    private static Network buildPlainCQN() {
        Network model = new Network("cqn");
        Delay delay = new Delay(model, "Think");
        Queue queue = new Queue(model, "Q1", SchedStrategy.PS);
        ClosedClass jobclass = new ClosedClass(model, "C1", 3, delay);
        delay.setService(jobclass, new Exp(1.0));
        queue.setService(jobclass, new Exp(2.0));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobclass, jobclass, model.serialRouting(delay, queue));
        model.link(P);
        return model;
    }

    private static double relErr(Matrix a, Matrix b, int rows) {
        double num = 0;
        double den = 0;
        for (int i = 0; i < rows; i++) {
            double d = a.get(i, 0) - b.get(i, 0);
            num += d * d;
            den += b.get(i, 0) * b.get(i, 0);
        }
        return Math.sqrt(num) / Math.sqrt(den);
    }

    @Test
    public void testMvaForkJoinValuesUnchanged() throws Exception {
        // Golden values recorded before the fixed point was moved out of
        // MVARunner into jline.solvers.fj.FJFixedPoint
        SolverMVA open = new SolverMVA(buildOpenFJ(), silent());
        Matrix Q = open.getAvgQLen();
        Matrix T = open.getAvgTput();
        // Refreshed 2026-07-22 together with the MATLAB and Python mirrors: with
        // MVA's own iter_tol (see silent()) the JAR reproduces the MATLAB values
        // to 1e-13, so all three codebases now carry the same golden.
        double[] Qgold = {0.0, 0.333332200208543, 0.199999847107973, 0.283332462218689};
        double[] Tgold = {0.5, 0.499999880790713, 0.499999880790713, 0.5};
        for (int i = 0; i < 4; i++) {
            assertEquals(Qgold[i], Q.get(i, 0), 1e-9, "QN row " + i);
            assertEquals(Tgold[i], T.get(i, 0), 1e-9, "TN row " + i);
        }

        SolverMVA closed = new SolverMVA(buildClosedFJ(4), silent());
        Matrix Q4 = closed.getAvgQLen();
        Matrix T4 = closed.getAvgTput();
        double[] Q4gold = {1.43153004164563, 2.18851387932703, 1.23958258137759, 1.84538508086124};
        double[] T4gold = {1.4315299485931, 1.43152994188856, 1.43152994188856, 1.4315299485931};
        for (int i = 0; i < 4; i++) {
            assertEquals(Q4gold[i], Q4.get(i, 0), 1e-8, "closed QN row " + i);
            assertEquals(T4gold[i], T4.get(i, 0), 1e-8, "closed TN row " + i);
        }
    }

    @Test
    public void testDriverIsTransparentWithoutForks() throws Exception {
        // With no Fork in the model the driver must run the inner solve exactly
        // once and return it untouched, for both solvers
        Matrix Qm = new SolverMVA(buildPlainCQN(), silent()).getAvgQLen();
        SolverOptions exact = silent();
        exact.method = "exact";
        Matrix Qe = new SolverMVA(buildPlainCQN(), exact).getAvgQLen();
        Matrix Qn = new SolverNC(buildPlainCQN(), silent()).getAvgQLen();
        for (int i = 0; i < Qe.getNumRows(); i++) {
            assertEquals(Qe.get(i, 0), Qm.get(i, 0), 1e-6 * Math.abs(Qe.get(i, 0)) + 1e-9);
            assertEquals(Qe.get(i, 0), Qn.get(i, 0), 1e-6 * Math.abs(Qe.get(i, 0)) + 1e-9);
        }
    }

    @Test
    public void testNcAcceptsForkJoin() {
        Network open = buildOpenFJ();
        Network closed = buildClosedFJ(3);
        assertTrue(new SolverNC(open, silent()).supports(open));
        assertTrue(new SolverNC(closed, silent()).supports(closed));
    }

    @Test
    public void testCtmcRejectsClassSwitchOnAForkOutputEdge() {
        // A class switch declared ON a fork output edge is encoded by link() as an
        // auto-inserted ClassSwitch node, which lives in the sn index space but not
        // in the user-supplied routing matrix that fjtag writes sibling routing
        // into. It must be refused with the documented message, not with a raw
        // out-of-bounds error. Model of matlab/examples/basic/forkJoin/fj_cs_postfork.m.
        Network model = new Network("cspostfork");
        Delay delay = new Delay(model, "Delay");
        Fork fork = new Fork(model, "Fork1");
        Join join = new Join(model, "Join1", fork);
        Queue q1 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue q2 = new Queue(model, "Queue2", SchedStrategy.PS);
        ClosedClass c1 = new ClosedClass(model, "class1", 1, delay);
        ClosedClass c2 = new ClosedClass(model, "class2", 1, delay);
        delay.setService(c1, new Exp(0.25));
        q1.setService(c1, new Exp(2.0));
        q2.setService(c1, new Exp(2.0));
        delay.setService(c2, new Exp(0.25));
        q1.setService(c2, new Exp(2.0));
        q2.setService(c2, new Exp(2.0));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(c1, c1, delay, fork, 1.0);
        P.set(c1, c1, fork, q1, 1.0);
        P.set(c1, c1, fork, q2, 1.0);
        P.set(c1, c1, q1, join, 1.0);
        P.set(c1, c1, q2, join, 1.0);
        P.set(c1, c1, join, delay, 1.0);
        P.set(c2, c2, delay, fork, 1.0);
        P.set(c2, c1, fork, q1, 1.0);
        P.set(c2, c1, fork, q2, 1.0);
        model.link(P);
        SolverOptions options = silentCtmc();
        Exception thrown = assertThrows(Exception.class,
            new org.junit.jupiter.api.function.Executable() {
                @Override public void execute() throws Throwable {
                    new SolverCTMC(model, options).getAvgQLen();
                }
            });
        assertTrue(thrown.getMessage() != null
                && thrown.getMessage().contains("Class switching between fork and join"),
            "expected the guarded refusal, got: " + thrown.getMessage());
    }

    @Test
    public void testNcOpenForkJoinMatchesMva() throws Exception {
        // On an open transformed model both solvers use exact open-network
        // formulas inside the same fixed point, so they must agree closely
        SolverMVA mva = new SolverMVA(buildOpenFJ(), silent());
        SolverNC nc = new SolverNC(buildOpenFJ(), silent());
        Matrix Qm = mva.getAvgQLen();
        Matrix Tm = mva.getAvgTput();
        Matrix Qn = nc.getAvgQLen();
        Matrix Tn = nc.getAvgTput();
        for (int i = 0; i < Qm.getNumRows(); i++) {
            assertEquals(Qm.get(i, 0), Qn.get(i, 0), 1e-3 * Math.abs(Qm.get(i, 0)) + 1e-8, "QN row " + i);
            assertEquals(Tm.get(i, 0), Tn.get(i, 0), 1e-3 * Math.abs(Tm.get(i, 0)) + 1e-8, "TN row " + i);
        }
    }

    @Test
    public void testNcClosedForkJoinAgreesWithMva() throws Exception {
        // Both routes drive the same fixed point and approximate the same
        // system, so they must land within a few percent of each other. The
        // auxiliary classes start at GlobalConstants.FineTol, so this is also
        // the test that the normalizing constant does not degenerate on them:
        // before the struct was recompiled after the fixed point, NC published
        // the first inner solve untouched and the Join queue length came out
        // exactly zero.
        SolverMVA mva = new SolverMVA(buildClosedFJ(3), silent());
        SolverNC nc = new SolverNC(buildClosedFJ(3), silent());
        Matrix Qm = mva.getAvgQLen();
        Matrix Tm = mva.getAvgTput();
        Matrix Qn = nc.getAvgQLen();
        Matrix Tn = nc.getAvgTput();

        assertTrue(relErr(Qn, Qm, 4) < 0.05,
            "NC and MVA queue lengths must agree within 5%, got " + relErr(Qn, Qm, 4));
        assertTrue(Math.abs(Tn.get(0, 0) - Tm.get(0, 0)) / Tm.get(0, 0) < 0.05,
            "NC and MVA throughputs must agree within 5%");
        // the Join must carry a synchronisation delay, not the zero that a
        // fixed point stuck at its first iteration would leave behind
        assertTrue(Qn.get(3, 0) > 0.1, "the Join queue length must be strictly positive, got " + Qn.get(3, 0));
        // flow balance around the fork-join subnetwork
        assertEquals(Tn.get(0, 0), Tn.get(3, 0), 1e-6 * Tn.get(0, 0),
            "the Join throughput must match the reference station throughput");
    }

    @Test
    public void testCtmcSolvesClosedForkJoinExactly() throws Exception {
        // SolverCTMC solves the same model natively on the tag-augmented copy
        // (ModelAdapter.fjtag), which is the exact reference the two fixed-point
        // routes are approximating. Its values are those of the MATLAB mirror to
        // 12 significant digits.
        SolverCTMC ctmc = new SolverCTMC(buildClosedFJ(3), silentCtmc());
        Matrix Qc = ctmc.getAvgQLen();
        Matrix Tc = ctmc.getAvgTput();
        double[] Qgold = {1.418899863876, 1.196465903699, 0.857919086961, 1.107815281588};
        for (int i = 0; i < 4; i++) {
            assertEquals(Qgold[i], Qc.get(i, 0), 1e-9, "CTMC QN row " + i);
            // every station of this single-chain model sits on the same cycle,
            // the Join included: its DEP is the join firing, which fires only
            // from the vanishing marking in which the sibling set is complete
            // and is therefore counted through the rate complement
            assertEquals(Qgold[0], Tc.get(i, 0), 1e-9, "CTMC TN row " + i);
        }
    }

    @Test
    public void testNcClosedForkJoinIsNoWorseThanMvaAgainstCtmc() throws Exception {
        // Mirrors the accuracy check of line-test.git/test/testsFJ/test_fj_driver_nc.m. Station
        // 4 is the Join, whose queue length is the synchronisation delay of the
        // transformed formulation and is not comparable with the exact one.
        Matrix Qc = new SolverCTMC(buildClosedFJ(3), silentCtmc()).getAvgQLen();
        Matrix Tc = new SolverCTMC(buildClosedFJ(3), silentCtmc()).getAvgTput();
        SolverMVA mva = new SolverMVA(buildClosedFJ(3), silent());
        SolverNC nc = new SolverNC(buildClosedFJ(3), silent());
        Matrix Qm = mva.getAvgQLen();
        Matrix Qn = nc.getAvgQLen();
        Matrix Tm = mva.getAvgTput();
        Matrix Tn = nc.getAvgTput();

        double errMVA = relErr(Qm, Qc, 3);
        double errNC = relErr(Qn, Qc, 3);
        assertTrue(errNC < 0.20, "NC queue lengths must be within 20% of the exact solution, got " + errNC);
        assertTrue(errNC <= errMVA + 1e-6,
            "NC must be no less accurate than MVA on queue lengths, got NC " + errNC + " vs MVA " + errMVA);

        double tputErrMVA = Math.abs(Tm.get(0, 0) - Tc.get(0, 0)) / Tc.get(0, 0);
        double tputErrNC = Math.abs(Tn.get(0, 0) - Tc.get(0, 0)) / Tc.get(0, 0);
        assertTrue(tputErrNC < 0.15, "NC throughput must be within 15% of the exact solution, got " + tputErrNC);
        assertTrue(tputErrNC <= tputErrMVA + 1e-6,
            "NC must be no less accurate than MVA on throughput, got NC " + tputErrNC + " vs MVA " + tputErrMVA);
    }

    @Test
    public void testNcForkJoinMetricsAreFinite() throws Exception {
        // A degenerating normalizing constant would surface as Inf/NaN rather
        // than as an error, so assert it explicitly
        SolverNC nc = new SolverNC(buildClosedFJ(5), silent());
        Matrix[] all = {nc.getAvgQLen(), nc.getAvgUtil(), nc.getAvgRespT(), nc.getAvgTput()};
        for (Matrix M : all) {
            for (int i = 0; i < M.getNumRows(); i++) {
                for (int j = 0; j < M.getNumCols(); j++) {
                    double v = M.get(i, j);
                    assertTrue(Double.isFinite(v), "metric must be finite, got " + v);
                }
            }
        }
        Matrix Q = nc.getAvgQLen();
        for (int i = 0; i < Q.getNumRows(); i++) {
            assertTrue(Q.get(i, 0) >= 0, "queue lengths must be non-negative");
        }
    }
}
