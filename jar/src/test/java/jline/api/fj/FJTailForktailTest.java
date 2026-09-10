/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.fj;

import jline.VerboseLevel;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Fork;
import jline.lang.nodes.Join;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.lang.processes.HyperExp;
import jline.solvers.Solver;
import jline.solvers.SolverOptions;
import jline.solvers.mva.SolverMVA;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Tests for the ForkTail tail-latency approximation and for the
 * NetworkSolver.getPerctRespT(percentiles, "forktail") entry point.
 *
 * <p>The golden values come from the MATLAB implementation
 * (matlab/src/api/fj/fj_tail_forktail.m), which native python reproduces to
 * 1e-6 on the same model.
 */
public class FJTailForktailTest {

    private static final double TOL = 1e-6;

    @Test
    public void testExponentialBranchIsExactFit() {
        // SCV = 1 gives alpha = 1, beta = E[T], so one branch reproduces the
        // exponential quantile exactly
        double ET = 2.5;
        FJ_tail_forktail.GEFit fit = FJ_tail_forktail.geFit(ET, ET * ET);
        assertEquals(1.0, fit.alpha, 1e-9);
        assertEquals(ET, fit.beta, 1e-9);
        assertEquals(-ET * Math.log(1 - 0.99), FJ_tail_forktail.fj_tail_forktail(ET, ET * ET, 1, 99), 1e-9);
    }

    @Test
    public void testHomogeneousClosedForm() {
        double ET = 1.3;
        int K = 5;
        double p = 0.95;
        assertEquals(-ET * Math.log(1 - Math.pow(p, 1.0 / K)),
            FJ_tail_forktail.fj_tail_forktail(ET, ET * ET, K, p), 1e-9);
    }

    @Test
    public void testVectorPathMatchesScalarPath() {
        double ET = 0.8;
        double VT = 1.7;
        int K = 4;
        double scalar = FJ_tail_forktail.fj_tail_forktail(ET, VT, K, 99);
        double[] ETv = new double[K];
        double[] VTv = new double[K];
        for (int i = 0; i < K; i++) {
            ETv[i] = ET;
            VTv[i] = VT;
        }
        assertEquals(scalar, FJ_tail_forktail.fj_tail_forktail(ETv, VTv, 99), 1e-6);
    }

    @Test
    public void testMg1MomentsMatchMM1() {
        // M/M/1: E[T] = 1/(mu-lambda) and the response time is exponential
        double mu = 1.0;
        double lambda = 0.7;
        FJ_tail_forktail.ResptMoments m =
            FJ_tail_forktail.fj_mg1_respt_moments(lambda, 1 / mu, 2 / (mu * mu), 6 / (mu * mu * mu));
        assertEquals(1 / (mu - lambda), m.mean, 1e-9);
        assertEquals(m.mean * m.mean, m.variance, 1e-9);
    }

    @Test
    public void testMg1RejectsUnstableBranch() {
        assertThrows(IllegalArgumentException.class,
            () -> FJ_tail_forktail.fj_mg1_respt_moments(1.2, 1.0, 2.0, 6.0));
    }

    @Test
    public void testRandomFanoutMixture() {
        double ET = 1.4;
        double VT = 2.6;
        double p = 99;
        // a degenerate mixture reproduces the fixed-fanout answer
        assertEquals(FJ_tail_forktail.fj_tail_forktail(ET, VT, 3, p),
            FJ_tail_forktail.fj_tail_forktail(ET, VT, new int[]{3, 5}, new double[]{1.0, 0.0}, p), 1e-6);
        // a genuine mixture sits between the extreme fanouts
        double x3 = FJ_tail_forktail.fj_tail_forktail(ET, VT, 3, p);
        double x9 = FJ_tail_forktail.fj_tail_forktail(ET, VT, 9, p);
        double xm = FJ_tail_forktail.fj_tail_forktail(ET, VT, new int[]{3, 9}, new double[]{0.5, 0.5}, p);
        assertTrue(xm > x3 && xm < x9, "mixture percentile must sit between the extreme fanouts");
    }

    @Test
    public void testGetPerctRespTForktailMatchesMatlab() {
        Network model = new Network("fj3");
        Source source = new Source(model, "Source");
        Fork fork = new Fork(model, "Fork");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.FCFS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.FCFS);
        Queue q3 = new Queue(model, "Q3", SchedStrategy.FCFS);
        Join join = new Join(model, "Join", fork);
        Sink sink = new Sink(model, "Sink");
        OpenClass jobclass = new OpenClass(model, "C1", 0);

        source.setArrival(jobclass, Exp.fitMean(1.0 / 0.8));
        q1.setService(jobclass, Exp.fitMean(1.0));
        q2.setService(jobclass, Erlang.fitMeanAndOrder(1.0, 2));
        q3.setService(jobclass, HyperExp.fitMeanAndSCV(1.0, 4.0));

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobclass, jobclass, source, fork, 1.0);
        P.set(jobclass, jobclass, fork, q1, 1.0);
        P.set(jobclass, jobclass, fork, q2, 1.0);
        P.set(jobclass, jobclass, fork, q3, 1.0);
        P.set(jobclass, jobclass, q1, join, 1.0);
        P.set(jobclass, jobclass, q2, join, 1.0);
        P.set(jobclass, jobclass, q3, join, 1.0);
        P.set(jobclass, jobclass, join, sink, 1.0);
        model.link(P);

        // SolverMVA.defaultOptions(), not Solver.defaultOptions(): the generic
        // defaults carry iter_tol = 1e-4, which stops the fork-join fixed point
        // with a residual of 4e-6 on the branch throughput. ForkTail reads that
        // throughput as the branch arrival rate and the M/G/1 moments amplify it
        // by rho/(1-rho)^2, so the percentile moves by ~1e-3, far more than the
        // golden band. MATLAB and python use the MVA defaults (iter_tol = 1e-6).
        SolverOptions options = SolverMVA.defaultOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverMVA solver = new SolverMVA(model, options);
        Matrix perct = solver.getPerctRespT(new double[]{95, 99}, "forktail");

        // MATLAB and native python both give 45.838545 and 80.521692
        assertEquals(45.838545, perct.get(0, 0), 1e-4);
        assertEquals(80.521692, perct.get(0, 1), 1e-4);
    }

    @Test
    public void testGetPerctRespTRejectsModelWithoutFork() {
        Network model = new Network("noFork");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Q1", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass jobclass = new OpenClass(model, "C1", 0);
        source.setArrival(jobclass, Exp.fitMean(2.0));
        queue.setService(jobclass, Exp.fitMean(1.0));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobclass, jobclass, model.serialRouting(source, queue, sink));
        model.link(P);

        SolverOptions options = Solver.defaultOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverMVA solver = new SolverMVA(model, options);
        assertThrows(RuntimeException.class, () -> solver.getPerctRespT(new double[]{99}, "forktail"));
    }
}
