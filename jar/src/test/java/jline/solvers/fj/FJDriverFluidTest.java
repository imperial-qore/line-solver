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
import jline.lang.constant.SolverType;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Fork;
import jline.lang.nodes.Join;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.solvers.SolverOptions;
import jline.solvers.fluid.SolverFluid;
import jline.solvers.mva.SolverMVA;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Regression tests for the SolverFluid fork-join route.
 *
 * <p>SolverFluid used to throw on any Fork or Join. It now drives the same
 * solver-agnostic fixed point MVA and NC drive ({@link FJFixedPoint}), with the
 * fluid analyzer as the inner solve: the MMT transformation emits only Source,
 * Delay, Queue, Router and ClassSwitch, every one of which the fluid drift
 * already carries, so no fork-join code was added to the fluid solver itself.
 *
 * <p>What is asserted is what the fixed point is responsible for rather than a
 * golden: flow balance across the fork (every station on the cycle carries the
 * same class throughput), a synchronisation delay that is actually charged at
 * the join, and agreement with SolverMVA to within the gap the two
 * approximations are expected to leave. A regression to the state before the
 * wiring shows up as an exception, and a regression in the auxiliary-class
 * refresh shows up as a join whose response time is zero.
 *
 * <p>Java mirror of the MATLAB @SolverFLD/fldDispatch.m route.
 */
public class FJDriverFluidTest {

    private static SolverOptions silentFluid() {
        SolverOptions options = new SolverOptions(SolverType.FLUID);
        options.verbose = VerboseLevel.SILENT;
        return options;
    }

    private static SolverOptions silentMva() {
        SolverOptions options = new SolverOptions(SolverType.MVA);
        options.verbose = VerboseLevel.SILENT;
        return options;
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
        q2.setService(jobclass, new Exp(3.0));
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

    private static Network buildOpenFJ() {
        Network model = new Network("fjopen");
        Source source = new Source(model, "Source");
        Fork fork = new Fork(model, "Fork");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.PS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.PS);
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

    /** SolverFluid used to refuse every fork-join model outright. */
    @Test
    public void fluidNoLongerRefusesForkJoin() {
        Network model = buildClosedFJ(6);
        SolverFluid solver = new SolverFluid(model, silentFluid());
        assertTrue(solver.supports(model),
                "SolverFluid must declare Fork/Join now that it drives the MMT fixed point");
        Matrix QN = solver.getAvgQLen();
        assertNotNull(QN);
        assertEquals(4, QN.getNumRows());
        for (int i = 0; i < QN.getNumRows(); i++) {
            assertTrue(QN.get(i, 0) >= 0 && Double.isFinite(QN.get(i, 0)),
                    "queue length at station " + i + " must be finite and non-negative");
        }
    }

    /**
     * Flow balance across the fork. Every station on the cycle carries the same
     * class throughput once the auxiliary open classes are merged back; a driver
     * that failed to merge them reports the delay at twice the queues' rate.
     */
    @Test
    public void closedForkJoinConservesFlow() {
        Network model = buildClosedFJ(6);
        SolverFluid solver = new SolverFluid(model, silentFluid());
        Matrix TN = solver.getAvgTput();
        double delay = TN.get(0, 0);
        assertTrue(delay > 0, "the reference station must carry throughput");
        for (int i = 1; i < TN.getNumRows(); i++) {
            assertEquals(delay, TN.get(i, 0), 1e-6 * delay,
                    "station " + i + " must carry the reference station's throughput");
        }
    }

    /**
     * The synchronisation delay must actually be charged. The auxiliary source
     * rate is written back with setRate, which leaves the SCV at 1, so a refresh
     * that only touches sn.rates leaves the fluid drift integrating the initial
     * GlobalConstants.FineTol rate and the join response time collapses to zero.
     */
    @Test
    public void joinChargesASynchronisationDelay() {
        Network model = buildClosedFJ(6);
        SolverFluid solver = new SolverFluid(model, silentFluid());
        Matrix RN = solver.getAvgRespT();
        double joinRespT = RN.get(3, 0);
        assertTrue(joinRespT > 1e-3,
                "the join must charge a synchronisation delay, got " + joinRespT);
        Matrix QN = solver.getAvgQLen();
        assertTrue(QN.get(3, 0) > 1e-3,
                "the join must hold the jobs waiting on their siblings, got " + QN.get(3, 0));
    }

    /**
     * The two approximations solve the same transformed model through the same
     * fixed point, so they may differ only by the gap between mean-value
     * analysis and the fluid limit, not by a structural error.
     */
    @Test
    public void fluidAgreesWithMvaOnTheClosedModel() {
        Network model = buildClosedFJ(6);
        Matrix TNfluid = new SolverFluid(model, silentFluid()).getAvgTput();
        Matrix TNmva = new SolverMVA(buildClosedFJ(6), silentMva()).getAvgTput();
        double relErr = Math.abs(TNfluid.get(0, 0) - TNmva.get(0, 0)) / TNmva.get(0, 0);
        assertTrue(relErr < 0.15,
                "fluid and MVA throughput must agree to 15%, got " + relErr);
    }

    /** The open route exercises the transformed model's own Source alongside the auxiliary one. */
    @Test
    public void openForkJoinSolves() {
        Network model = buildOpenFJ();
        SolverFluid solver = new SolverFluid(model, silentFluid());
        Matrix TN = solver.getAvgTput();
        for (int i = 0; i < TN.getNumRows(); i++) {
            assertTrue(Double.isFinite(TN.get(i, 0)),
                    "throughput at station " + i + " must be finite");
        }
    }
}
