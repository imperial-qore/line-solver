/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Exp;
import jline.solvers.mva.SolverMVA;
import jline.solvers.ssa.SolverSSA;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Regression tests for the SIMULATION path of the finite-difference branch of
 * {@link NetworkSolver#getSensitivityTable()}.
 *
 * <p>A simulator has no product-form recursion to differentiate, so it always takes
 * the finite-difference branch. Two behaviours are specific to it and are what these
 * tests pin: common random numbers (the base and the perturbed runs must share a
 * seed, otherwise the difference quotient measures Monte Carlo noise rather than a
 * derivative) and the coarser default step 1e-2, chosen so that the signal exceeds
 * the simulation error, against 1e-4 for the deterministic solvers.</p>
 *
 * <p>The model is a two-station closed network -- a delay with Exp(1) and a PS queue
 * with Exp(2), three jobs in a single class -- small enough that a modest sample
 * count resolves the derivatives, and in scope for SolverMVA, whose analytic branch
 * supplies the reference magnitudes.</p>
 */
public class SimulationSensitivityTest {

    /** Sample count for every simulated sweep here; keeps the class under ~1 min. */
    private static final int SAMPLES = 20000;

    /** Seed used wherever a test needs two sweeps to be comparable. */
    private static final int SEED = 4242;

    /**
     * Closed, one class of 3 jobs: Delay "D" with Exp(1), Queue "Q" under PS with
     * Exp(2), serial routing D -> Q -> D.
     */
    private static Network closedDelayQueue() {
        Network m = new Network("sim_sens");
        Delay d = new Delay(m, "D");
        Queue q = new Queue(m, "Q", SchedStrategy.PS);
        ClosedClass c = new ClosedClass(m, "C1", 3, d, 0);
        d.setService(c, new Exp(1.0));
        q.setService(c, new Exp(2.0));
        RoutingMatrix P = m.initRoutingMatrix();
        P.set(c, c, Network.serialRouting(d, q));
        m.link(P);
        return m;
    }

    private static SolverOptions simOptions(int seed) {
        SolverOptions o = new SolverOptions(SolverType.SSA);
        o.samples = SAMPLES;
        o.seed = seed;
        o.verbose = jline.VerboseLevel.SILENT;
        return o;
    }

    private static SolverSSA sim(Network model, int seed) {
        return new SolverSSA(model, simOptions(seed));
    }

    private static void assertColumnsEqual(NetworkSensitivityTable a, NetworkSensitivityTable b,
                                           double tol, String msg) {
        assertEquals(a.getDTput().size(), b.getDTput().size(), msg + ": row count");
        for (int k = 0; k < a.getDTput().size(); k++) {
            assertEquals(a.getDTput().get(k).doubleValue(), b.getDTput().get(k).doubleValue(), tol,
                    msg + ": dTput row " + k);
            assertEquals(a.getDRespT().get(k).doubleValue(), b.getDRespT().get(k).doubleValue(), tol,
                    msg + ": dRespT row " + k);
            assertEquals(a.getDQLen().get(k).doubleValue(), b.getDQLen().get(k).doubleValue(), tol,
                    msg + ": dQLen row " + k);
            assertEquals(a.getDUtil().get(k).doubleValue(), b.getDUtil().get(k).doubleValue(), tol,
                    msg + ": dUtil row " + k);
        }
    }

    @Test
    @DisplayName("a simulation solver takes the finite-difference branch")
    public void simulatorUsesFiniteDifferences() {
        NetworkSensitivityTable T = sim(closedDelayQueue(), SEED).getSensitivityTable();
        assertEquals("fd", T.getMethod());
        assertEquals(1, T.getDTput().size(), "one queue, one class: one row");
        assertEquals("Q", T.getStationNames().get(0));
        assertEquals("C1", T.getClassNames().get(0));
    }

    /**
     * The load-bearing test: it is what fails if the base and the perturbed runs
     * stop sharing a seed. Two separately constructed solvers carrying the same seed
     * must return bitwise-identical derivative columns.
     */
    @Test
    @DisplayName("common random numbers make the sweep reproducible")
    public void commonRandomNumbersAreReproducible() {
        NetworkSensitivityTable T1 = sim(closedDelayQueue(), SEED).getSensitivityTable();
        NetworkSensitivityTable T2 = sim(closedDelayQueue(), SEED).getSensitivityTable();
        assertColumnsEqual(T1, T2, 1e-12, "same seed, two solver instances");
    }

    /**
     * An unset seed is pinned to 23000 before the sweep. In the JAR a
     * {@link SolverOptions} is born with a random positive draw, so "unset" is
     * expressed as a non-positive seed, which is the guard the FD path tests.
     */
    @Test
    @DisplayName("an unset seed is pinned to 23000 for the sweep")
    public void unsetSeedIsPinned() {
        SolverOptions opt = simOptions(0);
        SolverSSA solver = new SolverSSA(closedDelayQueue(), opt);
        solver.getSensitivityTable();
        assertEquals(23000, solver.getOptions().seed,
                "the FD path must pin an unset seed so that base and perturbed runs pair up");
    }

    /**
     * The simulator default step is 1e-2, asserted behaviourally: at a fixed seed the
     * default-step table reproduces the explicit step 1e-2 exactly, and differs from
     * the table obtained at step 1e-3.
     */
    @Test
    @DisplayName("the simulator default finite-difference step is 1e-2")
    public void simulatorDefaultStepIsCoarse() {
        NetworkSensitivityTable dflt = sim(closedDelayQueue(), SEED)
                .getSensitivityTable("fd", Double.NaN, "forward");
        NetworkSensitivityTable h2 = sim(closedDelayQueue(), SEED)
                .getSensitivityTable("fd", 1e-2, "forward");
        NetworkSensitivityTable h3 = sim(closedDelayQueue(), SEED)
                .getSensitivityTable("fd", 1e-3, "forward");
        assertColumnsEqual(dflt, h2, 1e-12, "default step against explicit 1e-2");
        assertNotEquals(dflt.getDTput().get(0).doubleValue(), h3.getDTput().get(0).doubleValue(),
                "a 1e-3 step must not reproduce the 1e-2 default");
    }

    /**
     * Sign structure: a faster server raises throughput and shortens the queue, so
     * dTput is positive while dRespT, dQLen and dUtil are negative.
     */
    @Test
    @DisplayName("signs follow a faster server: dTput > 0, dRespT/dQLen/dUtil < 0")
    public void derivativeSignsAreConsistent() {
        NetworkSensitivityTable T = sim(closedDelayQueue(), SEED).getSensitivityTable();
        assertTrue(T.getDTput().get(0).doubleValue() > 0,
                "a faster server must raise throughput, got " + T.getDTput().get(0));
        assertTrue(T.getDRespT().get(0).doubleValue() < 0,
                "a faster server must shorten the response time, got " + T.getDRespT().get(0));
        assertTrue(T.getDQLen().get(0).doubleValue() < 0,
                "a faster server must shorten the queue, got " + T.getDQLen().get(0));
        assertTrue(T.getDUtil().get(0).doubleValue() < 0,
                "a faster server must lower the utilization, got " + T.getDUtil().get(0));
    }

    /**
     * Magnitude sanity against the analytic branch of SolverMVA on the same model.
     * A factor of two is the honest band at this sample count -- the Monte Carlo
     * error on a difference quotient is tens of percent -- and it still catches a
     * missing visit factor or an unpaired seed, both of which move the estimate by
     * an order of magnitude or destroy it outright.
     */
    @Test
    @DisplayName("simulated derivatives are within a factor of two of the analytic ones")
    public void magnitudesAgreeWithTheAnalyticBranch() {
        Network model = closedDelayQueue();
        NetworkSensitivityTable exact = new SolverMVA(model).getSensitivityTable("exact", Double.NaN, "forward");
        NetworkSensitivityTable T = sim(closedDelayQueue(), SEED).getSensitivityTable();
        assertEquals("exact", exact.getMethod());
        assertWithinFactorTwo(exact.getDTput(), T.getDTput(), "dTput");
        assertWithinFactorTwo(exact.getDRespT(), T.getDRespT(), "dRespT");
        assertWithinFactorTwo(exact.getDQLen(), T.getDQLen(), "dQLen");
        assertWithinFactorTwo(exact.getDUtil(), T.getDUtil(), "dUtil");
    }

    private static void assertWithinFactorTwo(List<Double> reference, List<Double> simulated,
                                              String column) {
        assertEquals(reference.size(), simulated.size(), column + ": row count");
        for (int k = 0; k < reference.size(); k++) {
            double ref = reference.get(k).doubleValue();
            double sim = simulated.get(k).doubleValue();
            double ratio = sim / ref;
            assertTrue(ratio > 0.5 && ratio < 2.0,
                    column + " row " + k + ": simulated " + sim + " is not within a factor of two of "
                            + "the analytic " + ref);
        }
    }
}
