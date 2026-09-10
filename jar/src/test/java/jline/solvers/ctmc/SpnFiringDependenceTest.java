package jline.solvers.ctmc;

import jline.lang.ClosedClass;
import jline.lang.Mode;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.nodes.Place;
import jline.lang.nodes.Transition;
import jline.lang.processes.Exp;
import jline.solvers.ssa.SolverSSA;
import jline.util.matrix.Matrix;
import jline.GlobalConstants;
import jline.VerboseLevel;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Marking-dependent transition firing rates (Transition.setFiringRateDependence).
 *
 * <p>A closed single-class SPN with N tokens on two places and mass-action rates on
 * both transitions (rate = k * tokens-at-input-place, single server) makes each token
 * an independent two-state CTMC, so the stationary marking of P1 is Binomial(N,
 * k2/(k1+k2)) and E[n1] = N*k2/(k1+k2). SolverCTMC evaluates the g(marking) multiplier
 * per enumerated state and must match that closed form exactly.
 *
 * <p>Neither SSA engine applies the multiplier, so SolverSSA must refuse the model
 * rather than silently simulate the nominal (unscaled) rate. Twin of the MATLAB guard
 * in solver_ssa_nrm.m and of python/tests/test_spn_firing_dependence.py.
 */
public class SpnFiringDependenceTest {

    private static final int N = 5;
    private static final double K1 = 2, K2 = 3;

    /** Mass-action loop P1 -(T1)-> P2 -(T2)-> P1 with N tokens initially on P1. */
    private static Network massActionLoop() {
        Network model = new Network("firingdep");
        Place p1 = new Place(model, "P1");
        Place p2 = new Place(model, "P2");
        Transition t1 = new Transition(model, "T1");
        Transition t2 = new Transition(model, "T2");
        ClosedClass jobClass = new ClosedClass(model, "C", N, p1, 0);
        final int i1 = model.getNodeIndex(p1);
        final int i2 = model.getNodeIndex(p2);

        Mode m1 = t1.addMode("m1");
        t1.setDistribution(m1, new Exp(K1));
        t1.setEnablingConditions(m1, jobClass, p1, 1);
        t1.setFiringOutcome(m1, jobClass, p2, 1);
        t1.setFiringRateDependence(m1, (Matrix marking) -> marking.get(i1, 0));

        Mode m2 = t2.addMode("m2");
        t2.setDistribution(m2, new Exp(K2));
        t2.setEnablingConditions(m2, jobClass, p2, 1);
        t2.setFiringOutcome(m2, jobClass, p1, 1);
        t2.setFiringRateDependence(m2, (Matrix marking) -> marking.get(i2, 0));

        RoutingMatrix routing = model.initRoutingMatrix();
        routing.set(jobClass, jobClass, p1, t1, 1.0);
        routing.set(jobClass, jobClass, p2, t2, 1.0);
        routing.set(jobClass, jobClass, t1, p2, 1.0);
        routing.set(jobClass, jobClass, t2, p1, 1.0);
        model.link(routing);

        Matrix s1 = new Matrix(1, 1);
        s1.set(0, 0, N);
        p1.setState(s1);
        Matrix s2 = new Matrix(1, 1);
        s2.set(0, 0, 0);
        p2.setState(s2);
        return model;
    }

    @Test
    public void testCtmcMatchesBinomialClosedForm() {
        Matrix qlen = new SolverCTMC(massActionLoop()).getAvgQLen();
        assertEquals(N * K2 / (K1 + K2), qlen.get(0, 0), 1e-6);
        assertEquals(N * K1 / (K1 + K2), qlen.get(1, 0), 1e-6);
    }

    @Test
    public void testSsaRejectsDependence() {
        // the refusal is the expected outcome here, so its SEVERE report is noise
        VerboseLevel savedVerbose = GlobalConstants.getVerbose();
        GlobalConstants.Verbose = VerboseLevel.SILENT;
        RuntimeException e;
        try {
            e = assertThrows(RuntimeException.class,
                    () -> new SolverSSA(massActionLoop()).getAvgTable());
        } finally {
            GlobalConstants.Verbose = savedVerbose;
        }
        String trace = "";
        for (Throwable t = e; t != null; t = t.getCause()) {
            trace = trace + t.getMessage();
        }
        assertTrue(trace.contains("marking-dependent firing rate"),
                "SolverSSA must refuse a marking-dependent firing rate, got: " + trace);
    }
}
