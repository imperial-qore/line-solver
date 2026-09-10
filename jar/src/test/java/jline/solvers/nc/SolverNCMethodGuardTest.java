/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.nc;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.api.pfqn.nc.Pfqn_nc;
import jline.io.Ret;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Exp;
import jline.solvers.SolverOptions;
import jline.solvers.mva.SolverMVA;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.util.Arrays;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Guards against the two silent-failure paths of the normalizing-constant solver:
 *
 * <ul>
 *   <li>method='comom' on a model with more than one queueing station. comom solves the
 *       repairman model, so it needs M=1. This case used to assign no normalizing
 *       constant at all, leaving lG null and surfacing downstream as an opaque
 *       NullPointerException (in MATLAB, as all-zero queue lengths reported as a
 *       completed analysis).</li>
 *   <li>an unadvertised method name (e.g. 'mom'), which used to fall through to the
 *       default analyzer and return ALL-ZERO queue lengths while still reporting a
 *       completed analysis.</li>
 * </ul>
 *
 * <p>Both assert the message, not merely that something was raised: the point of the
 * change is that the refusal is actionable, and a path that merely throws something is
 * indistinguishable from the NullPointerException it replaced.</p>
 */
public class SolverNCMethodGuardTest {

    private static final double TOL = 1e-6;

    @BeforeAll
    public static void setUp() {
        GlobalConstants.setVerbose(VerboseLevel.SILENT);
    }

    /** Delay + ONE queue, 2 classes: M=1, R=2, the repairman shape comom supports. */
    private static Network repairmanM1() {
        Network model = new Network("comom_M1");
        Delay d = new Delay(model, "Delay");
        Queue q = new Queue(model, "Queue1", SchedStrategy.PS);
        ClosedClass c1 = new ClosedClass(model, "Class1", 3, d, 0);
        ClosedClass c2 = new ClosedClass(model, "Class2", 2, d, 0);
        d.setService(c1, Exp.fitMean(1.0));
        d.setService(c2, Exp.fitMean(2.0));
        q.setService(c1, Exp.fitMean(1.5));
        q.setService(c2, Exp.fitMean(0.8));
        RoutingMatrix rm = model.initRoutingMatrix();
        rm.set(c1, c1, d, q, 1.0);
        rm.set(c1, c1, q, d, 1.0);
        rm.set(c2, c2, d, q, 1.0);
        rm.set(c2, c2, q, d, 1.0);
        model.link(rm);
        return model;
    }

    /** Delay + TWO queues, 2 classes: M=2, R=2, which comom cannot solve. */
    private static Network tandemM2() {
        Network model = new Network("comom_M2");
        Delay d = new Delay(model, "Delay");
        Queue q1 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue q2 = new Queue(model, "Queue2", SchedStrategy.PS);
        ClosedClass c1 = new ClosedClass(model, "Class1", 3, d, 0);
        ClosedClass c2 = new ClosedClass(model, "Class2", 2, d, 0);
        d.setService(c1, Exp.fitMean(1.0));
        d.setService(c2, Exp.fitMean(2.0));
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

    private static SolverNC nc(Network model, String method) {
        SolverOptions o = new SolverOptions();
        o.method = method;
        return new SolverNC(model, o);
    }

    /**
     * comom with M>1 must refuse, and the refusal must name the offending station count
     * and the way out. It must NOT return a null/zeroed normalizing constant.
     */
    @Test
    public void testComomRefusesMultipleQueueingStations() {
        RuntimeException e = assertThrows(RuntimeException.class,
                () -> nc(tandemM2(), "comom").getAvgQLen(),
                "comom on a 2-queue model must refuse rather than return a null/zero result");

        String msg = e.getMessage();
        assertNotNull(msg, "the refusal must carry a message");
        assertTrue(msg.contains("The 'comom' method supports a single queueing station, but this model has 2."),
                "message must name the actual station count, was: " + msg);
        assertTrue(msg.contains("Use 'default' or 'ca' for an exact normalizing constant"),
                "message must offer an actionable alternative, was: " + msg);
        // The old failure mode: an opaque NPE on the null lG.
        assertFalse(msg.contains("NullPointerException") || msg.contains("because \"lG\" is null"),
                "must not surface as a null-lG NullPointerException, was: " + msg);
    }

    /** The same guard at the API level, where the message is unwrapped. */
    @Test
    public void testPfqnNcComomRefusesMultipleQueueingStations() {
        // Two queueing stations (M=2), two classes (R=2), nonzero think time.
        Matrix L = new Matrix(new double[][]{{1.5, 0.8}, {0.9, 1.1}});
        Matrix N = new Matrix(new double[][]{{3.0, 2.0}});
        Matrix Z = new Matrix(new double[][]{{1.0, 2.0}});
        SolverOptions o = new SolverOptions();
        o.method = "comom";

        RuntimeException e = assertThrows(RuntimeException.class,
                () -> Pfqn_nc.compute_norm_const(L, N, Z, o));
        assertTrue(e.getMessage().contains(
                        "The 'comom' method supports a single queueing station, but this model has 2."),
                "was: " + e.getMessage());
    }

    /**
     * The M=1 repairman path is what SolverNC selects internally for 2-station
     * Delay+queue LN submodels, so it must keep working and stay exact.
     */
    @Test
    public void testComomSingleQueueingStationStillMatchesExact() {
        Matrix qComom = nc(repairmanM1(), "comom").getAvgQLen();
        Matrix qExact = nc(repairmanM1(), "exact").getAvgQLen();

        SolverOptions mvaOpt = new SolverOptions();
        mvaOpt.method = "exact";
        Matrix qMva = new SolverMVA(repairmanM1(), mvaOpt).getAvgQLen();

        assertEquals(2, qComom.getNumRows());
        assertEquals(2, qComom.getNumCols());

        double total = 0.0;
        for (int i = 0; i < qComom.getNumRows(); i++) {
            for (int r = 0; r < qComom.getNumCols(); r++) {
                assertEquals(qExact.get(i, r), qComom.get(i, r), TOL,
                        "comom must match NC 'exact' at station " + i + ", class " + r);
                assertEquals(qMva.get(i, r), qComom.get(i, r), TOL,
                        "comom must match exact MVA at station " + i + ", class " + r);
                total += qComom.get(i, r);
            }
        }
        // Guards against the all-zero regression: a zeroed Q would agree with nothing,
        // but pin the population explicitly so a doubly-zeroed comparison cannot pass.
        assertEquals(5.0, total, 1e-6, "total queue length must conserve the closed population");
    }

    /**
     * 'mom' names nothing in any codebase: the method-of-moments solver and the
     * method name were both removed, so no solver implements or advertises it. It
     * must be rejected rather than silently produce all-zero queue lengths via the
     * default analyzer.
     */
    @Test
    public void testUnadvertisedMomMethodIsRejected() {
        String[] valid = new SolverNC(repairmanM1()).listValidMethods();
        assertFalse(Arrays.asList(valid).contains("mom"),
                "'mom' must not be advertised by SolverNC (no codebase advertises it)");
        assertTrue(Arrays.asList(valid).contains("comom"), "'comom' is advertised");

        RuntimeException e = assertThrows(RuntimeException.class,
                () -> nc(repairmanM1(), "mom").getAvgQLen(),
                "an unadvertised method must be refused, not silently zeroed");
        assertTrue(e.getMessage().contains("The 'mom' method is unsupported by this solver."),
                "was: " + e.getMessage());
    }

    /** Every advertised method must survive its own validation gate. */
    @Test
    public void testAdvertisedComomIsNotRejectedByTheMethodGate() {
        // Regression guard: the unknown-method rejection must not catch a valid method.
        assertDoesNotThrow(() -> nc(repairmanM1(), "comom").getAvgQLen());
        assertDoesNotThrow(() -> nc(repairmanM1(), "default").getAvgQLen());
    }

    /**
     * 'comomld' used to be whitelisted at the gate without being advertised, which kept
     * it out of every consumer that enumerates listValidMethods -- the SanityQN harness
     * among them, so the exact load-dependent normalizing constant had no golden.
     */
    @Test
    public void testComomldIsAdvertisedAndExactOnAProductFormModel() {
        assertTrue(Arrays.asList(new SolverNC(repairmanM1()).listValidMethods()).contains("comomld"),
                "'comomld' is a user-selectable method and must be advertised");
        Matrix ld = nc(repairmanM1(), "comomld").getAvgQLen();
        Matrix exact = nc(repairmanM1(), "exact").getAvgQLen();
        assertEquals(exact.getNumRows() * exact.getNumCols(), ld.getNumRows() * ld.getNumCols(),
                "comomld must return a full queue-length table");
        for (int i = 0; i < exact.getNumRows(); i++) {
            for (int r = 0; r < exact.getNumCols(); r++) {
                assertEquals(exact.get(i, r), ld.get(i, r), 1e-9,
                        "comomld is exact on a product-form model, station " + i + " class " + r);
            }
        }
    }
}
