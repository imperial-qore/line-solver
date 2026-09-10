/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.ctmc;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.processes.Exp;
import jline.solvers.ba.SolverBA;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * The QRF reduction bounds moved out of SolverCTMC into SolverBA, and the
 * refusal has to carry that forwarding address.
 * <p>
 * Dropping the names from {@code listValidMethods} is what makes SolverCTMC
 * refuse them, and it is also what loses the address, so the two have to be
 * declared together. They were not: {@code NetworkSolver.checkDeclaredMethod}
 * sits above every dispatcher and reported the flat "the 'qrf.bas' method is
 * unsupported by this solver", while the redirect SolverCTMC did carry sat
 * downstream in {@code runAnalyzer} and never ran. The gate asks
 * {@code unsupportedMethodReason} first now.
 */
public class SolverCTMCQrfRedirectTest {

    /** A two-queue closed cycle: no delay station, which the QRF formulation refuses. */
    private static Network cqn() {
        Network model = new Network("cqn");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.FCFS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.FCFS);
        ClosedClass c = new ClosedClass(model, "C", 2, q1, 0);
        q1.setService(c, new Exp(1.0));
        q2.setService(c, new Exp(2.0));
        model.link(model.serialRouting(q1, q2));
        return model;
    }

    @Test
    public void testQrfMethodsRedirectToSolverBA() {
        for (String moved : new String[]{"qrf", "qrf.bas", "qrf.mmi", "qrf.rsrd"}) {
            RuntimeException e = assertThrows(RuntimeException.class,
                    () -> new SolverCTMC(cqn(), moved).getAvgTable());
            String msg = String.valueOf(e.getMessage());
            assertTrue(msg.contains("moved out of SolverCTMC into the dedicated SolverBA solver"),
                    "'" + moved + "' must name where it went, got: " + msg);
            assertTrue(msg.contains("SolverBA(model, \"" + moved + "\")"),
                    "'" + moved + "' must show the call that replaces it, got: " + msg);
        }
    }

    /** A name that never existed still gets the flat refusal, not a redirect. */
    @Test
    public void testUnknownMethodKeepsTheFlatRefusal() {
        RuntimeException e = assertThrows(RuntimeException.class,
                () -> new SolverCTMC(cqn(), "nosuchmethod").getAvgTable());
        String msg = String.valueOf(e.getMessage());
        assertTrue(msg.contains("unsupported by this solver"), msg);
        assertFalse(msg.contains("moved out of SolverCTMC"), msg);
    }

    /** The address is worth giving because it answers: SolverBA does serve them. */
    @Test
    public void testSolverBAServesTheMethodTheRedirectNames() {
        assertDoesNotThrow(() -> new SolverBA(cqn(), "qrf.bas").getAvgTable());
    }
}
