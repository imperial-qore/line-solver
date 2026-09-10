/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.nc;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.processes.Exp;
import jline.lang.nodes.Queue;
import jline.solvers.SolverResult;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Regression test for the specialized LCFS + LCFS-PR product-form solver
 * (solver_nc_lcfsqn): class multiplicities N_r > 1 must be handled by
 * expanding classes into exchangeable single-job copies. Expected values
 * are the exact CTMC solutions (cross-validated in MATLAB and Python).
 */
public class SolverNCLcfsqnTest {

    private static SolverResult solveCycle(int[] N) throws Exception {
        double[] L1 = {0.20, 0.45};
        double[] L2 = {0.30, 0.12};
        Network model = new Network("lcfsqn_cycle");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.LCFS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.LCFSPR);
        ClosedClass[] cl = new ClosedClass[2];
        for (int r = 0; r < 2; r++) {
            cl[r] = new ClosedClass(model, "C" + (r + 1), N[r], q1);
            q1.setService(cl[r], new Exp(1 / L1[r]));
            q2.setService(cl[r], new Exp(1 / L2[r]));
        }
        RoutingMatrix P = model.initRoutingMatrix();
        for (int r = 0; r < 2; r++) {
            P.set(cl[r], Network.serialRouting(q1, q2));
        }
        model.link(P);
        SolverNC nc = new SolverNC(model);
        return nc.getAvg();
    }

    @Test
    public void testSingleJobPerClass() throws Exception {
        // N = [1 1]: original one-job-per-class regime of the permanents
        SolverResult res = solveCycle(new int[]{1, 1});
        assertEquals(0.578947368421053, res.QN.get(0, 0), 1e-8);
        assertEquals(0.714285714285714, res.QN.get(0, 1), 1e-8);
        assertEquals(1.203007518796992, res.TN.get(0, 0), 1e-8);
        assertEquals(1.253132832080200, res.TN.get(0, 1), 1e-8);
    }

    @Test
    public void testClassMultiplicities() throws Exception {
        // N = [2 2]: silently returned all-zero metrics before the fix
        SolverResult res = solveCycle(new int[]{2, 2});
        assertEquals(1.457574, res.QN.get(0, 0), 1e-5);
        assertEquals(1.482294, res.QN.get(0, 1), 1e-5);
        assertEquals(1.113855, res.TN.get(0, 0), 1e-5);
        assertEquals(1.556857, res.TN.get(0, 1), 1e-5);
        // conservation at the LCFS-PR station
        assertEquals(2.0 - 1.457574, res.QN.get(1, 0), 1e-5);
        assertEquals(2.0 - 1.482294, res.QN.get(1, 1), 1e-5);
    }
}
