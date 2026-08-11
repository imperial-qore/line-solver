/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.ln;

import jline.VerboseLevel;
import jline.lang.constant.SchedStrategy;
import jline.lang.layered.Activity;
import jline.lang.layered.Entry;
import jline.lang.layered.LayeredNetwork;
import jline.lang.layered.Processor;
import jline.lang.layered.Task;
import jline.lang.processes.Exp;
import jline.solvers.LayeredNetworkSensitivityTable;
import jline.solvers.mva.SolverMVA;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Regression tests for {@link SolverLN#getSensitivityTable()}, the layer-wise view
 * over the rate sensitivities of a layered network.
 *
 * <p>The goldens are MATLAB's, from {@code @SolverLN/getSensitivityTable.m} on the
 * two-tier LQN built here, and are reproduced by native Python as well, so these
 * tests pin cross-codebase parity. Agreement is asserted to three significant
 * digits, the LN fixed point differing in its last digits across codebases.</p>
 *
 * <p>Every entry is a derivative WITHIN ITS LAYER, the other layers' fixed-point
 * parameters being held constant; it is a partial derivative of the layer submodel
 * and not the total derivative of the layered model.</p>
 */
public class LnSensitivityTableTest {

    /** Three significant digits of agreement against the MATLAB goldens. */
    private static final double REL_TOL = 1e-3;

    /**
     * Two-tier LQN: reference task T1 (multiplicity 5, think Exp(1/2)) on processor
     * P1 runs AS1 = Exp(10), which makes one synchronous call to entry E2 on task T2
     * (multiplicity 1, FCFS, think Exp(1/3)) on processor P2, where AS2 = Exp(20)
     * replies. Both processors are single-server PS.
     */
    private static LayeredNetwork twoTierLqn() {
        LayeredNetwork model = new LayeredNetwork("lqn_sens");
        Processor P1 = new Processor(model, "P1", 1, SchedStrategy.PS);
        Processor P2 = new Processor(model, "P2", 1, SchedStrategy.PS);
        Task T1 = new Task(model, "T1", 5, SchedStrategy.REF).on(P1);
        T1.setThinkTime(new Exp(1.0 / 2.0));
        Task T2 = new Task(model, "T2", 1, SchedStrategy.FCFS).on(P2);
        T2.setThinkTime(new Exp(1.0 / 3.0));
        Entry E1 = new Entry(model, "E1").on(T1);
        Entry E2 = new Entry(model, "E2").on(T2);
        new Activity(model, "AS1", new Exp(10)).on(T1).boundTo(E1).synchCall(E2, 1);
        new Activity(model, "AS2", new Exp(20)).on(T2).boundTo(E2).repliesTo(E2);
        return model;
    }

    private static SolverLN solver() {
        SolverLN s = new SolverLN(twoTierLqn(), new jline.solvers.ln.SolverFactory() {
            public jline.solvers.NetworkSolver at(jline.lang.Network m) {
                return new SolverMVA(m);
            }
        });
        s.options.verbose = VerboseLevel.SILENT;
        return s;
    }

    private static int rowOfLayer(LayeredNetworkSensitivityTable T, String layer) {
        for (int i = 0; i < T.getLayerNames().size(); i++) {
            if (T.getLayerNames().get(i).contains(layer)) {
                return i;
            }
        }
        throw new AssertionError("no row for layer " + layer + " in " + T.getLayerNames());
    }

    private static void assertRelEquals(double expected, double actual, String msg) {
        double denom = Math.max(Math.abs(expected), 1e-300);
        assertTrue(Math.abs(actual - expected) <= REL_TOL * denom,
                msg + ": expected " + expected + " but got " + actual
                        + " (rel err " + (Math.abs(actual - expected) / denom) + ")");
    }

    private static void assertGoldenLayer(LayeredNetworkSensitivityTable T, String layer,
                                          double dTput, double dRespT, double dQLen, double dUtil) {
        int k = rowOfLayer(T, layer);
        assertRelEquals(dTput, T.getDTput().get(k).doubleValue(), "dTput_dRate " + layer);
        assertRelEquals(dRespT, T.getDRespT().get(k).doubleValue(), "dRespT_dRate " + layer);
        assertRelEquals(dQLen, T.getDQLen().get(k).doubleValue(), "dQLen_dRate " + layer);
        assertRelEquals(dUtil, T.getDUtil().get(k).doubleValue(), "dUtil_dRate " + layer);
    }

    @Test
    @DisplayName("the layer table matches the MATLAB goldens, one row per layer")
    public void testLayerTableMatchesMatlab() {
        LayeredNetworkSensitivityTable T = solver().getSensitivityTable();
        assertEquals(3, T.getLayerNames().size(), "one row per layer: P1, P2 and T2");
        // A layer submodel is chain-based, which the analytic branch handles by
        // aggregating classes into chains, so every layer differentiates analytically.
        assertEquals("exact", T.getMethod(), "every layer uses the analytic branch");
        assertEquals(Arrays.asList("exact", "exact", "exact"), T.getLayerMethods(),
                "the per-layer branch labels must be reported individually");
        assertGoldenLayer(T, "P1", 0.015211, -0.014406, -0.031257, -0.021453);
        assertGoldenLayer(T, "P2", 0.00026872, -0.0024998, -0.00080616, -0.00080616);
        assertGoldenLayer(T, "T2", 0.0031701, -0.003003, -0.0067251, -0.0055844);
    }

    @Test
    @DisplayName("the column names carry the leading Layer column")
    public void testVariableNames() {
        LayeredNetworkSensitivityTable T = solver().getSensitivityTable();
        List<String> expected = Arrays.asList("Layer", "Station", "JobClass", "dTput_dRate",
                "dRespT_dRate", "dQLen_dRate", "dUtil_dRate");
        assertEquals(expected, T.getVariableNames());
        // get(int) must agree with the named accessors, i.e. the column order is real.
        assertEquals(T.getDTput(), T.get(0));
        assertEquals(T.getDRespT(), T.get(1));
        assertEquals(T.getDQLen(), T.get(2));
        assertEquals(T.getDUtil(), T.get(3));
    }

    @Test
    @DisplayName("the ensemble is solved on demand, and the explicit-option overload agrees")
    public void testSolvesOnDemandAndOverloadAgrees() {
        // getSensitivityTable is called on an unsolved solver: it must solve the
        // ensemble itself rather than differentiate an unconverged fixed point.
        LayeredNetworkSensitivityTable lazy = solver().getSensitivityTable();

        SolverLN eager = solver();
        eager.getAvgTable();
        LayeredNetworkSensitivityTable T = eager.getSensitivityTable("fd", Double.NaN, "forward");
        assertEquals(lazy.getLayerNames(), T.getLayerNames());
        for (int k = 0; k < T.getLayerNames().size(); k++) {
            assertRelEquals(lazy.getDQLen().get(k).doubleValue(), T.getDQLen().get(k).doubleValue(),
                    "dQLen_dRate row " + k);
        }
    }
}
