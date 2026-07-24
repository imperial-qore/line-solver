/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api;

import jline.VerboseLevel;
import jline.api.lqn.Lqn_boxbounds;
import jline.lang.constant.SchedStrategy;
import jline.lang.layered.Activity;
import jline.lang.layered.Entry;
import jline.lang.layered.LayeredNetwork;
import jline.lang.layered.Processor;
import jline.lang.layered.Task;
import jline.lang.processes.Exp;
import jline.solvers.AvgTable;
import jline.solvers.LayeredNetworkAvgTable;
import jline.solvers.mva.SolverMVA;
import jline.solvers.ln.SolverLN;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Validates the LQN Majumdar-Woodside robust box bounds ({@link Lqn_boxbounds}
 * and the SolverLN mwba.upper/mwba.lower methods) on a two-tier client-server
 * LQN: reference task T1 (mult 10, think 5.0) on processor P1 with a single
 * activity (host demand 1.0) that makes 2.5 synchronous calls to entry E2 on an
 * infinite DB task on processor P2 (host demand 0.8).
 *
 * <p>Per-chain processor demands: D_P1 = 1.0, D_P2 = 2.5*0.8 = 2.0, Z = 5,
 * N = 10. The processor-utilization upper bound gives X_ref+ = 1/D_P2 = 0.5
 * (P2 is the bottleneck) and the Majumdar-Woodside lower bound gives
 * X_ref- = N/(Z + N*(D_P1+D_P2)) = 10/35 = 0.2857. The exact LN throughput
 * (~0.499) lies inside [0.2857, 0.5].</p>
 */
public class LqnBoxBoundsTest {

    private static final double TOL = 1e-4;

    private LayeredNetwork buildModel() {
        LayeredNetwork model = new LayeredNetwork("cd");
        Processor P1 = new Processor(model, "P1", 1, SchedStrategy.PS);
        Processor P2 = new Processor(model, "P2", 1, SchedStrategy.PS);
        Task T1 = new Task(model, "T1", 10, SchedStrategy.REF).on(P1);
        T1.setThinkTime(new Exp(1.0 / 5.0));
        Task T2 = new Task(model, "T2", 1, SchedStrategy.INF).on(P2);
        Entry E1 = new Entry(model, "E1").on(T1);
        Entry E2 = new Entry(model, "E2").on(T2);
        Activity A1 = new Activity(model, "A1", new Exp(1.0)).on(T1).boundTo(E1).synchCall(E2, 2.5);
        Activity A2 = new Activity(model, "A2", new Exp(1.0 / 0.8)).on(T2).boundTo(E2).repliesTo(E2);
        return model;
    }

    @Test
    public void testApiThroughputBounds() {
        LayeredNetwork model = buildModel();
        Lqn_boxbounds.Result r = Lqn_boxbounds.compute(model.getStruct());
        assertEquals(1, r.Xup.length, "single reference chain");
        assertEquals(0.5, r.Xup[0], TOL, "upper (P2 bottleneck)");
        assertEquals(10.0 / 35.0, r.Xlo[0], TOL, "lower (Muntz-Wong)");
        assertTrue(r.Xlo[0] <= r.Xup[0] + TOL);
    }

    @Test
    public void testSolverLNBoundsBracketExact() {
        SolverLN exactSolver = new SolverLN(buildModel(), m -> new SolverMVA(m));
        exactSolver.options.verbose = VerboseLevel.SILENT;
        double exact = refTput((LayeredNetworkAvgTable) exactSolver.getAvgTable(), "T1");

        SolverLN up = new SolverLN(buildModel(), m -> new SolverMVA(m));
        up.options.method = "mwba.upper";
        up.options.verbose = VerboseLevel.SILENT;
        double xUp = refTput((LayeredNetworkAvgTable) up.getAvgTable(), "T1");

        SolverLN lo = new SolverLN(buildModel(), m -> new SolverMVA(m));
        lo.options.method = "mwba.lower";
        lo.options.verbose = VerboseLevel.SILENT;
        double xLo = refTput((LayeredNetworkAvgTable) lo.getAvgTable(), "T1");

        assertEquals(0.5, xUp, TOL, "SolverLN mwba.upper T1");
        assertEquals(10.0 / 35.0, xLo, TOL, "SolverLN mwba.lower T1");
        if (!Double.isNaN(exact)) {
            assertTrue(xLo <= exact + TOL && exact <= xUp + TOL,
                    "exact " + exact + " must lie within [" + xLo + ", " + xUp + "]");
        }
    }

    private static double refTput(LayeredNetworkAvgTable t, String node) {
        List<String> names = t.getNodeNames();
        List<Double> tput = t.getTput();
        for (int i = 0; i < names.size(); i++) {
            if (node.equals(names.get(i))) {
                return tput.get(i);
            }
        }
        return Double.NaN;
    }
}
