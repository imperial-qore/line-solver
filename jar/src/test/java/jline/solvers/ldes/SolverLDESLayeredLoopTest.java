/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.ldes;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;

import org.junit.jupiter.api.Test;

import jline.VerboseLevel;
import jline.lang.constant.SchedStrategy;
import jline.lang.layered.Activity;
import jline.lang.layered.ActivityPrecedence;
import jline.lang.layered.Entry;
import jline.lang.layered.LayeredNetwork;
import jline.lang.layered.Processor;
import jline.lang.layered.Task;
import jline.lang.processes.Exp;
import jline.lang.processes.Immediate;
import jline.solvers.LayeredNetworkAvgTable;

/**
 * Validates LDES simulation of LQN activity loops (POST_LOOP precedence).
 *
 * <p>The models use a single-job reference task on an infinite-server processor with an
 * immediate think time, so there is no queueing contention and the reference-task throughput
 * has the closed form {@code 1/E[work]}, where {@code E[work]} is the expected total host
 * demand along the activity graph. With unit host demands the activity chain
 * {@code AStart -> [ABody] -> AEnd} yields:</p>
 * <ul>
 *   <li>no loop: {@code E[work] = 3}, throughput {@code 1/3};</li>
 *   <li>loop mean 2: {@code E[work] = 1 + 2 + 1 = 4}, throughput {@code 1/4};</li>
 *   <li>loop mean 3: {@code E[work] = 1 + 3 + 1 = 5}, throughput {@code 1/5}.</li>
 * </ul>
 *
 * <p>The mean-2 case is the regression guard: the loop is encoded as a geometric branch with
 * probabilities {@code 1-1/2 = 0.5} (loop back) and {@code 1/2 = 0.5} (exit). These equal
 * weights previously fell through the generic uniform-fork fallback and the body executed a
 * fixed single time. The loop-body throughput must equal the loop count times the invocation
 * throughput.</p>
 */
public class SolverLDESLayeredLoopTest {

    private static final double REL_TOL = 5e-2;
    private static final int SAMPLES = 500000;
    private static final int SEED = 23000;

    private static LayeredNetwork buildModel(String name, boolean useLoop, double nloops) {
        LayeredNetwork m = new LayeredNetwork(name);
        Processor p1 = new Processor(m, "P1", 1, SchedStrategy.INF);
        Task t1 = new Task(m, "T1", 1, SchedStrategy.REF).on(p1).setThinkTime(new Immediate());
        Entry e1 = new Entry(m, "E1").on(t1);
        Activity aStart = new Activity(m, "AStart", new Exp(1.0)).on(t1);
        aStart.boundTo(e1);
        new Activity(m, "ABody", new Exp(1.0)).on(t1);
        new Activity(m, "AEnd", new Exp(1.0)).on(t1);
        if (useLoop) {
            t1.addPrecedence(ActivityPrecedence.Loop("AStart", "ABody", "AEnd", nloops));
        } else {
            t1.addPrecedence(ActivityPrecedence.Serial("AStart", "ABody"));
            t1.addPrecedence(ActivityPrecedence.Serial("ABody", "AEnd"));
        }
        return m;
    }

    private static LayeredNetworkAvgTable solve(LayeredNetwork m) {
        LDESOptions opt = new LDESOptions();
        opt.verbose = VerboseLevel.SILENT;
        opt.seed = SEED;
        opt.samples = SAMPLES;
        SolverLDES solver = new SolverLDES(m, opt);
        solver.getAvg();
        return solver.getLNAvgTable();
    }

    private static double tputOf(LayeredNetworkAvgTable table, String elementName) {
        List<String> names = table.getNodeNames();
        List<Double> tput = table.getTput();
        for (int i = 0; i < names.size(); i++) {
            if (names.get(i).equals(elementName)) {
                return tput.get(i);
            }
        }
        throw new IllegalArgumentException("Element not found in LN table: " + elementName);
    }

    @Test
    public void testNoLoopBaseline() {
        LayeredNetworkAvgTable table = solve(buildModel("noloop", false, 0));
        double t1 = tputOf(table, "T1");
        assertEquals(1.0 / 3.0, t1, REL_TOL * (1.0 / 3.0),
                "No-loop reference-task throughput should be 1/E[work] = 1/3");
    }

    @Test
    public void testLoopMeanTwo() {
        // Regression guard: loop-back and exit probabilities are both 0.5 here.
        LayeredNetworkAvgTable table = solve(buildModel("loop2", true, 2.0));
        double t1 = tputOf(table, "T1");
        double body = tputOf(table, "ABody");
        assertEquals(0.25, t1, REL_TOL * 0.25,
                "Loop mean 2 reference-task throughput should be 1/(1+2+1) = 0.25");
        assertTrue(Math.abs(body - 2.0 * t1) <= REL_TOL * (2.0 * t1),
                "Loop body should execute a mean of 2 times per invocation (ABody Tput = 2 x T1 Tput), got "
                        + body + " vs " + (2.0 * t1));
    }

    @Test
    public void testLoopMeanThree() {
        LayeredNetworkAvgTable table = solve(buildModel("loop3", true, 3.0));
        double t1 = tputOf(table, "T1");
        double body = tputOf(table, "ABody");
        assertEquals(0.20, t1, REL_TOL * 0.20,
                "Loop mean 3 reference-task throughput should be 1/(1+3+1) = 0.20");
        assertTrue(Math.abs(body - 3.0 * t1) <= REL_TOL * (3.0 * t1),
                "Loop body should execute a mean of 3 times per invocation (ABody Tput = 3 x T1 Tput), got "
                        + body + " vs " + (3.0 * t1));
    }
}
