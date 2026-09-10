/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.wrappers.lqns;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assumptions.assumeTrue;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.junit.jupiter.api.Test;

import jline.VerboseLevel;
import jline.lang.Network;
import jline.lang.constant.SchedStrategy;
import jline.lang.layered.Activity;
import jline.lang.layered.ActivityPrecedence;
import jline.lang.layered.Entry;
import jline.lang.layered.LayeredNetwork;
import jline.lang.layered.Processor;
import jline.lang.layered.Task;
import jline.lang.processes.Exp;
import jline.solvers.LayeredNetworkAvgTable;
import jline.solvers.NetworkSolver;
import jline.solvers.SolverOptions;
import jline.solvers.ln.SolverFactory;
import jline.solvers.ln.SolverLN;
import jline.solvers.mva.SolverMVA;
import jline.util.matrix.Matrix;

/**
 * An entry lqns never invoked has a service time of ZERO, not an absent one.
 *
 * <p>lqns omits {@code phase1-service-time} from {@code result-entry} exactly
 * when the entry's throughput is zero -- nothing was served, so there is no
 * per-invocation mean to report. Read verbatim that is a NaN, and it landed in
 * the RespT column where the LQN table says an entry HAS a response time and
 * every other solver reports one. The NaN mask is part of the answer (see
 * {@code _kb/06-solver-catalog.md}), and a tolerance-based parity comparison
 * cannot see a break in it: {@code compare_values} passes any cell where either
 * side is NaN.</p>
 *
 * <p>The value is derived from the activity rows and ONLY where they are
 * unanimous: if every activity reachable from the entry reports a zero service
 * time then every aggregation law agrees on zero -- the serial sum, the
 * branch-weighted mean of an OrFork, the order statistic of an AndFork. It is
 * deliberately not generalised the way {@link SolverLQNSEntryProcUtilTest}
 * generalises utilization, which DOES add over an activity graph; the last test
 * below is what pins that distinction down.</p>
 */
public class SolverLQNSEntryServiceTimeTest {

    /**
     * A working two-tier model plus a component no reference task reaches.
     *
     * <p>T3 is not a reference task and nobody calls E3, so lqns solves it at
     * zero throughput and writes a {@code result-entry} with no
     * {@code phase1-service-time}. This is the shape {@code lqn_ofbiz} carries
     * in its USAGE_DELAY component.</p>
     */
    private static LayeredNetwork unreachableEntry() {
        LayeredNetwork m = new LayeredNetwork("lqn_unreachable_entry");
        Processor p1 = new Processor(m, "P1", 1, SchedStrategy.INF);
        Task t1 = new Task(m, "T1", 1, SchedStrategy.REF).on(p1);
        Entry e1 = new Entry(m, "E1").on(t1);
        Processor p2 = new Processor(m, "P2", 1, SchedStrategy.INF);
        Task t2 = new Task(m, "T2", 1, SchedStrategy.INF).on(p2);
        Entry e2 = new Entry(m, "E2").on(t2);
        t1.setThinkTime(new Exp(1.0));
        Activity a1 = new Activity(m, "A1", new Exp(1.0)).on(t1);
        a1.boundTo(e1);
        a1.synchCall(e2, 1.0);
        Activity a2 = new Activity(m, "A2", new Exp(1.0)).on(t2);
        a2.boundTo(e2);
        a2.repliesTo(e2);

        Processor p3 = new Processor(m, "P3", 1, SchedStrategy.INF);
        Task t3 = new Task(m, "T3", 1, SchedStrategy.FCFS).on(p3);
        Entry e3 = new Entry(m, "E3").on(t3);
        // TWO activities in series, so the writer emits the activity-graph
        // (task-activities) form. A single bound activity is written as
        // entry-phase-activities instead, where the reader already has a
        // fallback of its own and the omission never surfaces.
        Activity a3 = new Activity(m, "A3", new Exp(1.0)).on(t3);
        a3.boundTo(e3);
        Activity a3b = new Activity(m, "A3b", new Exp(1.0)).on(t3);
        a3b.repliesTo(e3);
        t3.addPrecedence(ActivityPrecedence.Serial(a3, a3b));
        return m;
    }

    /**
     * An entry whose activity graph BRANCHES, so its service time is not a sum.
     *
     * <p>E2 runs A20 and then an OrFork to A21 or A22. lqns reports the entry's
     * {@code phase1-service-time} as the branch-weighted total; summing the
     * activity rows would over-count the arm not taken.</p>
     */
    private static LayeredNetwork branchingEntry() {
        LayeredNetwork m = new LayeredNetwork("lqn_branching_entry");
        Processor p1 = new Processor(m, "P1", 1, SchedStrategy.PS);
        Task t1 = new Task(m, "T1", 10, SchedStrategy.REF).on(p1);
        Entry e1 = new Entry(m, "E1").on(t1);
        Processor p2 = new Processor(m, "P2", 1, SchedStrategy.PS);
        Task t2 = new Task(m, "T2", 1, SchedStrategy.INF).on(p2);
        Entry e2 = new Entry(m, "E2").on(t2);

        t1.setThinkTime(new Exp(0.1));
        Activity a1 = new Activity(m, "A1", new Exp(1.0)).on(t1);
        a1.boundTo(e1);
        a1.synchCall(e2, 1.0);

        Activity a20 = new Activity(m, "A20", new Exp(1.0)).on(t2);
        a20.boundTo(e2);
        Activity a21 = new Activity(m, "A21", new Exp(1.0)).on(t2);
        a21.repliesTo(e2);
        Activity a22 = new Activity(m, "A22", new Exp(1.0)).on(t2);
        a22.repliesTo(e2);
        Matrix probs = new Matrix(1, 2);
        probs.set(0, 0, 0.5);
        probs.set(0, 1, 0.5);
        t2.addPrecedence(ActivityPrecedence.OrFork(a20, java.util.Arrays.asList(a21, a22), probs));
        return m;
    }

    private static Map<String, Integer> indexOf(LayeredNetworkAvgTable t) {
        Map<String, Integer> out = new HashMap<String, Integer>();
        List<String> names = t.getNodeNames();
        for (int i = 0; i < names.size(); i++) {
            out.put(names.get(i), i);
        }
        return out;
    }

    private static LayeredNetworkAvgTable lqnsTable(LayeredNetwork m) {
        SolverLQNS solver = new SolverLQNS(m);
        return (LayeredNetworkAvgTable) solver.getAvgTable();
    }

    @Test
    public void anEntryLqnsNeverInvokedReportsZeroNotNaN() {
        assumeTrue(SolverLQNS.isAvailable(), "no lqns binary on the PATH");
        LayeredNetworkAvgTable t = lqnsTable(unreachableEntry());
        Map<String, Integer> ix = indexOf(t);

        // the unreachable component solves at zero throughput
        assertEquals(0.0, t.getTput().get(ix.get("E3")), 1e-12);
        assertEquals(0.0, t.getTput().get(ix.get("A3")), 1e-12);
        assertEquals(0.0, t.getTput().get(ix.get("A3b")), 1e-12);

        // the regression itself: RespT was NaN, because lqns omits the attribute
        assertFalse(Double.isNaN(t.getRespT().get(ix.get("E3"))), "E3.RespT must not be NaN");
        assertEquals(0.0, t.getRespT().get(ix.get("E3")), 1e-12);

        // the reachable entries are untouched and still carry lqns' own numbers
        assertEquals(2.0, t.getRespT().get(ix.get("E1")), 2e-3);
        assertEquals(1.0, t.getRespT().get(ix.get("E2")), 1e-3);
    }

    /** SolverLN is the reference for which cells exist; LQNS must not differ. */
    @Test
    public void itAgreesWithLnOnTheSameModel() {
        assumeTrue(SolverLQNS.isAvailable(), "no lqns binary on the PATH");
        LayeredNetworkAvgTable lqns = lqnsTable(unreachableEntry());
        SolverOptions lnOpts = SolverLN.defaultOptions();
        lnOpts.verbose = VerboseLevel.SILENT;
        SolverLN ln = new SolverLN(unreachableEntry(), new SolverFactory() {
            @Override
            public NetworkSolver at(Network net) {
                SolverOptions o = new SolverOptions(jline.lang.constant.SolverType.MVA);
                o.verbose = VerboseLevel.SILENT;
                return new SolverMVA(net, o);
            }
        }, lnOpts);
        LayeredNetworkAvgTable lnt = (LayeredNetworkAvgTable) ln.getEnsembleAvg();

        Map<String, Integer> a = indexOf(lqns);
        Map<String, Integer> b = indexOf(lnt);
        String[] rows = {"P3", "T3", "E3", "A3", "A3b"};
        String[] metrics = {"QLen", "Util", "RespT", "ArvR", "Tput"};
        for (int r = 0; r < rows.length; r++) {
            for (int c = 0; c < metrics.length; c++) {
                // ResidT is excluded by decision: lqns reports no residence time at all.
                boolean x = Double.isNaN(column(lqns, metrics[c]).get(a.get(rows[r])));
                boolean y = Double.isNaN(column(lnt, metrics[c]).get(b.get(rows[r])));
                assertEquals(y, x, rows[r] + "." + metrics[c] + " differs between LQNS and LN");
            }
        }
    }

    /**
     * Where lqns DOES report the attribute, its value survives verbatim.
     *
     * <p>Utilizations add over an activity graph and response times do not, so
     * the fallback must never become a sum: on this OrFork the sum over the
     * entry's activities exceeds what lqns reports, and the reported value is
     * the one that must come back.</p>
     */
    @Test
    public void theDerivationDoesNotTouchABranchingEntry() {
        assumeTrue(SolverLQNS.isAvailable(), "no lqns binary on the PATH");
        LayeredNetworkAvgTable t = lqnsTable(branchingEntry());
        Map<String, Integer> ix = indexOf(t);

        double e2 = t.getRespT().get(ix.get("E2"));
        assertFalse(Double.isNaN(e2), "E2.RespT must not be NaN");
        double summed = t.getRespT().get(ix.get("A20")) + t.getRespT().get(ix.get("A21"))
                + t.getRespT().get(ix.get("A22"));
        // one arm of the fork is not taken, so the sum over-counts
        assertTrue(summed > e2 + 1e-6,
                "the fixture no longer branches: sum " + summed + " against reported " + e2);
    }

    private static List<Double> column(LayeredNetworkAvgTable t, String metric) {
        if ("QLen".equals(metric)) return t.getQLen();
        if ("Util".equals(metric)) return t.getUtil();
        if ("RespT".equals(metric)) return t.getRespT();
        if ("ArvR".equals(metric)) return t.getArvR();
        return t.getTput();
    }
}
