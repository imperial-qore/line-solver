/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.ln;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.junit.jupiter.api.Test;

import jline.VerboseLevel;
import jline.lang.constant.SchedStrategy;
import jline.lang.layered.Activity;
import jline.lang.layered.Entry;
import jline.lang.layered.LayeredNetwork;
import jline.lang.layered.Processor;
import jline.lang.layered.Task;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.solvers.LayeredNetworkAvgTable;
import jline.solvers.SolverOptions;
import jline.solvers.ldes.LDESOptions;
import jline.solvers.ldes.SolverLDES;
import jline.solvers.nc.SolverNC;

/**
 * Every LQN solver must agree about WHICH CELLS OF THE TABLE EXIST.
 *
 * <p>A NaN in an LQN average table is not a failed computation: it says the
 * quantity is not defined for that element kind. A processor has no queue
 * length, response time, residence or throughput of its own; an entry has no
 * residence; nothing reports an arrival rate. Two solvers that disagree about
 * the mask disagree about the MODEL rather than about arithmetic, and a
 * tolerance-based comparison cannot see it -- the goldens carry the string
 * {@code "NaN"} beside the numbers precisely so it is not lost.</p>
 *
 * <p>This is the standing form of the 2026-08-21 cross-check. LQNS is
 * deliberately NOT in the panel: it reports no residence time and no arrival
 * rate at all, and leaves both undefined by decision rather than by omission
 * (see {@code _kb/06-solver-catalog.md}, LQNS wrapper) -- and its binary is not
 * always present.</p>
 */
public class LqnNaNMaskParityTest {

    private static final String[] METRICS = {"QLen", "Util", "RespT", "ResidT", "ArvR", "Tput"};

    /** The model of {@code matlab/examples/basic/layeredModel/lqn_multi_solvers.m}. */
    private static LayeredNetwork multiSolvers() {
        LayeredNetwork m = new LayeredNetwork("LQN1");
        Processor p1 = new Processor(m, "P1", 1, SchedStrategy.INF);
        Task t1 = new Task(m, "T1", 1, SchedStrategy.REF).on(p1);
        Entry e1 = new Entry(m, "E1").on(t1);
        Processor p2 = new Processor(m, "P2", 1, SchedStrategy.INF);
        Task t2 = new Task(m, "T2", 1, SchedStrategy.INF).on(p2);
        Entry e2 = new Entry(m, "E2").on(t2);
        t1.setThinkTime(Erlang.fitMeanAndOrder(0.0001, 2));
        Activity a1 = new Activity(m, "A1", new Exp(1.0)).on(t1);
        a1.boundTo(e1);
        a1.synchCall(e2, 3.0);
        Activity a2 = new Activity(m, "A2", new Exp(1.0)).on(t2);
        a2.boundTo(e2);
        a2.repliesTo(e2);
        return m;
    }

    /**
     * The same model plus a component NOTHING REACHES: no reference task, no
     * caller. SolverLN marks such a component {@code ignore} and reports zero
     * for it -- which is right for the measures its elements HAVE, and was
     * being written over the ones they do not, flattening the mask. The
     * `lqn_ofbiz` golden carries exactly this shape (its {@code USAGE_DELAY}
     * component), which is why the defect survived the 2026-08-21 audit: no
     * golden pairs LQNS with an LN row on a model that has one.
     */
    private static LayeredNetwork ignoredComponent() {
        LayeredNetwork m = new LayeredNetwork("LQNignore");
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
        // unreachable: T3 is not a reference task and nobody calls E3
        Processor p3 = new Processor(m, "P3", 1, SchedStrategy.INF);
        Task t3 = new Task(m, "T3", 1, SchedStrategy.FCFS).on(p3);
        Entry e3 = new Entry(m, "E3").on(t3);
        Activity a3 = new Activity(m, "A3", new Exp(1.0)).on(t3);
        a3.boundTo(e3);
        a3.repliesTo(e3);
        return m;
    }

    private static List<Double> column(LayeredNetworkAvgTable t, String metric) {
        if ("QLen".equals(metric)) return t.getQLen();
        if ("Util".equals(metric)) return t.getUtil();
        if ("RespT".equals(metric)) return t.getRespT();
        if ("ResidT".equals(metric)) return t.getResidT();
        if ("ArvR".equals(metric)) return t.getArvR();
        return t.getTput();
    }

    /** "P1.QLen=NaN, P1.Util=num, ..." -- the mask, spelled so a failure reads. */
    private static Map<String, Boolean> maskOf(LayeredNetworkAvgTable t) {
        Map<String, Boolean> mask = new LinkedHashMap<String, Boolean>();
        List<String> names = t.getNodeNames();
        for (String metric : METRICS) {
            List<Double> col = column(t, metric);
            for (int i = 0; i < names.size() && i < col.size(); i++) {
                mask.put(names.get(i) + "." + metric, Double.isNaN(col.get(i)));
            }
        }
        return mask;
    }

    private static LayeredNetworkAvgTable lnTable(LayeredNetwork m, SolverFactory factory) {
        return lnTable(m, factory, null);
    }

    private static LayeredNetworkAvgTable lnTable(LayeredNetwork m, SolverFactory factory,
                                                  String method) {
        SolverOptions lnOpts = SolverLN.defaultOptions();
        lnOpts.verbose = VerboseLevel.SILENT;
        if (method != null) {
            lnOpts.method = method;
        }
        SolverLN solver = (factory == null)
                ? new SolverLN(m, lnOpts)
                : new SolverLN(m, factory, lnOpts);
        return (LayeredNetworkAvgTable) solver.getEnsembleAvg();
    }

    private static LayeredNetworkAvgTable ldesTable(LayeredNetwork m) {
        LDESOptions opt = new LDESOptions();
        opt.verbose = VerboseLevel.SILENT;
        opt.seed = 23000;
        opt.samples = 200000;
        SolverLDES solver = new SolverLDES(m, opt);
        solver.getAvg();
        return solver.getLNAvgTable();
    }

    private static void assertSameMask(String refName, Map<String, Boolean> ref,
                                       String otherName, Map<String, Boolean> other) {
        List<String> disagree = new ArrayList<String>();
        for (Map.Entry<String, Boolean> e : ref.entrySet()) {
            Boolean o = other.get(e.getKey());
            if (o == null) continue; // a solver that omits a column is not a mask defect
            if (!o.equals(e.getValue())) {
                disagree.add(e.getKey() + " (" + refName + (e.getValue() ? "=NaN" : "=value")
                        + ", " + otherName + (o ? "=NaN" : "=value") + ")");
            }
        }
        assertTrue(disagree.isEmpty(),
                "LQN table NaN mask differs between " + refName + " and " + otherName + ": "
                        + String.join("; ", disagree));
    }

    @Test
    public void everySolverAgreesOnWhichCellsExist() {
        LayeredNetwork m = multiSolvers();

        Map<String, Boolean> mva = maskOf(lnTable(m, null));
        Map<String, Boolean> nc = maskOf(lnTable(multiSolvers(), new SolverFactory() {
            @Override
            public jline.solvers.NetworkSolver at(jline.lang.Network net) {
                SolverOptions o = new SolverOptions(jline.lang.constant.SolverType.NC);
                o.verbose = VerboseLevel.SILENT;
                return new SolverNC(net, o);
            }
        }));
        Map<String, Boolean> ldes = maskOf(ldesTable(multiSolvers()));

        assertSameMask("LN(MVA)", mva, "LN(NC)", nc);
        assertSameMask("LN(MVA)", mva, "LDES", ldes);
    }

    /**
     * A component no reference task reaches is IDLE, not undefined.
     *
     * <p>Its elements report zero for every measure their kind has, and NaN for
     * the ones it does not -- the same mask a reachable element carries. Until
     * 2026-08-25 {@code getEnsembleAvg} wrote a flat zero across all six
     * columns, so an unreachable processor claimed a queue length of 0 and an
     * arrival rate of 0 where every solver, LQNS included, reports neither.</p>
     */
    @Test
    public void anIgnoredComponentKeepsTheMask() {
        Map<String, Boolean> mask = maskOf(lnTable(ignoredComponent(), null));
        Map<String, Boolean> reachable = maskOf(lnTable(multiSolvers(), null));

        // the unreachable rows carry the SAME mask as the reachable ones
        assertEquals(reachable.get("P1.QLen"), mask.get("P3.QLen"), "P3.QLen");
        assertEquals(reachable.get("P1.RespT"), mask.get("P3.RespT"), "P3.RespT");
        assertEquals(reachable.get("P1.ResidT"), mask.get("P3.ResidT"), "P3.ResidT");
        assertEquals(reachable.get("P1.Tput"), mask.get("P3.Tput"), "P3.Tput");
        assertEquals(reachable.get("P1.Util"), mask.get("P3.Util"), "P3.Util");
        assertEquals(reachable.get("T2.RespT"), mask.get("T3.RespT"), "T3.RespT");
        assertEquals(reachable.get("T2.ResidT"), mask.get("T3.ResidT"), "T3.ResidT");
        assertEquals(reachable.get("E2.RespT"), mask.get("E3.RespT"), "E3.RespT");
        assertEquals(reachable.get("E2.ResidT"), mask.get("E3.ResidT"), "E3.ResidT");
        assertEquals(reachable.get("A2.ResidT"), mask.get("A3.ResidT"), "A3.ResidT");

        // spelled out, so the row above cannot pass by both sides regressing
        assertTrue(mask.get("P3.QLen"), "P3.QLen must be NaN");
        assertTrue(mask.get("P3.RespT"), "P3.RespT must be NaN");
        assertTrue(mask.get("P3.Tput"), "P3.Tput must be NaN");
        assertEquals(Boolean.FALSE, mask.get("P3.Util"), "P3.Util must be a value");
        assertTrue(mask.get("T3.RespT"), "T3.RespT must be NaN");
        assertEquals(Boolean.FALSE, mask.get("T3.ResidT"), "T3.ResidT must be a value");
        assertTrue(mask.get("E3.ResidT"), "E3.ResidT must be NaN");
        for (String x : new String[]{"P3", "T3", "E3", "A3"}) {
            assertTrue(mask.get(x + ".ArvR"), x + ".ArvR must be NaN");
        }

        // and every solver still agrees about it
        Map<String, Boolean> ldes = maskOf(ldesTable(ignoredComponent()));
        assertSameMask("LN(MVA)", mask, "LDES", ldes);
    }

    /**
     * {@code srvn.ph} assembles the table in its own routine, and had its own
     * copy of the flat-zero branch.
     *
     * <p>The two encodings rebuild every figure differently -- one reads class
     * rows off the ensemble, the other composes a phase-type law -- so agreeing
     * on the mask is a claim about the table, not about shared code.</p>
     */
    @Test
    public void thePhEncodingMasksAnIgnoredComponentTheSameWay() {
        Map<String, Boolean> dflt = maskOf(lnTable(ignoredComponent(), null));
        Map<String, Boolean> ph = maskOf(lnTable(ignoredComponent(), null, "srvn.ph"));
        assertSameMask("LN(MVA)", dflt, "LN(MVA, srvn.ph)", ph);
        assertTrue(ph.get("P3.QLen"), "P3.QLen must be NaN under srvn.ph");
        assertTrue(ph.get("P3.ArvR"), "P3.ArvR must be NaN under srvn.ph");
        assertEquals(Boolean.FALSE, ph.get("P3.Util"), "P3.Util must be a value under srvn.ph");
    }

    /**
     * The mask itself, pinned. Without this the test above would still pass if
     * every solver regressed the same way -- LDES reporting a flat 0.0 for
     * ResidT, say, which is what it did until 2026-08-21.
     */
    @Test
    public void theMaskIsTheDocumentedOne() {
        Map<String, Boolean> mask = maskOf(lnTable(multiSolvers(), null));

        // A processor has no queue, no response time, no residence, no
        // throughput of its own -- only a utilization. TN is set NaN explicitly
        // in getEnsembleAvg "for consistency with LQNS".
        for (String p : new String[]{"P1", "P2"}) {
            assertTrue(mask.get(p + ".QLen"), p + ".QLen must be NaN");
            assertTrue(mask.get(p + ".RespT"), p + ".RespT must be NaN");
            assertTrue(mask.get(p + ".ResidT"), p + ".ResidT must be NaN");
            assertTrue(mask.get(p + ".Tput"), p + ".Tput must be NaN");
            assertEquals(Boolean.FALSE, mask.get(p + ".Util"), p + ".Util must be a value");
        }
        // An entry has a response time but no residence.
        for (String e : new String[]{"E1", "E2"}) {
            assertTrue(mask.get(e + ".ResidT"), e + ".ResidT must be NaN");
            assertEquals(Boolean.FALSE, mask.get(e + ".RespT"), e + ".RespT must be a value");
        }
        // Tasks and activities carry both.
        for (String x : new String[]{"T1", "T2", "A1", "A2"}) {
            assertEquals(Boolean.FALSE, mask.get(x + ".ResidT"), x + ".ResidT must be a value");
        }
        // Nobody reports an arrival rate on an LQN.
        for (String x : new String[]{"P1", "T1", "E1", "A1"}) {
            assertTrue(mask.get(x + ".ArvR"), x + ".ArvR must be NaN");
        }
    }
}
