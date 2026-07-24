/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers;

import jline.io.Ret;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Disabled;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.fluid.SolverFluid;
import jline.solvers.mva.SolverMVA;
import jline.solvers.nc.SolverNC;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.function.Executable;

import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Regression tests for {@link NetworkSolver#getSensitivityTable()}, the solver-level
 * view over the analytic rate sensitivities.
 *
 * <p>The goldens are MATLAB's, dumped at 17 significant digits from
 * {@code @NetworkSolver/getSensitivityTable.m} on the very models built here, so
 * these tests pin cross-codebase parity and not merely self-consistency. MATLAB is
 * the ground truth.</p>
 *
 * <p>The mathematics behind the closed branch ({@code pfqn_sens}) is validated
 * separately by the api-level harness {@link jline.api.PfqnSensComomTest}; these
 * tests check the solver-level plumbing, which can go wrong independently of a
 * correct algorithm: the demands and rates handed to the api layer, the chain rule
 * from demand to rate, the row-skipping rule, and the error contract.</p>
 */
public class SensitivityTableTest {

    /** Agreement required against the MATLAB goldens. */
    private static final double REL_TOL = 1e-9;

    // ---------- helpers ---------------------------------------------------

    /**
     * Closed, two classes, Q1 FCFS and Q2 PS, with a think-time delay. Identical to
     * the model the MATLAB goldens were dumped from.
     */
    private static Network closedTwoClass() {
        Network m = new Network("mt_closed");
        Delay d = new Delay(m, "Think");
        Queue q1 = new Queue(m, "Q1", SchedStrategy.FCFS);
        Queue q2 = new Queue(m, "Q2", SchedStrategy.PS);
        ClosedClass c1 = new ClosedClass(m, "C1", 3, d, 0);
        ClosedClass c2 = new ClosedClass(m, "C2", 2, d, 0);
        d.setService(c1, new Exp(1 / 1.0));
        d.setService(c2, new Exp(1 / 0.5));
        q1.setService(c1, new Exp(1 / 0.4));
        q1.setService(c2, new Exp(1 / 0.4));
        q2.setService(c1, new Exp(1 / 0.3));
        q2.setService(c2, new Exp(1 / 0.2));
        RoutingMatrix P = m.initRoutingMatrix();
        P.set(c1, c1, Network.serialRouting(d, q1, q2));
        P.set(c2, c2, Network.serialRouting(d, q1, q2));
        m.link(P);
        return m;
    }

    /** Purely open, two classes. Identical to the MATLAB golden model. */
    private static Network openTwoClass() {
        Network m = new Network("mt_open");
        Source src = new Source(m, "Src");
        Queue q1 = new Queue(m, "Q1", SchedStrategy.FCFS);
        Queue q2 = new Queue(m, "Q2", SchedStrategy.PS);
        Sink snk = new Sink(m, "Sink");
        OpenClass o1 = new OpenClass(m, "O1");
        OpenClass o2 = new OpenClass(m, "O2");
        src.setArrival(o1, new Exp(0.5));
        src.setArrival(o2, new Exp(0.3));
        q1.setService(o1, new Exp(1 / 0.4));
        q1.setService(o2, new Exp(1 / 0.4));
        q2.setService(o1, new Exp(1 / 0.3));
        q2.setService(o2, new Exp(1 / 0.2));
        RoutingMatrix P = m.initRoutingMatrix();
        P.set(o1, o1, Network.serialRouting(src, q1, q2, snk));
        P.set(o2, o2, Network.serialRouting(src, q1, q2, snk));
        m.link(P);
        return m;
    }

    /** Closed with a two-server queue, which is out of scope. */
    private static Network closedMultiserver() {
        Network m = new Network("mt_multi");
        Delay d = new Delay(m, "Think");
        Queue q = new Queue(m, "Q1", SchedStrategy.FCFS);
        q.setNumberOfServers(2);
        ClosedClass c = new ClosedClass(m, "C1", 3, d, 0);
        d.setService(c, new Exp(1 / 1.0));
        q.setService(c, new Exp(1 / 0.4));
        RoutingMatrix P = m.initRoutingMatrix();
        P.set(c, c, Network.serialRouting(d, q));
        m.link(P);
        return m;
    }

    /** Mixed open and closed, which is out of scope. */
    private static Network mixedOpenClosed() {
        Network m = new Network("mt_mixed");
        Source src = new Source(m, "Src");
        Delay d = new Delay(m, "Think");
        Queue q = new Queue(m, "Q1", SchedStrategy.PS);
        Sink snk = new Sink(m, "Sink");
        ClosedClass cc = new ClosedClass(m, "C1", 2, d, 0);
        OpenClass oc = new OpenClass(m, "O1");
        d.setService(cc, new Exp(1 / 1.0));
        q.setService(cc, new Exp(1 / 0.4));
        src.setArrival(oc, new Exp(0.3));
        q.setService(oc, new Exp(1 / 0.3));
        d.setService(oc, new Disabled());
        RoutingMatrix P = m.initRoutingMatrix();
        P.set(cc, cc, Network.serialRouting(d, q));
        P.set(oc, oc, Network.serialRouting(src, q, snk));
        m.link(P);
        return m;
    }

    /**
     * Closed, two classes, where C2 never visits Q2. The (Q2, C2) pair has no demand
     * and no rate, so MATLAB skips its row; this model pins that rule.
     */
    private static Network closedWithSkippedRow() {
        Network m = new Network("mt_skip");
        Delay d = new Delay(m, "Think");
        Queue q1 = new Queue(m, "Q1", SchedStrategy.PS);
        Queue q2 = new Queue(m, "Q2", SchedStrategy.PS);
        ClosedClass c1 = new ClosedClass(m, "C1", 2, d, 0);
        ClosedClass c2 = new ClosedClass(m, "C2", 2, d, 0);
        d.setService(c1, new Exp(1 / 1.0));
        d.setService(c2, new Exp(1 / 0.5));
        q1.setService(c1, new Exp(1 / 0.4));
        q1.setService(c2, new Exp(1 / 0.6));
        q2.setService(c1, new Exp(1 / 0.3));
        q2.setService(c2, new Disabled());
        RoutingMatrix P = m.initRoutingMatrix();
        P.set(c1, c1, Network.serialRouting(d, q1, q2));
        P.set(c2, c2, Network.serialRouting(d, q1));
        m.link(P);
        return m;
    }

    /**
     * Closed, one chain made of two classes that switch into one another. The chain
     * population sits on the reference class C1 and the switched class C2 carries
     * njobs = 0, so the analytic branch must aggregate to the chain before it
     * differentiates.
     */
    private static Network closedClassSwitching() {
        Network m = new Network("mt_cs");
        Delay d = new Delay(m, "Think");
        Queue q1 = new Queue(m, "Q1", SchedStrategy.PS);
        Queue q2 = new Queue(m, "Q2", SchedStrategy.PS);
        ClosedClass c1 = new ClosedClass(m, "C1", 2, d, 0);
        ClosedClass c2 = new ClosedClass(m, "C2", 0, d, 0);
        d.setService(c1, new Exp(1 / 1.0));
        d.setService(c2, new Exp(1 / 1.0));
        q1.setService(c1, new Exp(1 / 0.4));
        q1.setService(c2, new Exp(1 / 0.4));
        q2.setService(c1, new Exp(1 / 0.3));
        q2.setService(c2, new Exp(1 / 0.3));
        RoutingMatrix P = m.initRoutingMatrix();
        P.set(c1, c1, d, q1, 1.0);
        P.set(c1, c2, q1, q2, 1.0);   // C1 becomes C2 on leaving Q1
        P.set(c2, c1, q2, d, 1.0);    // and back to C1 on leaving Q2
        m.link(P);
        return m;
    }

    private static int rowOf(NetworkSensitivityTable T, String station, String jobclass) {
        int k = T.findRow(station, jobclass);
        assertTrue(k >= 0, "row " + station + "/" + jobclass + " missing from the sensitivity table");
        return k;
    }

    private static void assertRelEquals(double expected, double actual, String msg) {
        double denom = Math.max(Math.abs(expected), 1e-300);
        assertTrue(Math.abs(actual - expected) <= REL_TOL * denom,
                msg + ": expected " + expected + " but got " + actual
                        + " (rel err " + (Math.abs(actual - expected) / denom) + ")");
    }

    /**
     * Asserts one golden row: the four MATLAB values in table-column order.
     */
    private static void assertGoldenRow(NetworkSensitivityTable T, String station, String jobclass,
                                        double dTput, double dRespT, double dQLen, double dUtil) {
        int k = rowOf(T, station, jobclass);
        String tag = station + "/" + jobclass;
        assertRelEquals(dTput, T.getDTput().get(k), "dTput_dRate " + tag);
        assertRelEquals(dRespT, T.getDRespT().get(k), "dRespT_dRate " + tag);
        assertRelEquals(dQLen, T.getDQLen().get(k), "dQLen_dRate " + tag);
        assertRelEquals(dUtil, T.getDUtil().get(k), "dUtil_dRate " + tag);
    }

    // ---------- tests -----------------------------------------------------

    @Test
    @DisplayName("closed multiclass single-server matches the MATLAB goldens")
    public void testClosedMatchesMatlab() {
        // Goldens from MATLAB @NetworkSolver/getSensitivityTable.m on this model.
        NetworkSensitivityTable T = new SolverMVA(closedTwoClass()).getSensitivityTable();
        assertEquals(4, T.getStationNames().size(), "one row per visited (station, class)");
        assertGoldenRow(T, "Q1", "C1",
                0.22451658372356662, -0.55774567430504551,
                -0.38366074588642357, -0.091727874240920729);
        assertGoldenRow(T, "Q1", "C2",
                0.25555546097345389, -0.47054330708748149,
                -0.23565954411541706, -0.067045142260340432);
        assertGoldenRow(T, "Q2", "C1",
                0.071671442426249982, -0.20212629560055964,
                -0.19102956404839269, -0.080611727870445404);
        assertGoldenRow(T, "Q2", "C2",
                0.038597229825392973, -0.082536980672119842,
                -0.073239521727041554, -0.034597385697351907);
    }

    @Test
    @DisplayName("purely open matches the MATLAB goldens of the closed-form BCMP branch")
    public void testOpenMatchesMatlab() {
        NetworkSensitivityTable T = new SolverMVA(openTwoClass()).getSensitivityTable();
        assertEquals(4, T.getStationNames().size(), "one row per visited (station, class)");
        assertGoldenRow(T, "Q1", "O1",
                0.0, -0.30449826989619377,
                -0.15224913494809689, -0.080000000000000002);
        assertGoldenRow(T, "Q1", "O2",
                0.0, -0.27681660899653987,
                -0.083044982698961947, -0.048000000000000001);
        assertGoldenRow(T, "Q2", "O1",
                0.0, -0.13555519948726164,
                -0.067777599743630818, -0.044999999999999998);
        assertGoldenRow(T, "Q2", "O2",
                0.0, -0.054478448966511772,
                -0.016343534689953532, -0.012);
        // Open throughput is lambda*visits, fixed by the arrival rate; the service
        // rate cannot move it. MATLAB reports exact zeros, not near-zeros.
        for (int k = 0; k < T.getDTput().size(); k++) {
            assertEquals(0.0, T.getDTput().get(k).doubleValue(), 0.0,
                    "open dTput_dRate must be exactly zero");
        }
    }

    @Test
    @DisplayName("multiserver throws with MATLAB's message when 'exact' is asked for")
    public void testMultiserverThrows() {
        // MATLAB line_errors here; the JAR must not silently return null or an empty
        // table, which would let a caller read a wrong answer as a missing one. Under
        // the default 'auto' the model is out of the analytic scope and MATLAB falls
        // back to finite differences instead of erroring.
        final Network m = closedMultiserver();
        RuntimeException e = assertThrows(RuntimeException.class, new Executable() {
            public void execute() {
                new SolverMVA(m).getSensitivityTable("exact", Double.NaN, "forward");
            }
        });
        assertTrue(e.getMessage().contains("supports single-server stations only"),
                "unexpected message: " + e.getMessage());
        assertEquals("fd", new SolverMVA(m).getSensitivityTable().getMethod(),
                "an out-of-scope model must fall back to finite differences under 'auto'");
    }

    @Test
    @DisplayName("mixed open+closed throws with MATLAB's message when 'exact' is asked for")
    public void testMixedThrows() {
        final Network m = mixedOpenClosed();
        RuntimeException e = assertThrows(RuntimeException.class, new Executable() {
            public void execute() {
                new SolverMVA(m).getSensitivityTable("exact", Double.NaN, "forward");
            }
        });
        assertTrue(e.getMessage().contains("does not yet support mixed (open+closed) networks"),
                "unexpected message: " + e.getMessage());
        assertEquals("fd", new SolverMVA(m).getSensitivityTable().getMethod(),
                "an out-of-scope model must fall back to finite differences under 'auto'");
    }

    @Test
    @DisplayName("class switching is differentiated analytically, at chain level")
    public void testClassSwitchingUsesExact() {
        // Class switching makes the model chain-based: the population sits on the
        // reference class C1 and the switched class C2 carries njobs = 0. The analytic
        // branch aggregates the demands and the populations into chains, so it no
        // longer differentiates a model whose served classes are empty.
        NetworkSensitivityTable T = new SolverMVA(closedClassSwitching()).getSensitivityTable();
        assertEquals("exact", T.getMethod(),
                "a class-switching model is now in the analytic scope under 'auto'");
        assertNotNull(T.getSens(), "the closed analytic branch must expose its pfqn_sens result");
        assertTrue(T.getStationNames().size() > 0, "the exact table must carry rows");
        // The whole point of the fix: the derivatives of the served classes are not
        // silently zero any more.
        for (int k = 0; k < T.getStationNames().size(); k++) {
            assertTrue(isFinite(T.getDQLen().get(k)), "the exact table must carry finite values");
            assertTrue(Math.abs(T.getDTput().get(k).doubleValue()) > 1e-9,
                    "dTput_dRate must not be zero on a chain-aggregated model");
        }
        // Explicit 'exact' must be accepted rather than refused.
        NetworkSensitivityTable Te = new SolverMVA(closedClassSwitching())
                .getSensitivityTable("exact", Double.NaN, "forward");
        assertEquals("exact", Te.getMethod());
        // And the numerical branch must land on the same numbers.
        NetworkSensitivityTable fd = new SolverMVA(closedClassSwitching())
                .getSensitivityTable("fd", Double.NaN, "central");
        assertEquals(Te.getStationNames(), fd.getStationNames(), "station column");
        assertEquals(Te.getClassNames(), fd.getClassNames(), "job class column");
        for (int k = 0; k < Te.getStationNames().size(); k++) {
            String tag = Te.getStationNames().get(k) + "/" + Te.getClassNames().get(k);
            assertFdAgrees(Te.getDTput().get(k), fd.getDTput().get(k), "dTput_dRate " + tag);
            assertFdAgrees(Te.getDRespT().get(k), fd.getDRespT().get(k), "dRespT_dRate " + tag);
            assertFdAgrees(Te.getDQLen().get(k), fd.getDQLen().get(k), "dQLen_dRate " + tag);
            assertFdAgrees(Te.getDUtil().get(k), fd.getDUtil().get(k), "dUtil_dRate " + tag);
        }
    }

    /**
     * Closed, one class, with non-unit visit ratios: the delay routes to Q1, Q1 to Q2,
     * and Q2 back to the delay or to Q1 with equal probability, so v(Q1) = v(Q2) = 2
     * against v(Think) = 1. A branch that reports the chain throughput instead of the
     * class throughput is off by exactly that factor here.
     */
    private static Network closedVisitRatios() {
        Network m = new Network("mt_visits");
        Delay d = new Delay(m, "Think");
        Queue q1 = new Queue(m, "Q1", SchedStrategy.PS);
        Queue q2 = new Queue(m, "Q2", SchedStrategy.PS);
        ClosedClass c = new ClosedClass(m, "C1", 3, d, 0);
        d.setService(c, new Exp(1.0));
        q1.setService(c, new Exp(2.0));
        q2.setService(c, new Exp(3.0));
        RoutingMatrix P = m.initRoutingMatrix();
        P.set(c, c, d, q1, 1.0);
        P.set(c, c, q1, q2, 1.0);
        P.set(c, c, q2, d, 0.5);
        P.set(c, c, q2, q1, 0.5);
        m.link(P);
        return m;
    }

    @Test
    @DisplayName("non-unit visit ratios are carried by the analytic branch")
    public void testVisitRatiosMatchMatlab() {
        // Regression for the visit-ratio fix: the class throughput at a station is
        // X_c*v(i,r) and the response time is per visit, whereas the chain recursion
        // returns a per-chain throughput and a chain residence time. Goldens from
        // MATLAB on this model.
        NetworkSensitivityTable T = new SolverMVA(closedVisitRatios()).getSensitivityTable();
        assertEquals("exact", T.getMethod());
        assertEquals(2, T.getStationNames().size(), "two queues and one class give two rows");
        int k1 = rowOf(T, "Q1", "C1");
        assertEquals(0.451500, T.getDTput().get(k1).doubleValue(), 5e-6, "dTput_dRate Q1/C1");
        assertEquals(-0.611781, T.getDRespT().get(k1).doubleValue(), 5e-6, "dRespT_dRate Q1/C1");
        assertEquals(-0.529216, T.getDQLen().get(k1).doubleValue(), 5e-6, "dQLen_dRate Q1/C1");
        assertEquals(-0.161834, T.getDUtil().get(k1).doubleValue(), 5e-6, "dUtil_dRate Q1/C1");
        int k2 = rowOf(T, "Q2", "C1");
        assertEquals(0.137336, T.getDTput().get(k2).doubleValue(), 5e-6, "dTput_dRate Q2/C1");
        assertEquals(-0.219655, T.getDRespT().get(k2).doubleValue(), 5e-6, "dRespT_dRate Q2/C1");
        assertEquals(-0.270979, T.getDQLen().get(k2).doubleValue(), 5e-6, "dQLen_dRate Q2/C1");
        assertEquals(-0.126481, T.getDUtil().get(k2).doubleValue(), 5e-6, "dUtil_dRate Q2/C1");

        // The numerical branch must land on the same numbers: the fix is in the
        // composition back to the classes, not in a rescaling of the goldens.
        NetworkSensitivityTable fd = new SolverMVA(closedVisitRatios())
                .getSensitivityTable("fd", Double.NaN, "central");
        assertEquals(T.getStationNames(), fd.getStationNames(), "station column");
        for (int k = 0; k < T.getStationNames().size(); k++) {
            String tag = T.getStationNames().get(k) + "/" + T.getClassNames().get(k);
            assertFdAgrees(T.getDTput().get(k), fd.getDTput().get(k), "dTput_dRate " + tag);
            assertFdAgrees(T.getDRespT().get(k), fd.getDRespT().get(k), "dRespT_dRate " + tag);
            assertFdAgrees(T.getDQLen().get(k), fd.getDQLen().get(k), "dQLen_dRate " + tag);
            assertFdAgrees(T.getDUtil().get(k), fd.getDUtil().get(k), "dUtil_dRate " + tag);
        }
    }

    @Test
    @DisplayName("rows are skipped exactly where MATLAB skips them")
    public void testRowsSkippedAsInMatlab() {
        // C2 never visits Q2, so MATLAB emits three rows and not four. The absent row
        // must be absent, not present-and-zero.
        NetworkSensitivityTable T = new SolverMVA(closedWithSkippedRow()).getSensitivityTable();
        assertEquals(3, T.getStationNames().size(),
                "the (Q2, C2) pair has no demand and must not produce a row");
        assertEquals(-1, T.findRow("Q2", "C2"), "(Q2, C2) must be absent");
        assertGoldenRow(T, "Q1", "C1",
                0.1615495740410465, -0.49716867657865543,
                -0.22201875284504291, -0.065244855852990058);
        assertGoldenRow(T, "Q1", "C2",
                0.51893028048967538, -0.96104399901875692,
                -0.25946514024483763, -0.062752738383419548);
        assertGoldenRow(T, "Q2", "C1",
                0.034394902145189844, -0.11268320348120352,
                -0.079746786248187165, -0.062730414932985409);
        // The delay is not a Queue node, so it never contributes a row, exactly as in
        // MATLAB, which iterates over sn.nodetype == NodeType.Queue.
        assertEquals(-1, T.findRow("Think", "C1"), "the delay must not produce a row");
    }

    @Test
    @DisplayName("the sens second output is exposed and carries the full Jacobian")
    public void testSensSecondOutput() {
        // MATLAB returns [SensTable, sens]; Java has no multiple returns, so sens
        // hangs off the table. It must be the pfqn_sens result, not a summary.
        NetworkSensitivityTable T = new SolverMVA(closedTwoClass()).getSensitivityTable();
        Ret.pfqnSens sens = T.getSens();
        assertNotNull(sens, "sens must be present on a closed model");
        // Two queues and two classes give four L parameters; the two think-time
        // demands add two Z parameters, so the Jacobian is 2 x 6, as MATLAB reports.
        assertEquals(2, sens.dX.getNumRows(), "dX rows = classes");
        assertEquals(6, sens.dX.getNumCols(), "dX cols = L params + Z params");
        assertEquals(6, sens.dQ.length, "one dQ block per parameter");

        // The table is the diagonal of that Jacobian under the rate chain rule
        // d/d(rate) = -(D/rate) d/dL, with D(Q1,C1) = 0.4 and rate = 1/0.4.
        int R = 2;
        int ist = 0;   // Q1
        int r = 0;     // C1
        int p = ist * R + r;
        double chain = -0.4 / (1 / 0.4);
        int k = rowOf(T, "Q1", "C1");
        assertRelEquals(sens.dX.get(r, p) * chain, T.getDTput().get(k), "dTput is the dX diagonal");
        assertRelEquals(sens.dQ[p].get(ist, r) * chain, T.getDQLen().get(k), "dQLen is the dQ diagonal");
        assertRelEquals(sens.dU[p].get(ist, r) * chain, T.getDUtil().get(k), "dUtil is the dU diagonal");
        assertRelEquals(sens.dR[p].get(ist, r) * chain, T.getDRespT().get(k), "dRespT is the dR diagonal");
    }

    @Test
    @DisplayName("sens is empty on an open model, as in MATLAB")
    public void testSensEmptyOnOpen() {
        // The open branch is closed-form and runs no differentiated MVA; MATLAB
        // returns sens = [] there.
        NetworkSensitivityTable T = new SolverMVA(openTwoClass()).getSensitivityTable();
        assertNull(T.getSens(), "the open branch must not fabricate a pfqn_sens result");
    }

    @Test
    @DisplayName("the column names are MATLAB's, in MATLAB's order")
    public void testVariableNames() {
        NetworkSensitivityTable T = new SolverMVA(closedTwoClass()).getSensitivityTable();
        List<String> expected = Arrays.asList("Station", "JobClass", "dTput_dRate",
                "dRespT_dRate", "dQLen_dRate", "dUtil_dRate");
        assertEquals(expected, T.getVariableNames());
        // get(int) must agree with the named accessors, i.e. the column order is real.
        assertEquals(T.getDTput(), T.get(0));
        assertEquals(T.getDRespT(), T.get(1));
        assertEquals(T.getDQLen(), T.get(2));
        assertEquals(T.getDUtil(), T.get(3));
    }

    @Test
    @DisplayName("an unknown column is rejected rather than silently empty")
    public void testUnknownColumnThrows() {
        final NetworkSensitivityTable T = new SolverMVA(closedTwoClass()).getSensitivityTable();
        assertThrows(RuntimeException.class, new Executable() {
            public void execute() {
                T.getColumn("dNotAColumn");
            }
        });
    }

    // ---------- finite-difference branch ----------------------------------

    /**
     * Closed, one class, one delay and one PS queue. Small enough for SolverCTMC and
     * SolverFluid to run the 1+M*R (or 2*M*R) solves the finite-difference branch
     * needs, and in scope for the analytic branch, so that the two can be compared.
     */
    private static Network closedSingleClass() {
        Network m = new Network("mt_fd");
        Delay d = new Delay(m, "Think");
        Queue q = new Queue(m, "Q1", SchedStrategy.PS);
        ClosedClass c = new ClosedClass(m, "C1", 3, d, 0);
        d.setService(c, new Exp(1.0));
        q.setService(c, new Exp(2.0));
        RoutingMatrix P = m.initRoutingMatrix();
        P.set(c, c, Network.serialRouting(d, q));
        m.link(P);
        return m;
    }

    /** The same topology with Erlang-2 service, i.e. a non-exponential family. */
    private static Network closedSingleClassErlang() {
        Network m = new Network("mt_fd_erl");
        Delay d = new Delay(m, "Think");
        Queue q = new Queue(m, "Q1", SchedStrategy.PS);
        ClosedClass c = new ClosedClass(m, "C1", 2, d, 0);
        d.setService(c, new Exp(1.0));
        q.setService(c, new Erlang(4.0, 2));
        RoutingMatrix P = m.initRoutingMatrix();
        P.set(c, c, Network.serialRouting(d, q));
        m.link(P);
        return m;
    }

    @Test
    @DisplayName("the branch is chosen by the solver: exact on MVA and NC, fd elsewhere")
    public void testMethodDispatchPerSolver() {
        // supportsExactSensitivity() is the switch; it is a solver property, so a model
        // that is fully in the analytic scope still goes to finite differences under a
        // solver that does not evaluate a product-form recursion.
        assertEquals("exact", new SolverMVA(closedSingleClass()).getSensitivityTable().getMethod());
        assertEquals("exact", new SolverNC(closedSingleClass()).getSensitivityTable().getMethod());
        assertEquals("fd", new SolverCTMC(closedSingleClass()).getSensitivityTable().getMethod());
        assertEquals("fd", new SolverFluid(closedSingleClass()).getSensitivityTable().getMethod());
    }

    @Test
    @DisplayName("central differences reproduce the analytic branch to 4 significant digits")
    public void testFdAgreesWithExact() {
        Network m1 = closedSingleClass();
        NetworkSensitivityTable exact = new SolverMVA(m1).getSensitivityTable("exact", Double.NaN, "forward");
        Network m2 = closedSingleClass();
        NetworkSensitivityTable fd = new SolverMVA(m2).getSensitivityTable("fd", Double.NaN, "central");
        assertEquals("exact", exact.getMethod());
        assertEquals("fd", fd.getMethod());
        assertEquals(1, exact.getStationNames().size(), "one queue and one class give one row");
        assertEquals(exact.getStationNames(), fd.getStationNames());
        assertEquals(exact.getClassNames(), fd.getClassNames());

        int ke = rowOf(exact, "Q1", "C1");
        int kf = rowOf(fd, "Q1", "C1");
        // The analytic branch is the reference; the quotient must land on it to four
        // significant digits, i.e. a relative error below 1e-4.
        assertFdAgrees(exact.getDTput().get(ke), fd.getDTput().get(kf), "dTput_dRate");
        assertFdAgrees(exact.getDRespT().get(ke), fd.getDRespT().get(kf), "dRespT_dRate");
        assertFdAgrees(exact.getDQLen().get(ke), fd.getDQLen().get(kf), "dQLen_dRate");
        assertFdAgrees(exact.getDUtil().get(ke), fd.getDUtil().get(kf), "dUtil_dRate");

        // Absolute pinning of the analytic values, each to the digits quoted in the
        // reference: a common error in both branches would otherwise pass unnoticed.
        assertEquals(0.4903, exact.getDTput().get(ke).doubleValue(), 5e-5, "dTput_dRate");
        assertEquals(-0.59, exact.getDRespT().get(ke).doubleValue(), 5e-3, "dRespT_dRate");
        assertEquals(-0.4903, exact.getDQLen().get(ke).doubleValue(), 5e-5, "dQLen_dRate");
        assertEquals(-0.14958, exact.getDUtil().get(ke).doubleValue(), 5e-6, "dUtil_dRate");
    }

    /** Four significant digits of agreement between the two branches. */
    private static void assertFdAgrees(Double exact, Double fd, String col) {
        double e = exact.doubleValue();
        double denom = Math.max(Math.abs(e), 1e-300);
        assertTrue(Math.abs(fd.doubleValue() - e) <= 1e-4 * denom,
                col + ": finite differences gave " + fd + " against the analytic " + exact
                        + " (rel err " + (Math.abs(fd.doubleValue() - e) / denom) + ")");
    }

    @Test
    @DisplayName("'exact' is refused by a solver that has no product-form recursion")
    public void testExactRefusedByNonProductFormSolver() {
        final Network m = closedSingleClass();
        RuntimeException e = assertThrows(RuntimeException.class, new Executable() {
            public void execute() {
                new SolverCTMC(m).getSensitivityTable("exact", Double.NaN, "forward");
            }
        });
        assertTrue(e.getMessage().contains("SolverMVA and SolverNC only"),
                "unexpected message: " + e.getMessage());
        // The model itself is in scope, so the refusal is about the solver, not the
        // model: the same solver must still deliver a table by finite differences.
        assertEquals("fd", new SolverCTMC(closedSingleClass()).getSensitivityTable().getMethod());
    }

    @Test
    @DisplayName("'auto' falls back to fd exactly where 'exact' is refused for the model")
    public void testAutoFallsBackWhereExactIsRefused() {
        // The behaviour change introduced with the fd branch: an out-of-scope model no
        // longer errors under the default, it is differentiated numerically instead.
        final Network multi = closedMultiserver();
        assertThrows(RuntimeException.class, new Executable() {
            public void execute() {
                new SolverMVA(multi).getSensitivityTable("exact", Double.NaN, "forward");
            }
        });
        NetworkSensitivityTable Tmulti = new SolverMVA(closedMultiserver()).getSensitivityTable();
        assertEquals("fd", Tmulti.getMethod(), "a multiserver model must now yield an fd table");
        assertEquals(1, Tmulti.getStationNames().size(), "one queue and one class give one row");
        assertTrue(isFinite(Tmulti.getDQLen().get(0)), "the fd table must carry finite values");
        assertNull(Tmulti.getSens(), "the fd branch runs no differentiated MVA");

        final Network mixed = mixedOpenClosed();
        assertThrows(RuntimeException.class, new Executable() {
            public void execute() {
                new SolverMVA(mixed).getSensitivityTable("exact", Double.NaN, "forward");
            }
        });
        assertEquals("fd", new SolverMVA(mixedOpenClosed()).getSensitivityTable().getMethod(),
                "a mixed open+closed model must now yield an fd table");
    }

    @Test
    @DisplayName("invalid method, scheme and step are all rejected")
    public void testOptionValidation() {
        final Network m = closedSingleClass();
        RuntimeException eMethod = assertThrows(RuntimeException.class, new Executable() {
            public void execute() {
                new SolverMVA(m).getSensitivityTable("approximate", Double.NaN, "forward");
            }
        });
        assertTrue(eMethod.getMessage().contains("'auto', 'exact', 'fd'"),
                "unexpected message: " + eMethod.getMessage());
        RuntimeException eScheme = assertThrows(RuntimeException.class, new Executable() {
            public void execute() {
                new SolverMVA(m).getSensitivityTable("fd", Double.NaN, "backward");
            }
        });
        assertTrue(eScheme.getMessage().contains("'forward' or 'central'"),
                "unexpected message: " + eScheme.getMessage());
        // The step is only read by the fd branch, so it is validated there.
        RuntimeException eZero = assertThrows(RuntimeException.class, new Executable() {
            public void execute() {
                new SolverMVA(m).getSensitivityTable("fd", 0.0, "forward");
            }
        });
        assertTrue(eZero.getMessage().contains("scalar in (0,1)"),
                "unexpected message: " + eZero.getMessage());
        assertThrows(RuntimeException.class, new Executable() {
            public void execute() {
                new SolverMVA(m).getSensitivityTable("fd", -1e-4, "forward");
            }
        });
        assertThrows(RuntimeException.class, new Executable() {
            public void execute() {
                new SolverMVA(m).getSensitivityTable("fd", 1.0, "forward");
            }
        });
        assertThrows(RuntimeException.class, new Executable() {
            public void execute() {
                new SolverMVA(m).getSensitivityTable("fd", 2.5, "forward");
            }
        });
    }

    @Test
    @DisplayName("the fd branch emits exactly the rows the exact branch emits")
    public void testFdRowSelectionMatchesExact() {
        // The two branches build the mask differently, the exact branch from the
        // product-form demands and the fd branch from the visit ratios; on a model
        // where a class does not visit a queue they must still agree row for row.
        NetworkSensitivityTable exact = new SolverMVA(closedWithSkippedRow())
                .getSensitivityTable("exact", Double.NaN, "forward");
        NetworkSensitivityTable fd = new SolverMVA(closedWithSkippedRow())
                .getSensitivityTable("fd", Double.NaN, "central");
        assertEquals(exact.getStationNames(), fd.getStationNames(), "station column");
        assertEquals(exact.getClassNames(), fd.getClassNames(), "job class column");
        assertEquals(3, fd.getStationNames().size(), "the (Q2, C2) pair must not produce a row");
        assertEquals(-1, fd.findRow("Q2", "C2"), "(Q2, C2) must be absent from the fd table");
        // Agreement on the rows that are emitted, to the same four significant digits.
        String[] stations = {"Q1", "Q1", "Q2"};
        String[] jobclasses = {"C1", "C2", "C1"};
        for (int i = 0; i < stations.length; i++) {
            int ke = rowOf(exact, stations[i], jobclasses[i]);
            int kf = rowOf(fd, stations[i], jobclasses[i]);
            String tag = stations[i] + "/" + jobclasses[i];
            assertFdAgrees(exact.getDTput().get(ke), fd.getDTput().get(kf), "dTput_dRate " + tag);
            assertFdAgrees(exact.getDRespT().get(ke), fd.getDRespT().get(kf), "dRespT_dRate " + tag);
            assertFdAgrees(exact.getDQLen().get(ke), fd.getDQLen().get(kf), "dQLen_dRate " + tag);
            assertFdAgrees(exact.getDUtil().get(ke), fd.getDUtil().get(kf), "dUtil_dRate " + tag);
        }
    }

    @Test
    @DisplayName("the fd branch runs end to end on a non-exponential service process")
    public void testFdWithErlangService() {
        // Smoke coverage of the rate-scaling primitive on a family other than Exp: the
        // perturbed model must still be solvable and the quotient finite.
        NetworkSensitivityTable T = new SolverCTMC(closedSingleClassErlang()).getSensitivityTable();
        assertEquals("fd", T.getMethod());
        assertEquals(1, T.getStationNames().size(), "one queue and one class give one row");
        assertTrue(isFinite(T.getDTput().get(0)), "dTput_dRate must be finite");
        assertTrue(isFinite(T.getDRespT().get(0)), "dRespT_dRate must be finite");
        assertTrue(isFinite(T.getDQLen().get(0)), "dQLen_dRate must be finite");
        assertTrue(isFinite(T.getDUtil().get(0)), "dUtil_dRate must be finite");
        // Faster service can only raise throughput and lower the queue length.
        assertTrue(T.getDTput().get(0).doubleValue() > 0, "dTput_dRate must be positive");
        assertTrue(T.getDQLen().get(0).doubleValue() < 0, "dQLen_dRate must be negative");
    }

    private static boolean isFinite(Double v) {
        return v != null && !v.isNaN() && !v.isInfinite();
    }
}
