/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers;

import jline.api.pfqn.PfqnResptPsResult;
import jline.api.pfqn.sens.Pfqn_sens_mom;
import jline.api.qsys.QsysMm1PsResult;
import jline.api.sn.SnHasProductForm;
import jline.api.pfqn.sens.Pfqn_sens_mva;
import jline.io.Ret;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.solvers.mva.SolverMVA;
import jline.solvers.nc.SolverNC;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.function.Executable;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static jline.api.pfqn.Pfqn_respt_ps_moments.pfqn_respt_ps_moments;
import static jline.api.qsys.Qsys_mm1_ps.qsys_mm1_ps;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Regression tests for {@link NetworkSolver#getMomentTable(int[])} and
 * {@link NetworkSolver#getMomentStationTable(int[])}, the solver-level views
 * over the {@code Pfqn_sens_*} moment family. Mirrors the MATLAB harness
 * {@code test_moment_table.m} check for check.
 *
 * <p>The underlying algorithms are validated against brute-force enumeration, the
 * published tables of Strelen (1990) and simulation by the api-level harnesses
 * ({@link jline.api.PfqnSensMvaTest}, {@code PfqnSensMomTest},
 * {@code PfqnSensResptTest}, {@code PfqnSensLinearizerTest}). These tests therefore
 * do NOT re-check the mathematics. They check the solver-level plumbing, which is
 * where a table method can go wrong independently of correct algorithms:</p>
 *
 * <ul>
 *   <li>that the demands, visit ratios and service times handed to the api layer are
 *       the ones the model actually expresses (asserted by requiring the table's
 *       means to reproduce {@link NetworkSolver#getAvgTable()}, which is computed by
 *       a different code path);</li>
 *   <li>that the per-class and per-station-total tables are mutually consistent,
 *       Var[Q_i] = sum over class pairs of the per-class covariances;</li>
 *   <li>that RespTVar is populated at exactly the FCFS stations and NaN elsewhere,
 *       since the sojourn-time distribution is unknown at PS/LCFS centers and a
 *       number there would be wrong rather than missing;</li>
 *   <li>that the dispatch by model type (closed / mixed / purely open, single-server
 *       / multiserver) picks a path that runs and agrees.</li>
 * </ul>
 */
public class MomentTableTest {

    // ---------- helpers ---------------------------------------------------

    /**
     * Closed, two classes. Q1 is FCFS with class-independent rates (as BCMP requires
     * of FCFS); Q2 is PS with class-dependent rates, which is legal and which the
     * FCFS gate must tolerate.
     */
    private static Network closedMixedSched() {
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

    private static int rowOf(NetworkMomentTable T, String station, String jobclass) {
        int k = T.findRow(station, jobclass);
        assertTrue(k >= 0, "row " + station + "/" + jobclass + " missing from the moment table");
        return k;
    }

    /** Index of the (station, class) row of an avg table, or -1. */
    private static int avgRowOf(NetworkAvgTable A, String station, String jobclass) {
        List<String> st = A.getStationNames();
        List<String> cl = A.getClassNames();
        for (int j = 0; j < st.size(); j++) {
            if (st.get(j).equals(station) && cl.get(j).equals(jobclass)) {
                return j;
            }
        }
        return -1;
    }

    private static void assertRelEquals(double expected, double actual, double relTol, String msg) {
        double denom = Math.max(Math.abs(expected), 1e-300);
        assertTrue(Math.abs(actual - expected) <= relTol * denom,
                msg + ": expected " + expected + " but got " + actual
                        + " (rel err " + (Math.abs(actual - expected) / denom) + ")");
    }

    // ---------- tests -----------------------------------------------------

    @Test
    @DisplayName("means reproduce getAvgTable")
    public void testMeansReproduceGetAvgTable() {
        // The table is only meaningful if the parameters handed to the api layer are
        // the model's. getAvgTable reaches the same means by a different path, so
        // agreement pins the translation.
        Network m = closedMixedSched();
        SolverMVA s = new SolverMVA(m);
        NetworkMomentTable T = s.getMomentTable();
        NetworkAvgTable A = s.getAvgTable();
        assertTrue(T.getStationNames().size() > 0);
        for (int k = 0; k < T.getStationNames().size(); k++) {
            String st = T.getStationNames().get(k);
            String cl = T.getClassNames().get(k);
            int j = avgRowOf(A, st, cl);
            assertTrue(j >= 0, "row " + st + "/" + cl + " missing from getAvgTable");
            assertRelEquals(A.getQLen().get(j), T.getQLen().get(k), 1e-9, "QLen " + st + "/" + cl);
            assertRelEquals(A.getRespT().get(j), T.getRespT().get(k), 1e-9, "RespT " + st + "/" + cl);
        }
    }

    @Test
    @DisplayName("per-class variance is exactly what pfqn_sens_mva reports")
    public void testPerclassVarianceMatchesApi() {
        // The QLenVar column must be exactly what Pfqn_sens_mva reports; nothing in
        // the table method may rescale it.
        Network m = closedMixedSched();
        SolverMVA s = new SolverMVA(m);
        NetworkMomentTable T = s.getMomentTable();
        NetworkMomentResult mom = T.getMoments();

        Matrix L = new Matrix(2, 2);
        L.set(0, 0, 0.4);
        L.set(0, 1, 0.4);
        L.set(1, 0, 0.3);
        L.set(1, 1, 0.2);
        Matrix N = new Matrix(1, 2);
        N.set(0, 0, 3);
        N.set(0, 1, 2);
        Matrix Z = new Matrix(1, 2);
        Z.set(0, 0, 1.0);
        Z.set(0, 1, 0.5);
        Ret.pfqnSensMva ref = Pfqn_sens_mva.pfqn_sens_mva(L, N, Z);

        for (int i = 0; i < 2; i++) {
            for (int r = 0; r < 2; r++) {
                assertRelEquals(ref.QVar.get(i, r), mom.QVar.get(i, r), 1e-9, "QVar(" + i + "," + r + ")");
            }
        }
        int k = rowOf(T, "Q1", "C1");
        assertRelEquals(ref.QVar.get(0, 0), T.getQLenVar().get(k), 1e-9, "Q1/C1 QLenVar");
        assertRelEquals(ref.QVar.get(0, 0) / (ref.Q.get(0, 0) * ref.Q.get(0, 0)),
                T.getQLenSCV().get(k), 1e-9, "Q1/C1 QLenSCV");
    }

    @Test
    @DisplayName("station-total Var equals the sum of the per-class covariance block")
    public void testStationTotalConsistentWithPerclass() {
        // Var[Q_i] of the station table must equal the sum of the per-class
        // covariance block of the per-class table. The two come from different
        // recursions (Strelen's column scaling vs de Souza e Silva-Muntz's per-class
        // one), so this is a real cross-check, not a tautology.
        Network m = closedMixedSched();
        SolverMVA s = new SolverMVA(m);
        NetworkMomentResult mom = s.getMomentTable().getMoments();
        NetworkMomentStationTable S = s.getMomentStationTable();
        for (int ist = 0; ist < 2; ist++) {
            double tot = 0.0;
            for (int r = 0; r < 2; r++) {
                for (int t = 0; t < 2; t++) {
                    tot += mom.QCov[ist].get(r, t);
                }
            }
            assertRelEquals(tot, S.getQLenVar().get(ist), 1e-8, "station " + ist + " total Var");
            double qtot = 0.0;
            for (int r = 0; r < 2; r++) {
                qtot += mom.Q.get(ist, r);
            }
            assertRelEquals(qtot, S.getQLen().get(ist), 1e-9, "station " + ist + " total QLen");
        }
    }

    @Test
    @DisplayName("RespTVar only at FCFS stations")
    public void testResptVarianceOnlyAtFcfs() {
        // RespTVar must be a number at the FCFS station and NaN at the PS station.
        // The PS moments of Mitra and Morrison cover the terminal-driven system, one
        // PS station and delays; this model has a second queueing station, which is
        // outside that scope, so a value at Q2 would be fabricated.
        Network m = closedMixedSched();
        NetworkMomentTable T = new SolverMVA(m).getMomentTable();
        String[] classes = new String[]{"C1", "C2"};
        for (int c = 0; c < classes.length; c++) {
            int kf = rowOf(T, "Q1", classes[c]);
            assertFalse(Double.isNaN(T.getRespTVar().get(kf)), "FCFS station must report RespTVar");
            assertTrue(T.getRespTVar().get(kf) > 0, "FCFS RespTVar must be positive");
            assertFalse(Double.isNaN(T.getRespTSCV().get(kf)), "FCFS station must report RespTSCV");
            int kp = rowOf(T, "Q2", classes[c]);
            assertTrue(Double.isNaN(T.getRespTVar().get(kp)), "PS station must NOT report RespTVar");
            assertTrue(Double.isNaN(T.getRespTSCV().get(kp)), "PS station must NOT report RespTSCV");
        }
    }

    @Test
    @DisplayName("open PS RespTVar is the Mitra-Morrison moment")
    public void testOpenPsResptVariance() {
        // A single PS station fed by Poisson streams is the open system of Mitra and
        // Morrison, where the sojourn-time moments are exact in closed form.
        double l1 = 0.3;
        double l2 = 0.4;
        double mu1 = 1.0;
        double mu2 = 3.0;
        Network m = new Network("mt_open_ps");
        Source src = new Source(m, "S");
        Sink snk = new Sink(m, "K");
        Queue q = new Queue(m, "CPU", SchedStrategy.PS);
        OpenClass o1 = new OpenClass(m, "C1", 0);
        OpenClass o2 = new OpenClass(m, "C2", 0);
        src.setArrival(o1, new Exp(l1));
        src.setArrival(o2, new Exp(l2));
        q.setService(o1, new Exp(mu1));
        q.setService(o2, new Exp(mu2));
        RoutingMatrix P = m.initRoutingMatrix();
        P.set(o1, o1, Network.serialRouting(src, q, snk));
        P.set(o2, o2, Network.serialRouting(src, q, snk));
        m.link(P);
        NetworkMomentTable T = new SolverMVA(m).getMomentTable();
        QsysMm1PsResult ref = qsys_mm1_ps(new double[]{l1, l2}, new double[]{mu1, mu2});
        String[] classes = new String[]{"C1", "C2"};
        for (int c = 0; c < classes.length; c++) {
            int k = rowOf(T, "CPU", classes[c]);
            assertRelEquals(ref.W[c], T.getRespT().get(k), 1e-9, "open PS mean RespT");
            assertRelEquals(ref.W2[c] - ref.W[c] * ref.W[c], T.getRespTVar().get(k),
                    1e-9, "open PS RespTVar");
        }
    }

    @Test
    @DisplayName("a PS station on a routing cycle reports no RespTVar")
    public void testOpenPsWithFeedbackIsBlank() {
        // With feedback the arrival stream at the station is no longer Poisson, so
        // the open formula does not apply and a blank is the honest answer.
        Network m = new Network("mt_open_ps_fb");
        Source src = new Source(m, "S");
        Sink snk = new Sink(m, "K");
        Queue q = new Queue(m, "CPU", SchedStrategy.PS);
        OpenClass o1 = new OpenClass(m, "C1", 0);
        src.setArrival(o1, new Exp(0.3));
        q.setService(o1, new Exp(1.0));
        RoutingMatrix P = m.initRoutingMatrix();
        P.set(o1, o1, src, q, 1.0);
        P.set(o1, o1, q, q, 0.5);
        P.set(o1, o1, q, snk, 0.5);
        m.link(P);
        NetworkMomentTable T = new SolverMVA(m).getMomentTable();
        int k = rowOf(T, "CPU", "C1");
        assertTrue(Double.isNaN(T.getRespTVar().get(k)),
                "a PS station with feedback must NOT report RespTVar");
    }

    @Test
    @DisplayName("closed terminal-driven PS RespTVar is the Mitra-Morrison moment")
    public void testClosedTerminalDrivenPsResptVariance() {
        // Terminals in series with one PS CPU: the closed system the paper analyses.
        Network m = new Network("mt_closed_ps");
        Delay d = new Delay(m, "Terminals");
        Queue cpu = new Queue(m, "CPU", SchedStrategy.PS);
        ClosedClass c1 = new ClosedClass(m, "C1", 4, d, 0);
        ClosedClass c2 = new ClosedClass(m, "C2", 2, d, 0);
        d.setService(c1, new Exp(0.02));
        d.setService(c2, new Exp(0.01));
        cpu.setService(c1, new Exp(1.0));
        cpu.setService(c2, new Exp(2.0));
        RoutingMatrix P = m.initRoutingMatrix();
        P.set(c1, c1, Network.serialRouting(d, cpu));
        P.set(c2, c2, Network.serialRouting(d, cpu));
        m.link(P);
        NetworkMomentTable T = new SolverMVA(m).getMomentTable();
        PfqnResptPsResult ref = pfqn_respt_ps_moments(new double[]{1.0, 0.5},
                new double[]{4, 2}, new double[]{50.0, 100.0});
        String[] classes = new String[]{"C1", "C2"};
        for (int c = 0; c < classes.length; c++) {
            int k = rowOf(T, "CPU", classes[c]);
            assertRelEquals(ref.W[c], T.getRespT().get(k), 1e-7,
                    "closed PS mean RespT must agree with Little's law");
            assertRelEquals(ref.W2[c] - ref.W[c] * ref.W[c], T.getRespTVar().get(k),
                    1e-9, "closed PS RespTVar");
        }
    }

    @Test
    @DisplayName("FCFS mean RespT from Theorem 4.1 matches MVA")
    public void testResptMeanFromTheorem41MatchesMva() {
        // At the FCFS station the mean RespT is overwritten by the t=1 case of
        // Strelen Theorem 4.1, which is a different expression from Little's law.
        // They must agree, which also validates the service-time/visit-ratio
        // factorization the table hands to Pfqn_sens_respt.
        Network m = closedMixedSched();
        SolverMVA s = new SolverMVA(m);
        NetworkMomentTable T = s.getMomentTable();
        NetworkAvgTable A = s.getAvgTable();
        String[] classes = new String[]{"C1", "C2"};
        for (int c = 0; c < classes.length; c++) {
            int k = rowOf(T, "Q1", classes[c]);
            int j = avgRowOf(A, "Q1", classes[c]);
            assertTrue(j >= 0);
            assertRelEquals(A.getRespT().get(j), T.getRespT().get(k), 1e-8, "Q1/" + classes[c] + " RespT");
        }
    }

    @Test
    @DisplayName("purely open M/M/1 takes the exact BCMP closed form")
    public void testOpenMm1ClosedForm() {
        // A purely open M/M/1 has no population lattice, so the table takes the exact
        // BCMP branch. Every moment is known in closed form: n ~ Geom(rho).
        double lambda = 0.5;
        double svc = 1.0;
        double rho = lambda * svc;
        Network m = new Network("mt_open");
        Source src = new Source(m, "S");
        Sink snk = new Sink(m, "K");
        Queue q = new Queue(m, "O1", SchedStrategy.FCFS);
        OpenClass oc = new OpenClass(m, "O", 0);
        src.setArrival(oc, new Exp(lambda));
        q.setService(oc, new Exp(1 / svc));
        m.link(Network.serialRouting(src, q, snk));
        NetworkMomentTable T = new SolverMVA(m).getMomentTable();
        int k = rowOf(T, "O1", "O");
        assertRelEquals(rho / (1 - rho), T.getQLen().get(k), 1e-10, "open QLen");
        assertRelEquals(rho / ((1 - rho) * (1 - rho)), T.getQLenVar().get(k), 1e-10, "open QLenVar");
        assertRelEquals(svc / (1 - rho), T.getRespT().get(k), 1e-10, "open RespT");
    }

    @Test
    @DisplayName("purely open multiclass M/M/1 takes the multinomial closed form")
    public void testOpenMulticlassMm1ClosedForm() {
        // Two open classes at one station: the total is geometric in the aggregate
        // utilization and the classes are multinomial given the total, so
        // Var[n_r] = E[n](p_r - p_r^2) + p_r^2 Var[n].
        double l1 = 0.3;
        double l2 = 0.2;
        double s1 = 1.0;
        double s2 = 1.5;
        Network m = new Network("mt_open2");
        Source src = new Source(m, "S");
        Sink snk = new Sink(m, "K");
        Queue q = new Queue(m, "O1", SchedStrategy.PS);
        OpenClass o1 = new OpenClass(m, "O1c", 0);
        OpenClass o2 = new OpenClass(m, "O2c", 0);
        src.setArrival(o1, new Exp(l1));
        src.setArrival(o2, new Exp(l2));
        q.setService(o1, new Exp(1 / s1));
        q.setService(o2, new Exp(1 / s2));
        RoutingMatrix P = m.initRoutingMatrix();
        P.set(o1, o1, Network.serialRouting(src, q, snk));
        P.set(o2, o2, Network.serialRouting(src, q, snk));
        m.link(P);
        NetworkMomentTable T = new SolverMVA(m).getMomentTable();
        double rho1 = l1 * s1;
        double rho2 = l2 * s2;
        double rho = rho1 + rho2;
        double En = rho / (1 - rho);
        double Vn = rho / ((1 - rho) * (1 - rho));
        double p1 = rho1 / rho;
        int k = rowOf(T, "O1", "O1c");
        assertRelEquals(p1 * En, T.getQLen().get(k), 1e-9, "open multiclass QLen");
        assertRelEquals(En * (p1 - p1 * p1) + p1 * p1 * Vn, T.getQLenVar().get(k), 1e-9,
                "open multiclass QLenVar");
    }

    @Test
    @DisplayName("Linearizer tracks the exact station moments within the published bands")
    public void testLinearizerTracksExact() {
        // The Linearizer is an approximation, so it is banded, not equated. The bands
        // are the accuracy the reference claims.
        Network m = closedMixedSched();
        SolverMVA s = new SolverMVA(m);
        // The algorithm is a solver property, not a table argument: the Linearizer is
        // selected by constructing the solver with method 'lin'.
        NetworkMomentStationTable E = s.getMomentStationTable(3);
        NetworkMomentStationTable L =
                new SolverMVA(m, "method", "lin").getMomentStationTable(3);
        assertEquals(E.getStationNames().size(), L.getStationNames().size());
        for (int k = 0; k < E.getStationNames().size(); k++) {
            assertRelEquals(E.getQLen().get(k), L.getQLen().get(k), 0.021, "linearizer QLen " + k);
            assertRelEquals(E.getQLenM3().get(k), L.getQLenM3().get(k), 0.062, "linearizer QLenM3 " + k);
        }
    }

    @Test
    @DisplayName("closed multiserver dispatches through pfqn_sens_mvaldmx")
    public void testMultiserverClosedDispatches() {
        // A multiserver closed model must route through Pfqn_sens_mvaldmx rather than
        // Pfqn_sens_mva, and still reproduce getAvgTable's means.
        Network m = new Network("mt_ms");
        Delay d = new Delay(m, "Think");
        Queue q1 = new Queue(m, "Q1", SchedStrategy.PS);
        q1.setNumberOfServers(2);
        Queue q2 = new Queue(m, "Q2", SchedStrategy.PS);
        ClosedClass c1 = new ClosedClass(m, "C1", 4, d, 0);
        d.setService(c1, new Exp(1 / 1.0));
        q1.setService(c1, new Exp(1 / 0.4));
        q2.setService(c1, new Exp(1 / 0.3));
        m.link(Network.serialRouting(d, q1, q2));
        SolverMVA s = new SolverMVA(m);
        NetworkMomentTable T = s.getMomentTable();
        NetworkAvgTable A = s.getAvgTable();
        // the multiserver branch must have been taken
        assertTrue(T.getMoments().qlenLdmx != null, "multiserver closed model must use pfqn_sens_mvaldmx");
        assertTrue(T.getMoments().qlenMva == null);
        assertTrue(T.getStationNames().size() > 0);
        for (int k = 0; k < T.getStationNames().size(); k++) {
            int j = avgRowOf(A, T.getStationNames().get(k), T.getClassNames().get(k));
            if (j >= 0) {
                assertRelEquals(A.getQLen().get(j), T.getQLen().get(k), 1e-8,
                        "multiserver QLen " + T.getStationNames().get(k));
            }
            assertTrue(T.getQLenVar().get(k) > 0, "multiserver QLenVar must be positive");
        }
    }

    @Test
    @DisplayName("order selects columns; a scalar k means 1..k, a set is literal")
    public void testOrderSelectsColumns() {
        // ORDER is a set of moment orders. A scalar k means 1:k; a vector is literal.
        Network m = closedMixedSched();
        SolverMVA s = new SolverMVA(m);

        NetworkMomentTable T1 = s.getMomentTable(1);
        assertEquals(Arrays.asList("Station", "JobClass", "QLen", "RespT"), T1.getVariableNames());

        NetworkMomentTable T2 = s.getMomentTable();                    // default is 2
        assertEquals(Arrays.asList("Station", "JobClass", "QLen", "QLenVar", "QLenSCV",
                "RespT", "RespTVar", "RespTSCV"), T2.getVariableNames());
        assertEquals(Arrays.asList("Station", "JobClass", "QLen", "QLenVar", "QLenSCV",
                "QLenSkew", "RespT", "RespTVar", "RespTSCV", "RespTSkew"),
                s.getMomentTable(3).getVariableNames(), "QLenSkew sits after QLenSCV");
        assertEquals(T2.getVariableNames(), s.getMomentTable(new int[]{1, 2}).getVariableNames(),
                "scalar 2 must equal the set [1 2]");
        // MATLAB's isscalar([2]) is true, so a one-element set takes the scalar path
        // and means 1:2 there; the JAR must agree or the same argument would mean
        // different things in different codebases.
        assertEquals(T2.getVariableNames(), s.getMomentTable(new int[]{2}).getVariableNames(),
                "the one-element set [2] must take the scalar path, i.e. mean 1:2");

        NetworkMomentTable T3 = s.getMomentTable(3);
        assertTrue(T3.hasColumn("RespTSkew"));
        // order 3 DOES add a per-class queue-length skewness: scaling L(i,r) alone is
        // Akyildiz-Strelen Theorem 1 with T={r}, so the quantity exists
        assertTrue(T3.hasColumn("QLenSkew"));

        NetworkMomentTable T23 = s.getMomentTable(new int[]{2, 3});    // no means
        assertFalse(T23.hasColumn("QLen"));
        assertTrue(T23.hasColumn("QLenVar"));
        assertTrue(T23.hasColumn("RespTSkew"));

        NetworkMomentStationTable S1 = s.getMomentStationTable(1);
        assertEquals(Arrays.asList("Station", "QLen"), S1.getVariableNames());
        NetworkMomentStationTable S2 = s.getMomentStationTable();
        assertEquals(Arrays.asList("Station", "QLen", "QLenVar", "QLenSCV"), S2.getVariableNames());
        NetworkMomentStationTable S3 = s.getMomentStationTable(3);
        assertTrue(S3.hasColumn("QLenM3"));
        assertTrue(S3.hasColumn("QLenSkew"));
        NetworkMomentStationTable S13 = s.getMomentStationTable(new int[]{1, 3});
        assertEquals(Arrays.asList("Station", "QLen", "QLenM3", "QLenSkew"), S13.getVariableNames());

        // the values must not depend on which columns were asked for
        for (int k = 0; k < S2.getStationNames().size(); k++) {
            assertRelEquals(S2.getQLen().get(k), S3.getQLen().get(k), 1e-12, "S3 vs S2 QLen " + k);
        }
        for (int k = 0; k < T2.getStationNames().size(); k++) {
            assertRelEquals(T2.getQLen().get(k), T3.getQLen().get(k), 1e-12, "T3 vs T2 QLen " + k);
            double v2 = T2.getRespTVar().get(k);
            double v3 = T3.getRespTVar().get(k);
            if (Double.isNaN(v2)) {
                assertTrue(Double.isNaN(v3), "RespTVar NaN-ness must not depend on order " + k);
            } else {
                assertRelEquals(v2, v3, 1e-12, "T3 vs T2 RespTVar " + k);
            }
        }
    }

    @Test
    @DisplayName("order rejects values outside 1..3")
    public void testOrderRejectsBadValues() {
        Network m = closedMixedSched();
        final SolverMVA s = new SolverMVA(m);
        // The MATLAB non-integer case (order 2.5) is rejected by the Java type
        // system: the order parameters are int / int[], so no such call compiles.
        assertThrows(RuntimeException.class, new Executable() {
            public void execute() {
                s.getMomentTable(0);
            }
        });
        assertThrows(RuntimeException.class, new Executable() {
            public void execute() {
                s.getMomentTable(4);
            }
        });
        assertThrows(RuntimeException.class, new Executable() {
            public void execute() {
                s.getMomentTable(new int[]{1, 5});
            }
        });
        assertThrows(RuntimeException.class, new Executable() {
            public void execute() {
                s.getMomentStationTable(new int[]{1, 5});
            }
        });
        assertThrows(RuntimeException.class, new Executable() {
            public void execute() {
                s.getMomentStationTable(0);
            }
        });
    }

    @Test
    @DisplayName("chain variance equals the summed per-class covariance block of that chain")
    public void testChainTableConsistentWithPerclass() {
        // Each chain's variance must equal the sum of the per-class covariance block
        // over the classes of that chain, since Var[sum_{r in c} n_ir] = sum_{r,s in c}
        // Cov. The two come from different groupings of the same recursion.
        Network m = closedMixedSched();
        SolverMVA s = new SolverMVA(m);
        NetworkMomentResult mom = s.getMomentTable().getMoments();
        NetworkStruct sn = m.getStruct();
        NetworkMomentChainTable Ch = s.getMomentChainTable(3);
        String[] stations = new String[]{"Q1", "Q2"};
        // both classes are closed and each forms its own chain in this model
        for (int k = 0; k < Ch.getStationNames().size(); k++) {
            int ist = -1;
            for (int i = 0; i < stations.length; i++) {
                if (stations[i].equals(Ch.getStationNames().get(k))) {
                    ist = i;
                }
            }
            assertTrue(ist >= 0);
            int c = Integer.parseInt(Ch.getChainNames().get(k).substring("Chain".length())) - 1;
            List<Integer> cls = new ArrayList<Integer>();
            for (int r = 0; r < sn.nclasses; r++) {
                if (sn.chains.get(c, r) != 0.0) {
                    cls.add(Integer.valueOf(r));
                }
            }
            assertTrue(cls.size() > 0);
            double tot = 0.0;
            double qtot = 0.0;
            for (int a = 0; a < cls.size(); a++) {
                qtot += mom.Q.get(ist, cls.get(a).intValue());
                for (int b = 0; b < cls.size(); b++) {
                    tot += mom.QCov[ist].get(cls.get(a).intValue(), cls.get(b).intValue());
                }
            }
            assertRelEquals(tot, Ch.getQLenVar().get(k), 1e-8, "chain Var " + k);
            assertRelEquals(qtot, Ch.getQLen().get(k), 1e-9, "chain QLen " + k);
        }
        assertTrue(Ch.hasColumn("QLenM3"));
        assertEquals(Arrays.asList("Station", "Chain", "QLen"),
                s.getMomentChainTable(1).getVariableNames());
    }

    @Test
    @DisplayName("a single chain makes the chain table equal the station table")
    public void testSingleChainEqualsStationTotal() {
        // When every class sits in one chain, the chain grouping IS the station total,
        // so getMomentChainTable must reproduce getMomentStationTable exactly.
        Network m = new Network("mt_1chain");
        Delay d = new Delay(m, "Think");
        Queue q1 = new Queue(m, "Q1", SchedStrategy.PS);
        Queue q2 = new Queue(m, "Q2", SchedStrategy.PS);
        ClosedClass c1 = new ClosedClass(m, "C1", 3, d, 0);
        ClosedClass c2 = new ClosedClass(m, "C2", 2, d, 0);
        d.setService(c1, new Exp(1 / 1.0));
        d.setService(c2, new Exp(1 / 0.5));
        q1.setService(c1, new Exp(1 / 0.4));
        q1.setService(c2, new Exp(1 / 0.6));
        q2.setService(c1, new Exp(1 / 0.3));
        q2.setService(c2, new Exp(1 / 0.2));
        // a class switch on the way back to the delay merges the two classes into a
        // single chain
        RoutingMatrix P = m.initRoutingMatrix();
        P.set(c1, c1, d, q1, 1.0);
        P.set(c1, c1, q1, q2, 1.0);
        P.set(c1, c2, q2, d, 1.0);
        P.set(c2, c2, d, q1, 1.0);
        P.set(c2, c2, q1, q2, 1.0);
        P.set(c2, c1, q2, d, 1.0);
        m.link(P);
        SolverMVA s = new SolverMVA(m);
        NetworkStruct sn = m.getStruct();
        assertEquals(1, sn.chains.getNumRows(), "the model must have a single chain");
        NetworkMomentChainTable Ch = s.getMomentChainTable(3);
        NetworkMomentStationTable St = s.getMomentStationTable(3);
        assertEquals(St.getStationNames().size(), Ch.getStationNames().size());
        for (int k = 0; k < St.getStationNames().size(); k++) {
            assertRelEquals(St.getQLen().get(k), Ch.getQLen().get(k), 1e-9, "QLen " + k);
            assertRelEquals(St.getQLenVar().get(k), Ch.getQLenVar().get(k), 1e-9, "QLenVar " + k);
            assertRelEquals(St.getQLenM3().get(k), Ch.getQLenM3().get(k), 1e-9, "QLenM3 " + k);
            assertRelEquals(St.getQLenSkew().get(k), Ch.getQLenSkew().get(k), 1e-9, "QLenSkew " + k);
        }
    }

    @Test
    @DisplayName("per-class QLenSkew exists and matches pfqn_sens_mom with groups=1:R")
    public void testPerclassSkewExistsAndMatchesApi() {
        // Per-class queue-length skewness IS defined: scaling L(i,r) alone is the
        // class-subset parameter T={r} of Akyildiz-Strelen Theorem 1. Order 3 must
        // therefore expose QLenSkew, and it must equal Pfqn_sens_mom with groups=1:R.
        Network m = closedMixedSched();
        SolverMVA s = new SolverMVA(m);
        NetworkMomentTable T = s.getMomentTable(3);
        assertTrue(T.hasColumn("QLenSkew"));

        Matrix L = new Matrix(2, 2);
        L.set(0, 0, 0.4);
        L.set(0, 1, 0.4);
        L.set(1, 0, 0.3);
        L.set(1, 1, 0.2);
        Matrix N = new Matrix(1, 2);
        N.set(0, 0, 3);
        N.set(0, 1, 2);
        Matrix Z = new Matrix(1, 2);
        Z.set(0, 0, 1.0);
        Z.set(0, 1, 0.5);
        Ret.pfqnSensMom ref = Pfqn_sens_mom.pfqn_sens_mom(L, N, Z, Matrix.ones(1, 2),
                Pfqn_sens_mom.perClassGroups(2));
        int k = rowOf(T, "Q1", "C1");
        assertFalse(Double.isNaN(T.getQLenSkew().get(k)));
        assertRelEquals(ref.Skew.get(0, 0), T.getQLenSkew().get(k), 1e-9, "Q1/C1 QLenSkew");
        // the per-class grouping must also reproduce Pfqn_sens_mva's second moments,
        // which is what pins groups=1:R to the per-class quantity rather than a total
        for (int i = 0; i < 2; i++) {
            for (int r = 0; r < 2; r++) {
                assertRelEquals(T.getMoments().QVar.get(i, r), ref.Var.get(i, r), 1e-9,
                        "groups=1:R Var vs pfqn_sens_mva QVar (" + i + "," + r + ")");
            }
        }
    }

    @Test
    @DisplayName("the method comes from the solver, not from the call")
    public void testMethodComesFromTheSolverNotTheCall() {
        // The algorithm is a property of the solver object, set at construction, not
        // an argument of the table call. Two solvers over the same model must
        // therefore disagree exactly as exact-vs-approximate.
        Network m = closedMixedSched();
        NetworkMomentStationTable ex = new SolverMVA(m).getMomentStationTable(3);
        NetworkMomentStationTable lin =
                new SolverMVA(m, "method", "lin").getMomentStationTable(3);
        assertEquals(ex.getStationNames().size(), lin.getStationNames().size());
        // same model, different solver method -> same quantity, different algorithm
        boolean identical = true;
        for (int k = 0; k < ex.getStationNames().size(); k++) {
            assertRelEquals(ex.getQLen().get(k), lin.getQLen().get(k), 0.021,
                    "solver-selected linearizer QLen " + k);
            assertRelEquals(ex.getQLenM3().get(k), lin.getQLenM3().get(k), 0.062,
                    "solver-selected linearizer QLenM3 " + k);
            if (ex.getQLenM3().get(k).doubleValue() != lin.getQLenM3().get(k).doubleValue()) {
                identical = false;
            }
        }
        // and they must not be bit-identical, else the method is being ignored
        assertFalse(identical, "the solver method is being ignored: the Linearizer "
                + "reproduced the exact recursion bit for bit");
    }

    @Test
    @DisplayName("the Linearizer token family all select the approximate path")
    public void testLinearizerMethodFamily() {
        // MATLAB's isLinearizerMethod accepts the whole Linearizer family, case
        // insensitively; every other method takes the exact path.
        Network m = closedMixedSched();
        NetworkMomentStationTable ex = new SolverMVA(m).getMomentStationTable(3);
        NetworkMomentStationTable lin =
                new SolverMVA(m, "method", "lin").getMomentStationTable(3);
        String[] family = new String[]{"lin", "amva.lin", "egflin", "gflin"};
        for (int f = 0; f < family.length; f++) {
            NetworkMomentStationTable T =
                    new SolverMVA(m, "method", family[f]).getMomentStationTable(3);
            assertTrue(T.getStationNames().size() > 0, family[f] + " produced no rows");
            for (int k = 0; k < T.getStationNames().size(); k++) {
                assertEquals(lin.getQLenM3().get(k), T.getQLenM3().get(k),
                        family[f] + " must take the same Linearizer path as 'lin'");
            }
        }
        // an exact-family method must NOT take the Linearizer path
        NetworkMomentStationTable T = new SolverMVA(m, "method", "exact").getMomentStationTable(3);
        for (int k = 0; k < T.getStationNames().size(); k++) {
            assertEquals(ex.getQLenM3().get(k), T.getQLenM3().get(k),
                    "'exact' must take the exact path");
        }
    }

    @Test
    @DisplayName("chain table rejects the linearizer unless there is a single chain")
    public void testChainTableRejectsLinearizerWithManyChains() {
        // pfqn_sens_linearizer approximates the per-station totals, not a per-chain
        // grouping, so it may only stand in when the two coincide.
        Network m = closedMixedSched();
        final SolverMVA s = new SolverMVA(m, "method", "lin");
        assertEquals(2, m.getStruct().chains.getNumRows(), "this model must have two chains");
        assertThrows(RuntimeException.class, new Executable() {
            public void execute() {
                s.getMomentChainTable(2);
            }
        });
    }

    @Test
    @DisplayName("station table rejects open models")
    public void testStationTableRejectsOpen() {
        // The higher-moment recursion of the reference is stated for closed networks;
        // the method must say so rather than return something.
        Network m = new Network("mt_open3");
        Source src = new Source(m, "S");
        Sink snk = new Sink(m, "K");
        Queue q = new Queue(m, "O1", SchedStrategy.FCFS);
        OpenClass oc = new OpenClass(m, "O", 0);
        src.setArrival(oc, new Exp(0.5));
        q.setService(oc, new Exp(1 / 1.0));
        m.link(Network.serialRouting(src, q, snk));
        final SolverMVA s = new SolverMVA(m);
        assertThrows(RuntimeException.class, new Executable() {
            public void execute() {
                s.getMomentStationTable();
            }
        });
    }

    @Test
    @DisplayName("both tables refuse a non-product-form model")
    public void testNonProductFormIsRefused() {
        // The moment identity Cov = L dQ/dL is a theorem about the product form.
        // Outside it, L dQ/dL is still computable and is simply NOT a covariance, so
        // returning a number would be a confident wrong answer. Both tables must
        // refuse. A heterogeneous-FCFS model (class-dependent rates at an FCFS
        // station) is the canonical non-product-form case.
        Network m = new Network("mt_npf");
        Delay d = new Delay(m, "Think");
        Queue f = new Queue(m, "F", SchedStrategy.FCFS);
        ClosedClass k1 = new ClosedClass(m, "C1", 2, d, 0);
        ClosedClass k2 = new ClosedClass(m, "C2", 1, d, 0);
        d.setService(k1, new Exp(1 / 1.0));
        d.setService(k2, new Exp(1 / 1.0));
        f.setService(k1, new Exp(1 / 0.4));
        f.setService(k2, new Exp(1 / 0.9));   // heterogeneous FCFS
        RoutingMatrix P = m.initRoutingMatrix();
        P.set(k1, k1, Network.serialRouting(d, f));
        P.set(k2, k2, Network.serialRouting(d, f));
        m.link(P);
        assertFalse(SnHasProductForm.snHasProductForm(m.getStruct()),
                "the fixture must be non-product-form");
        final SolverMVA s = new SolverMVA(m);
        // the message is asserted so that the test cannot pass on some unrelated
        // rejection and leave the product-form gate itself unreached
        Throwable e1 = assertThrows(RuntimeException.class, new Executable() {
            public void execute() {
                s.getMomentStationTable(2);
            }
        });
        assertTrue(e1.getMessage().contains("requires a product-form model"),
                "must be refused by the product-form gate, but got: " + e1.getMessage());
        final SolverMVA s2 = new SolverMVA(m);
        Throwable e2 = assertThrows(RuntimeException.class, new Executable() {
            public void execute() {
                s2.getMomentChainTable(2);
            }
        });
        assertTrue(e2.getMessage().contains("requires a product-form model"),
                "must be refused by the product-form gate, but got: " + e2.getMessage());
    }

    @Test
    @DisplayName("the numerical-derivative path tracks the exact one")
    public void testFdFallbackTracksExact() {
        // A product-form model solved by an approximate MVA method has no
        // hand-differentiated moment implementation, so the solver's own means are
        // differentiated numerically and the identity applied. The moments must track
        // the exact ones, and must inherit the accuracy of THAT method's means:
        // 'amva' has near-exact means here, so its moments must be near-exact too.
        Network m = closedMixedSched();
        NetworkMomentStationTable ex = new SolverMVA(m).getMomentStationTable(3);
        NetworkMomentStationTable fd =
                new SolverMVA(m, "method", "amva").getMomentStationTable(3);
        List<String> st = ex.getStationNames();
        assertEquals(st.size(), fd.getStationNames().size(),
                "the two paths must report the same stations");
        boolean anyDifferent = false;
        for (int i = 0; i < st.size(); i++) {
            assertEquals(st.get(i), fd.getStationNames().get(i), "station order must agree");
            assertRelEquals(ex.getQLen().get(i).doubleValue(),
                    fd.getQLen().get(i).doubleValue(), 0.02, "QLen at " + st.get(i));
            assertRelEquals(ex.getQLenVar().get(i).doubleValue(),
                    fd.getQLenVar().get(i).doubleValue(), 0.05, "QLenVar at " + st.get(i));
            assertRelEquals(ex.getQLenM3().get(i).doubleValue(),
                    fd.getQLenM3().get(i).doubleValue(), 0.05, "QLenM3 at " + st.get(i));
            if (ex.getQLenM3().get(i).doubleValue() != fd.getQLenM3().get(i).doubleValue()) {
                anyDifferent = true;
            }
        }
        // the numerical path must not be silently returning the analytic answer
        assertTrue(anyDifferent,
                "the numerical path returned the analytic answer bit-for-bit, so it did "
                        + "not actually run");
        // and the chain table must take the same path
        NetworkMomentChainTable fdc =
                new SolverMVA(m, "method", "amva").getMomentChainTable(3);
        assertEquals(new SolverMVA(m).getMomentChainTable(3).getStationNames().size(),
                fdc.getStationNames().size(),
                "the chain table must report the same rows on both paths");
    }

    @Test
    @DisplayName("the oracle covers every solver and method")
    public void testOracleCoversEverySolverAndMethod() {
        // The identity needs only the mean queue lengths as a function of the demands,
        // so the numerical-derivative oracle re-runs THIS solver, whatever it is. That
        // makes every method of every product-form solver usable, not just those with
        // a hand-differentiated implementation. Restricting it to one analyzer would
        // have restricted the moments to that analyzer's methods for no mathematical
        // reason.
        Network m = closedMixedSched();
        NetworkMomentStationTable ex = new SolverMVA(m).getMomentStationTable(3);

        // 'sum' is a SolverMVA method that no hand-differentiated implementation
        // reaches; it must now produce moments rather than an error.
        NetworkMomentStationTable T =
                new SolverMVA(m, "method", "sum").getMomentStationTable(3);
        assertEquals(ex.getStationNames().size(), T.getStationNames().size(),
                "the 'sum' method must report the same stations as the exact path");
        for (int i = 0; i < T.getQLenVar().size(); i++) {
            assertTrue(T.getQLenVar().get(i).doubleValue() > 0,
                    "'sum' must produce a positive variance at "
                            + T.getStationNames().get(i));
        }

        // SolverNC's convolution algorithm is EXACT, so differentiating its means must
        // reproduce the exact moments. This is the strongest evidence the oracle is
        // sound: a different solver, a different algorithm, the same answer.
        NetworkMomentStationTable NCca =
                new SolverNC(m, "method", "ca").getMomentStationTable(3);
        List<String> st = ex.getStationNames();
        assertEquals(st.size(), NCca.getStationNames().size(),
                "SolverNC 'ca' must report the same stations as the exact path");
        for (int i = 0; i < st.size(); i++) {
            assertEquals(st.get(i), NCca.getStationNames().get(i), "station order must agree");
            assertRelEquals(ex.getQLen().get(i).doubleValue(),
                    NCca.getQLen().get(i).doubleValue(), 1e-6, "QLen at " + st.get(i));
            assertRelEquals(ex.getQLenVar().get(i).doubleValue(),
                    NCca.getQLenVar().get(i).doubleValue(), 1e-4, "QLenVar at " + st.get(i));
            assertRelEquals(ex.getQLenM3().get(i).doubleValue(),
                    NCca.getQLenM3().get(i).doubleValue(), 1e-4, "QLenM3 at " + st.get(i));
        }
    }
}
