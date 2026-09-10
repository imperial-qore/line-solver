/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.mva;

import jline.api.pfqn.mva.Pfqn_mvaoi;
import jline.api.pfqn.mva.Pfqn_mvaoi_marg;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Exp;
import jline.solvers.NetworkAvgTable;
import jline.solvers.SolverOptions;
import jline.solvers.ctmc.SolverCTMC;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.ToDoubleFunction;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;

/**
 * Regression tests for the exact mean-value MVA path of order-independent (OI)
 * closed queueing networks ({@link SolverMVAOIAnalyzer} -> {@link Pfqn_mvaoi}).
 *
 * <p>An OI station is a class-dependent load-dependent server whose total service
 * rate mu(c) is a permutation-invariant function of the per-class occupancy. A
 * closed network of infinite-server (delay) and load-independent (single-server,
 * product-form) stations plus any number of OI stations is product-form, so the
 * mean-value CMVA is EXACT: its per-class throughput and queue-lengths coincide
 * with the CTMC stationary solution to machine precision. These tests lock that
 * agreement to guard against future drift; the model family mirrors the MATLAB
 * research script test_oi_6 (multiserver OI stations with server-compatibility
 * rates and an IS delay).
 *
 * <p>Utilization is compared too: at an OI station LINE defines U_r = E[sir_r]/c,
 * the mean number of class-r jobs receiving a strictly positive rank rate over
 * the number of servers, and the MVA path evaluates it exactly from the OI count
 * marginal (see {@link jline.api.pfqn.ld.Pfqn_oi_insvc}), matching the CTMC. It
 * previously reported the offered load X_r/rate_r, which coincides only when a
 * job engages a single server.
 *
 * @see SolverMVAOIAnalyzer
 * @see jline.api.pfqn.mva.Pfqn_mvaoi
 */
public class SolverMVAOITest {

    /** Exact solvers: lock agreement tightly to catch any numerical drift. */
    private static final double TOL = 1e-9;

    // ---- model builders (closed OI networks: IS delay + OI/LI queues) -------

    /** Single OI station + IS delay, composition-dependent rate. */
    private static Network singleOI() {
        Network model = new Network("OI1");
        Delay d = new Delay(model, "IS");
        Queue q = new Queue(model, "OI", SchedStrategy.OI);
        ClosedClass c1 = new ClosedClass(model, "C1", 2, d);
        ClosedClass c2 = new ClosedClass(model, "C2", 1, d);
        d.setService(c1, new Exp(1.0));
        d.setService(c2, new Exp(2.0));
        q.setService((Matrix c) -> 1.0 + 0.5 * count(c, 0) + 1.0 * count(c, 1));
        q.setNumberOfServers(1);
        q.setCap(3);
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(c1, c1, Network.serialRouting(d, q));
        P.set(c2, c2, Network.serialRouting(d, q));
        model.link(P);
        return model;
    }

    /** Two OI stations + IS delay. */
    private static Network twoOI() {
        Network model = new Network("OI2");
        Delay d = new Delay(model, "IS");
        Queue q1 = new Queue(model, "OIa", SchedStrategy.OI);
        Queue q2 = new Queue(model, "OIb", SchedStrategy.OI);
        ClosedClass c1 = new ClosedClass(model, "C1", 2, d);
        ClosedClass c2 = new ClosedClass(model, "C2", 1, d);
        d.setService(c1, new Exp(1.0));
        d.setService(c2, new Exp(2.0));
        q1.setService((Matrix c) -> 1.0 + 0.5 * count(c, 0) + 1.0 * count(c, 1));
        q2.setService((Matrix c) -> 1.5 + 1.0 * count(c, 0) + 0.3 * count(c, 1));
        q1.setNumberOfServers(1); q1.setCap(3);
        q2.setNumberOfServers(1); q2.setCap(3);
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(c1, c1, Network.serialRouting(d, q1, q2));
        P.set(c2, c2, Network.serialRouting(d, q1, q2));
        model.link(P);
        return model;
    }

    /** OI station + load-independent PS queue + IS delay (mixed product form). */
    private static Network oiPlusPS() {
        Network model = new Network("OIPS");
        Delay d = new Delay(model, "IS");
        Queue q1 = new Queue(model, "OI", SchedStrategy.OI);
        Queue q2 = new Queue(model, "PS", SchedStrategy.PS);
        ClosedClass c1 = new ClosedClass(model, "C1", 2, d);
        ClosedClass c2 = new ClosedClass(model, "C2", 1, d);
        d.setService(c1, new Exp(1.0));
        d.setService(c2, new Exp(2.0));
        q1.setService((Matrix c) -> 1.0 + 0.5 * count(c, 0) + 1.0 * count(c, 1));
        q2.setService(c1, new Exp(1.4));
        q2.setService(c2, new Exp(1.1));
        q1.setNumberOfServers(1); q1.setCap(3);
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(c1, c1, Network.serialRouting(d, q1, q2));
        P.set(c2, c2, Network.serialRouting(d, q1, q2));
        model.link(P);
        return model;
    }

    /**
     * OI station + MULTISERVER (c = 2) PS queue + IS delay. Pfqn_mvaoi models an LI
     * queue as a single server, so the multiserver station must be promoted to the
     * OI representation mu(n) = (min(|n|,c)/|n|) sum_r n_r/D_r; otherwise MVA
     * silently returns the c = 1 answer.
     */
    private static Network oiPlusMultiserverPS() {
        Network model = new Network("OIPSms");
        Delay d = new Delay(model, "IS");
        Queue q1 = new Queue(model, "OI", SchedStrategy.OI);
        Queue q2 = new Queue(model, "PS", SchedStrategy.PS);
        ClosedClass c1 = new ClosedClass(model, "C1", 2, d);
        ClosedClass c2 = new ClosedClass(model, "C2", 1, d);
        d.setService(c1, new Exp(1.0));
        d.setService(c2, new Exp(2.0));
        q1.setService((Matrix c) -> 1.0 + 0.5 * count(c, 0) + 1.0 * count(c, 1));
        q2.setService(c1, new Exp(1.4));
        q2.setService(c2, new Exp(1.1));
        q1.setNumberOfServers(1); q1.setCap(3);
        q2.setNumberOfServers(2);
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(c1, c1, Network.serialRouting(d, q1, q2));
        P.set(c2, c2, Network.serialRouting(d, q1, q2));
        model.link(P);
        return model;
    }

    /**
     * OI station + class-DEPENDENT single-server FCFS queue + IS delay: NOT product
     * form, hence outside the exact OI path.
     */
    private static Network oiPlusClassDepFCFS() {
        Network model = new Network("OInonPF");
        Delay d = new Delay(model, "IS");
        Queue q1 = new Queue(model, "OI", SchedStrategy.OI);
        Queue q2 = new Queue(model, "Q", SchedStrategy.FCFS);
        ClosedClass c1 = new ClosedClass(model, "C1", 2, d);
        ClosedClass c2 = new ClosedClass(model, "C2", 1, d);
        d.setService(c1, new Exp(1.0));
        d.setService(c2, new Exp(2.0));
        q1.setService((Matrix c) -> 1.0 + 0.5 * count(c, 0) + 1.0 * count(c, 1));
        q2.setService(c1, new Exp(1.4));   // class-dependent FCFS rates
        q2.setService(c2, new Exp(0.6));
        q1.setNumberOfServers(1); q1.setCap(3);
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(c1, c1, Network.serialRouting(d, q1, q2));
        P.set(c2, c2, Network.serialRouting(d, q1, q2));
        model.link(P);
        return model;
    }

    /**
     * Three multiserver OI stations + IS delay with server-compatibility rates,
     * mirroring the MATLAB research script test_oi_6. Each station's rate is the
     * sum of the capacities of the servers compatible with at least one class
     * currently present.
     */
    private static Network oi6Compat() {
        final int[][][] compat = {
            {{1, 0}, {1, 1}, {0, 1}},
            {{1, 2}, {1, 1}, {1, 1}},
            {{1, 1}, {0, 1}, {1, 0}}
        };
        final double[][] muCap = {{1.0, 2.0, 3.0}, {1.0, 2.0, 3.0}, {1.5, 2.5, 3.5}};
        double[] sigma = {5.0, 1.0};
        int[] N = {2, 2};
        Network model = new Network("OI6");
        Delay d = new Delay(model, "IS");
        Queue[] q = new Queue[3];
        for (int m = 0; m < 3; m++) q[m] = new Queue(model, "OI" + (m + 1), SchedStrategy.OI);
        ClosedClass c1 = new ClosedClass(model, "C1", N[0], d);
        ClosedClass c2 = new ClosedClass(model, "C2", N[1], d);
        d.setService(c1, new Exp(sigma[0]));
        d.setService(c2, new Exp(sigma[1]));
        for (int m = 0; m < 3; m++) {
            final int[][] cm = compat[m];
            final double[] mc = muCap[m];
            q[m].setService((Matrix c) -> compatRate(c, cm, mc));
            q[m].setNumberOfServers(3);
            q[m].setCap(N[0] + N[1]);
        }
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(c1, c1, Network.serialRouting(d, q[0], q[1], q[2]));
        P.set(c2, c2, Network.serialRouting(d, q[0], q[1], q[2]));
        model.link(P);
        return model;
    }

    // ---- OI rate helpers (class ids in the ordered list are 0-based) --------

    /** Number of class-r jobs in the ordered OI job list c. */
    private static int count(Matrix c, int r) {
        int n = 0;
        for (int k = 0; k < c.getNumCols(); k++) {
            if ((int) Math.round(c.get(0, k)) == r) n++;
        }
        return n;
    }

    /** Sum of capacities of servers compatible with a class present in c. */
    private static double compatRate(Matrix c, int[][] comp, double[] muCap) {
        Set<Integer> present = new HashSet<Integer>();
        for (int k = 0; k < c.getNumCols(); k++) {
            int v = (int) Math.round(c.get(0, k));
            if (v >= 0) present.add(v);
        }
        double sum = 0;
        for (int srv = 0; srv < comp.length; srv++) {
            for (int cls : present) {
                if (comp[srv][cls] != 0) { sum += muCap[srv]; break; }
            }
        }
        return sum;
    }

    // ---- solver helpers -----------------------------------------------------

    private static NetworkAvgTable ctmc(Network model, int cutoff) {
        SolverOptions opt = SolverCTMC.defaultOptions();
        opt.cutoff = Matrix.singleton(cutoff);
        return new SolverCTMC(model, opt).getAvgTable();
    }

    private static NetworkAvgTable mva(Network model) {
        return new SolverMVA(model).getAvgTable();
    }

    /** Assert MVA == CTMC on the exact invariants (QLen, Util, RespT, Tput). */
    private void assertExactVsCTMC(Network mvaModel, Network ctmcModel, int cutoff) {
        NetworkAvgTable ref = ctmc(ctmcModel, cutoff);
        NetworkAvgTable got = mva(mvaModel);
        assertVector("QLen", ref.getQLen(), got.getQLen());
        assertVector("Util", ref.getUtil(), got.getUtil());
        assertVector("RespT", ref.getRespT(), got.getRespT());
        assertVector("Tput", ref.getTput(), got.getTput());
    }

    private void assertVector(String metric, List<Double> ref, List<Double> got) {
        assertEquals(ref.size(), got.size(), metric + ": row count mismatch");
        for (int i = 0; i < ref.size(); i++) {
            assertEquals(ref.get(i), got.get(i), TOL,
                    "OI MVA vs CTMC " + metric + "[" + i + "] mismatch");
        }
    }

    // ---- solver-level exactness tests --------------------------------------

    @Test
    @DisplayName("MVA OI: single OI + delay matches CTMC")
    public void singleOIMatchesCTMC() {
        assertExactVsCTMC(singleOI(), singleOI(), 3);
    }

    @Test
    @DisplayName("MVA OI: two OI stations + delay matches CTMC")
    public void twoOIMatchesCTMC() {
        assertExactVsCTMC(twoOI(), twoOI(), 3);
    }

    @Test
    @DisplayName("MVA OI: OI + PS queue + delay matches CTMC")
    public void oiPlusPSMatchesCTMC() {
        assertExactVsCTMC(oiPlusPS(), oiPlusPS(), 3);
    }

    @Test
    @DisplayName("MVA OI: multiserver compatibility model (test_oi_6) matches CTMC")
    public void oi6CompatMatchesCTMC() {
        assertExactVsCTMC(oi6Compat(), oi6Compat(), 4);
    }

    @Test
    @DisplayName("MVA OI: OI + multiserver (c=2) PS queue + delay matches CTMC")
    public void oiPlusMultiserverPSMatchesCTMC() {
        assertExactVsCTMC(oiPlusMultiserverPS(), oiPlusMultiserverPS(), 3);
    }

    @Test
    @DisplayName("MVA OI: rejects a non-product-form OI model instead of returning zeros")
    public void mvaRejectsNonProductForm() {
        // AMVA cannot represent an OI rank rate mu(n) (it only sees sn.rates), so
        // SolverMVA must refuse rather than silently return a zero queue-length at
        // the OI station.
        Network model = oiPlusClassDepFCFS();
        assertFalse(jline.solvers.nc.handlers.Solver_nc_oi.nc_is_oi_model(model.getStruct()));
        assertThrows(RuntimeException.class, () -> new SolverMVA(model).getAvgTable());
    }

    // ---- API-level cross-check: mean-value CMVA == marginal LD-MVA ----------

    @Test
    @DisplayName("Pfqn_mvaoi (CMVA) equals Pfqn_mvaoi_marg (marginal) incl. LI queues")
    public void mvaoiEqualsMarg() {
        // Z (delay), K OI rate handles over count vectors, J LI-queue demands.
        double[] Z = {0.5, 0.6};
        int[] N = {3, 3};
        double[][] Dli = {{0.4, 0.7}};
        List<ToDoubleFunction<int[]>> mu = new ArrayList<ToDoubleFunction<int[]>>();
        mu.add(n -> 1.0 + 0.5 * n[0] + 1.0 * n[1]);
        mu.add(n -> 1.5 + 1.0 * n[0] + 0.2 * n[1]);
        int R = 2, K = mu.size(), J = Dli.length, M = 1 + J + K;

        Pfqn_mvaoi.Result r = Pfqn_mvaoi.pfqn_mvaoi(Z, N, mu, Dli);

        // Marginal form over the full station list: row 0 delay, then LI, then OI.
        double[][] D = new double[M][R];
        boolean[] isDelay = new boolean[M];
        List<ToDoubleFunction<int[]>> mm = new ArrayList<ToDoubleFunction<int[]>>();
        D[0] = Z.clone(); isDelay[0] = true; mm.add(null);
        for (int j = 0; j < J; j++) { D[1 + j] = Dli[j].clone(); mm.add(null); }
        for (int i = 0; i < K; i++) mm.add(mu.get(i));
        Pfqn_mvaoi_marg.Result rm = Pfqn_mvaoi_marg.pfqn_mvaoi_marg(D, N, isDelay, mm);

        for (int x = 0; x < R; x++) {
            assertEquals(rm.XN[x], r.X[x], TOL, "throughput X[" + x + "]");
            assertEquals(rm.QN[0][x], r.Qdelay[x], TOL, "delay QLen[" + x + "]");
        }
        for (int j = 0; j < J; j++) {
            for (int x = 0; x < R; x++) {
                assertEquals(rm.QN[1 + j][x], r.Qli[j][x], TOL, "LI QLen[" + j + "][" + x + "]");
            }
        }
        for (int i = 0; i < K; i++) {
            for (int x = 0; x < R; x++) {
                assertEquals(rm.QN[1 + J + i][x], r.Qoi[i][x], TOL, "OI QLen[" + i + "][" + x + "]");
            }
        }
    }
}
