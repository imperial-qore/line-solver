package jline.solvers.ba;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Arrays;
import java.util.List;

import org.junit.jupiter.api.Test;

import jline.api.sn.SnToQrfAlpha;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.solvers.NetworkAvgTable;
import jline.solvers.ctmc.SolverCTMC;
import jline.util.matrix.Matrix;

/**
 * The load-dependent QRF arms on delay, multiserver and load-dependent models.
 *
 * <p>WHAT CHANGED AND WHY IT IS NOT AN APPROXIMATION. {@code qrf.mmi.ld} and
 * {@code qrf.mmi.linear} carry a scaling alpha(i,n) that multiplies every rate
 * out of station i at population n, completions and background phase changes
 * alike. That is exactly the rate law of an infinite server (alpha = n), of a
 * c-server station (alpha = min(n,c)) and of limited load dependence, so
 * deriving alpha from the model makes the relaxed chain the model's OWN chain
 * rather than an approximation of it. The two arms therefore serve models the
 * rest of the QRF family still refuses.
 *
 * <p>THE ORACLE IS EXACTNESS, NOT CLOSENESS. At M = 2 the pairwise joint of a
 * closed chain is fully determined by the marginal, so the QRF polytope is tight
 * and the answer must EQUAL the exact CTMC one. Every two-station fixture here
 * is asserted at 1e-9 against SolverCTMC, which pins the alpha derivation, the
 * BN readout and the utilization normalizer at once; a looser tolerance would
 * let a wrong alpha through.
 */
public class SolverBaQrfLdTest {

    private static final double TOL = 1e-9;

    /** Two stations in a cycle; station 0 is a Delay when delay0. */
    private static Network cqn(int N, double c0, boolean delay0, double[] lld0, int phases0) {
        Network model = new Network("qrfLd");
        jline.lang.nodes.Station s0;
        if (delay0) {
            s0 = new Delay(model, "D0");
        } else {
            Queue q = new Queue(model, "Q0", SchedStrategy.FCFS);
            q.setNumberOfServers((int) c0);
            if (lld0 != null) {
                q.setLoadDependence(new Matrix(new double[][]{lld0}));
            }
            s0 = q;
        }
        Queue q1 = new Queue(model, "Q1", SchedStrategy.FCFS);
        ClosedClass c = new ClosedClass(model, "C", N, s0);
        if (phases0 > 1) {
            s0.setService(c, Erlang.fitMeanAndOrder(1.0, phases0));
        } else {
            s0.setService(c, new Exp(1.0));
        }
        q1.setService(c, new Exp(1.5));
        model.link(model.serialRouting(s0, q1));
        return model;
    }

    private static Network plain(int N, double c0) {
        return cqn(N, c0, false, null, 1);
    }

    private static double[] col(List<Double> v) {
        double[] out = new double[v.size()];
        for (int i = 0; i < v.size(); i++) out[i] = v.get(i);
        return out;
    }

    private static void assertSameAsCtmc(Network model, String method, String label) {
        NetworkAvgTable q = new SolverBA(model, method).getAvgTable();
        NetworkAvgTable e = new SolverCTMC(model).getAvgTable();
        double[] qq = col(q.getQLen()), qe = col(e.getQLen());
        double[] uq = col(q.getUtil()), ue = col(e.getUtil());
        double[] tq = col(q.getTput()), te = col(e.getTput());
        for (int i = 0; i < qq.length; i++) {
            assertEquals(qe[i], qq[i], TOL, label + " QLen[" + i + "]");
            assertEquals(ue[i], uq[i], TOL, label + " Util[" + i + "]");
            assertEquals(te[i], tq[i], TOL, label + " Tput[" + i + "]");
        }
    }

    // ---- the alpha derivation ----

    @Test
    public void alphaIsMinNCAtAMultiserverStation() {
        SnToQrfAlpha.Result r = SnToQrfAlpha.snToQrfAlpha(plain(4, 3).getStruct());
        assertTrue(r.msg.isEmpty(), r.msg);
        assertTrue(r.ld);
        assertArrayEquals3(new double[]{1, 2, 3, 3}, r.alpha[0]);
        assertArrayEquals3(new double[]{1, 1, 1, 1}, r.alpha[1]);
        assertEquals(3.0, r.peak[0], TOL);
        assertEquals(1.0, r.peak[1], TOL);
    }

    @Test
    public void alphaPeakIsTheDeclaredServersNotTheReachableMaximum() {
        // c = 3 with N = 2: alpha reaches 2, but the station still has three
        // servers and LINE reports U = T*S/c. Normalizing by max(alpha) would
        // overstate the utilization by 3/2.
        SnToQrfAlpha.Result r = SnToQrfAlpha.snToQrfAlpha(plain(2, 3).getStruct());
        assertTrue(r.msg.isEmpty(), r.msg);
        assertArrayEquals3(new double[]{1, 2}, r.alpha[0]);
        assertEquals(3.0, r.peak[0], TOL);
    }

    @Test
    public void alphaIsNAtADelay() {
        SnToQrfAlpha.Result r = SnToQrfAlpha.snToQrfAlpha(cqn(3, 1, true, null, 1).getStruct());
        assertTrue(r.msg.isEmpty(), r.msg);
        assertTrue(r.ld);
        assertArrayEquals3(new double[]{1, 2, 3}, r.alpha[0]);
        assertTrue(Double.isInfinite(r.peak[0]));
        assertEquals(1.0, r.peak[1], TOL);
    }

    @Test
    public void alphaComposesLldWithTheServerCount() {
        SnToQrfAlpha.Result r = SnToQrfAlpha.snToQrfAlpha(
                cqn(3, 2, false, new double[]{1.0, 1.5, 2.0}, 1).getStruct());
        assertTrue(r.msg.isEmpty(), r.msg);
        assertTrue(r.ld);
        assertArrayEquals3(new double[]{1 * 1.0, 2 * 1.5, 2 * 2.0}, r.alpha[0]);
        // The peak is the PRODUCT of the server count and the reachable lld peak.
        assertEquals(4.0, r.peak[0], TOL);
    }

    @Test
    public void alphaIsAllOnesOnASingleServerModel() {
        // The load-independent model must reach the ld arms unchanged, which is
        // what keeps this change a no-op there.
        SnToQrfAlpha.Result r = SnToQrfAlpha.snToQrfAlpha(plain(3, 1).getStruct());
        assertTrue(r.msg.isEmpty(), r.msg);
        assertFalse(r.ld);
        for (int i = 0; i < r.alpha.length; i++)
            for (int n = 0; n < r.alpha[i].length; n++) assertEquals(1.0, r.alpha[i][n], TOL);
        assertEquals(1.0, r.peak[0], TOL);
    }

    @Test
    public void phasetypeWhereSeveralJobsAreServedAtOnceIsRefused() {
        // The QRF local state carries ONE phase per station, which describes one
        // job in service and no more, so a PH multiserver would answer a
        // different chain and the number would bound nothing.
        for (Network model : Arrays.asList(cqn(2, 2, false, null, 2), cqn(2, 1, true, null, 2))) {
            SnToQrfAlpha.Result r = SnToQrfAlpha.snToQrfAlpha(model.getStruct());
            assertTrue(r.msg.contains("one phase per station"), r.msg);
            assertTrue(r.ld, "ld must stay set through the refusal");
            assertThrows(RuntimeException.class,
                    () -> new SolverBA(model, "qrf.mmi.ld").getAvgTable());
        }
    }

    // ---- dispatch ----

    @Test
    public void alphaFreeArmsRefuseAndNameTheArmsThatServe() {
        for (String method : new String[]{"qr", "qrf.mmi", "qrf.mem", "qrf.bethe"}) {
            for (Network model : Arrays.asList(plain(2, 2), cqn(2, 1, true, null, 1),
                    cqn(2, 1, false, new double[]{1.0, 2.0}, 1))) {
                assertThrows(RuntimeException.class,
                        () -> new SolverBA(model, method).getAvgTable(),
                        method + " must refuse a load-dependent model");
            }
        }
    }

    @Test
    public void listValidMethodsKeepsTheLdArmsOnADelayModel() {
        // A caller enumerating the list must still see the only two bound
        // methods the model has; dropping them with the rest would hide them.
        List<String> methods =
                Arrays.asList(new SolverBA(cqn(2, 1, true, null, 1)).listValidMethods());
        assertTrue(methods.contains("qrf.mmi.ld"), methods.toString());
        assertTrue(methods.contains("qrf.mmi.linear"), methods.toString());
        assertFalse(methods.contains("qrf.mmi"), methods.toString());
        assertFalse(methods.contains("qrf.bas"), methods.toString());
    }

    // ---- the numerical oracle ----

    @Test
    public void twoStationAnswersAreExact() {
        // M = 2 makes the polytope tight, so these are equalities, not bounds.
        for (String method : new String[]{"qrf.mmi.ld", "qrf.mmi.linear"}) {
            assertSameAsCtmc(cqn(3, 1, true, null, 1), method, method + " delay");
            assertSameAsCtmc(plain(3, 2), method, method + " c=2");
            assertSameAsCtmc(plain(4, 3), method, method + " c=3");
            assertSameAsCtmc(cqn(3, 1, false, new double[]{1.0, 1.5, 2.0}, 1), method,
                    method + " lld");
            assertSameAsCtmc(plain(3, 1), method, method + " single-server");
        }
    }

    @Test
    public void multiserverAtPopulationOneIsTheSingleServerModel() {
        // min(1,c) = 1, so the two chains are identical and only the utilization
        // normalizer differs. No CTMC needed: the oracle is the other run.
        for (String method : new String[]{"qrf.mmi.ld", "qrf.mmi.linear"}) {
            NetworkAvgTable one = new SolverBA(plain(1, 1), method).getAvgTable();
            NetworkAvgTable three = new SolverBA(plain(1, 3), method).getAvgTable();
            double[] q1 = col(one.getQLen()), q3 = col(three.getQLen());
            double[] t1 = col(one.getTput()), t3 = col(three.getTput());
            double[] u1 = col(one.getUtil()), u3 = col(three.getUtil());
            for (int i = 0; i < q1.length; i++) {
                assertEquals(q1[i], q3[i], TOL, method + " QLen[" + i + "]");
                assertEquals(t1[i], t3[i], TOL, method + " Tput[" + i + "]");
            }
            assertEquals(u1[0], u3[0] * 3.0, TOL, method + " Util[0]");
        }
    }

    private static void assertArrayEquals3(double[] expected, double[] actual) {
        assertEquals(expected.length, actual.length, "length");
        for (int i = 0; i < expected.length; i++) assertEquals(expected[i], actual[i], TOL, "[" + i + "]");
    }
}
