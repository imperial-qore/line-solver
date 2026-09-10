package jline.api;

import jline.api.nc.MeGegecnResult;
import jline.api.nc.MeOqnBlkResult;
import jline.api.nc.Me_gegecn;
import jline.api.nc.Me_oqn;
import jline.api.nc.Me_oqn_blk;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Validates the censored GE/GE/c/K;N building block of Kouvatsos (1994),
 * Section 4.1, against the exact Markovian queue-length distributions, and
 * the transfer-blocking network algorithm of Tahilramani, Manjunath and Bose
 * (1999) against the results published in their Table 1.
 */
public class MeOqnBlkTest {

    private static final double TOL = 1e-10;

    private static Matrix col(double... v) {
        Matrix m = new Matrix(v.length, 1);
        for (int i = 0; i < v.length; i++) {
            m.set(i, 0, v[i]);
        }
        return m;
    }

    /** Censored GE/GE/1/0;N reduces to M/M/1/N when the streams are Markovian. */
    @Test
    public void testCensoredSingleServerIsExact() {
        double[] rhos = {0.3, 0.8, 1.0, 1.5};
        int[] caps = {1, 2, 5, 10};
        for (double rho : rhos) {
            for (int n : caps) {
                MeGegecnResult r = Me_gegecn.me_gegecn(rho, 1.0, 1.0, 1.0, 1, 0, n);
                double[] un = new double[n + 1];
                double sum = 0.0;
                for (int k = 0; k <= n; k++) {
                    un[k] = Math.pow(rho, k);
                    sum += un[k];
                }
                double lex = 0.0;
                for (int k = 0; k <= n; k++) {
                    un[k] /= sum;
                    lex += k * un[k];
                    assertEquals(un[k], r.getP()[k], TOL,
                            "state probability rho=" + rho + " N=" + n + " k=" + k);
                }
                assertEquals(lex, r.getL(), TOL, "mean queue length rho=" + rho + " N=" + n);
                // A Poisson stream is blocked exactly when the buffer is full (PASTA)
                assertEquals(un[n], r.getPB(), TOL, "blocking probability rho=" + rho + " N=" + n);
                assertEquals(1 - un[0], r.getU(), TOL, "utilization rho=" + rho + " N=" + n);
            }
        }
    }

    /** Censored GE/GE/c/0;N reduces to M/M/c/N when the streams are Markovian. */
    @Test
    public void testCensoredMultiserverIsExact() {
        int[] servers = {2, 3, 5};
        for (int c : servers) {
            for (int n : new int[]{c, c + 3, c + 8}) {
                double lam = 0.7 * c;
                MeGegecnResult r = Me_gegecn.me_gegecn(lam, 1.0, 1.0, 1.0, c, 0, n);
                double[] un = new double[n + 1];
                double sum = 0.0;
                double fact = 1.0;
                for (int k = 0; k <= n; k++) {
                    if (k > 0) {
                        fact *= k;
                    }
                    if (k <= c) {
                        un[k] = Math.pow(lam, k) / fact;
                    } else {
                        double factc = 1.0;
                        for (int t = 1; t <= c; t++) {
                            factc *= t;
                        }
                        un[k] = Math.pow(lam, c) / factc * Math.pow(lam / c, k - c);
                    }
                    sum += un[k];
                }
                double lex = 0.0;
                double ubusy = 0.0;
                for (int k = 0; k <= n; k++) {
                    un[k] /= sum;
                    lex += k * un[k];
                    ubusy += Math.min(k, c) * un[k];
                }
                assertEquals(lex, r.getL(), TOL, "mean queue length c=" + c + " N=" + n);
                assertEquals(ubusy / c, r.getU(), TOL, "utilization c=" + c + " N=" + n);
                assertEquals(un[n], r.getPB(), TOL, "blocking probability c=" + c + " N=" + n);
            }
        }
    }

    /** The network algorithm on a single M/M/1/N with loss is exact. */
    @Test
    public void testSingleStationWithLoss() {
        for (double rho : new double[]{0.5, 0.9, 1.2}) {
            int n = 4;
            MeOqnBlkResult r = Me_oqn_blk.me_oqn_blk(1, col(rho), col(1.0), col(1.0), col(1.0),
                    new Matrix(1, 1), col(1.0), col(n), new int[]{Me_oqn_blk.RULE_LOSS});
            double sum = 0.0;
            double lex = 0.0;
            for (int k = 0; k <= n; k++) {
                sum += Math.pow(rho, k);
                lex += k * Math.pow(rho, k);
            }
            lex /= sum;
            double pn = Math.pow(rho, n) / sum;
            assertEquals(lex, r.getQ().get(0, 0), 1e-6, "queue length rho=" + rho);
            assertEquals(rho * (1 - pn), r.getT().get(0, 0), 1e-5, "carried throughput rho=" + rho);
            assertEquals(pn, r.getPBa().get(0, 0), 1e-5, "blocking probability rho=" + rho);
        }
    }

    /** With every buffer unbounded the algorithm reproduces me_oqn. */
    @Test
    public void testUnboundedReducesToMeOqn() {
        int m = 3;
        Matrix p = new Matrix(m, m);
        p.set(0, 1, 0.6);
        p.set(0, 2, 0.4);
        p.set(1, 2, 1.0);
        Matrix inf = col(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);
        MeOqnBlkResult r = Me_oqn_blk.me_oqn_blk(m, col(1.0, 0, 0), col(1.0, 1.0, 1.0),
                col(3.0, 2.5, 2.0), col(1.0, 2.0, 0.5), p, col(1.0, 1.0, 1.0), inf,
                new int[]{Me_oqn_blk.RULE_LOSS, Me_oqn_blk.RULE_LOSS, Me_oqn_blk.RULE_LOSS});
        Matrix[][] pr = new Matrix[m][m];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                pr[i][j] = new Matrix(1, 1);
                pr[i][j].set(0, 0, p.get(i, j));
            }
        }
        Matrix lam0 = new Matrix(m, 1);
        lam0.set(0, 0, 1.0);
        Matrix ca0 = new Matrix(m, 1);
        for (int i = 0; i < m; i++) {
            ca0.set(i, 0, 1.0);
        }
        Matrix mu = new Matrix(m, 1);
        mu.set(0, 0, 3.0);
        mu.set(1, 0, 2.5);
        mu.set(2, 0, 2.0);
        Matrix cs = new Matrix(m, 1);
        cs.set(0, 0, 1.0);
        cs.set(1, 0, 2.0);
        cs.set(2, 0, 0.5);
        jline.api.nc.MeOqnResult ref = Me_oqn.me_oqn(m, 1, lam0, ca0, mu, cs, pr, col(1.0, 1.0, 1.0),
                new boolean[m], new jline.api.nc.MeOqnOptions());
        for (int i = 0; i < m; i++) {
            assertEquals(ref.getL().get(i, 0), r.getQ().get(i, 0), 1e-10, "station " + i);
        }
    }

    /**
     * Feed-forward network of Tahilramani, Manjunath and Bose (1999), Table 1.
     * The published values come from an independent implementation of the same
     * algorithm, so agreement is asserted at the level of the algorithm rather
     * than bitwise: 5% on the queue lengths, 2% on the throughputs.
     */
    @Test
    public void testTransferBlockingAgainstPublishedTable() {
        int m = 3;
        Matrix p = new Matrix(m, m);
        p.set(0, 1, 0.4);
        p.set(0, 2, 0.4);
        p.set(1, 2, 0.5);
        MeOqnBlkResult r = Me_oqn_blk.me_oqn_blk(m, col(1.5, 0, 0), col(2.0, 1.0, 1.0),
                col(2.0, 2.0, 2.0), col(1.0, 1.0, 1.0), p, col(1.0, 3.0, 1.0), col(5, 4, 3),
                new int[]{Me_oqn_blk.RULE_BAS, Me_oqn_blk.RULE_BAS, Me_oqn_blk.RULE_BAS});
        double[] kPaper = {1.702, 0.261, 0.578};
        double[] tPaper = {1.284, 0.514, 0.771};
        for (int i = 0; i < m; i++) {
            assertTrue(Math.abs(r.getQ().get(i, 0) - kPaper[i]) / kPaper[i] < 0.05,
                    "queue length at station " + i + " is " + r.getQ().get(i, 0)
                            + ", published " + kPaper[i]);
            assertTrue(Math.abs(r.getT().get(i, 0) - tPaper[i]) / tPaper[i] < 0.02,
                    "throughput at station " + i + " is " + r.getT().get(i, 0)
                            + ", published " + tPaper[i]);
        }
    }

    /** The GE domain is enforced rather than approximated. */
    @Test
    public void testHypoExponentialIsRejected() {
        assertThrows(RuntimeException.class, () -> Me_gegecn.me_gegecn(0.5, 1.0, 1.0, 0.5, 1, 0, 4),
                "a service scv below 1 must be refused, the GE distribution being undefined there");
    }
}
