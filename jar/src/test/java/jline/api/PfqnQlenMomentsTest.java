package jline.api;

import java.util.ArrayList;
import java.util.List;

import jline.io.Ret;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static jline.api.moment.Moment_binomial_from_tail.moment_binomial_from_tail;
import static jline.api.moment.Moment_joint_binomial_from_tail.moment_joint_binomial_from_tail;
import static jline.api.moment.Moment_joint_central_from_tail.moment_joint_central_from_tail;
import static jline.api.moment.Moment_joint_tail_from_binomial.moment_joint_tail_from_binomial;
import static jline.api.moment.Moment_tail_from_binomial.moment_tail_from_binomial;
import static jline.api.pfqn.nc.Pfqn_qlen_joint_moments.pfqn_qlen_joint_moments;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Unit tests for the tail vertex of the house of moments and for
 * Pfqn_qlen_joint_moments, which turns normalizing constants into the joint
 * moments of the queue-length vector of a closed product-form network.
 *
 * <p>Both routes of the routine are refereed by the same independent oracle:
 * the exact product-form state distribution, enumerated over the whole state
 * space and summed directly. The oracle knows nothing about normalizing
 * constants, survival identities or moment conversions, so an agreement checks
 * the whole chain rather than a round trip. The identities behind the two
 * routes are proved symbolically in sage/proofs/qlen_tail_moments.py, and the
 * same tests exist in python/tests/test_pfqn_qlen_moments.py and in
 * test_pfqn_qlen_moments.m in line-test.git.
 */
public class PfqnQlenMomentsTest {

    private static final double TOL = 1e-9;

    // ---------- product-form oracle --------------------------------------

    /** One state of the enumeration: per-class populations at every station. */
    private static class State {
        int[][] n;   // (M+1) x R, the last row being the delay
        double p;

        State(int[][] n) {
            this.n = n;
        }
    }

    /**
     * Enumerates the exact product-form state distribution.
     *
     * @param L demand matrix (M x R)
     * @param N population vector
     * @param Z think time vector
     * @return the states with their probabilities
     */
    private static List<State> oracleStates(double[][] L, int[] N, double[] Z) {
        int M = L.length;
        int R = L[0].length;
        List<int[][]> partial = new ArrayList<int[][]>();
        partial.add(new int[M + 1][R]);
        for (int r = 0; r < R; r++) {
            List<int[][]> next = new ArrayList<int[][]>();
            for (int[][] base : partial) {
                for (int[] alloc : compositions(N[r], M + 1)) {
                    int[][] copy = new int[M + 1][R];
                    for (int i = 0; i <= M; i++) {
                        System.arraycopy(base[i], 0, copy[i], 0, R);
                    }
                    for (int i = 0; i <= M; i++) {
                        copy[i][r] = alloc[i];
                    }
                    next.add(copy);
                }
            }
            partial = next;
        }
        List<State> out = new ArrayList<State>();
        double tot = 0.0;
        for (int[][] n : partial) {
            State s = new State(n);
            double w = 1.0;
            for (int i = 0; i < M; i++) {
                int sum = 0;
                for (int r = 0; r < R; r++) {
                    sum += n[i][r];
                }
                w *= fact(sum);
                for (int r = 0; r < R; r++) {
                    w *= Math.pow(L[i][r], n[i][r]) / fact(n[i][r]);
                }
            }
            for (int r = 0; r < R; r++) {
                w *= Math.pow(Z[r], n[M][r]) / fact(n[M][r]);
            }
            s.p = w;
            tot += w;
            out.add(s);
        }
        for (State s : out) {
            s.p /= tot;
        }
        return out;
    }

    /**
     * All ways of splitting n items into k nonnegative parts.
     *
     * @param n total
     * @param k parts
     * @return the compositions
     */
    private static List<int[]> compositions(int n, int k) {
        List<int[]> out = new ArrayList<int[]>();
        int[] cur = new int[k];
        comp(n, k, 0, cur, out);
        return out;
    }

    private static void comp(int n, int k, int pos, int[] cur, List<int[]> out) {
        if (pos == k - 1) {
            cur[pos] = n;
            out.add(cur.clone());
            return;
        }
        for (int v = 0; v <= n; v++) {
            cur[pos] = v;
            comp(n - v, k, pos + 1, cur, out);
        }
    }

    private static double fact(int n) {
        double f = 1.0;
        for (int t = 2; t <= n; t++) {
            f *= t;
        }
        return f;
    }

    private static Matrix mat(double[][] v) {
        Matrix m = new Matrix(v.length, v[0].length);
        for (int i = 0; i < v.length; i++) {
            for (int j = 0; j < v[0].length; j++) {
                m.set(i, j, v[i][j]);
            }
        }
        return m;
    }

    private static Matrix row(double... v) {
        return mat(new double[][]{v});
    }

    /**
     * Checks the mean vector and the covariance matrix of the selected
     * coordinates against the enumerated distribution.
     *
     * @param L demand matrix
     * @param N population vector
     * @param Z think time vector
     * @param pairs selected 0-based (station,class) coordinates
     * @param out the routine's result
     */
    private static void assertAgainstOracle(double[][] L, int[] N, double[] Z, int[][] pairs,
                                            Ret.pfqnQlenMoments out) {
        List<State> states = oracleStates(L, N, Z);
        int d = pairs.length;
        double[] mean = new double[d];
        for (int j = 0; j < d; j++) {
            for (State s : states) {
                mean[j] += s.p * s.n[pairs[j][0]][pairs[j][1]];
            }
            assertTrue(Math.abs(out.mean.get(j) - mean[j]) < TOL,
                    "mean " + j + ": " + out.mean.get(j) + " vs " + mean[j]);
        }
        for (int j = 0; j < d; j++) {
            for (int l = 0; l < d; l++) {
                double cov = 0.0;
                for (State s : states) {
                    cov += s.p * (s.n[pairs[j][0]][pairs[j][1]] - mean[j])
                            * (s.n[pairs[l][0]][pairs[l][1]] - mean[l]);
                }
                assertTrue(Math.abs(out.cov.get(j, l) - cov) < TOL,
                        "cov " + j + "," + l + ": " + out.cov.get(j, l) + " vs " + cov);
            }
        }
    }

    // ---------- the tail vertex ------------------------------------------

    @Test
    public void testTailEdgeOnADeterministicVariable() {
        // N = 3 with probability one: t = (1,1,1,1,0) and b_j = C(3,j)
        double[] t = {1, 1, 1, 1, 0};
        Matrix b = moment_binomial_from_tail(row(t).transpose());
        double[] exp = {1, 3, 3, 1, 0};
        for (int j = 0; j < 5; j++) {
            assertTrue(Math.abs(b.get(j) - exp[j]) < TOL, "order " + j);
        }
        Matrix back = moment_tail_from_binomial(b);
        for (int j = 0; j < 5; j++) {
            assertTrue(Math.abs(back.get(j) - t[j]) < TOL, "inverse " + j);
        }
    }

    @Test
    public void testJointTailEdgeAndCentralMoments() {
        // a dependent bivariate law on {0,..,3}^2 given by its pmf
        int n = 4;
        double[] pmf = new double[n * n];
        double tot = 0.0;
        for (int u = 0; u < n; u++) {
            for (int v = 0; v < n; v++) {
                pmf[u * n + v] = 1.0 + u + 2.0 * v + u * v;
                tot += pmf[u * n + v];
            }
        }
        for (int i = 0; i < pmf.length; i++) {
            pmf[i] /= tot;
        }
        double[] tail = new double[n * n];
        for (int a = 0; a < n; a++) {
            for (int b = 0; b < n; b++) {
                double acc = 0.0;
                for (int u = a; u < n; u++) {
                    for (int v = b; v < n; v++) {
                        acc += pmf[u * n + v];
                    }
                }
                tail[a * n + b] = acc;
            }
        }
        int[] dims = {n, n};
        double[] bin = moment_joint_binomial_from_tail(tail, dims);
        for (int a = 0; a < n; a++) {
            for (int b = 0; b < n; b++) {
                double exp = 0.0;
                for (int u = 0; u < n; u++) {
                    for (int v = 0; v < n; v++) {
                        exp += binom(u, a) * binom(v, b) * pmf[u * n + v];
                    }
                }
                assertTrue(Math.abs(bin[a * n + b] - exp) < TOL, "multi-order " + a + "," + b);
            }
        }
        double[] back = moment_joint_tail_from_binomial(bin, dims);
        for (int i = 0; i < tail.length; i++) {
            assertTrue(Math.abs(back[i] - tail[i]) < TOL, "inverse " + i);
        }
        double[] mc = moment_joint_central_from_tail(tail, dims);
        double mu0 = 0.0;
        double mu1 = 0.0;
        for (int u = 0; u < n; u++) {
            for (int v = 0; v < n; v++) {
                mu0 += u * pmf[u * n + v];
                mu1 += v * pmf[u * n + v];
            }
        }
        double cov = 0.0;
        for (int u = 0; u < n; u++) {
            for (int v = 0; v < n; v++) {
                cov += pmf[u * n + v] * (u - mu0) * (v - mu1);
            }
        }
        assertTrue(Math.abs(mc[1 * n + 1] - cov) < TOL, "covariance");
    }

    private static double binom(int n, int k) {
        if (k < 0 || k > n) {
            return 0.0;
        }
        double c = 1.0;
        for (int i = 0; i < k; i++) {
            c = c * (n - i) / (i + 1);
        }
        return c;
    }

    // ---------- the queueing routine -------------------------------------

    @Test
    public void testSingleClassTailRoute() {
        double[][] L = {{2.0}, {1.0}};
        int[] N = {5};
        double[] Z = {0.7};
        int[][] pairs = {{0, 0}, {1, 0}};
        Ret.pfqnQlenMoments out = pfqn_qlen_joint_moments(mat(L), row(5), row(Z), pairs);
        assertTrue("tail".equals(out.route), "route");
        assertAgainstOracle(L, N, Z, pairs, out);
        // the third central moment of the first coordinate
        List<State> states = oracleStates(L, N, Z);
        double m3 = 0.0;
        double mu = out.mean.get(0);
        for (State s : states) {
            m3 += s.p * Math.pow(s.n[0][0] - mu, 3);
        }
        assertTrue(Math.abs(out.central[3 * out.dims[1]] - m3) < 1e-8, "third central moment");
    }

    @Test
    public void testMulticlassPmfRoute() {
        // one station, two classes: the shape a class-oriented method of
        // moments produces. The naive survival formula fails here, so this
        // exercises the complementary-network route
        double[][] L = {{2.0, 1.0}};
        int[] N = {4, 3};
        double[] Z = {0.5, 0.8};
        int[][] pairs = {{0, 0}, {0, 1}};
        Ret.pfqnQlenMoments out = pfqn_qlen_joint_moments(mat(L), row(4, 3), row(Z), pairs);
        assertTrue("pmf".equals(out.route), "route");
        assertAgainstOracle(L, N, Z, pairs, out);
    }

    @Test
    public void testCrossStationAndCrossClassPairs() {
        double[][] L = {{2.0, 1.0}, {1.0, 3.0}};
        int[] N = {3, 2};
        double[] Z = {0.5, 0.0};
        int[][] cross = {{0, 0}, {1, 1}};
        assertAgainstOracle(L, N, Z, cross,
                pfqn_qlen_joint_moments(mat(L), row(3, 2), row(Z), cross));
        int[][] all = {{0, 0}, {0, 1}, {1, 0}, {1, 1}};
        assertAgainstOracle(L, N, Z, all,
                pfqn_qlen_joint_moments(mat(L), row(3, 2), row(Z), all));
    }

    @Test
    public void testBothRoutesAgreeWhenBothApply() {
        // with one class the two routes are both valid and touch different
        // networks, so agreement is a genuine cross-check
        double[][] L = {{2.0}, {1.0}, {0.5}};
        int[][] pairs = {{0, 0}, {2, 0}};
        Ret.pfqnQlenMoments a = pfqn_qlen_joint_moments(mat(L), row(4), row(0.3), pairs,
                "tail", null, null);
        Ret.pfqnQlenMoments b = pfqn_qlen_joint_moments(mat(L), row(4), row(0.3), pairs,
                "pmf", null, null);
        for (int i = 0; i < a.tail.length; i++) {
            assertTrue(Math.abs(a.tail[i] - b.tail[i]) < TOL, "tail " + i);
            assertTrue(Math.abs(a.central[i] - b.central[i]) < 1e-8, "central " + i);
        }
    }

    @Test
    public void testInjectedSourceIsUsedAndCounted() {
        double[][] L = {{2.0}, {1.0}};
        int[][] pairs = {{0, 0}, {1, 0}};
        final int[] calls = new int[1];
        Ret.pfqnQlenMoments ref = pfqn_qlen_joint_moments(mat(L), row(4), row(0.6), pairs);
        Ret.pfqnQlenMoments got = pfqn_qlen_joint_moments(mat(L), row(4), row(0.6), pairs,
                "tail", (Lsub, pops) -> {
                    calls[0]++;
                    Matrix lg = new Matrix(pops.getNumRows(), 1);
                    for (int p = 0; p < pops.getNumRows(); p++) {
                        // serve the even populations only, the rest must fall back
                        if (((int) Math.round(pops.get(p, 0))) % 2 == 0) {
                            Matrix n = new Matrix(1, 1);
                            n.set(0, 0, pops.get(p, 0));
                            // the same method the routine uses for the rows it
                            // has to fall back on, otherwise the two halves of
                            // the array would come from different algorithms
                            jline.solvers.nc.NCOptions opt = new jline.solvers.nc.NCOptions();
                            opt.method = "exact";
                            lg.set(p, 0, jline.api.pfqn.nc.Pfqn_nc.pfqn_nc(new Matrix(1, 1),
                                    Lsub, n, row(0.6), opt).lG);
                        } else {
                            lg.set(p, 0, Double.NaN);
                        }
                    }
                    return lg;
                }, null);
        assertTrue(calls[0] == 1, "the source must be called exactly once");
        assertTrue(got.served > 0 && got.served < got.points, "partial service");
        assertTrue(got.evals == got.points - got.served, "the rest must fall back");
        for (int i = 0; i < ref.central.length; i++) {
            assertTrue(Math.abs(ref.central[i] - got.central[i]) < 1e-8, "central " + i);
        }
    }

    @Test
    public void testErrorPaths() {
        Matrix L = mat(new double[][]{{2.0, 1.0}});
        Matrix N = row(3, 2);
        Matrix Z = row(0.1, 0.1);
        assertThrows(IllegalArgumentException.class,
                () -> pfqn_qlen_joint_moments(L, N, Z, null, "tail", null, null));
        assertThrows(IllegalArgumentException.class,
                () -> pfqn_qlen_joint_moments(L, N, Z, new int[][]{{0, 0}, {0, 0}}));
        assertThrows(IllegalArgumentException.class,
                () -> pfqn_qlen_joint_moments(L, N, Z, new int[][]{{1, 0}}));
        assertThrows(IllegalArgumentException.class,
                () -> pfqn_qlen_joint_moments(L, row(3), Z, new int[][]{{0, 0}}));
        assertThrows(IllegalArgumentException.class,
                () -> pfqn_qlen_joint_moments(L, N, Z, null, "nosuch", null, null));
    }
}
