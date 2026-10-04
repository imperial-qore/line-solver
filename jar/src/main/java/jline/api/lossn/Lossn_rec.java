package jline.api.lossn;

import java.util.ArrayList;
import java.util.List;

import jline.api.mdd.MDD;
import jline.api.mdd.MddNextState;
import jline.api.mdd.MddStruct;
import jline.api.mdd.Mdd_rec;
import jline.api.mdd.Mdd_reachset;
import jline.util.matrix.Matrix;

/**
 * Exact analysis of a loss network by MDD-rec: the normalising constant is the
 * sum of a product form over the admissible set {n &gt;= 0 : A n &lt;= C}, which
 * is what a decision diagram holding that set computes in one memoised walk.
 *
 * <p>A Kelly loss network carries offered load nu_r on route r and admits a call
 * only while the resource constraint A n &lt;= C still holds after it. The
 * stationary law is the truncation of independent Poisson counts to that set,
 * P(n) = (1/G) prod_r nu_r^(n_r)/n_r!, so g_r(k) = nu_r^k/k! and
 * {@link Mdd_rec} returns G. By PASTA the acceptance probability of a class-r
 * call is the ratio of two such constants, 1 - B_r = G(C - A e_r) / G(C), which
 * is one further diagram per class.</p>
 *
 * <p><b>Why this exists alongside {@link Lossn_manjunath}.</b> The
 * Manjunath-Sikdar transform evaluates G exactly as a multidimensional residue,
 * and the residue argument counts WHOLE UNITS: it needs an integral A and C. On
 * a region declaring a fractional class size or capacity the analyzer had no
 * exact route at all and fell back to the Erlang fixed point, an approximation.
 * MDD-rec needs only that the admissible set be finite and bounded coordinate by
 * coordinate, which a fractional constraint still is, so it is exact there too.
 * It is also an exact alternative to the Monte Carlo summation
 * {@link Lossn_mci} estimates.</p>
 *
 * <p>References: F. P. Kelly, "Loss networks", Annals of Applied Probability
 * 1(3), 1991. S. Balsamo, A. Marin, I. Stojic, "Computation of the normalising
 * constant for product-form models of distributed systems with synchronisation",
 * Future Generation Computer Systems 111 (2020) 475-490.</p>
 */
public class Lossn_rec {

    private Lossn_rec() {
    }

    /** Carried load, blocking, log normalising constant and walk count. */
    public static class LossnRecResult {
        /** Mean number of class-r calls in progress, the carried load. */
        public double[] QLen;
        /** Blocking probability per class. */
        public double[] Loss;
        /** log of the normalising constant G(C). */
        public double lG;
        /** Number of diagram walks performed, K + 1. */
        public int niter;
    }

    /**
     * Exact loss-network analysis by MDD-rec.
     *
     * @param nu offered load per class, length K
     * @param A  J x K non-negative resource requirement matrix
     * @param C  capacity vector, length J
     */
    public static LossnRecResult lossn_rec(Matrix nu, Matrix A, Matrix C) {
        int K = nu.getNumElements();
        int J = A.getNumRows();
        if (A.getNumCols() != K) {
            throw new RuntimeException("lossn_rec: A has " + A.getNumCols() + " columns but there "
                    + "are " + K + " classes");
        }
        if (C.getNumElements() != J) {
            throw new RuntimeException("lossn_rec: A has " + J + " rows but C has "
                    + C.getNumElements() + " entries");
        }
        double[][] Aa = new double[J][K];
        for (int j = 0; j < J; j++) {
            for (int r = 0; r < K; r++) {
                Aa[j][r] = A.get(j, r);
                if (Aa[j][r] < 0) {
                    throw new RuntimeException("lossn_rec: the resource matrix A must be "
                            + "non-negative");
                }
            }
        }
        double[] Cv = new double[J];
        for (int j = 0; j < J; j++) {
            Cv[j] = C.get(j);
        }

        // ---- per-class bound: the most calls the tightest constraint alone admits
        int[] bound = new int[K];
        for (int r = 0; r < K; r++) {
            double b = Double.POSITIVE_INFINITY;
            for (int j = 0; j < J; j++) {
                if (Aa[j][r] > 0) {
                    b = Math.min(b, Math.floor(Cv[j] / Aa[j][r]));
                }
            }
            if (Double.isInfinite(b)) {
                throw new RuntimeException("lossn_rec: class " + (r + 1) + " consumes no resource, "
                        + "so the admissible set is unbounded in that coordinate and its "
                        + "normalising constant diverges");
            }
            bound[r] = (int) Math.max(0, b);
        }

        double[][] g = new double[K][];
        for (int r = 0; r < K; r++) {
            g[r] = new double[bound[r] + 1];
            double fact = 1;
            for (int k = 0; k <= bound[r]; k++) {
                if (k > 0) {
                    fact *= k;
                }
                g[r][k] = Math.pow(nu.get(r), k) / fact;
            }
        }

        double lG = logG(Aa, Cv, bound, g);
        if (Double.isInfinite(lG)) {
            throw new RuntimeException("lossn_rec: the admissible set is empty: no call of any "
                    + "class fits within C");
        }

        // ---- carried load per class, from the marginals of the same diagram
        MddStruct mdds = diagram(Aa, Cv, bound);
        double G = Math.exp(lG);
        double[] qlen = new double[K];
        for (int r = 0; r < K; r++) {
            double[] pk = Mdd_rec.mdd_rec_marginal(mdds, g, r);
            double s = 0;
            for (int k = 0; k < pk.length; k++) {
                s += k * pk[k] / G;
            }
            qlen[r] = s;
        }

        // ---- blocking: 1 - B_r = G(C - A e_r)/G(C), Kelly's ratio, by PASTA
        double[] loss = new double[K];
        for (int r = 0; r < K; r++) {
            double[] Cr = new double[J];
            boolean fits = true;
            for (int j = 0; j < J; j++) {
                Cr[j] = Cv[j] - Aa[j][r];
                if (Cr[j] < 0) {
                    fits = false;
                }
            }
            if (!fits) {
                loss[r] = 1;                               // the call never fits
                continue;
            }
            double lGr = logG(Aa, Cr, bound, g);
            loss[r] = Double.isInfinite(lGr) ? 1 : 1 - Math.exp(lGr - lG);
            loss[r] = Math.min(1, Math.max(0, loss[r]));
        }

        LossnRecResult out = new LossnRecResult();
        out.QLen = qlen;
        out.Loss = loss;
        out.lG = lG;
        out.niter = K + 1;
        return out;
    }

    /**
     * The admissible set, generated one call at a time from the empty network.
     * Adding a call is the only move, so the breadth-first closure visits exactly
     * the admissible vectors.
     */
    private static MddStruct diagram(final double[][] A, final double[] C, final int[] bound) {
        final int K = bound.length;
        final int J = C.length;
        int[] domain = new int[K];
        for (int r = 0; r < K; r++) {
            domain[r] = bound[r] + 1;
        }
        MddNextState nextfun = new MddNextState() {
            public int[][] next(int[] state) {
                List<int[]> out = new ArrayList<int[]>();
                for (int r = 0; r < K; r++) {
                    if (state[r] >= bound[r]) {
                        continue;
                    }
                    int[] t = new int[K];
                    System.arraycopy(state, 0, t, 0, K);
                    t[r]++;
                    boolean ok = true;
                    for (int j = 0; j < J && ok; j++) {
                        double s = 0;
                        for (int q = 0; q < K; q++) {
                            s += A[j][q] * t[q];
                        }
                        ok = s <= C[j] + 1e-12;
                    }
                    if (ok) {
                        out.add(t);
                    }
                }
                return out.toArray(new int[out.size()][]);
            }
        };
        MDD mdd = Mdd_reachset.mdd_reachset(domain, new int[K], nextfun);
        return mdd.toStruct();
    }

    /**
     * log G over the admissible set at capacity C, keeping the per-class domains
     * of the FULL problem so that one set of factors g serves every reduced
     * capacity.
     */
    private static double logG(double[][] A, double[] C, int[] bound, double[][] g) {
        for (int j = 0; j < C.length; j++) {
            if (C[j] < 0) {
                return Double.NEGATIVE_INFINITY;
            }
        }
        double G = Mdd_rec.mdd_rec(diagram(A, C, bound), g);
        return G > 0 ? Math.log(G) : Double.NEGATIVE_INFINITY;
    }
}
