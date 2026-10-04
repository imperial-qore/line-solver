/**
 * @file Distinct-load Normalizing Constant (DNC) at a nonintegral population
 *
 * Normalizing constant and throughput of a single-class closed product-form network at a
 * REAL-VALUED population, by partial-fraction inversion of the network generating
 * function (Dowdy and Gordon 1984). Ported at parity from MATLAB pfqn_dnc.m.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.FastMath;

import jline.util.matrix.Matrix;

public final class Pfqn_dnc {
    private Pfqn_dnc() {}

    /** Result of the DNC evaluation at a real population. */
    public static final class Result {
        /** Throughput G(N-1)/G(N); NaN when the population is not positive. */
        public final double X;
        /** Normalizing constant at population N. */
        public final double G;
        /** Logarithm of the normalizing constant. */
        public final double lG;

        public Result(double X, double G, double lG) {
            this.X = X;
            this.G = G;
            this.lG = lG;
        }
    }

    /**
     * Normalizing constant and throughput at a real-valued population.
     *
     * <p>With distinct loads x_1..x_G of multiplicities m_1..m_G the generating function
     * prod_g (1-x_g u)^{-m_g} expands as
     * G(n) = sum_g sum_{j=1..m_g} A_gj C(n+j-1,j-1) x_g^n, every term of which is analytic
     * in n, so evaluating at a real n interpolates the integral normalizing constants
     * exactly and gives a smooth throughput curve through the integral points. For
     * all-distinct loads A_g = prod_{l!=g} x_g/(x_g - x_l) is used directly; with repeated
     * loads the coefficients are recovered from G(0..M-1).</p>
     *
     * <p>Only the queueing part admits this continuation: the delay sequence Z^n/n! is
     * entire and has no partial-fraction expansion, so a think time is not accepted here.
     * Use {@link jline.api.pfqn.mva.Pfqn_nintmva} for nonintegral populations with a
     * delay.</p>
     *
     * @param L service demand vector (M x 1) of the queueing stations
     * @param N population (real, nonnegative; may be fractional)
     * @return throughput, normalizing constant and its logarithm
     */
    public static Result pfqn_dnc(Matrix L, double N) {
        if (N < 0) {
            throw new IllegalArgumentException("pfqn_dnc requires a nonnegative population.");
        }
        List<Double> pos = new ArrayList<Double>();
        for (int i = 0; i < L.length(); i++) {
            if (L.get(i) > 0) {
                pos.add(L.get(i));
            }
        }
        if (pos.isEmpty()) {
            throw new IllegalArgumentException("pfqn_dnc requires at least one station with positive demand.");
        }
        int M = pos.size();
        double[] y = new double[M];
        double xmax = 0.0;
        for (int i = 0; i < M; i++) {
            xmax = FastMath.max(xmax, pos.get(i));
        }
        for (int i = 0; i < M; i++) {
            y[i] = pos.get(i) / xmax;
        }

        // Distinct loads merged under a relative tolerance, so numerically coincident
        // loads go to the multiplicity branch rather than to a near-singular denominator.
        double[] ys = y.clone();
        Arrays.sort(ys);
        List<Double> uL = new ArrayList<Double>();
        List<Integer> multL = new ArrayList<Integer>();
        uL.add(ys[M - 1]);
        multL.add(1);
        for (int i = M - 2; i >= 0; i--) {
            double last = uL.get(uL.size() - 1);
            if (ys[i] > last * (1 - 1e-9)) {
                multL.set(multL.size() - 1, multL.get(multL.size() - 1) + 1);
            } else {
                uL.add(ys[i]);
                multL.add(1);
            }
        }
        int Gd = uL.size();
        double[] u = new double[Gd];
        int[] mult = new int[Gd];
        boolean allSimple = true;
        for (int g = 0; g < Gd; g++) {
            u[g] = uL.get(g);
            mult[g] = multL.get(g);
            if (mult[g] > 1) {
                allSimple = false;
            }
        }

        double[] A;
        int[] node;
        double[] j;
        if (allSimple) {
            A = new double[Gd];
            node = new int[Gd];
            j = new double[Gd];
            for (int g = 0; g < Gd; g++) {
                double prod = 1.0;
                for (int l = 0; l < Gd; l++) {
                    if (l != g) {
                        prod *= u[g] / (u[g] - u[l]);
                    }
                }
                A[g] = prod;
                node[g] = g;
                j[g] = 1.0;
            }
        } else {
            // Recover the coefficients from G(0..M-1), computed by convolution on the
            // scaled loads (bounded by construction since max u = 1).
            double[] gint = new double[M];
            gint[0] = 1.0;
            for (int i = 1; i < M; i++) {
                gint[i] = 0.0;
            }
            for (int i = 0; i < M; i++) {
                double[] next = new double[M];
                double pw = 1.0;
                for (int k = 0; k < M; k++) {
                    for (int q = 0; k + q < M; q++) {
                        next[k + q] += gint[q] * pw;
                    }
                    pw *= y[i];
                }
                gint = next;
            }
            node = new int[M];
            j = new double[M];
            int c = 0;
            for (int g = 0; g < Gd; g++) {
                for (int jj = 1; jj <= mult[g]; jj++) {
                    node[c] = g;
                    j[c] = jj;
                    c++;
                }
            }
            Matrix F = new Matrix(M, M);
            for (int k = 0; k < M; k++) {
                for (int q = 0; q < M; q++) {
                    F.set(k, q, FastMath.exp(Gamma.logGamma(k + j[q]) - Gamma.logGamma(j[q])
                            - Gamma.logGamma(k + 1.0) + k * FastMath.log(u[node[q]])));
                }
            }
            Matrix b = new Matrix(M, 1);
            for (int k = 0; k < M; k++) {
                b.set(k, 0, gint[k]);
            }
            Matrix sol = new Matrix(M, 1);
            Matrix.solve(F, b, sol);
            A = new double[M];
            for (int k = 0; k < M; k++) {
                A[k] = sol.get(k, 0);
            }
        }

        double GN = dncEval(N, A, u, node, j);
        double GN1 = dncEval(N - 1, A, u, node, j);

        double lG = FastMath.log(GN) + N * FastMath.log(xmax);
        double X;
        if (N <= 0 || GN <= 0 || Double.isNaN(GN1)) {
            X = Double.NaN;
        } else {
            X = (GN1 / GN) / xmax;
        }
        return new Result(X, FastMath.exp(lG), lG);
    }

    /**
     * Partial-fraction series evaluated at a real population n. The continuation is
     * analytic for n &gt; -1; below that the binomial factor changes sign and the
     * log-domain evaluation would lose it, so it is not extended there.
     */
    private static double dncEval(double n, double[] A, double[] u, int[] node, double[] j) {
        if (n <= -1) {
            return Double.NaN;
        }
        double acc = 0.0;
        for (int q = 0; q < A.length; q++) {
            acc += A[q] * FastMath.exp(Gamma.logGamma(n + j[q]) - Gamma.logGamma(j[q])
                    - Gamma.logGamma(n + 1.0) + n * FastMath.log(u[node[q]]));
        }
        return acc;
    }
}
