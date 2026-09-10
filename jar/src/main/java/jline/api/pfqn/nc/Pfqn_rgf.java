/**
 * @file Recursion by Generating Functions (RGF) for product-form normalizing constants
 *
 * Exact normalizing constant of a single-class closed product-form network obtained by
 * convolving the per-node generating-function sequences of Coury and Harrison (1997),
 * Property 1, instead of the per-station Buzen recursion. Ported at parity from MATLAB
 * pfqn_rgf.m.
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

public final class Pfqn_rgf {
    private Pfqn_rgf() {}

    /**
     * Result of the RGF recursion: the whole log-normalizing-constant sequence.
     */
    public static final class Result {
        /** Normalizing constant G(N). */
        public final double G;
        /** Logarithm of G(N). */
        public final double lG;
        /** Logarithms of g(0), g(1), ..., g(N). */
        public final double[] lg;

        public Result(double G, double lG, double[] lg) {
            this.G = G;
            this.lG = lG;
            this.lg = lg;
        }
    }

    /**
     * Exact normalizing constant by convolution of per-node generating-function sequences.
     *
     * <p>A GROUP of m stations sharing the same demand p collapses into the single
     * negative-binomial sequence r(k) = C(k+m-1,k) p^k, so the whole group costs one
     * sequence rather than m convolution passes; the delay contributes the Poisson
     * sequence Z^k/k!. Cost O(G N^2) against Buzen's O(M N), so RGF is the cheaper route
     * on heavily replicated models with moderate populations (G N &lt; M). The recursion
     * runs entirely in the log domain, so no intermediate overflow or underflow is
     * possible.</p>
     *
     * @param L service demand vector (M x 1) of the queueing stations
     * @param N population (nonnegative integer)
     * @param Z think time
     * @return the normalizing constant, its logarithm, and the full log g(0..N) sequence
     */
    public static Result pfqn_rgf(Matrix L, double N, double Z) {
        if (N < 0 || N != FastMath.rint(N)) {
            throw new IllegalArgumentException("pfqn_rgf requires a nonnegative integer population.");
        }
        if (Z < 0) {
            throw new IllegalArgumentException("pfqn_rgf requires a nonnegative think time.");
        }
        int n = (int) FastMath.rint(N);

        double[] lg = new double[n + 1];
        Arrays.fill(lg, Double.NEGATIVE_INFINITY);
        lg[0] = 0.0;

        double[] kk = new double[n + 1];
        for (int k = 0; k <= n; k++) {
            kk[k] = k;
        }

        if (Z > 0) {
            double[] lr = new double[n + 1];
            for (int k = 0; k <= n; k++) {
                lr[k] = kk[k] * FastMath.log(Z) - Gamma.logGamma(kk[k] + 1);
            }
            lg = logconv(lg, lr);
        }

        List<Double> pos = new ArrayList<Double>();
        for (int i = 0; i < L.length(); i++) {
            double d = L.get(i);
            if (d < 0) {
                throw new IllegalArgumentException("pfqn_rgf requires nonnegative demands.");
            }
            if (d > 0) {
                pos.add(d);
            }
        }
        double[] p = new double[pos.size()];
        for (int i = 0; i < pos.size(); i++) {
            p[i] = pos.get(i);
        }
        Arrays.sort(p);

        int i = 0;
        while (i < p.length) {
            int m = 1;
            while (i + m < p.length && p[i + m] == p[i]) {
                m++;
            }
            double[] lr = new double[n + 1];
            if (m == 1) {
                for (int k = 0; k <= n; k++) {
                    lr[k] = kk[k] * FastMath.log(p[i]);
                }
            } else {
                for (int k = 0; k <= n; k++) {
                    lr[k] = Gamma.logGamma(kk[k] + m) - Gamma.logGamma(kk[k] + 1)
                            - Gamma.logGamma(m) + kk[k] * FastMath.log(p[i]);
                }
            }
            lg = logconv(lg, lr);
            i += m;
        }

        double lG = lg[n];
        return new Result(FastMath.exp(lG), lG, lg);
    }

    /** Log-domain linear convolution truncated at the common length. */
    private static double[] logconv(double[] u, double[] v) {
        int n = u.length;
        double[] c = new double[n];
        for (int k = 0; k < n; k++) {
            double vm = Double.NEGATIVE_INFINITY;
            for (int j = 0; j <= k; j++) {
                double t = u[j] + v[k - j];
                if (t > vm) {
                    vm = t;
                }
            }
            if (Double.isInfinite(vm)) {
                c[k] = vm;
            } else {
                double acc = 0.0;
                for (int j = 0; j <= k; j++) {
                    acc += FastMath.exp(u[j] + v[k - j] - vm);
                }
                c[k] = vm + FastMath.log(acc);
            }
        }
        return c;
    }
}
