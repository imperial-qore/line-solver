/**
 * @file Single-class LIMITED load-dependent normalizing constant auxiliary computation
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.ld;

import java.util.Arrays;

import org.apache.commons.math3.util.FastMath;

import jline.io.Ret;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;

/**
 * Same recursion, same arithmetic and bit-identical results to
 * {@link Pfqn_gldsingle}, but with the rate-offset axis truncated at the
 * limited load-dependence threshold instead of at the population.
 *
 * <p>Unrolling the recursion of Pfqn_gldsingle,
 *
 * <pre>
 *   g(m,n,t) = g(m-1,n,1) + L_m/mu(m,t) * g(m,n-1,t+1)
 * </pre>
 *
 * shows that its third index is an offset into station m's rate function,
 *
 * <pre>
 *   g(m,n,t) = sum_{j=0..n} prod_{i=0..j-1} L_m/alpha_m(t+i) * g(m-1,n-j,1)
 * </pre>
 *
 * so once t &gt;= s_m, where s_m is the population past which alpha_m stays
 * constant, every factor is alpha_m(s_m), the product collapses to
 * (L_m/alpha_m(s_m))^j and
 *
 * <pre>
 *   g(m,n,t) = g(m,n,s_m)   for all t &gt;= s_m
 * </pre>
 *
 * The N-s_m upper slices that Pfqn_gldsingle computes are duplicates of one
 * another. Capping the offset at s_m and reading g(m,n-1,min(t+1,s_m)) keeps
 * every value the answer needs.
 *
 * <p>COST. O(N * sum_k s_k) time against O(M N^2), and O(N * max_k s_k) space
 * against the dense map Pfqn_gldsingle builds, the levels being rolled. On a
 * multiserver model, where s_k is the server count, this is LINEAR in the
 * population rather than quadratic. The log-domain branch remains a sum of
 * nonnegative terms and loses no digits to cancellation, unlike the closed form
 * of Pfqn_explicit_ld, which reaches the same asymptotics through an
 * alternating sum.
 *
 * <p>There is no gain on a station whose rates never settle, an infinite server
 * alpha(n)=n being the usual case: it gets s_k = N and costs what it costs in
 * Pfqn_gldsingle. The saving is over the OTHER stations, so a model carrying
 * one delay among M queues drops from O(M N^2) to O(N^2 + N sum_k s_k).
 *
 * <p>The threshold is detected per station rather than declared, so an
 * arbitrary rate matrix is accepted and simply yields s_k = N, at which point
 * this is Pfqn_gldsingle with its slices rolled.
 */
public final class Pfqn_lldsingle {
    private Pfqn_lldsingle() {}

    /**
     * Auxiliary function used by Pfqn_ncld and Pfqn_nre to compute the normalizing
     * constant in a single-class limited load-dependent model.
     *
     * @param L       demands at all stations
     * @param N       number of jobs for each class
     * @param mu      load-dependent scaling factors
     * @param options solver options
     * @return normalizing constant (G) and its logarithm (lG)
     */
    public static Ret.pfqnNc pfqn_lldsingle(Matrix L, Matrix N, Matrix mu, SolverOptions options) {
        int M = L.getNumRows();
        int R = L.getNumCols();

        if (R > 1) {
            throw new RuntimeException("pfqn_lldsingle: multiclass model detected. pfqn_lldsingle is for single class models.");
        }

        int Nt = (int) N.get(0);
        if (Nt <= 0) {
            // the empty product, as Pfqn_gldsingle returns from its unexecuted loops
            return new Ret.pfqnNc(1.0, 0.0);
        }

        int ncol = Math.min(mu.getNumCols(), Nt);
        int[] s = thresholds(mu, M, Nt, ncol);

        // see _kb/03-api-layer.md for rationale
        boolean useLog = true;
        for (int m = 0; m < M && useLog; m++) {
            if (!(L.get(m) >= 0)) useLog = false;
        }
        for (int m = 0; m < mu.getNumRows() && useLog; m++) {
            for (int c = 0; c < mu.getNumCols() && useLog; c++) {
                if (!(mu.get(m, c) > 0)) useLog = false;
            }
        }

        if (useLog) {
            // see _kb/03-api-layer.md for rationale
            double[] lgprev = new double[Nt + 1];
            Arrays.fill(lgprev, Double.NEGATIVE_INFINITY);
            lgprev[0] = 0.0; // g(0,0,1)=1, and g(0,n,1)=0 for n>=1: no station holds them
            for (int m = 1; m <= M; m++) {
                int sm = s[m - 1];
                double lL = FastMath.log(L.get(m - 1)); // -Inf where the demand is zero
                double[] cur = new double[(Nt + 1) * sm];
                Arrays.fill(cur, Double.NEGATIVE_INFINITY);
                for (int tm = 0; tm < sm; tm++) cur[tm] = 0.0; // log(1): zero jobs
                for (int n = 1; n <= Nt; n++) {
                    // offsets above Nt-n+1 are never read back, exactly as in
                    // Pfqn_gldsingle, so the triangle is kept
                    int tmax = Math.min(sm, Nt - n + 1);
                    for (int tm = 1; tm <= tmax; tm++) {
                        double a = lgprev[n];
                        int tsrc = Math.min(tm + 1, sm);
                        // +Inf rate zeroes the term, matching a division by Inf
                        double b = lL + cur[(n - 1) * sm + (tsrc - 1)] - FastMath.log(rate(mu, m - 1, tm - 1, ncol));
                        cur[n * sm + (tm - 1)] = logSumExp2(a, b);
                    }
                }
                for (int n = 0; n <= Nt; n++) lgprev[n] = cur[n * sm];
            }
            double lG = lgprev[Nt];
            return new Ret.pfqnNc(FastMath.exp(lG), lG);
        }

        double[] gprev = new double[Nt + 1];
        gprev[0] = 1.0; // g(0,0,1)=1, and g(0,n,1)=0 for n>=1
        for (int m = 1; m <= M; m++) {
            int sm = s[m - 1];
            double[] cur = new double[(Nt + 1) * sm];
            for (int tm = 0; tm < sm; tm++) cur[tm] = 1.0; // g(m,0,t)=1
            for (int n = 1; n <= Nt; n++) {
                int tmax = Math.min(sm, Nt - n + 1);
                for (int tm = 1; tm <= tmax; tm++) {
                    int tsrc = Math.min(tm + 1, sm);
                    cur[n * sm + (tm - 1)] = gprev[n]
                            + L.get(m - 1) * cur[(n - 1) * sm + (tsrc - 1)] / rate(mu, m - 1, tm - 1, ncol);
                }
            }
            for (int n = 0; n <= Nt; n++) gprev[n] = cur[n * sm];
        }
        double G = gprev[Nt];
        double lG = FastMath.log(G);
        return new Ret.pfqnNc(G, lG);
    }

    /**
     * Per-station limited load-dependence threshold: the smallest offset past
     * which the rate row is constant, so that alpha_k(n) = alpha_k(s_k) for
     * every n &gt;= s_k. Equality is tested first so that an infinite rate,
     * which the recursion admits and zeroes through log(mu) = +Inf, ties with
     * itself instead of producing Inf-Inf. The tolerance is then confined to
     * FINITE pairs: at tail = Inf the bound ulp*max(abs(tail),1) is itself Inf
     * and abs(prev-Inf) &lt;= Inf would tie every finite rate to it, collapsing
     * the row on a false tie. A MISSED tie only costs time; a FALSE tie would
     * be a wrong answer, hence the strict tolerance, which follows
     * Pfqn_explicit_ld's scan.
     */
    private static int[] thresholds(Matrix mu, int M, int Nt, int ncol) {
        int[] s = new int[M];
        for (int m = 0; m < M; m++) {
            double tail = rate(mu, m, Nt - 1, ncol);
            s[m] = Nt;
            for (int n = Nt - 1; n >= 1; n--) {
                double prev = rate(mu, m, n - 1, ncol);
                if (prev == tail || (Double.isFinite(tail) && Double.isFinite(prev)
                        && Math.abs(prev - tail) <= Math.ulp(1.0) * Math.max(Math.abs(tail), 1.0))) {
                    s[m] = n;
                } else {
                    break;
                }
            }
            if (s[m] < 1) s[m] = 1;
        }
        return s;
    }

    /** Rates past the last column of mu are 1, as in the scalar recursion. */
    private static double rate(Matrix mu, int m, int c, int ncol) {
        return c < ncol ? mu.get(m, c) : 1.0;
    }

    /**
     * Pairwise log-sum-exp, stable when either argument is -Inf.
     */
    private static double logSumExp2(double a, double b) {
        if (a > b) {
            if (b == Double.NEGATIVE_INFINITY) return a;
            return a + FastMath.log1p(FastMath.exp(b - a));
        }
        if (a == Double.NEGATIVE_INFINITY) return b;
        return b + FastMath.log1p(FastMath.exp(a - b));
    }
}
