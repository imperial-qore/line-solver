/**
 * @file Single-class load-dependent normalizing constant auxiliary computation
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.ld;

import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math3.util.FastMath;

import jline.io.Ret;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;

public final class Pfqn_gldsingle {
    private Pfqn_gldsingle() {}

    /**
     * Auxiliary function used by pfqn_gld to compute the normalizing constant in a
     * single-class load-dependent model.
     *
     * @param L       demands at all stations
     * @param N       number of jobs for each class
     * @param mu      load-dependent scaling factors
     * @param options solver options
     * @return normalizing constant (G) and its logarithm (lG)
     */
    public static Ret.pfqnNc pfqn_gldsingle(Matrix L, Matrix N, Matrix mu, SolverOptions options) {
        int M = L.getNumRows();
        int R = L.getNumCols();

        if (R > 1) {
            throw new RuntimeException("pfqn_gldsingle: multiclass model detected. pfqn_gldsingle is for single class models.");
        }

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
            Map<Ret.pfqnGldIndex, Double> lg = new HashMap<Ret.pfqnGldIndex, Double>();
            // lg(0+1,n+1,1+1) stays -Inf for n>=1: no station can hold n>=1 jobs.
            int n = 1;
            while (n <= N.get(0)) {
                lg.put(new Ret.pfqnGldIndex(0 + 1, n + 1, 1 + 1), Double.NEGATIVE_INFINITY);
                n++;
            }
            for (int m = 1; m <= M; m++) {
                int tm = 1;
                while (tm <= N.get(0) + 1) {
                    lg.put(new Ret.pfqnGldIndex(m + 1, 0 + 1, tm + 1), 0.0); // log(1): zero jobs
                    tm++;
                }
                double lL = FastMath.log(L.get(m - 1)); // -Inf where the demand is zero
                int n2 = 1;
                while (n2 <= N.get(0)) {
                    int tm2 = 1;
                    while (tm2 <= N.get(0) - n2 + 1) {
                        Double lgPrev = lg.get(new Ret.pfqnGldIndex(m - 1 + 1, n2 + 1, 1 + 1));
                        if (lgPrev == null) lgPrev = Double.NEGATIVE_INFINITY;
                        Double lgCurr = lg.get(new Ret.pfqnGldIndex(m + 1, n2 - 1 + 1, tm2 + 1 + 1));
                        if (lgCurr == null) lgCurr = Double.NEGATIVE_INFINITY;
                        // +Inf rate zeroes the term, matching a division by Inf
                        double a = lgPrev;
                        double b = lL + lgCurr - FastMath.log(mu.get(m - 1, tm2 - 1));
                        lg.put(new Ret.pfqnGldIndex(m + 1, n2 + 1, tm2 + 1), logSumExp2(a, b));
                        tm2++;
                    }
                    n2++;
                }
            }
            Double lGv = lg.get(new Ret.pfqnGldIndex(M + 1, (int) N.get(0) + 1, 1 + 1));
            if (lGv == null) lGv = Double.NEGATIVE_INFINITY;
            double lG = lGv;
            return new Ret.pfqnNc(FastMath.exp(lG), lG);
        }

        Map<Ret.pfqnGldIndex, Double> g = new HashMap<Ret.pfqnGldIndex, Double>();

        // Initialize boundary conditions: g(0+1, n+1, 1+1) = 0 for n=1:N
        int n = 1;
        while (n <= N.get(0)) {
            g.put(new Ret.pfqnGldIndex(0 + 1, n + 1, 1 + 1), 0.0);
            n++;
        }

        for (int m = 1; m <= M; m++) {
            // Initialize boundary conditions: g(m+1, 0+1, tm+1) = 1 for tm=1:(N+1)
            int tm = 1;
            while (tm <= N.get(0) + 1) {
                g.put(new Ret.pfqnGldIndex(m + 1, 0 + 1, tm + 1), 1.0);
                tm++;
            }
            int n2 = 1;
            while (n2 <= N.get(0)) {
                int tm2 = 1;
                while (tm2 <= N.get(0) - n2 + 1) {
                    Double gPrev = g.get(new Ret.pfqnGldIndex(m - 1 + 1, n2 + 1, 1 + 1));
                    if (gPrev == null) gPrev = 0.0;
                    Double gCurr = g.get(new Ret.pfqnGldIndex(m + 1, n2 - 1 + 1, tm2 + 1 + 1));
                    if (gCurr == null) gCurr = 0.0;
                    g.put(new Ret.pfqnGldIndex(m + 1, n2 + 1, tm2 + 1),
                            gPrev + L.get(m - 1) * gCurr / mu.get(m - 1, tm2 - 1));
                    tm2++;
                }
                n2++;
            }
        }
        double G = g.get(new Ret.pfqnGldIndex(M + 1, (int) N.get(0) + 1, 1 + 1));
        double lG = FastMath.log(G);
        return new Ret.pfqnNc(G, lG);
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
