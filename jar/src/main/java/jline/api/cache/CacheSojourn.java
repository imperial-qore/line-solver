/**
 * @file Grid-independent sojourn-weighted average of a cache mean-field transient.
 *
 * @since LINE 3.0
 */
package jline.api.cache;

import jline.util.matrix.Matrix;
import odesolver.LSODA;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

import static jline.api.mam.Map_normalize.map_normalize;
import static jline.api.mam.Map_pie.map_pie;

/**
 * The average of a cache mean-field trajectory against a phase-type clock,
 * integrated WITH the drift. Java twin of MATLAB {@code cache_sojourn_ode.m}.
 *
 * <p>For a drift {@code dx/dt = f(x)} and a clock with row phase vector
 * {@code phi}, {@code dphi/dt = phi A}, density {@code g(t) = phi(t) c}, it
 * returns {@code xbar = int_0^T x g dt / int_0^T g dt} and {@code wtot = int_0^T
 * g dt}, by augmenting the state with {@code phi}, {@code int g x} and
 * {@code int g}. The value is accurate to the ODE tolerance and does NOT depend
 * on any output grid: a Riemann sum over the grid did, and that is what made the
 * adaptive (MATLAB) and fixed-grid (JAR/Python/C++) ENV cache mean fields differ
 * in the third digit.</p>
 */
public final class CacheSojourn {
    private CacheSojourn() {}

    /** The drift of the cache mean field. */
    public interface Drift {
        double[] apply(double[] x);
    }

    /** A phase-type clock: {@code dphi/dt = phi A}, {@code phi(0) = phi0}, density {@code phi c}. */
    public static final class Clock {
        public final double[][] A;
        public final double[] phi0;
        public final double[] c;

        public Clock(double[][] A, double[] phi0, double[] c) {
            this.A = A;
            this.phi0 = phi0;
            this.c = c;
        }
    }

    /** A sojourn average: occupancy {@code xbar} (null when {@code wtot} is 0), clock mass, per-user miss rate. */
    public static final class Result {
        public final double[] xbar;
        public final double wtot;
        /** Per-user miss rate of {@link #xbar}, filled by the policy that knows its out-of-cache map. */
        public double[] MU;

        public Result(double[] xbar, double wtot) {
            this.xbar = xbar;
            this.wtot = wtot;
        }
    }

    /**
     * The holding-time clock of a random-environment stage whose holding time is
     * the MAP {D0,D1}, in a drift time unit that is {@code lam} times real time:
     * {@code A = D0/lam}, {@code c = -D0 1/lam}, {@code phi0 = map_pie}, so that
     * {@code int_0^T g = F(T/lam)} with F the holding-time CDF.
     */
    public static Clock holdingClock(Matrix D0in, Matrix D1in, double lam) {
        Matrix D0 = D0in.copy();
        Matrix D1 = D1in.copy();
        map_normalize(D0, D1);
        int nph = D0.getNumRows();
        Matrix pie = map_pie(D0, D1);
        double[][] A = new double[nph][nph];
        double[] c = new double[nph];
        double[] phi0 = new double[nph];
        for (int i = 0; i < nph; i++) {
            double rs = 0.0;
            for (int j = 0; j < nph; j++) {
                A[i][j] = D0.get(i, j) / lam;
                rs += D0.get(i, j);
            }
            c[i] = -rs / lam;
            phi0[i] = pie.get(i);
        }
        return new Clock(A, phi0, c);
    }

    /** Integrates the augmented system over [0, t1] from {@code x0} and returns the clock-weighted average. */
    public static Result average(final Drift f, double[] x0, double t1, final Clock clk) {
        final int n = x0.length;
        final int nph = clk.phi0.length;
        final int dim = 2 * n + nph + 1;
        FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
            @Override
            public int getDimension() {
                return dim;
            }

            @Override
            public void computeDerivatives(double t, double[] z, double[] dz) {
                double[] x = new double[n];
                System.arraycopy(z, 0, x, 0, n);
                double[] d = f.apply(x);
                System.arraycopy(d, 0, dz, 0, n);
                double g = 0.0;
                for (int j = 0; j < nph; j++) {
                    double s = 0.0;
                    for (int i = 0; i < nph; i++) {
                        s += z[n + i] * clk.A[i][j];
                    }
                    dz[n + j] = s;
                    g += z[n + j] * clk.c[j];
                }
                for (int i = 0; i < n; i++) {
                    dz[n + nph + i] = g * x[i];
                }
                dz[dim - 1] = g;
            }
        };
        double[] z0 = new double[dim];
        System.arraycopy(x0, 0, z0, 0, n);
        System.arraycopy(clk.phi0, 0, z0, n, nph);
        double[] z1 = new double[dim];
        LSODA lsoda = new LSODA(1e-12, 1.0, 1e-8, 1e-10, 12, 5);
        lsoda.integrate(ode, 0.0, z0, t1, z1);
        double wtot = z1[dim - 1];
        double[] xbar = null;
        if (wtot > 0) {
            xbar = new double[n];
            for (int i = 0; i < n; i++) {
                xbar[i] = z1[n + nph + i] / wtot;
            }
        }
        return new Result(xbar, wtot);
    }

    /**
     * Per-user miss rate {@code MU(v) = sum_i lambda_v(i) pi0(i)} of a per-item
     * out-of-cache probability {@code pi0}, the same functional as the transient
     * {@code MU_t}, which is affine in the occupancy.
     */
    public static double[] missRate(Matrix[] lambdaCache, double[] pi0) {
        int u = lambdaCache.length;
        double[] mu = new double[u];
        for (int v = 0; v < u; v++) {
            double s = 0.0;
            for (int i = 0; i < pi0.length; i++) {
                double val = lambdaCache[v].get(i, 0);
                if (Double.isFinite(val)) {
                    s += val * pi0[i];
                }
            }
            mu[v] = s;
        }
        return mu;
    }
}
