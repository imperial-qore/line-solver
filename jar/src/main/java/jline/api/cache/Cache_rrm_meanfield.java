/**
 * @file Steady state of the RANDOM(m) multi-list mean field
 *
 * @since LINE 3.0
 */
package jline.api.cache;

import jline.util.matrix.Matrix;
import odesolver.LSODA;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

/**
 * Steady state of the RANDOM(m) multi-list mean field.
 *
 * <p>Twin of matlab/src/api/cache/cache_rrm_meanfield.m, which is a SCRIPT
 * rather than a function: it fixes the item count, the capacities, the
 * popularity vector and the initial condition, integrates the drift of
 * cache_rrm_meanfield_ode.m with a stiff integrator over [0, 10000] and prints
 * the terminal occupancy, the miss rate lambda'*x(:,1) and the miss ratio.
 * This class is that computation with the data as arguments, matching the
 * Python and C++ twins.</p>
 *
 * <p>The initial condition is the reference's: every item outside the cache,
 * x(k,0) = 1. The state is flattened COLUMN-MAJOR, x[k + s*n], which is what
 * the MATLAB reshape of the same vector means.</p>
 *
 * <p>Reference:
 *   N. Gast, B. Van Houdt, "Transient and steady-state regime of a family of
 *   list-based cache replacement algorithms", Queueing Syst. 83, 2016.</p>
 */
public final class Cache_rrm_meanfield {
    private Cache_rrm_meanfield() {}

    /** Integration horizon of the reference script. */
    private static final double DEFAULT_TMAX = 10000.0;

    /**
     * Steady state of the RANDOM(m) mean field with the reference horizon.
     *
     * @param lambda per-item request rates (n entries)
     * @param m      list capacities (h entries)
     * @return terminal occupancy, miss rate and miss ratio
     */
    public static CacheRrmMeanfieldResult cache_rrm_meanfield(Matrix lambda, Matrix m) {
        return cache_rrm_meanfield(lambda, m, DEFAULT_TMAX);
    }

    /**
     * Steady state of the RANDOM(m) mean field.
     *
     * @param lambda per-item request rates (n entries)
     * @param m      list capacities (h entries)
     * @param tmax   integration horizon
     * @return terminal occupancy, miss rate and miss ratio
     */
    public static CacheRrmMeanfieldResult cache_rrm_meanfield(Matrix lambda, Matrix m, double tmax) {
        final int n = (int) lambda.getNumElements();
        final int h = (int) m.getNumElements();
        if (n == 0) {
            throw new IllegalArgumentException("cache_rrm_meanfield: no items");
        }
        if (h == 0) {
            throw new IllegalArgumentException("cache_rrm_meanfield: no cache lists");
        }
        final int dim = n * (1 + h);
        final Matrix lambdaLocal = lambda;
        final Matrix mLocal = m;

        FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
            @Override
            public int getDimension() {
                return dim;
            }

            @Override
            public void computeDerivatives(double t, double[] y, double[] yDot) {
                Matrix x = new Matrix(n, 1 + h);
                for (int s = 0; s <= h; s++) {
                    for (int k = 0; k < n; k++) {
                        x.set(k, s, y[k + s * n]);
                    }
                }
                Matrix dxdt = Cache_rrm_meanfield_ode.cache_rrm_meanfield_ode(x, lambdaLocal, mLocal, n, h);
                for (int s = 0; s <= h; s++) {
                    for (int k = 0; k < n; k++) {
                        yDot[k + s * n] = dxdt.get(k, s);
                    }
                }
            }
        };

        double[] x0 = new double[dim];
        for (int k = 0; k < n; k++) {
            x0[k] = 1.0;
        }
        double[] xf = new double[dim];
        LSODA lsoda = new LSODA(1e-12, 1.0, 1e-8, 1e-10, 12, 5);
        lsoda.integrate(ode, 0.0, x0, tmax, xf);

        Matrix x = new Matrix(n, 1 + h);
        for (int s = 0; s <= h; s++) {
            for (int k = 0; k < n; k++) {
                x.set(k, s, xf[k + s * n]);
            }
        }

        double missrate = 0.0;
        double tot = 0.0;
        for (int k = 0; k < n; k++) {
            missrate += lambda.get(k) * xf[k];
            tot += lambda.get(k);
        }
        double missratio = tot > 0 ? missrate / tot : 0.0;
        return new CacheRrmMeanfieldResult(x, missrate, missratio);
    }
}
