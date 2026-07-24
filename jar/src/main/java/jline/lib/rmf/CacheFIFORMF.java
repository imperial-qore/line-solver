/**
 * @file Position-resolved mean field for FIFO(m) caches
 *
 * @since LINE 3.0
 */
package jline.lib.rmf;

import odesolver.LSODA;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

/**
 * Position-resolved density-dependent population process (DDPP) mean field for
 * FIFO(m) cache replacement.
 *
 * <p>FIFO(m) and RANDOM(m) share the exact stationary distribution (Gast and
 * Van Houdt, SIGMETRICS 2015, Thm 1: {@code pi_FIFO(m) = pi_RAND(m)}), so their
 * steady-state hit ratios coincide. Their mean-field TRANSIENTS differ: FIFO
 * evicts the deterministic tail (fixed residence of m insertions) whereas
 * RANDOM evicts a uniformly random victim (geometric residence), so H(t) from a
 * cold cache ramps differently even though H(inf) agrees. This model provides
 * that dedicated FIFO transient.</p>
 *
 * <p>FIFO(m) differs from strict FIFO(m) ({@link CacheSFIFORMF}) only in the
 * reinsertion position on a hit: the demoted tail of list i+1 lands at the
 * vacated position j of list i (in place, no within-list shift), whereas strict
 * FIFO reinserts it at position 1.</p>
 *
 * <p>Reference: N. Gast and B. Van Houdt, "Transient and Steady-state Regime of
 * a Family of List-based Cache Replacement Algorithms", ACM SIGMETRICS 2015.</p>
 */
public class CacheFIFORMF {

    private final int n;
    private final int h;
    private final int[] m;
    private final double[] p;
    private final int[][] slotOf;
    private final int[][] sidx;
    private final int slots;
    private final int dim;
    private final double[] x0;

    public CacheFIFORMF(double[] popularity, int[] capacities) {
        this.n = popularity.length;
        this.h = capacities.length;
        this.m = capacities.clone();
        this.p = popularity.clone();

        int total = 0;
        int maxM = 0;
        for (int i = 0; i < h; i++) {
            total += m[i];
            if (m[i] > maxM) {
                maxM = m[i];
            }
        }
        this.slots = total;
        this.dim = n * total;
        this.slotOf = new int[total][2];
        this.sidx = new int[h + 1][maxM + 1];
        int s = 0;
        for (int i = 1; i <= h; i++) {
            for (int j = 1; j <= m[i - 1]; j++) {
                slotOf[s][0] = i;
                slotOf[s][1] = j;
                sidx[i][j] = s;
                s++;
            }
        }

        Integer[] order = new Integer[n];
        for (int i = 0; i < n; i++) {
            order[i] = i;
        }
        java.util.Arrays.sort(order, new java.util.Comparator<Integer>() {
            @Override
            public int compare(Integer a, Integer b) {
                return Double.compare(p[b], p[a]);
            }
        });
        this.x0 = new double[dim];
        int pos = 0;
        for (s = 0; s < slots; s++) {
            if (pos < n) {
                x0[order[pos] * slots + s] = 1.0;
                pos++;
            }
        }
    }

    public int getDimension() {
        return dim;
    }

    private int flat(int k, int i, int j) {
        return k * slots + sidx[i][j];
    }

    private double outOf(double[] x, int k) {
        double acc = 0.0;
        for (int s = 0; s < slots; s++) {
            acc += x[k * slots + s];
        }
        double o = 1.0 - acc;
        if (o < 0.0) {
            o = 0.0;
        } else if (o > 1.0) {
            o = 1.0;
        }
        return o;
    }

    private double occList(double[] x, int k, int i) {
        double acc = 0.0;
        for (int j = 1; j <= m[i - 1]; j++) {
            acc += x[flat(k, i, j)];
        }
        return acc;
    }

    /**
     * Mean-field drift F(x) for FIFO(m).
     *
     * @param xin state vector of dimension {@link #getDimension()}.
     * @return dX of the same dimension.
     */
    public double[] drift(double[] xin) {
        double[] x = new double[dim];
        for (int a = 0; a < dim; a++) {
            double v = xin[a];
            x[a] = v < 0.0 ? 0.0 : (v > 1.0 ? 1.0 : v);
        }

        double[][] hpos = new double[h + 1][];
        double[] hi = new double[h + 1];
        for (int i = 1; i <= h; i++) {
            hpos[i] = new double[m[i - 1] + 1];
        }
        for (int s = 0; s < slots; s++) {
            int i = slotOf[s][0];
            int j = slotOf[s][1];
            double acc = 0.0;
            for (int k = 0; k < n; k++) {
                acc += p[k] * x[flat(k, i, j)];
            }
            hpos[i][j] = acc;
            hi[i] += acc;
        }
        double miss = 0.0;
        for (int k = 0; k < n; k++) {
            miss += p[k] * outOf(x, k);
        }

        double[] sfull = new double[h + 1];
        sfull[1] = miss;
        for (int i = 2; i <= h; i++) {
            sfull[i] = hi[i - 1];
        }

        double[] dX = new double[dim];
        for (int k = 0; k < n; k++) {
            for (int s = 0; s < slots; s++) {
                int i = slotOf[s][0];
                int j = slotOf[s][1];
                double xk = x[flat(k, i, j)];
                // outflow: full shift toward j+1 (tail leaves); promote up if i<h
                double o = sfull[i] * xk;
                if (i < h) {
                    o += p[k] * xk;
                }
                dX[flat(k, i, j)] -= o;
                // inflow: full shift from j-1, or front insertion at j == 1
                if (j >= 2) {
                    dX[flat(k, i, j)] += sfull[i] * x[flat(k, i, j - 1)];
                } else {
                    if (i == 1) {
                        dX[flat(k, 1, 1)] += p[k] * outOf(x, k);           // miss inserts item k
                    } else {
                        dX[flat(k, i, 1)] += p[k] * occList(x, k, i - 1);  // promotion from i-1
                    }
                }
                // FIFO demotion: tail of i+1 lands in place at the same position j
                if (i < h) {
                    dX[flat(k, i, j)] += hpos[i][j] * x[flat(k, i + 1, m[i])];
                }
            }
        }
        return dX;
    }

    /**
     * Mean-field fixed point by integrating dx/dt = F(x) to steady state.
     *
     * @param tmax integration horizon.
     * @return fixed-point state vector of dimension {@link #getDimension()}.
     */
    public double[] fixedPoint(double tmax) {
        FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
            @Override
            public int getDimension() {
                return dim;
            }

            @Override
            public void computeDerivatives(double t, double[] y, double[] yDot) {
                double[] d = drift(y);
                System.arraycopy(d, 0, yDot, 0, d.length);
            }
        };
        double[] result = new double[dim];
        LSODA lsoda = new LSODA(1e-12, 1.0, 1e-8, 1e-10, 12, 5);
        lsoda.integrate(ode, 0.0, x0, tmax, result);
        return result;
    }

    public double[] fixedPoint() {
        return fixedPoint(20000.0);
    }

    /**
     * Integrate the drift over a finite window on a uniform time grid from a
     * supplied (or default) initial occupancy. Transient counterpart of
     * {@link #fixedPoint()}, mirroring {@code CacheRMF.driftTrajectory}.
     *
     * @param time    end time of the window (start is 0).
     * @param nPoints number of uniform grid points (>= 2).
     * @param xinit   initial occupancy (dim,), or null for the default warm start.
     * @return {@code Object[]{ T (double[nPoints]), X (double[nPoints][dim]) }}.
     */
    public Object[] driftTrajectory(double time, int nPoints, double[] xinit) {
        double[] T = new double[nPoints];
        for (int i = 0; i < nPoints; i++) {
            T[i] = time * i / (nPoints - 1);
        }
        FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
            @Override
            public int getDimension() {
                return dim;
            }

            @Override
            public void computeDerivatives(double t, double[] y, double[] yDot) {
                double[] d = drift(y);
                System.arraycopy(d, 0, yDot, 0, d.length);
            }
        };
        double[][] X = new double[nPoints][dim];
        double[] state = new double[dim];
        if (xinit != null) {
            System.arraycopy(xinit, 0, state, 0, dim);
        } else {
            System.arraycopy(x0, 0, state, 0, dim);
        }
        System.arraycopy(state, 0, X[0], 0, dim);
        for (int i = 1; i < nPoints; i++) {
            double[] result = new double[dim];
            LSODA lsoda = new LSODA(1e-12, 1.0, 1e-8, 1e-10, 12, 5);
            lsoda.integrate(ode, T[i - 1], state, T[i], result);
            System.arraycopy(result, 0, X[i], 0, dim);
            System.arraycopy(result, 0, state, 0, dim);
        }
        return new Object[]{T, X};
    }

    /**
     * Per-item out-of-cache (miss) probability of a state vector.
     *
     * @param x state vector.
     * @return pi0[k] = P(item k out of cache), length n.
     */
    public double[] missProb(double[] x) {
        double[] pi0 = new double[n];
        for (int k = 0; k < n; k++) {
            pi0[k] = outOf(x, k);
        }
        return pi0;
    }
}
