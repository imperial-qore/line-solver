/**
 * @file Position-resolved mean field for strict FIFO(m) caches (SFIFO)
 *
 * @since LINE 3.0
 */
package jline.lib.rmf;

import odesolver.LSODA;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

/**
 * Position-resolved density-dependent population process (DDPP) mean field for
 * strict FIFO(m) cache replacement.
 *
 * <p>Strict FIFO(m) is NOT equivalent to RANDOM(m)/FIFO(m). Gast and Van Houdt
 * (SIGMETRICS 2015) prove {@code pi_FIFO(m) = pi_RAND(m)} exactly but show
 * strict FIFO(m) differs and give it no mean-field model. The difference is the
 * within-list age ordering: on a hit in list {@code i < h} the demoted tail of
 * list {@code i+1} is reinserted at position 1 of list {@code i} (positions
 * {@code 1..j-1} shift back), which the per-item per-list occupancy of
 * RANDOM(m) cannot represent. This model tracks {@code x[k,i,j] = P(item k in
 * position j of list i)} with deterministic (age-based) demotion/eviction and
 * returns the plain mean-field fixed point; it reduces to RANDOM(m)/FIFO(m)
 * when {@code m_1 = ... = m_{h-1} = 1}.</p>
 *
 * <p>Reference: N. Gast and B. Van Houdt, "Transient and Steady-state Regime of
 * a Family of List-based Cache Replacement Algorithms", ACM SIGMETRICS 2015.</p>
 */
public class CacheSFIFORMF {

    private final int n;               // number of items
    private final int h;               // number of lists
    private final int[] m;             // list capacities m_1..m_h (length h)
    private final double[] p;          // aggregate popularity (sums to 1)
    private final int[][] slotOf;      // slotOf[s] = {list i (1-based), position j (1-based)}
    private final int[][] sidx;        // sidx[i][j] = flat slot index (i in 1..h, j in 1..m_i)
    private final int slots;           // S = sum(m)
    private final int dim;             // n * S
    private final double[] x0;

    /**
     * @param popularity aggregate per-item request probabilities (length n).
     * @param capacities list capacities m_1..m_h (length h).
     */
    public CacheSFIFORMF(double[] popularity, int[] capacities) {
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

        // Popularity-ordered warm start: most popular items fill the slots.
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
     * Mean-field drift F(x) for strict FIFO(m).
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

        // per-position and per-list hit rates, and the miss rate
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

        // full-shift rate of each list (list 1 on a miss; list i on a hit in i-1)
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
                // outflow: shift toward j+1 (or leave list at tail); promote up if i<h
                double o = (sfull[i] + gi(hpos, i, j)) * xk;
                if (i < h) {
                    o += p[k] * xk;
                }
                dX[flat(k, i, j)] -= o;
                // inflow
                if (j >= 2) {
                    dX[flat(k, i, j)] += (sfull[i] + gi(hpos, i, j - 1)) * x[flat(k, i, j - 1)];
                } else {
                    if (i == 1) {
                        dX[flat(k, 1, 1)] += p[k] * outOf(x, k);            // miss inserts item k
                    } else {
                        dX[flat(k, i, 1)] += p[k] * occList(x, k, i - 1);   // promotion from i-1
                    }
                    if (i < h) {
                        dX[flat(k, i, 1)] += hi[i] * x[flat(k, i + 1, m[i])]; // demotion from i+1 tail
                    }
                }
            }
        }
        return dX;
    }

    /**
     * Partial-shift rate of a slot at position jp of list i: aggregate hit rate
     * at deeper positions of list i. The top list h never moves on a hit.
     */
    private double gi(double[][] hpos, int i, int jp) {
        if (i == h) {
            return 0.0;
        }
        double g = 0.0;
        for (int jj = jp + 1; jj <= m[i - 1]; jj++) {
            g += hpos[i][jj];
        }
        return g;
    }

    /**
     * Compute the mean-field fixed point by integrating dx/dt = F(x) to steady
     * state with LSODA.
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
     * Per-item out-of-cache (miss) probability at the fixed point.
     *
     * @param xss fixed-point state vector.
     * @return pi0[k] = P(item k out of cache), length n.
     */
    public double[] missProb(double[] xss) {
        double[] pi0 = new double[n];
        for (int k = 0; k < n; k++) {
            pi0[k] = outOf(xss, k);
        }
        return pi0;
    }
}
