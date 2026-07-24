/**
 * @file General position-resolved mean field with a per-item access graph
 *
 * @since LINE 3.0
 */
package jline.lib.rmf;

import odesolver.LSODA;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

/**
 * General position-resolved density-dependent population process (DDPP) mean
 * field for FIFO(m) / strict FIFO(m) caches honouring a per-item access graph.
 *
 * <p>{@code G[k]} is the (h+1)x(h+1) access graph of item k: row 0 is miss
 * admission (col 0 = reject, col l = admit to list l), row 1+i is a hit in list
 * i (col b = promote to list b>=i; b==i means STAY in place, the FIFO/SFIFO
 * convention). A miss admits at the head of the target list (its tail evicted);
 * a hit at position j of list i promotes to the head of target b>i, and the tail
 * of b is demoted to list i -- to the vacated position j for FIFO
 * ({@code reinsertHead=false}) or to the head with a 1..j-1 shift for strict
 * FIFO ({@code reinsertHead=true}). Reduces exactly to the linear drift when G
 * is the standard chain; the fixed point is integrated from a COLD (empty)
 * cache so non-admissible items drain.</p>
 */
public class CachePosGraphRMF {

    private final int n;
    private final int h;
    private final int[] m;
    private final double[] p;
    private final boolean reinsertHead;
    private final int[][] slotOf;   // slotOf[s] = {list i (1-based), position j (1-based)}
    private final int[][] sidx;     // sidx[i][j] = flat slot index
    private final int slots;
    private final int dim;
    private double[][][] G;

    public CachePosGraphRMF(double[] popularity, int[] capacities, boolean reinsertHead) {
        this.n = popularity.length;
        this.h = capacities.length;
        this.m = capacities.clone();
        this.p = popularity.clone();
        this.reinsertHead = reinsertHead;
        int total = 0, maxM = 0;
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
    }

    public void setItemGraph(double[][][] g) {
        this.G = g;
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
        return o < 0.0 ? 0.0 : (o > 1.0 ? 1.0 : o);
    }

    private double occList(double[] x, int k, int i) {
        double acc = 0.0;
        for (int j = 1; j <= m[i - 1]; j++) {
            acc += x[flat(k, i, j)];
        }
        return acc;
    }

    /** Partial-shift rate (SFIFO): promote-outs at positions deeper than jp. */
    private double gg(double[][] pop, int i, int jp) {
        double g = 0.0;
        for (int jj = jp + 1; jj <= m[i - 1]; jj++) {
            g += pop[i][jj];
        }
        return g;
    }

    public double[] driftGraph(double[] xin) {
        double[] x = new double[dim];
        for (int a = 0; a < dim; a++) {
            double v = xin[a];
            x[a] = v < 0.0 ? 0.0 : (v > 1.0 ? 1.0 : v);
        }
        double[] MI = new double[h + 1];         // MI[l] miss admission into list l
        double[][] HP = new double[h + 1][h + 1]; // HP[i][b] promotion i->b (b>i)
        for (int k = 0; k < n; k++) {
            double ok = outOf(x, k);
            double[][] gk = G[k];
            for (int l = 1; l <= h; l++) {
                MI[l] += p[k] * ok * gk[0][l];
            }
            for (int i = 1; i <= h; i++) {
                double oc = occList(x, k, i);
                for (int b = i + 1; b <= h; b++) {
                    HP[i][b] += p[k] * oc * gk[i][b];
                }
            }
        }
        double[] Sin = new double[h + 1];
        for (int l = 1; l <= h; l++) {
            Sin[l] = MI[l];
            for (int s = 1; s < l; s++) {
                Sin[l] += HP[s][l];
            }
        }
        double[][] pop = new double[h + 1][];
        for (int i = 1; i <= h; i++) {
            pop[i] = new double[m[i - 1] + 1];
            for (int j = 1; j <= m[i - 1]; j++) {
                double acc = 0.0;
                for (int k = 0; k < n; k++) {
                    acc += p[k] * x[flat(k, i, j)] * (1.0 - G[k][i][i]);
                }
                pop[i][j] = acc;
            }
        }
        double[] dX = new double[dim];
        for (int k = 0; k < n; k++) {
            double[][] gk = G[k];
            double ok = outOf(x, k);
            for (int i = 1; i <= h; i++) {
                for (int j = 1; j <= m[i - 1]; j++) {
                    double xk = x[flat(k, i, j)];
                    double o = p[k] * xk * (1.0 - gk[i][i]);
                    o += reinsertHead ? (Sin[i] + gg(pop, i, j)) * xk : Sin[i] * xk;
                    dX[flat(k, i, j)] -= o;
                    if (j >= 2) {
                        double sh = reinsertHead ? (Sin[i] + gg(pop, i, j - 1)) : Sin[i];
                        dX[flat(k, i, j)] += sh * x[flat(k, i, j - 1)];
                    } else {
                        dX[flat(k, i, 1)] += p[k] * ok * gk[0][i];
                        for (int s = 1; s < i; s++) {
                            dX[flat(k, i, 1)] += p[k] * occList(x, k, s) * gk[s][i];
                        }
                    }
                    for (int b = i + 1; b <= h; b++) {
                        if (reinsertHead) {
                            if (j == 1) {
                                dX[flat(k, i, 1)] += HP[i][b] * x[flat(k, b, m[b - 1])];
                            }
                        } else {
                            double poj = 0.0;
                            for (int kk = 0; kk < n; kk++) {
                                poj += p[kk] * x[flat(kk, i, j)] * G[kk][i][b];
                            }
                            dX[flat(k, i, j)] += poj * x[flat(k, b, m[b - 1])];
                        }
                    }
                }
            }
        }
        return dX;
    }

    /** Fixed point of the general drift from a cold (empty) cache. */
    public double[] fixedPointGraph(double tmax) {
        FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
            @Override
            public int getDimension() {
                return dim;
            }

            @Override
            public void computeDerivatives(double t, double[] y, double[] yDot) {
                double[] d = driftGraph(y);
                System.arraycopy(d, 0, yDot, 0, d.length);
            }
        };
        double[] result = new double[dim];
        double[] cold = new double[dim];
        LSODA lsoda = new LSODA(1e-12, 1.0, 1e-8, 1e-10, 12, 5);
        lsoda.integrate(ode, 0.0, cold, tmax, result);
        return result;
    }

    public double[] fixedPointGraph() {
        return fixedPointGraph(20000.0);
    }

    /**
     * Integrate the general drift over a finite window on a uniform time grid.
     * The start is a cold (empty) cache unless {@code xinit} is supplied.
     *
     * @param time    end time of the window (start is 0).
     * @param nPoints number of uniform grid points (>= 2).
     * @param xinit   initial occupancy (dim,), or null for a cold start.
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
                double[] d = driftGraph(y);
                System.arraycopy(d, 0, yDot, 0, d.length);
            }
        };
        double[][] X = new double[nPoints][dim];
        double[] state = new double[dim];
        if (xinit != null) {
            System.arraycopy(xinit, 0, state, 0, dim);
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

    public double[] missProb(double[] x) {
        double[] pi0 = new double[n];
        for (int k = 0; k < n; k++) {
            pi0[k] = outOf(x, k);
        }
        return pi0;
    }
}
