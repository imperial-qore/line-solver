/**
 * @file Cache Miss Analysis for FIFO(m) via position-resolved mean field
 *
 * @since LINE 3.0
 */
package jline.api.cache;

import jline.lib.rmf.CacheFIFORMF;
import jline.lib.rmf.CachePosGraphRMF;
import jline.util.matrix.Matrix;

/**
 * Position-resolved mean-field miss rates for FIFO(m) caches.
 *
 * <p>FIFO(m) and RANDOM(m) share the exact stationary distribution (Gast15
 * Thm 1), so FLD serves FIFO steady state from {@link Cache_miss_rmf}; the
 * dedicated value of this class is the FIFO mean-field TRANSIENT
 * ({@link #cache_miss_fifo_rmf_tran}), which differs from RANDOM(m) even though
 * the fixed points agree. Reduces to RANDOM(m) when
 * {@code m_1 = ... = m_{h-1} = 1}.</p>
 *
 * <p>Reference: N. Gast and B. Van Houdt, "Transient and Steady-state Regime of
 * a Family of List-based Cache Replacement Algorithms", ACM SIGMETRICS 2015.</p>
 */
public final class Cache_miss_fifo_rmf {
    private Cache_miss_fifo_rmf() {}

    private static double[] popularity(Matrix[] lambdaCache, int u, int n) {
        double[] lamI = new double[n];
        for (int v = 0; v < u; v++) {
            for (int i = 0; i < n; i++) {
                double val = lambdaCache[v].get(i, 0);
                if (Double.isFinite(val)) {
                    lamI[i] += val;
                }
            }
        }
        double sumLam = 0.0;
        for (int i = 0; i < n; i++) {
            sumLam += lamI[i];
        }
        double[] p = new double[n];
        for (int i = 0; i < n; i++) {
            p[i] = sumLam > 0 ? lamI[i] / sumLam : 1.0 / n;
        }
        return p;
    }

    private static int[] capacities(Matrix m, int h) {
        int[] mi = new int[h];
        for (int k = 0; k < h; k++) {
            mi[k] = (int) Math.round(m.get(k));
        }
        return mi;
    }

    /**
     * Steady-state FIFO(m) miss rates via the position-resolved mean field.
     * Note FIFO(m) equals RANDOM(m) at steady state (Gast15 Thm 1); prefer
     * {@link Cache_miss_rmf#cache_miss_rmf} when only the fixed point is needed.
     */
    public static CacheMissFpiResult cache_miss_fifo_rmf(Matrix gamma, Matrix m, Matrix[] lambdaCache) {
        return cache_miss_fifo_rmf(gamma, m, lambdaCache, null);
    }

    /**
     * FIFO(m) miss rates honouring a custom access graph. A non-linear
     * {@code accost} uses the general position-resolved drift from a cold cache
     * (see {@link CachePosGraphRMF}); the linear default keeps the pre-filled path.
     */
    public static CacheMissFpiResult cache_miss_fifo_rmf(Matrix gamma, Matrix m, Matrix[] lambdaCache, Matrix[][] accost) {
        int u = lambdaCache.length;
        int n = lambdaCache[0].getNumRows();
        int h = m.length();
        double[] p = popularity(lambdaCache, u, n);

        double[][][] G = Cache_miss_rmf.buildItemGraphs(accost, lambdaCache, n, h);
        double[] pi0arr;
        if (G != null) {
            CachePosGraphRMF gen = new CachePosGraphRMF(p, capacities(m, h), false);
            gen.setItemGraph(G);
            pi0arr = gen.missProb(gen.fixedPointGraph());
        } else {
            CacheFIFORMF fifo = new CacheFIFORMF(p, capacities(m, h));
            pi0arr = fifo.missProb(fifo.fixedPoint());
        }

        double[] lamI = new double[n];
        for (int v = 0; v < u; v++) {
            for (int i = 0; i < n; i++) {
                double val = lambdaCache[v].get(i, 0);
                if (Double.isFinite(val)) {
                    lamI[i] += val;
                }
            }
        }

        Matrix pi0 = new Matrix(n, 1);
        Matrix MI = new Matrix(n, 1);
        for (int i = 0; i < n; i++) {
            double v = Math.min(1.0, Math.max(0.0, pi0arr[i]));
            pi0.set(i, 0, v);
            MI.set(i, 0, lamI[i] * v);
        }
        Matrix MU = new Matrix(u, 1);
        for (int v = 0; v < u; v++) {
            double s = 0.0;
            for (int i = 0; i < n; i++) {
                double val = lambdaCache[v].get(i, 0);
                if (Double.isFinite(val)) {
                    s += val * pi0.get(i, 0);
                }
            }
            MU.set(v, 0, s);
        }
        double M = MI.elementSum();
        return new CacheMissFpiResult(M, MU, MI, pi0);
    }

    /**
     * Transient FIFO(m) cache trajectory via the position-resolved mean field.
     * Mirrors {@link Cache_miss_rmf#cache_miss_rmf_tran}; the FIFO transient
     * differs from RANDOM(m) even though the steady states agree.
     *
     * @param time    end time of the transient window (start is 0).
     * @param nPoints number of uniform grid points (>= 2).
     * @param xinit   initial occupancy (dim,), or null for the default warm start.
     */
    public static CacheMissRmfTranResult cache_miss_fifo_rmf_tran(Matrix gamma, Matrix m, Matrix[] lambdaCache,
                                                                  double time, int nPoints, double[] xinit) {
        return cache_miss_fifo_rmf_tran(gamma, m, lambdaCache, time, nPoints, xinit, null);
    }

    /**
     * Transient FIFO(m) cache trajectory honouring a custom access graph. A
     * non-linear {@code accost} uses the general position-resolved drift (cold
     * start when xinit is null); the linear default keeps the FIFO drift.
     */
    public static CacheMissRmfTranResult cache_miss_fifo_rmf_tran(Matrix gamma, Matrix m, Matrix[] lambdaCache,
                                                                  double time, int nPoints, double[] xinit,
                                                                  Matrix[][] accost) {
        int u = lambdaCache.length;
        int n = lambdaCache[0].getNumRows();
        int h = m.length();
        double[] p = popularity(lambdaCache, u, n);
        int[] mi = capacities(m, h);

        double[][][] G = Cache_miss_rmf.buildItemGraphs(accost, lambdaCache, n, h);
        Object[] traj;
        if (G != null) {
            CachePosGraphRMF gen = new CachePosGraphRMF(p, mi, false);
            gen.setItemGraph(G);
            traj = gen.driftTrajectory(time, nPoints, xinit);
        } else {
            traj = new CacheFIFORMF(p, mi).driftTrajectory(time, nPoints, xinit);
        }
        double[] T = (double[]) traj[0];
        double[][] X = (double[][]) traj[1]; // nPoints x dim
        int nt = T.length;
        int slots = 0;
        for (int i = 0; i < h; i++) {
            slots += mi[i];
        }
        int dim = n * slots;

        double[][] xtraj = new double[dim][nt];
        for (int t = 0; t < nt; t++) {
            for (int d = 0; d < dim; d++) {
                xtraj[d][t] = X[t][d];
            }
        }
        double[][] pi0_t = new double[n][nt];
        for (int t = 0; t < nt; t++) {
            for (int i = 0; i < n; i++) {
                double acc = 0.0;
                for (int s = 0; s < slots; s++) {
                    acc += X[t][i * slots + s];
                }
                double o = 1.0 - acc;
                pi0_t[i][t] = o < 0.0 ? 0.0 : (o > 1.0 ? 1.0 : o);
            }
        }
        double[][] MU_t = new double[u][nt];
        for (int v = 0; v < u; v++) {
            for (int t = 0; t < nt; t++) {
                double s = 0.0;
                for (int i = 0; i < n; i++) {
                    double val = lambdaCache[v].get(i, 0);
                    if (Double.isFinite(val)) {
                        s += val * pi0_t[i][t];
                    }
                }
                MU_t[v][t] = s;
            }
        }
        return new CacheMissRmfTranResult(T, pi0_t, MU_t, xtraj);
    }
}
