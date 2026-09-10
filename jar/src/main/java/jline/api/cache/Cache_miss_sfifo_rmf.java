/**
 * @file Cache Miss Analysis for strict FIFO(m) via position-resolved mean field
 *
 * @since LINE 3.0
 */
package jline.api.cache;

import jline.lib.rmf.CacheSFIFORMF;
import jline.lib.rmf.CachePosGraphRMF;
import jline.util.matrix.Matrix;

/**
 * Position-resolved mean-field miss rates for strict FIFO(m) caches.
 *
 * <p>Java port of the MATLAB {@code cache_miss_sfifo_rmf.m} / Python
 * {@code cache_miss_sfifo_rmf}. Mirrors the {@link Cache_miss_fpi} /
 * {@link Cache_miss_rmf} contract: the popularity is recovered from the
 * per-user per-item arrival rates (list-0 column), and the per-item miss
 * probability is the out-of-cache occupancy of the position-resolved DDPP
 * fixed point (see {@link CacheSFIFORMF}).</p>
 *
 * <p>Unlike RANDOM(m)/FIFO(m), strict FIFO(m) has no closed-form or per-list
 * mean field (Gast15). This routine reduces to RANDOM(m)/FIFO(m) when
 * {@code m_1 = ... = m_{h-1} = 1}.</p>
 *
 * <p>Reference: N. Gast and B. Van Houdt, "Transient and Steady-state Regime of
 * a Family of List-based Cache Replacement Algorithms", ACM SIGMETRICS 2015.</p>
 */
public final class Cache_miss_sfifo_rmf {
    private Cache_miss_sfifo_rmf() {}

    /**
     * Compute strict FIFO(m) cache miss rates via the position-resolved mean field.
     *
     * @param gamma       item access factors (accepted for interface parity with
     *                    {@link Cache_miss_fpi}; unused).
     * @param m           cache capacity vector (h,).
     * @param lambdaCache per-user per-item per-list arrival rates: {@code lambdaCache[v]}
     *                    is an (n x (h+1)) matrix; column 0 carries the request rates.
     * @return {@link CacheMissFpiResult} with global miss rate M, per-user MU,
     *         per-item MI, and per-item miss probabilities pi0.
     */
    public static CacheMissFpiResult cache_miss_sfifo_rmf(Matrix gamma, Matrix m, Matrix[] lambdaCache) {
        return cache_miss_sfifo_rmf(gamma, m, lambdaCache, null);
    }

    /**
     * Strict FIFO(m) miss rates honouring a custom access graph. A non-linear
     * {@code accost} uses the general position-resolved drift from a cold cache
     * (see {@link CachePosGraphRMF}); the linear default keeps the pre-filled path.
     */
    public static CacheMissFpiResult cache_miss_sfifo_rmf(Matrix gamma, Matrix m, Matrix[] lambdaCache, Matrix[][] accost) {
        int u = lambdaCache.length;
        int n = lambdaCache[0].getNumRows();
        int h = m.length();

        // Aggregate per-item request rates over users (list-0 column).
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

        int[] mi = new int[h];
        for (int k = 0; k < h; k++) {
            mi[k] = (int) Math.round(m.get(k));
        }

        double[][][] G = Cache_miss_rmf.buildItemGraphs(accost, lambdaCache, n, h);
        double[] pi0arr;
        if (G != null) {
            CachePosGraphRMF gen = new CachePosGraphRMF(p, mi, true);
            gen.setItemGraph(G);
            pi0arr = gen.missProb(gen.fixedPointGraph());
        } else {
            CacheSFIFORMF sfifo = new CacheSFIFORMF(p, mi);
            pi0arr = sfifo.missProb(sfifo.fixedPoint());
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
     * Transient strict FIFO(m) cache trajectory via the position-resolved mean
     * field. Mirrors {@link Cache_miss_rmf#cache_miss_rmf_tran}.
     *
     * @param time    end time of the transient window (start is 0).
     * @param nPoints number of uniform grid points (>= 2).
     * @param xinit   initial occupancy (dim,), or null for the default warm start.
     */
    public static CacheMissRmfTranResult cache_miss_sfifo_rmf_tran(Matrix gamma, Matrix m, Matrix[] lambdaCache,
                                                                   double time, int nPoints, double[] xinit) {
        return cache_miss_sfifo_rmf_tran(gamma, m, lambdaCache, time, nPoints, xinit, null);
    }

    /**
     * Transient strict FIFO(m) cache trajectory honouring a custom access graph.
     * A non-linear {@code accost} uses the general position-resolved drift (cold
     * start when xinit is null); the linear default keeps the strict-FIFO drift.
     */
    public static CacheMissRmfTranResult cache_miss_sfifo_rmf_tran(Matrix gamma, Matrix m, Matrix[] lambdaCache,
                                                                   double time, int nPoints, double[] xinit,
                                                                   Matrix[][] accost) {
        int u = lambdaCache.length;
        int n = lambdaCache[0].getNumRows();
        int h = m.length();

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
        int[] mi = new int[h];
        for (int k = 0; k < h; k++) {
            mi[k] = (int) Math.round(m.get(k));
        }

        double[][][] G = Cache_miss_rmf.buildItemGraphs(accost, lambdaCache, n, h);
        Object[] traj;
        if (G != null) {
            CachePosGraphRMF gen = new CachePosGraphRMF(p, mi, true);
            gen.setItemGraph(G);
            traj = gen.driftTrajectory(time, nPoints, xinit);
        } else {
            traj = new CacheSFIFORMF(p, mi).driftTrajectory(time, nPoints, xinit);
        }
        double[] T = (double[]) traj[0];
        double[][] X = (double[][]) traj[1];
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
