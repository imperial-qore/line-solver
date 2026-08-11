/**
 * @file Cache Miss Analysis via Refined Mean Field (RMF)
 *
 * @since LINE 3.0
 */
package jline.api.cache;

import jline.lib.rmf.CacheRMF;
import jline.util.matrix.Matrix;

/**
 * Refined mean-field (1/N-accurate) miss rates for RANDOM(m) caches.
 *
 * <p>Java port of the MATLAB {@code cache_miss_rmf.m} / Python
 * {@code cache_miss_rmf}. Mirrors the {@link Cache_miss_fpi} contract: the
 * popularity is recovered from the per-user per-item arrival rates (list-0
 * column), the mean-field fixed point is refined with its 1/N correction via
 * {@link CacheRMF}, and the per-item miss probability is the steady-state
 * occupancy of list 0.</p>
 *
 * <p>Reference: N. Gast, "Expected Values Estimated via Mean-Field
 * Approximation are 1/N-Accurate", Proc. ACM Meas. Anal. Comput. Syst., 2017.</p>
 */
public final class Cache_miss_rmf {
    private Cache_miss_rmf() {}

    /**
     * Compute cache miss rates using the refined mean-field method.
     *
     * @param gamma       item access factors (accepted for interface parity with
     *                    {@link Cache_miss_fpi}; used only for sizing fallbacks).
     * @param m           cache capacity vector (h,).
     * @param lambdaCache per-user per-item per-list arrival rates: {@code lambdaCache[v]}
     *                    is an (n x (h+1)) matrix; column 0 carries the request rates.
     * @return {@link CacheMissFpiResult} with global miss rate M, per-user MU,
     *         per-item MI, and per-item miss probabilities pi0.
     */
    public static CacheMissFpiResult cache_miss_rmf(Matrix gamma, Matrix m, Matrix[] lambdaCache) {
        return cache_miss_rmf(gamma, m, lambdaCache, null);
    }

    /**
     * Compute cache miss rates, honouring a custom access graph. {@code accost}
     * is the per-(user,item) access graph ({@code accost[v][k]} an (h+1)x(h+1)
     * Matrix: row 0 = miss admission per list, row 1+i = hit-in-list-i promotion
     * target). A non-linear graph uses the general RANDOM(m) drift with the
     * plain mean-field fixed point; the standard linear chain (or {@code null})
     * keeps the 1/N-refined path.
     */
    public static CacheMissFpiResult cache_miss_rmf(Matrix gamma, Matrix m, Matrix[] lambdaCache, Matrix[][] accost) {
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
            p[i] = sumLam > 0 ? lamI[i] / sumLam : 0.0;
        }

        int[] mi = new int[h];
        for (int k = 0; k < h; k++) {
            mi[k] = (int) Math.round(m.get(k));
        }

        CacheRMF rmf = new CacheRMF(p, mi);

        // A non-default access graph modulates admission (row 0) and promotion
        // (row 1+i) per item; the general drift honours it with the plain
        // mean-field fixed point, while the linear chain keeps the 1/N-refined
        // path unchanged. see _kb/03-api-layer.md for rationale.
        double[][][] G = buildItemGraphs(accost, lambdaCache, n, h);
        double[] xss;
        if (G != null) {
            rmf.setItemGraph(G);
            xss = rmf.fixedPointGraph();
        } else {
            xss = rmf.fixedPoint();
            try {
                Object[] res = rmf.meanFieldExpansionSteadyState();
                double[] pi = (double[]) res[0];
                double[] V = (double[]) res[1];
                double[] xref = new double[pi.length];
                boolean allFinite = true;
                for (int i = 0; i < pi.length; i++) {
                    xref[i] = pi[i] + V[i] / n;
                    if (!Double.isFinite(xref[i])) {
                        allFinite = false;
                        break;
                    }
                }
                if (allFinite) {
                    xss = xref;
                }
            } catch (RuntimeException e) {
                // keep plain mean field
            }
        }

        // Per-item miss probability = occupancy of list 0, clipped to [0,1].
        Matrix pi0 = new Matrix(n, 1);
        for (int i = 0; i < n; i++) {
            double val = xss[rmf.index(i, 0)];
            pi0.set(i, 0, Math.min(1.0, Math.max(0.0, val)));
        }

        Matrix MI = new Matrix(n, 1);
        for (int i = 0; i < n; i++) {
            MI.set(i, 0, lamI[i] * pi0.get(i, 0));
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
     * Per-item (h+1)x(h+1) access graph aggregated over users by request rate.
     * Returns null when accost is absent or the standard linear chain (so the
     * caller keeps the refined linear path). {@code accost} is [user][item].
     */
    static double[][][] buildItemGraphs(Matrix[][] accost, Matrix[] lambdaCache, int n, int h) {
        if (accost == null) {
            return null;
        }
        int u = accost.length;
        double[][] lin = linearGraph(h);
        double[][][] G = new double[n][h + 1][h + 1];
        boolean isLinear = true;
        for (int k = 0; k < n; k++) {
            double[][] num = new double[h + 1][h + 1];
            double den = 0.0;
            for (int v = 0; v < u; v++) {
                if (accost[v] == null || k >= accost[v].length || accost[v][k] == null) {
                    continue;
                }
                double wv = 0.0;
                if (v < lambdaCache.length) {
                    double val = lambdaCache[v].get(k, 0);
                    wv = Double.isFinite(val) ? val : 0.0;
                }
                Matrix gvk = accost[v][k];
                for (int a = 0; a <= h; a++) {
                    for (int b = 0; b <= h; b++) {
                        num[a][b] += wv * gvk.get(a, b);
                    }
                }
                den += wv;
            }
            double[][] gk = new double[h + 1][h + 1];
            if (den > 0) {
                for (int a = 0; a <= h; a++) {
                    for (int b = 0; b <= h; b++) {
                        gk[a][b] = num[a][b] / den;
                    }
                }
            } else if (accost[0] != null && k < accost[0].length && accost[0][k] != null) {
                for (int a = 0; a <= h; a++) {
                    for (int b = 0; b <= h; b++) {
                        gk[a][b] = accost[0][k].get(a, b);
                    }
                }
            } else {
                gk = lin;
            }
            // normalize rows
            for (int a = 0; a <= h; a++) {
                double rs = 0.0;
                for (int b = 0; b <= h; b++) {
                    rs += gk[a][b];
                }
                if (rs > 0) {
                    for (int b = 0; b <= h; b++) {
                        gk[a][b] /= rs;
                    }
                }
            }
            G[k] = gk;
            for (int a = 0; a <= h && isLinear; a++) {
                for (int b = 0; b <= h; b++) {
                    if (Math.abs(gk[a][b] - lin[a][b]) > 1e-9) {
                        isLinear = false;
                        break;
                    }
                }
            }
        }
        return isLinear ? null : G;
    }

    static double[][] linearGraph(int h) {
        double[][] g = new double[h + 1][h + 1];
        g[0][1] = 1.0;
        for (int a = 1; a < h; a++) {
            g[a][a + 1] = 1.0;
        }
        g[h][h] = 1.0;
        return g;
    }

    /**
     * Transient refined mean-field cache trajectory. Integrates the same plain
     * mean-field drift that {@link #cache_miss_rmf} drives to steady state over
     * a finite window, from a supplied (or default) initial occupancy.
     *
     * <p>Mirrors the MATLAB {@code cache_miss_rmf.m} tspan path (order-0, no 1/N
     * correction in the transient).</p>
     *
     * @param gamma       item access factors (accepted for interface parity;
     *                    used only for sizing).
     * @param m           cache capacity vector (h,).
     * @param lambdaCache per-user per-item per-list arrival rates.
     * @param time        end time of the transient window (start is 0).
     * @param nPoints     number of uniform grid points (>= 2).
     * @param xinit       initial occupancy (modelDimension,), or null for the
     *                    default first-m-in-list initial state.
     * @return {@link CacheMissRmfTranResult} with the time grid, per-item list-0
     *         occupancy, per-user miss-rate trajectory, and full occupancy
     *         trajectory.
     */
    public static CacheMissRmfTranResult cache_miss_rmf_tran(Matrix gamma, Matrix m, Matrix[] lambdaCache,
                                                             double time, int nPoints, double[] xinit) {
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
            p[i] = sumLam > 0 ? lamI[i] / sumLam : 0.0;
        }

        int[] mi = new int[h];
        for (int k = 0; k < h; k++) {
            mi[k] = (int) Math.round(m.get(k));
        }

        CacheRMF rmf = new CacheRMF(p, mi);
        Object[] traj = rmf.driftTrajectory(time, nPoints, xinit);
        double[] T = (double[]) traj[0];
        double[][] X = (double[][]) traj[1]; // nPoints x dim
        int nt = T.length;
        int dim = n * (h + 1);

        double[][] xtraj = new double[dim][nt];
        for (int t = 0; t < nt; t++) {
            for (int d = 0; d < dim; d++) {
                xtraj[d][t] = X[t][d];
            }
        }

        double[][] pi0_t = new double[n][nt];
        for (int i = 0; i < n; i++) {
            for (int t = 0; t < nt; t++) {
                double val = X[t][rmf.index(i, 0)];
                pi0_t[i][t] = Math.min(1.0, Math.max(0.0, val));
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
