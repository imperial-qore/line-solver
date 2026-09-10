/**
 * @file TTL approximation for LRU caches with arrival routing (access graphs)
 *
 * Characteristic-time (TTL) approximation for LRU caches whose lists form an
 * arbitrary access graph (linear chains and trees), mirroring the MATLAB
 * reference cache_ttl_lrua.m. Each item is modeled by an embedded DTMC over
 * the lists, with transition probabilities given by exponential timer races,
 * converted to time-stationary probabilities via mean holding times. The
 * characteristic times are fixed so that the expected occupancy of each list
 * equals its capacity; the joint fixed point is computed by a damped Newton
 * iteration on the log characteristic times with a finite-difference
 * Jacobian (no external nonlinear solver required).
 *
 * @since LINE 3.0
 */
package jline.api.cache;

import org.apache.commons.math3.util.FastMath;

import jline.api.mc.Dtmc_solve;
import jline.util.matrix.Matrix;

public final class Cache_ttl_lrua {
    private Cache_ttl_lrua() {}

    private static final int MAX_NEWTON = 200;
    private static final double TOL = 1e-10;

    /**
     * Solve LRU caches with access graphs using the TTL tree approximation.
     *
     * @param lambda Per-user arrival rate matrices, lambda[v] is (n x h+1)
     *               with column 0 the out-of-cache request rate and column l
     *               the request rate seen while the item is in list l.
     * @param R      Access graph, R[v][k] is the (h+1 x h+1) routing matrix
     *               of item k for user v (row/col 0 is the out-of-cache node).
     * @param m      Cache capacity vector (1 x h).
     * @return Matrix (n x h+1) of per-item probabilities; column 0 is the
     *         out-of-cache (miss) probability, column l the probability of
     *         residing in list l.
     */
    public static Matrix cache_ttl_lrua(Matrix[] lambda, Matrix[][] R, Matrix m) {
        if (lambda == null || R == null || m == null) {
            throw new IllegalArgumentException("Lambda, R, and m parameters cannot be null");
        }

        int n = lambda[0].getNumRows();
        int h = lambda[0].getNumCols() - 1;

        // Damped Newton on y = log(T): solve occ(exp(y)) - m = 0 jointly.
        // The log parameterization keeps the characteristic times positive.
        double[] y = new double[h];
        double[] res = residual(lambda, R, y, m, n, h);
        double nrm = norm(res);
        for (int it = 0; it < MAX_NEWTON && nrm > TOL; it++) {
            // finite-difference Jacobian
            double[][] J = new double[h][h];
            double eps = 1e-6;
            for (int c = 0; c < h; c++) {
                double[] yp = new double[h];
                System.arraycopy(y, 0, yp, 0, h);
                yp[c] += eps;
                double[] rp = residual(lambda, R, yp, m, n, h);
                for (int r = 0; r < h; r++) {
                    J[r][c] = (rp[r] - res[r]) / eps;
                }
            }
            double[] step = solveLinear(J, res, h);
            if (step == null) {
                break;
            }
            // backtracking line search on the residual norm
            double alpha = 1.0;
            double[] ynew = new double[h];
            double[] rnew = null;
            double nrmNew = Double.POSITIVE_INFINITY;
            for (int ls = 0; ls < 40; ls++) {
                for (int l = 0; l < h; l++) {
                    ynew[l] = y[l] - alpha * step[l];
                }
                rnew = residual(lambda, R, ynew, m, n, h);
                nrmNew = norm(rnew);
                if (nrmNew < nrm) {
                    break;
                }
                alpha *= 0.5;
            }
            if (nrmNew >= nrm) {
                break;
            }
            System.arraycopy(ynew, 0, y, 0, h);
            res = rnew;
            nrm = nrmNew;
        }

        double[] x = new double[h];
        for (int l = 0; l < h; l++) {
            x[l] = FastMath.exp(y[l]);
        }
        return randProb(lambda, R, x, n, h);
    }

    /* Occupancy residuals occ_l(exp(y)) - m_l. */
    private static double[] residual(Matrix[] lambda, Matrix[][] R, double[] y,
                                     Matrix m, int n, int h) {
        double[] x = new double[h];
        for (int l = 0; l < h; l++) {
            x[l] = FastMath.exp(y[l]);
        }
        Matrix prob = randProb(lambda, R, x, n, h);
        double[] res = new double[h];
        for (int l = 0; l < h; l++) {
            double capa = 0.0;
            for (int i = 0; i < n; i++) {
                capa += prob.get(i, l + 1);
            }
            res[l] = capa - m.get(l);
        }
        return res;
    }

    private static double norm(double[] v) {
        double s = 0.0;
        for (int i = 0; i < v.length; i++) {
            s += v[i] * v[i];
        }
        return FastMath.sqrt(s);
    }

    /* Solve J s = r by Gaussian elimination with partial pivoting. */
    private static double[] solveLinear(double[][] J, double[] r, int h) {
        double[][] a = new double[h][h + 1];
        for (int i = 0; i < h; i++) {
            System.arraycopy(J[i], 0, a[i], 0, h);
            a[i][h] = r[i];
        }
        for (int c = 0; c < h; c++) {
            int piv = c;
            for (int i = c + 1; i < h; i++) {
                if (FastMath.abs(a[i][c]) > FastMath.abs(a[piv][c])) {
                    piv = i;
                }
            }
            if (FastMath.abs(a[piv][c]) < 1e-300) {
                return null;
            }
            double[] tmp = a[c];
            a[c] = a[piv];
            a[piv] = tmp;
            for (int i = c + 1; i < h; i++) {
                double f = a[i][c] / a[c][c];
                for (int j = c; j <= h; j++) {
                    a[i][j] -= f * a[c][j];
                }
            }
        }
        double[] s = new double[h];
        for (int i = h - 1; i >= 0; i--) {
            double v = a[i][h];
            for (int j = i + 1; j < h; j++) {
                v -= a[i][j] * s[j];
            }
            s[i] = v / a[i][i];
        }
        return s;
    }

    // see _kb/03-api-layer.md for rationale
    private static Matrix randProb(Matrix[] lambda, Matrix[][] R, double[] x, int n, int h) {
        Matrix randprob = new Matrix(n, h + 1);
        for (int i = 0; i < n; i++) {
            Matrix Ri = R[0][i];
            double[][] trans = new double[h + 1][h + 1];
            for (int j = 0; j <= h; j++) {
                for (int k = 0; k <= h; k++) {
                    if (Ri.get(j, k) > 0) {
                        if (j == 0) {
                            trans[j][k] = Ri.get(j, k);
                        } else {
                            trans[j][k] = (1.0 - FastMath.exp(-lambda[0].get(i, j) * x[j - 1]))
                                    * Ri.get(j, k);
                        }
                        if (j != k && k > 0) {
                            trans[k][j] = FastMath.exp(-lambda[0].get(i, k) * x[k - 1]);
                        }
                    }
                }
            }

            // see _kb/03-api-layer.md for rationale
            boolean[] conn = new boolean[h + 1];
            int nconn = 0;
            for (int j = 0; j <= h; j++) {
                for (int k = 0; k <= h; k++) {
                    if (trans[j][k] != 0.0) {
                        conn[k] = true;
                    }
                }
            }
            int[] map = new int[h + 1];
            for (int j = 0; j <= h; j++) {
                if (conn[j]) {
                    map[j] = nconn;
                    nconn++;
                } else {
                    map[j] = -1;
                }
            }
            if (nconn == 0) {
                continue;
            }
            Matrix sub = new Matrix(nconn, nconn);
            for (int j = 0; j <= h; j++) {
                if (!conn[j]) {
                    continue;
                }
                for (int k = 0; k <= h; k++) {
                    if (conn[k]) {
                        sub.set(map[j], map[k], trans[j][k]);
                    }
                }
            }
            Matrix dtmcprob = Dtmc_solve.dtmc_solve(sub);

            double denom = 0.0;
            double[] avgtime = new double[h + 1];
            double[] ssprob = new double[h + 1];
            for (int node = 0; node <= h; node++) {
                if (!conn[node]) {
                    continue;
                }
                ssprob[node] = dtmcprob.get(map[node]);
                double rate = lambda[0].get(i, node);
                if (node > 0) {
                    avgtime[node] = rate > 0
                            ? (1.0 - FastMath.exp(-rate * x[node - 1])) / rate : 0.0;
                } else {
                    avgtime[node] = rate > 0 ? 1.0 / rate : 1.0;
                }
                denom += ssprob[node] * avgtime[node];
            }
            if (denom <= 0) {
                continue;
            }
            for (int node = 0; node <= h; node++) {
                if (conn[node]) {
                    randprob.set(i, node, ssprob[node] * avgtime[node] / denom);
                }
            }
        }
        return randprob;
    }
}
