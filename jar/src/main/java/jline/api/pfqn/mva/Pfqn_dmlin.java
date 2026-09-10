/**
 * @file de Souza e Silva-Muntz Improved Linearizer (IL)
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva;

import org.apache.commons.math3.util.FastMath;

import jline.io.Ret;
import jline.util.matrix.Matrix;

/**
 * de Souza e Silva-Muntz Improved Linearizer (IL).
 *
 * <p>E. de Souza e Silva, R. R. Muntz, "A note on the computational cost of the
 * Linearizer algorithm for queueing networks", IEEE Trans. Computers 39(6),
 * 1990. Linearizer evaluates the arrival-instant queue length as</p>
 *
 * <pre>A_k^(c)(n) = sum_i (n_i - delta_c^(i)) [Q_ik(n)/n_i + Delta^(i)_ck],</pre>
 *
 * <p>re-summing the C Delta-terms at every Core iteration, at every one of the
 * C+1 populations: O(K C^3) per refresh pass. IL splits that sum into the part
 * that moves with the Core iterate and the part that does not,</p>
 *
 * <pre>
 * A_k^(c)(n)     = sum_i (n_i - delta_c^(i)) Q_ik(n)/n_i + xi_ck(n),
 * xi_ck(N)       = sum_i (N_i - delta_c^(i)) Delta^(i)_ck,
 * xi_ck(N - 1_j) = xi_ck(N) - Delta^(j)_ck,
 * </pre>
 *
 * <p>so the C K aggregates xi are computed ONCE per refresh pass and each Core
 * iteration then costs O(K C) instead of O(K C^2). Time drops to O(K C^2) with
 * the space unchanged at O(K C^2), and, because the split is an identity and
 * not an approximation, the fixed point is the one Linearizer reaches:
 * pfqn_dmlin and pfqn_linearizer agree to round-off.</p>
 */
public final class Pfqn_dmlin {
    private Pfqn_dmlin() {}

    public static Ret.pfqnAMVA pfqn_dmlin(Matrix L, Matrix N) {
        return pfqn_dmlin(L, N, new Matrix(1, L.getNumCols()), 1e-8, 1000, null, 3);
    }

    public static Ret.pfqnAMVA pfqn_dmlin(Matrix L, Matrix N, Matrix Z) {
        return pfqn_dmlin(L, N, Z, 1e-8, 1000, null, 3);
    }

    public static Ret.pfqnAMVA pfqn_dmlin(Matrix L, Matrix N, Matrix Z, double tol, int maxiter) {
        return pfqn_dmlin(L, N, Z, tol, maxiter, null, 3);
    }

    public static Ret.pfqnAMVA pfqn_dmlin(Matrix L, Matrix N, Matrix Z, double tol, int maxiter, Matrix QN0) {
        return pfqn_dmlin(L, N, Z, tol, maxiter, QN0, 3);
    }

    public static Ret.pfqnAMVA pfqn_dmlin(Matrix L, Matrix N, Matrix Zin, double tol, int maxiter,
                                          Matrix QN0, int npasses) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        Matrix Z = (Zin == null || Zin.isEmpty()) ? new Matrix(1, R) : Zin;

        boolean allZero = true;
        for (int i = 0; i < M && allZero; i++) {
            for (int r = 0; r < R; r++) {
                if (L.get(i, r) != 0.0) {
                    allZero = false;
                    break;
                }
            }
        }
        if (M == 0 || allZero) {
            Matrix X0 = new Matrix(1, R);
            Matrix U0 = new Matrix(M, R);
            for (int r = 0; r < R; r++) {
                X0.set(r, Z.get(r) > 0 ? N.get(r) / Z.get(r) : 0.0);
                for (int i = 0; i < M; i++) {
                    U0.set(i, r, X0.get(r) * L.get(i, r));
                }
            }
            return new Ret.pfqnAMVA(new Matrix(M, R), U0, new Matrix(M, R), null,
                    new Matrix(1, R), X0, 0);
        }

        // Initialize, as Linearizer does, from Bard-Schweitzer at every population
        Matrix[] Qs = new Matrix[R + 1];
        for (int s = 0; s <= R; s++) {
            Matrix N1 = oner(N, s);
            Ret.pfqnAMVA seed = (QN0 == null || QN0.isEmpty())
                    ? Pfqn_bs.pfqn_bs(L, N1, Z)
                    : Pfqn_bs.pfqn_bs(L, N1, Z, tol, maxiter, QN0.copy());
            Qs[s] = seed.Q.copy();
        }
        // Delta[r][c] holds the column of Delta^(r)_c over the stations
        double[][][] Delta = new double[M][R][R];
        Matrix xi = new Matrix(M, R);

        int totiter = 0;
        for (int pass = 0; pass < npasses; pass++) {
            for (int s = 0; s <= R; s++) {
                Matrix N1 = oner(N, s);
                Matrix xis = new Matrix(M, R);
                for (int i = 0; i < M; i++) {
                    for (int c = 0; c < R; c++) {
                        // xi at population N - 1_s, exactly; s == 0 leaves xi at N
                        double v = xi.get(i, c);
                        if (s > 0) {
                            v -= Delta[i][s - 1][c];
                        }
                        xis.set(i, c, v);
                    }
                }
                CoreResult cr = core(L, M, R, N1, Z, Qs[s], xis, tol, maxiter - totiter);
                Qs[s] = cr.Q;
                totiter += cr.iter;
            }
            // Refresh the Delta-terms, then aggregate them into xi once per pass
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < R; r++) {
                    if (N.get(r) == 1.0) {
                        Qs[r + 1].set(i, r, 0.0);
                    }
                    for (int s = 1; s <= R; s++) {
                        double ns = N.get(r) - (r == s - 1 ? 1.0 : 0.0);
                        if (N.get(r) > 0 && ns > 0) {
                            Delta[i][r][s - 1] = Qs[s].get(i, r) / ns - Qs[0].get(i, r) / N.get(r);
                        } else if (N.get(r) > 0) {
                            Delta[i][r][s - 1] = -Qs[0].get(i, r) / N.get(r);
                        } else {
                            Delta[i][r][s - 1] = 0.0;
                        }
                    }
                }
            }
            for (int i = 0; i < M; i++) {
                for (int c = 0; c < R; c++) {
                    double acc = 0.0;
                    for (int r = 0; r < R; r++) {
                        double w = N.get(r) - (r == c ? 1.0 : 0.0);
                        if (w > 0) {
                            acc += w * Delta[i][r][c];
                        }
                    }
                    xi.set(i, c, acc);
                }
            }
        }

        CoreResult fin = core(L, M, R, N, Z, Qs[0], xi, tol, maxiter - totiter);
        totiter += fin.iter;
        Matrix U = new Matrix(M, R);
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                U.set(i, r, fin.X.get(r) * L.get(i, r));
            }
        }
        Matrix C = new Matrix(1, R);
        for (int r = 0; r < R; r++) {
            C.set(r, fin.X.get(r) > 0 ? N.get(r) / fin.X.get(r) - Z.get(r) : 0.0);
        }
        return new Ret.pfqnAMVA(fin.Q, U, fin.W, null, C, fin.X, totiter);
    }

    private static Matrix oner(Matrix N, int s) {
        Matrix out = N.copy();
        if (s > 0) {
            out.set(s - 1, N.get(s - 1) - 1);
        }
        return out;
    }

    private static final class CoreResult {
        Matrix Q;
        Matrix W;
        Matrix X;
        int iter;
    }

    /** Fixed point of the aggregated arrival-instant estimate with the MVA equations. */
    private static CoreResult core(Matrix L, int M, int R, Matrix N1, Matrix Z,
                                   Matrix Qin, Matrix xi, double tol, int maxiter) {
        Matrix Q = Qin.copy();
        Matrix W = new Matrix(M, R);
        Matrix T = new Matrix(1, R);
        int iter = 0;
        boolean hasConverged = false;
        while (!hasConverged) {
            Matrix Qlast = Q.copy();
            for (int c = 0; c < R; c++) {
                for (int i = 0; i < M; i++) {
                    double acc = 0.0;
                    for (int r = 0; r < R; r++) {
                        if (N1.get(r) > 0) {
                            double nr = N1.get(r) - (r == c ? 1.0 : 0.0);
                            if (nr > 0) {
                                acc += nr * Q.get(i, r) / N1.get(r);
                            }
                        }
                    }
                    W.set(i, c, L.get(i, c) * (1 + acc + xi.get(i, c)));
                }
            }
            for (int r = 0; r < R; r++) {
                double wsum = 0.0;
                for (int i = 0; i < M; i++) {
                    wsum += W.get(i, r);
                }
                T.set(r, N1.get(r) > 0 ? N1.get(r) / (Z.get(r) + wsum) : 0.0);
                for (int i = 0; i < M; i++) {
                    Q.set(i, r, T.get(r) * W.get(i, r));
                }
            }
            double nrm = 0.0;
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < R; r++) {
                    double d = Q.get(i, r) - Qlast.get(i, r);
                    nrm += d * d;
                }
            }
            if (FastMath.sqrt(nrm) < tol || iter > maxiter) {
                hasConverged = true;
            }
            iter++;
        }
        CoreResult out = new CoreResult();
        out.Q = Q;
        out.W = W;
        out.X = T;
        out.iter = iter;
        return out;
    }
}
