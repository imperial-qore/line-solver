/**
 * @file Linearizer++ approximate MVA with higher-order moment corrections
 *
 * Implements the Linearizer++ algorithm for closed queueing networks with enhanced accuracy
 * through higher-order moment corrections. Supports multiple approximation levels with
 * configurable precision-performance trade-offs for complex multi-class systems.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva;

import jline.io.Ret;
import jline.util.matrix.Matrix;

/**
 * Linearizer++ algorithm for closed networks without think times.
 */
public final class Pfqn_linearizerpp {
    private Pfqn_linearizerpp() {}

    public static Ret.pfqnAMVA pfqn_linearizerpp(Matrix L, Matrix N, int level) {
        return pfqn_linearizerpp(L, N, new Matrix(1, L.getNumCols()), level, 1e-4, 1000, 0);
    }

    public static Ret.pfqnAMVA pfqn_linearizerpp(Matrix L, Matrix N, Matrix Z, int level) {
        return pfqn_linearizerpp(L, N, Z, level, 1e-4, 1000, 0);
    }

    public static Ret.pfqnAMVA pfqn_linearizerpp(Matrix L, Matrix N, Matrix Z, int level, double tol) {
        return pfqn_linearizerpp(L, N, Z, level, tol, 1000, 0);
    }

    public static Ret.pfqnAMVA pfqn_linearizerpp(Matrix L, Matrix N, Matrix Z, int level, double tol, int maxiter) {
        return pfqn_linearizerpp(L, N, Z, level, tol, maxiter, 0);
    }

    public static Ret.pfqnAMVA pfqn_linearizerpp(Matrix L, Matrix N, Matrix Z, int level,
                                                 double tol, int maxiter, int flag) {
        int M = L.getNumRows();
        int R = L.getNumCols();

        if (Z == null || Z.isEmpty()) {
            Z = new Matrix(1, R);
        }

        double[] Nvec = N.toArray1D();
        double[] Zvec = Z.toArray1D();
        double[][] Lmat = L.toArray2D();
        int J;
        if (level == 1) {
            J = 1;
        } else if (level == 2) {
            J = 1 + R;
        } else {
            J = (R + 1) * (R + 2) / 2;
        }

        double[][][] X = new double[R][M][J];
        double[][][] Xn = new double[R][M][J];
        double[][][] Y = new double[R][J][J];
        double[][][] W = new double[R][J][J];

        double Ntot = 0.0;
        for (double v : Nvec) {
            Ntot += v;
        }

        for (int r = 0; r < R; r++) {
            for (int i = 0; i < M; i++) {
                for (int v = 0; v < J; v++) {
                    X[r][i][v] = nvr(Nvec, v + 1, r + 1, R) / M;
                }
            }
        }

        for (int r = 0; r < R; r++) {
            for (int u = 0; u < J; u++) {
                for (int v = 0; v < J; v++) {
                    Y[r][u][v] = ff(Nvec, r + 1, u + 1, v + 1, R);
                }
            }
        }

        for (int r = 0; r < R; r++) {
            double[][] Ymat = new double[J][J];
            for (int i = 0; i < J; i++) {
                for (int j = 0; j < J; j++) {
                    Ymat[i][j] = Y[r][i][j];
                }
            }
            Matrix Ymatrix = new Matrix(Ymat);
            Matrix Winv = Matrix.inv(Ymatrix);
            for (int i = 0; i < J; i++) {
                for (int j = 0; j < J; j++) {
                    W[r][i][j] = Winv.get(i, j);
                }
            }
        }

        double err = tol + 1;
        int iter = 0;
        while (err > tol && iter < maxiter) {
            for (int r = 0; r < R; r++) {
                double[] N1 = Nvec.clone();
                N1[r] -= 1.0;
                for (int v = 0; v < J; v++) {
                    double den = Zvec[r];
                    for (int j = 0; j < M; j++) {
                        double tmp = 1.0;
                        for (int s = 0; s < R; s++) {
                            double[] vettX = X[s][j];
                            for (int u = 0; u < J; u++) {
                                double f = ff(N1, s + 1, u + 1, v + 1, R);
                                double[] vettW = new double[J];
                                for (int k = 0; k < J; k++) {
                                    vettW[k] = W[s][k][u];
                                }
                                double inner = 0.0;
                                for (int k = 0; k < J; k++) {
                                    inner += vettW[k] * vettX[k];
                                }
                                tmp += inner * f;
                            }
                        }
                        den += tmp * Lmat[j][r];
                    }
                    for (int i = 0; i < M; i++) {
                        double tmp = 1.0;
                        for (int s = 0; s < R; s++) {
                            double[] vettX = X[s][i];
                            for (int u = 0; u < J; u++) {
                                double f = ff(N1, s + 1, u + 1, v + 1, R);
                                double[] vettW = new double[J];
                                for (int k = 0; k < J; k++) {
                                    vettW[k] = W[s][k][u];
                                }
                                double inner = 0.0;
                                for (int k = 0; k < J; k++) {
                                    inner += vettW[k] * vettX[k];
                                }
                                tmp += inner * f;
                            }
                        }
                        Xn[r][i][v] = tmp * Lmat[i][r] * nvr(Nvec, v + 1, r + 1, R) / den;
                    }
                }
            }

            err = 0.0;
            for (int r = 0; r < R; r++) {
                for (int i = 0; i < M; i++) {
                    err = Math.max(err, Math.abs(X[r][i][0] - Xn[r][i][0]));
                }
            }
            err /= Ntot;
            iter++;
            for (int r = 0; r < R; r++) {
                for (int i = 0; i < M; i++) {
                    System.arraycopy(Xn[r][i], 0, X[r][i], 0, J);
                }
            }
        }

        Matrix Q = new Matrix(M, R);
        for (int r = 0; r < R; r++) {
            for (int i = 0; i < M; i++) {
                Q.set(i, r, X[r][i][0]);
            }
        }

        double Zsum = 0.0;
        for (double z : Zvec) {
            Zsum += z;
        }
        Matrix Xres = new Matrix(1, R);
        if (Zsum == 0.0) {
            double Nsum = 0.0;
            for (double val : Nvec) {
                Nsum += val;
            }
            for (int r = 0; r < R; r++) {
                double sumQdivL = 0.0;
                for (int i = 0; i < M; i++) {
                    sumQdivL += Q.get(i, r) / Lmat[i][r];
                }
                Xres.set(0, r, sumQdivL / (Nsum + M - 1));
            }
        } else {
            throw new RuntimeException("pfqn_linearizerpp does not support think times");
        }

        Matrix U = new Matrix(M, R);
        Matrix Rmat = new Matrix(M, R);
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                U.set(i, r, Xres.get(0, r) * L.get(i, r));
                if (Xres.get(0, r) > 0) {
                    Rmat.set(i, r, Q.get(i, r) / Xres.get(0, r));
                }
            }
        }

        return new Ret.pfqnAMVA(Q, U, Rmat, null, null, Xres, iter);
    }

    private static double[] nv(double[] N, int v, int R) {
        double[] res = N.clone();
        if (v - 1 != 0) {
            v = v - 1;
            if (v <= R) {
                res[v - 1] -= 1.0;
            } else {
                v = v - R;
                outer:
                for (int i = 1; i <= R; i++) {
                    for (int j = i; j <= R; j++) {
                        v -= 1;
                        if (v == 0) {
                            res[i - 1] -= 1.0;
                            res[j - 1] -= 1.0;
                            break outer;
                        }
                    }
                }
            }
        }
        return res;
    }

    private static double nvr(double[] N, int v, int r, int R) {
        return nv(N, v, R)[r - 1];
    }

    private static double ff(double[] N, int r, int u, int v, int R) {
        double[] Nv = nv(N, v, R);
        if (u == 1) {
            return Nv[r - 1];
        } else {
            u -= 1;
            if (u <= R) {
                return Nv[r - 1] * Nv[u - 1];
            } else {
                u -= R;
                for (int i = 1; i <= R; i++) {
                    for (int j = i; j <= R; j++) {
                        u -= 1;
                        if (u == 0) {
                            return Nv[r - 1] * Nv[i - 1] * Nv[j - 1];
                        }
                    }
                }
            }
        }
        return 0.0;
    }
}
