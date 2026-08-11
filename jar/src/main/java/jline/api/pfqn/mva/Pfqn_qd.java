/**
 * @file Queue-Dependent (QD) approximate MVA solver
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva;

import org.apache.commons.math3.util.FastMath;

import jline.io.Ret;
import jline.util.matrix.Matrix;

public final class Pfqn_qd {
    private Pfqn_qd() {}

    public interface GammaFunction {
        Matrix apply(Matrix A);
    }

    public interface BetaFunction {
        Matrix apply(Matrix Akr);
    }

    public static Ret.pfqnQd pfqn_qd(Matrix L, Matrix N) {
        return pfqn_qd(L, N, null, null, null);
    }

    public static Ret.pfqnQd pfqn_qd(Matrix L, Matrix N, GammaFunction ga) {
        return pfqn_qd(L, N, ga, null, null);
    }

    public static Ret.pfqnQd pfqn_qd(Matrix L, Matrix N, GammaFunction ga, BetaFunction be) {
        return pfqn_qd(L, N, ga, be, null);
    }

    /**
     * Queue-Dependent (QD) approximate MVA solver
     */
    public static Ret.pfqnQd pfqn_qd(Matrix L, Matrix N, GammaFunction ga, BetaFunction be, Matrix Q0) {
        final int M = L.getNumRows();
        final int R = L.getNumCols();

        GammaFunction gaFunc = (ga != null) ? ga : new GammaFunction() {
            @Override
            public Matrix apply(Matrix A) {
                return Matrix.ones(M, 1);
            }
        };

        BetaFunction beFunc = (be != null) ? be : new BetaFunction() {
            @Override
            public Matrix apply(Matrix Akr) {
                return Matrix.ones(M, R);
            }
        };

        Matrix Q;
        if (Q0 != null) {
            Q = Q0.copy();
        } else {
            Q = new Matrix(M, R);
            for (int r = 0; r < R; r++) {
                double sumLr = 0.0;
                for (int i = 0; i < M; i++) {
                    sumLr += L.get(i, r);
                }
                if (sumLr > 0) {
                    for (int i = 0; i < M; i++) {
                        Q.set(i, r, L.get(i, r) / sumLr * N.get(r));
                    }
                }
            }
        }

        double Ntot = N.elementSum();
        double delta = (Ntot - 1) / Ntot;
        Matrix deltar = new Matrix(1, R);
        for (int r = 0; r < R; r++) {
            deltar.set(r, (N.get(r) - 1) / N.get(r));
        }

        Matrix Q_1 = Q.copy();
        for (int i = 0; i < Q_1.getNumRows(); i++) {
            for (int j = 0; j < Q_1.getNumCols(); j++) {
                Q_1.set(i, j, Q_1.get(i, j) * 10);
            }
        }

        double tol = 1e-6;
        int iter = 0;
        Matrix C = new Matrix(M, R);
        Matrix X = new Matrix(1, R);
        Matrix U = new Matrix(M, R);

        double maxDiff = Double.MAX_VALUE;
        while (maxDiff > tol) {
            iter++;

            for (int i = 0; i < M; i++) {
                for (int j = 0; j < R; j++) {
                    Q_1.set(i, j, Q.get(i, j));
                }
            }

            Matrix[] Ak = new Matrix[R];
            for (int r = 0; r < R; r++) {
                Ak[r] = new Matrix(M, 1);
            }
            Matrix Akr = new Matrix(M, R);
            for (int k = 0; k < M; k++) {
                double sumQk = 0.0;
                for (int s = 0; s < R; s++) {
                    sumQk += Q.get(k, s);
                }
                for (int r = 0; r < R; r++) {
                    Ak[r].set(k, 0, 1 + delta * sumQk);
                    Akr.set(k, r, 1 + deltar.get(r) * Q.get(k, r));
                }
            }

            for (int r = 0; r < R; r++) {
                Matrix g = gaFunc.apply(Ak[r]);
                Matrix b = beFunc.apply(Akr);

                double sumCr = 0.0;
                for (int k = 0; k < M; k++) {
                    double sumQk = 0.0;
                    for (int s = 0; s < R; s++) {
                        sumQk += Q.get(k, s);
                    }
                    C.set(k, r, L.get(k, r) * g.get(k, 0) * b.get(k, r) * (1 + delta * sumQk));
                    sumCr += C.get(k, r);
                }

                X.set(r, N.get(r) / sumCr);

                for (int k = 0; k < M; k++) {
                    Q.set(k, r, X.get(r) * C.get(k, r));
                    U.set(k, r, L.get(k, r) * g.get(k, 0) * b.get(k, r) * X.get(r));
                }
            }

            maxDiff = 0.0;
            for (int i = 0; i < M; i++) {
                for (int j = 0; j < R; j++) {
                    double diff = FastMath.abs(Q.get(i, j) - Q_1.get(i, j));
                    if (diff > maxDiff) {
                        maxDiff = diff;
                    }
                }
            }
        }

        return new Ret.pfqnQd(Q, X, U, iter);
    }
}
