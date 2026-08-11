/**
 * @file Extended General-Form linearizer approximate MVA with class-specific parameters
 */
package jline.api.pfqn.mva;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import org.apache.commons.math3.util.FastMath;

import jline.io.Ret;
import jline.lang.constant.SchedStrategy;
import jline.util.matrix.Matrix;

public final class Pfqn_egflinearizer {
    private Pfqn_egflinearizer() {}

    public static Ret.pfqnAMVA pfqn_egflinearizer(Matrix L,
                                                   Matrix N,
                                                   Matrix Z,
                                                   SchedStrategy[] type,
                                                   double tol,
                                                   int maxiter,
                                                   Matrix alpha) {
        return pfqn_egflinearizer(L, N, Z, type, tol, maxiter, alpha, null);
    }

    public static Ret.pfqnAMVA pfqn_egflinearizer(Matrix L,
                                                   Matrix N,
                                                   Matrix Z,
                                                   SchedStrategy[] type,
                                                   double tol,
                                                   int maxiter,
                                                   Matrix alpha,
                                                   Matrix QN0) {
        int M = L.getNumRows();
        int R = L.getNumCols();

        if (Z == null || Z.isEmpty()) {
            Z = new Matrix(1, R);
        }
        Z = Z.sumCols();
        boolean Lmaxcols = true;
        for (int j = 0; j < L.getNumCols(); j++) {
            double maxcol = 0.0;
            for (int i = 0; i < L.getNumRows(); i++) {
                if (L.get(i, j) > maxcol) {
                    maxcol = L.get(i, j);
                }
            }
            if (maxcol != 0.0) {
                Lmaxcols = false;
                break;
            }
        }
        if (L.isEmpty() || Lmaxcols) {
            Matrix X = new Matrix(N.getNumRows(), N.getNumCols());
            for (int i = 0; i < X.getNumRows(); i++) {
                for (int j = 0; j < X.getNumCols(); j++) {
                    X.set(i, j, N.get(i, j) / Z.get(i, j));
                }
            }
            Matrix Q = new Matrix(M, R);
            Matrix U = new Matrix(M, R);
            Matrix W = new Matrix(M, R);
            Matrix T = new Matrix(M, R);
            Matrix C = new Matrix(1, R);
            for (int r = 0; r < R; r++) {
                for (int i = 0; i < M; i++) {
                    U.set(i, r, X.get(r) * L.get(i, r));
                }
            }
            for (int r = 0; r < R; r++) {
                for (int i = 0; i < M; i++) {
                    T.set(i, r, Q.get(i, r) / W.get(i, r));
                }
            }
            int totiter = 0;
            return new Ret.pfqnAMVA(Q, U, W, T, C, X, totiter);
        }

        // Initialise
        Matrix[] Q = new Matrix[M];
        Matrix[] Delta = new Matrix[M];
        for (int i = 0; i < M; i++) {
            Q[i] = new Matrix(R, 1 + R);
            Delta[i] = new Matrix(R, R);
        }

        for (int s = -1; s < R; s++) {
            ArrayList<Integer> sList = new ArrayList<Integer>(Collections.singletonList(Integer.valueOf(s)));
            Matrix N_1 = Matrix.oner(N, sList);
            Ret.pfqnAMVA init0 = (QN0 == null || QN0.isEmpty())
                    ? Pfqn_bs.pfqn_bs(L, N_1, Z)
                    : Pfqn_bs.pfqn_bs(L, N_1, Z, tol, maxiter, QN0);
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < R; r++) {
                    Q[i].set(r, 1 + s, init0.Q.get(i, r));
                }
            }
        }

        int totiter = 0;

        for (int I = 0; I < 3; I++) {
            for (int s = -1; s < R; s++) {
                ArrayList<Integer> sList = new ArrayList<Integer>(Collections.singletonList(Integer.valueOf(s)));
                Matrix N_1 = Matrix.oner(N, sList);
                Matrix Q1 = new Matrix(M, R);
                for (int i = 0; i < M; i++) {
                    for (int j = 0; j < R; j++) {
                        Q1.set(i, j, Q[i].get(j, 1 + s));
                    }
                }
                Ret.LinearizerResult ret1 = egflinearizer_core(L, M, R, N_1, Z, Q1, Delta, type, tol, maxiter - totiter, alpha);
                for (int i = 0; i < M; i++) {
                    for (int j = 0; j < R; j++) {
                        Q[i].set(j, 1 + s, ret1.Q.get(i, j));
                    }
                }
                totiter += ret1.iter;
            }
            // Upgrade delta
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < R; r++) {
                    if (N.get(r) == 1.0) {
                        // At population N-e_r only class r itself vanishes; the other
                        // classes keep the queue lengths Core just computed.
                        Q[i].set(r, 1 + r, 0);
                    }
                    for (int s = 0; s < R; s++) {
                        ArrayList<Integer> sList = new ArrayList<Integer>(Collections.singletonList(Integer.valueOf(s)));
                        Matrix Ns = Matrix.oner(N, sList);
                        if (Ns.get(r) > 0) {
                            Delta[i].set(r, s,
                                    Q[i].get(r, 1 + s) / FastMath.pow(Ns.get(r), alpha.get(r))
                                            - Q[i].get(r, 0) / FastMath.pow(N.get(r), alpha.get(r)));
                        } else {
                            // see _kb/03-api-layer.md for rationale
                            Delta[i].set(r, s, -Q[i].get(r, 0) / FastMath.pow(N.get(r), alpha.get(r)));
                        }
                    }
                }
            }
        }

        Matrix Q1 = new Matrix(M, R);
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < R; j++) {
                Q1.set(i, j, Q[i].get(j, 0));
            }
        }
        Ret.LinearizerResult ret1 = egflinearizer_core(L, M, R, N, Z, Q1, Delta, type, tol, maxiter - totiter, alpha);
        Matrix newQ = ret1.Q;
        Matrix W = ret1.W;
        Matrix X = ret1.T;
        totiter += ret1.iter;

        Matrix U = new Matrix(M, R);
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                U.set(i, r, X.get(r) * L.get(i, r));
            }
        }
        Matrix C = new Matrix(1, R);
        for (int i = 0; i < R; i++) {
            C.set(0, i, N.get(i) / X.get(i) - Z.get(i));
        }
        return new Ret.pfqnAMVA(newQ, U, W, null, C, X, totiter);
    }

    static Ret.LinearizerResult egflinearizer_core(Matrix L,
                                                    int M,
                                                    int R,
                                                    Matrix N_1,
                                                    Matrix Z,
                                                    Matrix Q,
                                                    Matrix[] Delta,
                                                    SchedStrategy[] type,
                                                    double tol,
                                                    int maxiter,
                                                    Matrix alpha) {
        boolean hasConverged = false;
        Matrix W = new Matrix(L);
        int iter = 0;
        Matrix T = null;
        while (!hasConverged) {
            Matrix Qlast = new Matrix(Q);
            Ret.pfqnLinearizerEstimate ret1 = egflinearizer_estimate(L, M, R, N_1, Z, Q, Delta, W, alpha);
            Matrix[] Q_1 = ret1.Q_1;
            Ret.LinearizerResult ret2 = egflinearizer_forwardMVA(L, M, R, type, N_1, Z, Q_1);
            Q = ret2.Q;
            W = ret2.W;
            T = ret2.T;
            if (Q.sub(Qlast).norm() < tol || iter > maxiter) {
                hasConverged = true;
            }
            iter++;
        }
        return new Ret.LinearizerResult(Q, W, T, iter);
    }

    static Ret.pfqnLinearizerEstimate egflinearizer_estimate(Matrix L,
                                                              int M,
                                                              int R,
                                                              Matrix N_1,
                                                              Matrix Z,
                                                              Matrix Q,
                                                              Matrix[] Delta,
                                                              Matrix W,
                                                              Matrix alpha) {
        Matrix[] Q_1 = new Matrix[M];
        for (int i = 0; i < M; i++) {
            Q_1[i] = new Matrix(R, 1 + R);
        }
        Matrix T_1 = new Matrix(R, 1 + R);
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                for (int s = 0; s < R; s++) {
                    ArrayList<Integer> sList = new ArrayList<Integer>(Collections.singletonList(Integer.valueOf(s)));
                    Matrix Ns = Matrix.oner(N_1, sList);
                    // see _kb/03-api-layer.md for rationale
                    if (N_1.get(r) <= 0 || Ns.get(r) <= 0) {
                        Q_1[i].set(r, 1 + s, 0);
                    } else {
                        Q_1[i].set(r, 1 + s,
                                FastMath.pow(Ns.get(r), alpha.get(r))
                                        * (Q.get(i, r) / FastMath.pow(N_1.get(r), alpha.get(r))
                                                + Delta[i].get(r, s)));
                    }
                }
            }
        }
        return new Ret.pfqnLinearizerEstimate(Q_1, T_1);
    }

    static Ret.LinearizerResult egflinearizer_forwardMVA(Matrix L,
                                                          int M,
                                                          int R,
                                                          SchedStrategy[] type,
                                                          Matrix N_1,
                                                          Matrix Z,
                                                          Matrix[] Q_1) {
        Matrix W = new Matrix(M, R);
        Matrix T = new Matrix(1, R);
        Matrix Q = new Matrix(M, R);

        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                double sumQ = 0.0;
                for (int k = 0; k < Q_1[i].getNumRows(); k++) {
                    sumQ += Q_1[i].get(k, 1 + r);
                }
                W.set(i, r, L.get(i, r) * (1 + sumQ));
            }
        }

        for (int r = 0; r < R; r++) {
            T.set(r, N_1.get(r) / (Z.get(r) + Matrix.extractColumn(W, r, null).elementSum()));
            for (int i = 0; i < M; i++) {
                Q.set(i, r, T.get(r) * W.get(i, r));
            }
        }
        return new Ret.LinearizerResult(Q, W, T);
    }
}
