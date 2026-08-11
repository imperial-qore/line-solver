/**
 * @file Bard-Schweitzer approximate MVA for FCFS scheduling with weighted priorities
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva;

import org.apache.commons.math3.util.FastMath;

import jline.io.Ret;
import jline.util.Maths;
import jline.util.matrix.Matrix;

public final class Pfqn_bsfcfs {
    private Pfqn_bsfcfs() {}

    /**
     * Bard-Schweitzer approximate MVA for FCFS scheduling with weighted priorities.
     */
    public static Ret.pfqnAMVA pfqn_bsfcfs(Matrix L,
                                           Matrix N,
                                           Matrix Z,
                                           double tol,
                                           int maxiter,
                                           Matrix QN0,
                                           Matrix weight) {
        int M = L.getNumRows();
        int R = L.getNumCols();

        Matrix Zmat = (Z != null) ? Z : new Matrix(1, R);

        Matrix CN = new Matrix(M, R);
        Matrix QN;
        if (QN0 == null || QN0.isEmpty()) {
            QN = N.repmat(M, 1);
            for (int i = 0; i < QN.getNumRows(); i++) {
                for (int j = 0; j < QN.getNumCols(); j++) {
                    QN.set(i, j, QN.get(i, j) / M);
                }
            }
        } else {
            QN = QN0.copy();
            for (int i = 0; i < QN.getNumRows(); i++) {
                for (int j = 0; j < QN.getNumCols(); j++) {
                    QN.set(i, j, QN.get(i, j) + FastMath.ulp(1.0));
                }
            }
        }

        Matrix weightMatrix = (weight != null) ? weight : Matrix.ones(M, R);

        Matrix XN = new Matrix(1, R);
        Matrix UN = new Matrix(M, R);
        Matrix relprio = new Matrix(M, R);

        int it = 1;
        while (it <= maxiter) {
            Matrix QN_1 = new Matrix(QN);

            for (int ist = 0; ist < M; ist++) {
                for (int r = 0; r < R; r++) {
                    relprio.set(ist, r, QN.get(ist, r) * weightMatrix.get(ist, r));
                }
            }

            for (int r = 0; r < R; r++) {
                for (int ist = 0; ist < M; ist++) {
                    CN.set(ist, r, L.get(ist, r));
                    for (int s = 0; s < R; s++) {
                        if (s != r) {
                            CN.set(ist, r,
                                CN.get(ist, r) + L.get(ist, s) * QN.get(ist, s) * relprio.get(ist, s) / relprio.get(ist, r));
                        } else {
                            CN.set(ist, r,
                                CN.get(ist, r) + L.get(ist, r) * QN.get(ist, r) * (N.get(r) - 1) / N.get(r) * relprio.get(ist, s) / relprio.get(ist, r));
                        }
                    }
                }
                XN.set(r, N.get(r) / (Zmat.get(r) + Matrix.extractColumn(CN, r, null).elementSum()));
            }
            for (int r = 0; r < R; r++) {
                for (int ist = 0; ist < M; ist++) {
                    QN.set(ist, r, XN.get(r) * CN.get(ist, r));
                }
            }
            for (int r = 0; r < R; r++) {
                for (int ist = 0; ist < M; ist++) {
                    UN.set(ist, r, XN.get(r) * L.get(ist, r));
                }
            }

            double maxabs = Double.MIN_VALUE;
            for (int i = 0; i < QN.getNumRows(); i++) {
                for (int j = 0; j < QN.getNumCols(); j++) {
                    double absValue = FastMath.abs(1 - QN.get(i, j) / QN_1.get(i, j));
                    maxabs = Maths.max(maxabs, absValue);
                }
            }
            if (maxabs < tol) {
                break;
            }
            it++;
        }

        Matrix RN = XN.repmat(M, 1);
        for (int i = 0; i < RN.getNumRows(); i++) {
            for (int j = 0; j < RN.getNumCols(); j++) {
                RN.set(i, j, QN.get(i, j) / RN.get(i, j));
            }
        }
        return new Ret.pfqnAMVA(QN, UN, RN, null, CN, XN, it);
    }

    public static Ret.pfqnAMVA pfqn_bsfcfs(Matrix L, Matrix N, Matrix Z, double tol, int maxiter, Matrix QN0) {
        return pfqn_bsfcfs(L, N, Z, tol, maxiter, QN0, null);
    }

    public static Ret.pfqnAMVA pfqn_bsfcfs(Matrix L, Matrix N, Matrix Z, double tol, int maxiter) {
        return pfqn_bsfcfs(L, N, Z, tol, maxiter, null, null);
    }

    public static Ret.pfqnAMVA pfqn_bsfcfs(Matrix L, Matrix N, Matrix Z, double tol) {
        return pfqn_bsfcfs(L, N, Z, tol, 1000, null, null);
    }

    public static Ret.pfqnAMVA pfqn_bsfcfs(Matrix L, Matrix N, Matrix Z) {
        return pfqn_bsfcfs(L, N, Z, 1.0e-6, 1000, null, null);
    }

    public static Ret.pfqnAMVA pfqn_bsfcfs(Matrix L, Matrix N) {
        return pfqn_bsfcfs(L, N, null, 1.0e-6, 1000, null, null);
    }
}
