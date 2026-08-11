/**
 * @file Bard-Schweitzer approximate Mean Value Analysis with priority support
 *
 * Implements the classic Bard-Schweitzer approximate MVA algorithm for closed queueing networks
 * with optional weighted priority extensions.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva;

import java.util.Arrays;

import org.apache.commons.math3.util.FastMath;

import jline.io.Ret;
import jline.lang.constant.SchedStrategy;
import jline.util.Maths;
import jline.util.matrix.Matrix;

public final class Pfqn_bs {
    private Pfqn_bs() {}

    public static Ret.pfqnAMVA pfqn_bs(Matrix L, Matrix N) {
        Matrix Z = new Matrix(N.getNumRows(), N.getNumCols());
        return pfqn_bs(L, N, Z);
    }

    public static Ret.pfqnAMVA pfqn_bs(Matrix L, Matrix N, Matrix Z) {
        return pfqn_bs(L, N, Z, 1.0e-6, 1000, null);
    }

    public static Ret.pfqnAMVA pfqn_bs(Matrix L, Matrix N, Matrix Z, double tol, int maxiter) {
        return pfqn_bs(L, N, Z, tol, maxiter, null);
    }

    public static Ret.pfqnAMVA pfqn_bs(Matrix L, Matrix N, Matrix Z, double tol, int maxiter, Matrix QN0) {
        int M = L.getNumRows();
        if (QN0 == null || QN0.isEmpty()) {
            QN0 = N.repmat(M, 1);
            for (int i = 0; i < QN0.getNumRows(); i++) {
                for (int j = 0; j < QN0.getNumCols(); j++) {
                    QN0.set(i, j, QN0.get(i, j) / M);
                }
            }
        }
        SchedStrategy[] type = new SchedStrategy[M];
        Arrays.fill(type, SchedStrategy.PS);
        return pfqn_bs(L, N, Z, tol, maxiter, QN0, type);
    }

    /**
     * Bard-Schweitzer approximate mean value analysis algorithm with weighted priorities.
     */
    public static Ret.pfqnAMVA pfqn_bs(Matrix L, Matrix N, Matrix Z, double tol, int maxiter,
                                        Matrix QN0, Matrix weight) {
        int M = L.getNumRows();
        SchedStrategy[] type = new SchedStrategy[M];
        Arrays.fill(type, SchedStrategy.FCFS);
        return pfqn_bs(L, N, Z, tol, maxiter, QN0, type, weight);
    }

    public static Ret.pfqnAMVA pfqn_bs(Matrix L, Matrix N, Matrix Z, double tol, int maxiter,
                                        Matrix QN0, SchedStrategy[] type) {
        return pfqn_bs(L, N, Z, tol, maxiter, QN0, type, null);
    }

    /**
     * Bard-Schweitzer approximate mean value analysis algorithm with optional weighted priorities.
     */
    public static Ret.pfqnAMVA pfqn_bs(Matrix L, Matrix N, Matrix Z, double tol, int maxiter,
                                        Matrix QN0, SchedStrategy[] type, Matrix weight) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        if (QN0 == null || QN0.isEmpty()) {
            QN0 = N.repmat(M, 1);
            for (int i = 0; i < QN0.getNumRows(); i++) {
                for (int j = 0; j < QN0.getNumCols(); j++) {
                    QN0.set(i, j, QN0.get(i, j) / M);
                }
            }
        } else {
            // Add small epsilon to avoid zero problems as in MATLAB
            for (int i = 0; i < QN0.getNumRows(); i++) {
                for (int j = 0; j < QN0.getNumCols(); j++) {
                    QN0.set(i, j, QN0.get(i, j) + FastMath.ulp(1.0));
                }
            }
        }

        // Initialize weight matrix if not provided
        Matrix weightMatrix = (weight != null) ? weight : Matrix.ones(M, R);

        Matrix CN = new Matrix(M, R);
        Matrix QN = QN0;
        Matrix XN = new Matrix(1, R);
        Matrix UN = new Matrix(M, R);
        Matrix relprio = new Matrix(M, R);

        int it = 1;
        while (it <= maxiter) {
            Matrix QN_1 = new Matrix(QN);

            // Calculate relative priorities if using weighted FCFS
            boolean useWeightedFCFS = weight != null;
            if (useWeightedFCFS) {
                for (int ist = 0; ist < M; ist++) {
                    for (int r = 0; r < R; r++) {
                        relprio.set(ist, r, QN.get(ist, r) * weightMatrix.get(ist, r));
                    }
                }
            }

            for (int r = 0; r < R; r++) {
                if (N.get(r) == 0.0) {
                    // see _kb/03-api-layer.md for rationale
                    XN.set(r, 0.0);
                    for (int ist = 0; ist < M; ist++) {
                        CN.set(ist, r, 0.0);
                        QN.set(ist, r, 0.0);
                        UN.set(ist, r, 0.0);
                    }
                    continue;
                }
                for (int ist = 0; ist < M; ist++) {
                    CN.set(ist, r, L.get(ist, r));
                    if (L.get(ist, r) == 0.0) {
                        // 0 service demand at this station => this class does not visit the current node
                        continue;
                    }
                    for (int s = 0; s < R; s++) {
                        if (s != r) {
                            if (type[ist] == SchedStrategy.FCFS && useWeightedFCFS) {
                                // Weighted FCFS approximation
                                CN.set(ist, r, CN.get(ist, r) + L.get(ist, s) * QN.get(ist, s) * relprio.get(ist, s) / relprio.get(ist, r));
                            } else if (type[ist] == SchedStrategy.FCFS) {
                                // Standard FCFS approximation
                                CN.set(ist, r, CN.get(ist, r) + L.get(ist, s) * QN.get(ist, s));
                            } else {
                                // PS approximation
                                CN.set(ist, r, CN.get(ist, r) + L.get(ist, r) * QN.get(ist, s));
                            }
                        } else {
                            if (type[ist] == SchedStrategy.FCFS && useWeightedFCFS) {
                                CN.set(ist, r, CN.get(ist, r) + L.get(ist, r) * QN.get(ist, r) * (N.get(r) - 1) / N.get(r));
                            } else {
                                CN.set(ist, r, CN.get(ist, r) + L.get(ist, r) * QN.get(ist, r) * (N.get(r) - 1) / N.get(r));
                            }
                        }
                    }
                }
                XN.set(r, N.get(r) / (Z.get(r) + Matrix.extractColumn(CN, r, null).elementSum()));
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
            // Convergence is measured on the non-empty classes only: an empty class
            // has QN = QN_1 = 0, and 0/0 = NaN would make the test never fire.
            double maxabs = Double.MIN_VALUE;
            for (int i = 0; i < QN.getNumRows(); i++) {
                for (int j = 0; j < QN.getNumCols(); j++) {
                    if (N.get(j) == 0.0) {
                        continue;
                    }
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
                if (N.get(j) == 0.0) {
                    // 0/0 for an empty class; its residence time is 0, not NaN
                    RN.set(i, j, 0.0);
                } else {
                    RN.set(i, j, QN.get(i, j) / RN.get(i, j));
                }
            }
        }
        return new Ret.pfqnAMVA(QN, UN, RN, null, CN, XN, it);
    }
}
