/**
 * @file Aggregate Queue Length approximate MVA for closed queueing networks
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva;

import java.util.ArrayList;
import java.util.List;

import jline.io.Ret;
import jline.util.matrix.Matrix;

public final class Pfqn_aql {
    private Pfqn_aql() {}

    public static Ret.pfqnAMVA pfqn_aql(Matrix L, Matrix N, Matrix Z, int maxiter) {
        return pfqn_aql(L, N, Z, 1e-7, maxiter);
    }

    public static Ret.pfqnAMVA pfqn_aql(Matrix L, Matrix N) {
        return pfqn_aql(L, N, new Matrix(1, L.getNumCols()), 1e-7, 1000);
    }

    public static Ret.pfqnAMVA pfqn_aql(Matrix L, Matrix N, Matrix Z) {
        return pfqn_aql(L, N, Z, 1e-7, 1000);
    }

    public static Ret.pfqnAMVA pfqn_aql(Matrix Lin, Matrix N, Matrix Zin, double tol, int maxiter) {
        return pfqn_aql(Lin, N, Zin, tol, maxiter, null);
    }

    public static Ret.pfqnAMVA pfqn_aql(Matrix Lin, Matrix N, Matrix Zin, double tol, int maxiter, Matrix QN0ext) {
        Matrix L = Lin;
        Matrix Z = Zin;
        int M = L.getNumRows();
        int K = L.getNumCols();

        if (Z.isEmpty()) {
            Z = new Matrix(1, K);
        }

        Matrix QN0;
        if (QN0ext != null && !QN0ext.isEmpty()) {
            QN0 = QN0ext;   // warm start from supplied queue lengths
        } else {
            QN0 = new Matrix(M, K);
            double value = N.sumRows(0) / M;
            for (int i = 0; i < M; i++) {
                for (int j = 0; j < K; j++) {
                    QN0.set(i, j, value);
                }
            }
        }

        List<Matrix> Q = new ArrayList<Matrix>(K + 1);
        List<Matrix> R = new ArrayList<Matrix>(K + 1);
        List<Matrix> X = new ArrayList<Matrix>(K + 1);
        for (int i = 0; i <= K; i++) {
            Q.add(new Matrix(M, 1));
            R.add(new Matrix(M, K));
            X.add(new Matrix(1, K));
        }
        Matrix gamma = new Matrix(M, K);

        for (int t = 0; t <= K; t++) {
            for (int k = 0; k < M; k++) {
                Q.get(t).set(k, 0, QN0.get(k, 0));
            }
        }

        int it = 0;
        while (true) {
            List<Matrix> Q_olditer = new ArrayList<Matrix>();
            for (Matrix q : Q) {
                Q_olditer.add(q.copy());
            }
            it++;

            for (int t = 0; t <= K; t++) {
                Matrix n;
                if (t > 0) {
                    n = Matrix.oner(N, Integer.valueOf(t - 1));
                } else {
                    n = N;
                }
                for (int k = 0; k < M; k++) {
                    for (int s = 0; s < K; s++) {
                        double RValue = L.get(k, s) * (1 + (n.elementSum() - 1) * (Q.get(t).get(k, 0) / n.elementSum() - gamma.get(k, s)));
                        R.get(t).set(k, s, RValue);
                    }
                }

                for (int s = 0; s < K; s++) {
                    double sumR = R.get(t).sumCols(s);
                    double XValue = n.get(0, s) / (Z.get(0, s) + sumR);
                    X.get(t).set(0, s, XValue);
                }

                for (int k = 0; k < M; k++) {
                    double QValue = 0.0;
                    for (int s = 0; s < K; s++) {
                        QValue += X.get(t).get(0, s) * R.get(t).get(k, s);
                    }
                    Q.get(t).set(k, 0, QValue);
                }
            }

            for (int k = 0; k < M; k++) {
                for (int s = 0; s < K; s++) {
                    double gammaValue = (Q.get(0).get(k, 0) / N.elementSum()) - (Q.get(s + 1).get(k, 0) / (N.elementSum() - 1));
                    gamma.set(k, s, gammaValue);
                }
            }

            if (Matrix.maxAbsDiff(Q_olditer.get(0), Q.get(0)) < tol || it == maxiter) {
                break;
            }
        }

        Matrix XN = X.get(0);
        Matrix RN = R.get(0);
        Matrix UN = new Matrix(M, K);
        Matrix QN = new Matrix(M, K);
        Matrix AN = new Matrix(M, K);
        Matrix TN = new Matrix(M, K);

        for (int k = 0; k < M; k++) {
            for (int s = 0; s < K; s++) {
                double UNValue = XN.get(0, s) * L.get(k, s);
                UN.set(k, s, UNValue);
                double QNValue = UNValue * (1 + Q.get(s + 1).get(k, 0));
                QN.set(k, s, QNValue);
                AN.set(k, s, Q.get(s + 1).get(k, 0));
                TN.set(k, s, AN.get(k, s));
            }
        }

        Matrix CN = new Matrix(1, K);
        for (int r = 0; r < K; r++) {
            CN.set(0, r, N.get(0, r) / XN.get(0, r));
        }

        return new Ret.pfqnAMVA(QN, UN, RN, TN, CN, XN, it);
    }
}
