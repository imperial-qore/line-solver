/**
 * @file Bard Large Customer Population (LCP) approximate Mean Value Analysis
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

/**
 * Bard Large Customer Population (LCP) approximate MVA.
 *
 * <p>Y. Bard, "Some extensions to multiclass queueing network analysis", in
 * Performance of Computer Systems, North-Holland, 1979. The first approximate
 * MVA algorithm: it estimates the arrival-instant queue length by the
 * time-averaged one WITHOUT removing the arriving customer,</p>
 *
 * <pre>A_k^(c)(N) = Q_k(N - 1_c) ~= Q_k(N) = sum_s Q_ks(N),</pre>
 *
 * <p>since with a large population one customer less cannot change the mean
 * queue lengths appreciably. Setting the Bard-Schweitzer proportional term
 * Q_kc(N)/N_c to zero recovers this algorithm, so LCP is uniformly more
 * pessimistic than {@link Pfqn_bs} and is inaccurate at small populations.</p>
 */
public final class Pfqn_lcp {
    private Pfqn_lcp() {}

    public static Ret.pfqnAMVA pfqn_lcp(Matrix L, Matrix N) {
        return pfqn_lcp(L, N, new Matrix(1, L.getNumCols()));
    }

    public static Ret.pfqnAMVA pfqn_lcp(Matrix L, Matrix N, Matrix Z) {
        return pfqn_lcp(L, N, Z, 1.0e-6, 1000, null);
    }

    public static Ret.pfqnAMVA pfqn_lcp(Matrix L, Matrix N, Matrix Z, double tol, int maxiter, Matrix QN0) {
        int M = L.getNumRows();
        SchedStrategy[] type = new SchedStrategy[M];
        Arrays.fill(type, SchedStrategy.PS);
        return pfqn_lcp(L, N, Z, tol, maxiter, QN0, type);
    }

    public static Ret.pfqnAMVA pfqn_lcp(Matrix L, Matrix N, Matrix Z, double tol, int maxiter,
                                        Matrix QN0, SchedStrategy[] type) {
        int M = L.getNumRows();
        int R = L.getNumCols();
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
        }
        Matrix CN = new Matrix(M, R);
        Matrix XN = new Matrix(1, R);
        Matrix UN = new Matrix(M, R);

        int it = 1;
        while (it <= maxiter) {
            Matrix QN_1 = new Matrix(QN);
            for (int r = 0; r < R; r++) {
                if (N.get(r) == 0.0) {
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
                        continue;
                    }
                    for (int s = 0; s < R; s++) {
                        // LCP: the arriving customer is NOT removed, so the class-r
                        // term carries no (N(r)-1)/N(r) factor
                        if (type[ist] == SchedStrategy.FCFS && s != r) {
                            CN.set(ist, r, CN.get(ist, r) + L.get(ist, s) * QN.get(ist, s));
                        } else {
                            CN.set(ist, r, CN.get(ist, r) + L.get(ist, r) * QN.get(ist, s));
                        }
                    }
                }
                XN.set(r, N.get(r) / (Z.get(r) + Matrix.extractColumn(CN, r, null).elementSum()));
            }
            for (int r = 0; r < R; r++) {
                for (int ist = 0; ist < M; ist++) {
                    QN.set(ist, r, XN.get(r) * CN.get(ist, r));
                    UN.set(ist, r, XN.get(r) * L.get(ist, r));
                }
            }
            double maxabs = Double.MIN_VALUE;
            for (int i = 0; i < QN.getNumRows(); i++) {
                for (int j = 0; j < QN.getNumCols(); j++) {
                    if (N.get(j) == 0.0) {
                        continue;
                    }
                    maxabs = Maths.max(maxabs, FastMath.abs(1 - QN.get(i, j) / QN_1.get(i, j)));
                }
            }
            if (maxabs < tol) {
                break;
            }
            it++;
        }
        Matrix RN = new Matrix(M, R);
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < R; j++) {
                RN.set(i, j, N.get(j) == 0.0 ? 0.0 : QN.get(i, j) / XN.get(j));
            }
        }
        return new Ret.pfqnAMVA(QN, UN, RN, null, CN, XN, it);
    }
}
