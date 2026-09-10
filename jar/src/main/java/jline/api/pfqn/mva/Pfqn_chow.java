/**
 * @file Chow Second Approximation (SA) approximate Mean Value Analysis
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
 * Chow Second Approximation (SA) approximate MVA.
 *
 * <p>W.-M. Chow, "Approximations for large scale closed queueing networks",
 * Perform. Eval. 3(1), 1983. The arrival-instant queue length is written
 * exactly as</p>
 *
 * <pre>
 * A_k^(c)(N) = Q_k(N - 1_c) = Q_k(N) (1 + theta_ck),
 * theta_ck   = [Q_k(N - 1_c) - Q_k(N)] / Q_k(N),
 * </pre>
 *
 * <p>and the theta-terms are estimated ONCE, off the Bard LCP solution, before
 * the fixed point is run. Two estimators are given: the BACKWARD one uses
 * Qhat(N - 1_c), the FORWARD one Qhat(N + 1_c). Chow reports the forward form
 * to be the more accurate of the two, so it is the default here. Setting every
 * theta to zero recovers {@link Pfqn_lcp}.</p>
 */
public final class Pfqn_chow {
    private Pfqn_chow() {}

    /** Estimator of the theta-terms. */
    public static final String FORWARD = "forward";
    /** Estimator of the theta-terms. */
    public static final String BACKWARD = "backward";

    public static Ret.pfqnAMVA pfqn_chow(Matrix L, Matrix N) {
        return pfqn_chow(L, N, new Matrix(1, L.getNumCols()));
    }

    public static Ret.pfqnAMVA pfqn_chow(Matrix L, Matrix N, Matrix Z) {
        return pfqn_chow(L, N, Z, 1.0e-6, 1000, null);
    }

    public static Ret.pfqnAMVA pfqn_chow(Matrix L, Matrix N, Matrix Z, double tol, int maxiter, Matrix QN0) {
        int M = L.getNumRows();
        SchedStrategy[] type = new SchedStrategy[M];
        Arrays.fill(type, SchedStrategy.PS);
        return pfqn_chow(L, N, Z, tol, maxiter, QN0, type, FORWARD);
    }

    public static Ret.pfqnAMVA pfqn_chow(Matrix L, Matrix N, Matrix Z, double tol, int maxiter,
                                         Matrix QN0, SchedStrategy[] type) {
        return pfqn_chow(L, N, Z, tol, maxiter, QN0, type, FORWARD);
    }

    public static Ret.pfqnAMVA pfqn_chow(Matrix L, Matrix N, Matrix Z, double tol, int maxiter,
                                         Matrix QN0, SchedStrategy[] type, String variant) {
        int M = L.getNumRows();
        int R = L.getNumCols();

        // theta-terms from the LCP solution
        Ret.pfqnAMVA base = Pfqn_lcp.pfqn_lcp(L, N, Z, tol, maxiter, QN0, type);
        double[] Qtot = new double[M];
        for (int i = 0; i < M; i++) {
            double acc = 0.0;
            for (int r = 0; r < R; r++) {
                acc += base.Q.get(i, r);
            }
            Qtot[i] = acc;
        }
        Matrix theta = new Matrix(M, R);
        for (int r = 0; r < R; r++) {
            if (N.get(r) == 0.0) {
                continue;
            }
            Matrix Np = N.copy();
            double[] ref = new double[M];
            double[] delta = new double[M];
            if (BACKWARD.equals(variant)) {
                Np.set(r, N.get(r) - 1);
                Ret.pfqnAMVA alt = Pfqn_lcp.pfqn_lcp(L, Np, Z, tol, maxiter, QN0, type);
                for (int i = 0; i < M; i++) {
                    double acc = 0.0;
                    for (int s = 0; s < R; s++) {
                        acc += alt.Q.get(i, s);
                    }
                    ref[i] = Qtot[i];
                    delta[i] = acc - Qtot[i];
                }
            } else {
                Np.set(r, N.get(r) + 1);
                Ret.pfqnAMVA alt = Pfqn_lcp.pfqn_lcp(L, Np, Z, tol, maxiter, QN0, type);
                for (int i = 0; i < M; i++) {
                    double acc = 0.0;
                    for (int s = 0; s < R; s++) {
                        acc += alt.Q.get(i, s);
                    }
                    ref[i] = acc;
                    delta[i] = Qtot[i] - acc;
                }
            }
            for (int i = 0; i < M; i++) {
                if (ref[i] > 0) {
                    theta.set(i, r, delta[i] / ref[i]);
                }
            }
        }

        // fixed point with A_k^(c) = Q_k (1 + theta_ck)
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
                        if (type[ist] == SchedStrategy.FCFS && s != r) {
                            CN.set(ist, r, CN.get(ist, r)
                                    + L.get(ist, s) * QN.get(ist, s) * (1 + theta.get(ist, r)));
                        } else {
                            CN.set(ist, r, CN.get(ist, r)
                                    + L.get(ist, r) * QN.get(ist, s) * (1 + theta.get(ist, r)));
                        }
                    }
                    // a theta below -1 would make the arrival-instant queue negative
                    if (CN.get(ist, r) < L.get(ist, r)) {
                        CN.set(ist, r, L.get(ist, r));
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
