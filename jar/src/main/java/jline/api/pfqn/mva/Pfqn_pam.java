/**
 * @file Hsieh-Lam Proportional Approximation Methods (PAMB/PAMI/PAMT)
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva;

import jline.io.Ret;
import jline.util.matrix.Matrix;

/**
 * Hsieh-Lam Proportional Approximation Methods (PAMB/PAMI/PAMT).
 *
 * <p>C. T. Hsieh, S. S. Lam, "PAM - A noniterative approximate solution method
 * for closed multichain queueing networks", ACM SIGMETRICS Perform. Eval. Rev.
 * 16(1), 1988. The three variants are NONITERATIVE: the queue lengths are
 * seeded by the proportion of a class demand that falls at each centre,</p>
 *
 * <pre>E_ck = D_ck / sum_i D_ci,   Q_ck(N) = E_ck N_c,</pre>
 *
 * <p>and the MVA equations are then unrolled a fixed number of times. PAMB
 * applies the last MVA step; PAMI additionally scales a class down wherever it
 * would drive a centre past full utilization; PAMT seeds at N - 1_i - 1_j and
 * applies the last TWO MVA steps before that capping. The seed spreads the
 * whole class population over the queueing centres and ignores Z, exactly as
 * published: PAM buys speed, not accuracy.</p>
 */
public final class Pfqn_pam {
    private Pfqn_pam() {}

    public static final String PAMB = "pamb";
    public static final String PAMI = "pami";
    public static final String PAMT = "pamt";

    public static Ret.pfqnAMVA pfqn_pam(Matrix L, Matrix N) {
        return pfqn_pam(L, N, new Matrix(1, L.getNumCols()), PAMB);
    }

    public static Ret.pfqnAMVA pfqn_pam(Matrix L, Matrix N, Matrix Z) {
        return pfqn_pam(L, N, Z, PAMB);
    }

    public static Ret.pfqnAMVA pfqn_pam(Matrix L, Matrix N, Matrix Zin, String variant) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        Matrix Z = (Zin == null || Zin.isEmpty()) ? new Matrix(1, R) : Zin;
        String var = (variant == null) ? PAMB : variant.toLowerCase();

        // E_ck, the share of the class-r demand served at station ist
        Matrix E = new Matrix(M, R);
        for (int r = 0; r < R; r++) {
            double tot = 0.0;
            for (int i = 0; i < M; i++) {
                tot += L.get(i, r);
            }
            if (tot > 0) {
                for (int i = 0; i < M; i++) {
                    E.set(i, r, L.get(i, r) / tot);
                }
            }
        }
        Matrix Q = new Matrix(M, R);   // Q_ck(N) = E_ck N_c
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                Q.set(i, r, E.get(i, r) * N.get(r));
            }
        }

        Matrix RN = new Matrix(M, R);
        Matrix XN = new Matrix(1, R);
        if (PAMT.equals(var)) {
            for (int i = 0; i < R; i++) {
                Matrix Qmi = new Matrix(M, R);   // Q_jk(N - 1_i)
                for (int j = 0; j < R; j++) {
                    // Q_ck(N - 1_i - 1_j) = Q_ck(N) - E_ck [(c==i) + (c==j)]
                    double[] agg = new double[M];
                    for (int k = 0; k < M; k++) {
                        double acc = 0.0;
                        for (int c = 0; c < R; c++) {
                            double q = Q.get(k, c);
                            if (c == i) q -= E.get(k, c);
                            if (c == j) q -= E.get(k, c);
                            acc += q;
                        }
                        agg[k] = acc;
                    }
                    double[] Rj = new double[M];
                    double Rjtot = 0.0;
                    for (int k = 0; k < M; k++) {
                        Rj[k] = L.get(k, j) * (1 + agg[k]);
                        Rjtot += Rj[k];
                    }
                    double nj = N.get(j) - (i == j ? 1.0 : 0.0);
                    double Xj = (nj > 0) ? nj / (Rjtot + Z.get(j)) : 0.0;
                    for (int k = 0; k < M; k++) {
                        Qmi.set(k, j, Xj * Rj[k]);
                    }
                }
                double Ritot = 0.0;
                for (int k = 0; k < M; k++) {
                    double acc = 0.0;
                    for (int c = 0; c < R; c++) {
                        acc += Qmi.get(k, c);
                    }
                    RN.set(k, i, L.get(k, i) * (1 + acc));
                    Ritot += RN.get(k, i);
                }
                if (N.get(i) > 0) {
                    XN.set(i, N.get(i) / (Ritot + Z.get(i)));
                }
            }
        } else {
            for (int r = 0; r < R; r++) {
                // Q_jk(N - 1_r) = Q_jk(N) - E_jk [j == r]
                double Rtot = 0.0;
                for (int k = 0; k < M; k++) {
                    double acc = 0.0;
                    for (int c = 0; c < R; c++) {
                        acc += Q.get(k, c) - (c == r ? E.get(k, c) : 0.0);
                    }
                    RN.set(k, r, L.get(k, r) * (1 + acc));
                    Rtot += RN.get(k, r);
                }
                if (N.get(r) > 0) {
                    XN.set(r, N.get(r) / (Rtot + Z.get(r)));
                }
            }
        }

        if (PAMI.equals(var) || PAMT.equals(var)) {
            // scale a class down when it would drive a centre it visits past U = 1
            double[] U = new double[M];
            for (int k = 0; k < M; k++) {
                double acc = 0.0;
                for (int c = 0; c < R; c++) {
                    acc += L.get(k, c) * XN.get(c);
                }
                U[k] = acc;
            }
            for (int r = 0; r < R; r++) {
                double S = Double.NEGATIVE_INFINITY;
                for (int k = 0; k < M; k++) {
                    if (L.get(k, r) != 0.0 && U[k] > S) {
                        S = U[k];
                    }
                }
                if (S > 1.0) {
                    XN.set(r, XN.get(r) / S);
                }
            }
        }

        Matrix QN = new Matrix(M, R);
        Matrix UN = new Matrix(M, R);
        for (int k = 0; k < M; k++) {
            for (int r = 0; r < R; r++) {
                QN.set(k, r, XN.get(r) * RN.get(k, r));
                UN.set(k, r, XN.get(r) * L.get(k, r));
            }
        }
        return new Ret.pfqnAMVA(QN, UN, RN, null, RN, XN, 1);
    }
}
