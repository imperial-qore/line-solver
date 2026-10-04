/**
 * @file M3PP(2,m) exact count-based parameter fitting
 *
 * @since LINE 3.0
 */
package jline.api.mam.m3pp;

import jline.api.mam.Mmpp2_fitc;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class M3pp2m_fitc {
    private M3pp2m_fitc() {}

    // ---------------------------------------------------------------------
    // NOTE: the former (av, btv, binfv) convenience wrappers were removed. They
    // invented the statistics the fit needs but does not receive -- a third
    // moment 6a^3, time scales 1/2/3 and a per-class delta-variance of 0.1 --
    // so they fitted made-up data, and the second one silently aliased the
    // approx_ag_multiclass algorithm to this one. Call the full signature
    // below (mirrors m3pp2m_fitc.m), or M3pp2m_fitc_approx_ag_multiclass for
    // the approximate variant.
    // ---------------------------------------------------------------------

    /**
     * Fits a second-order Marked MMPP using exact count statistics.
     */
    public static Matrix[] m3pp2m_fitc(double a, double bt1, double bt2, double binf,
                                       double m3t2, double t1, double t2,
                                       double[] ai, double[] dvt3, double t3) {
        int m = ai.length;

        Matrix[] mmpp = Mmpp2_fitc.mmpp2_fitc(a, bt1, bt2, binf, m3t2, t1, t2);

        if (mmpp[0].getNumRows() == 1) {
            Matrix[] mmmpp = new Matrix[2 + m];
            mmmpp[0] = mmpp[0];
            mmmpp[1] = mmpp[1];
            for (int i = 0; i < m; i++) {
                mmmpp[2 + i] = new Matrix(new double[] { ai[i] });
            }
            return mmmpp;
        }

        Matrix[] mmmpp = new Matrix[2 + m];
        for (int idx = 0; idx < mmmpp.length; idx++) {
            mmmpp[idx] = new Matrix(2, 2);
        }
        mmmpp[0] = mmpp[0];
        mmmpp[1] = mmpp[1];

        double l1 = mmpp[1].get(0, 0);
        double l2 = mmpp[1].get(1, 1);
        double r1 = mmpp[0].get(0, 1);
        double r2 = mmpp[0].get(1, 0);

        double[][] q = new double[2][m];

        for (int i = 0; i < m - 1; i++) {
            double a_1 = ai[i];
            double dv_1 = dvt3[i];
            double t = t3;

            double sinhTerm = Math.sinh((r1 * t) / 2 + (r2 * t) / 2);
            double expTerm = Math.exp(-(r1 * t) / 2 - (r2 * t) / 2);
            double sinhExp = sinhTerm * expTerm;

            double numerator1 = -(dv_1 * Math.pow(r1, 4) + dv_1 * Math.pow(r2, 4) - 2 * a_1 * Math.pow(r1, 4) * t -
                    2 * a_1 * Math.pow(r2, 4) * t + 4 * dv_1 * r1 * Math.pow(r2, 3) + 4 * dv_1 * Math.pow(r1, 3) * r2 +
                    l1 * Math.pow(r2, 4) * t + l2 * Math.pow(r1, 4) * t + 6 * dv_1 * Math.pow(r1, 2) * Math.pow(r2, 2) +
                    4 * a_1 * l1 * Math.pow(r2, 3) * t - 4 * a_1 * l2 * Math.pow(r2, 3) * t -
                    8 * a_1 * r1 * Math.pow(r2, 3) * t - 8 * a_1 * Math.pow(r1, 3) * r2 * t +
                    3 * l1 * r1 * Math.pow(r2, 3) * t + l1 * Math.pow(r1, 3) * r2 * t +
                    l2 * r1 * Math.pow(r2, 3) * t + 3 * l2 * Math.pow(r1, 3) * r2 * t -
                    12 * a_1 * Math.pow(r1, 2) * Math.pow(r2, 2) * t + 3 * l1 * Math.pow(r1, 2) * Math.pow(r2, 2) * t +
                    2 * Math.pow(l1, 2) * r1 * Math.pow(r2, 2) * t + 2 * Math.pow(l1, 2) * Math.pow(r1, 2) * r2 * t +
                    3 * l2 * Math.pow(r1, 2) * Math.pow(r2, 2) * t + 2 * Math.pow(l2, 2) * r1 * Math.pow(r2, 2) * t +
                    2 * Math.pow(l2, 2) * Math.pow(r1, 2) * r2 * t - 8 * a_1 * l1 * Math.pow(r2, 2) * sinhExp +
                    8 * a_1 * l2 * Math.pow(r2, 2) * sinhExp - 4 * Math.pow(l1, 2) * r1 * r2 * sinhExp -
                    4 * Math.pow(l2, 2) * r1 * r2 * sinhExp + 8 * a_1 * l1 * r1 * Math.pow(r2, 2) * t +
                    4 * a_1 * l1 * Math.pow(r1, 2) * r2 * t - 8 * a_1 * l2 * r1 * Math.pow(r2, 2) * t -
                    4 * a_1 * l2 * Math.pow(r1, 2) * r2 * t - 4 * l1 * l2 * r1 * Math.pow(r2, 2) * t -
                    4 * l1 * l2 * Math.pow(r1, 2) * r2 * t - 8 * a_1 * l1 * r1 * r2 * sinhExp +
                    8 * a_1 * l2 * r1 * r2 * sinhExp + 8 * l1 * l2 * r1 * r2 * sinhExp);

            double denominator1 = 4 * l1 * r2 * (r1 + r2) * (2 * l1 * sinhExp -
                    2 * l2 * sinhExp - l1 * r1 * t - l1 * r2 * t + l2 * r1 * t + l2 * r2 * t);

            q[0][i] = numerator1 / denominator1;

            double numerator2 = dv_1 * Math.pow(r1, 4) + dv_1 * Math.pow(r2, 4) - 2 * a_1 * Math.pow(r1, 4) * t -
                    2 * a_1 * Math.pow(r2, 4) * t + 4 * dv_1 * r1 * Math.pow(r2, 3) + 4 * dv_1 * Math.pow(r1, 3) * r2 +
                    l1 * Math.pow(r2, 4) * t + l2 * Math.pow(r1, 4) * t + 6 * dv_1 * Math.pow(r1, 2) * Math.pow(r2, 2) -
                    4 * a_1 * l1 * Math.pow(r1, 3) * t + 4 * a_1 * l2 * Math.pow(r1, 3) * t -
                    8 * a_1 * r1 * Math.pow(r2, 3) * t - 8 * a_1 * Math.pow(r1, 3) * r2 * t +
                    3 * l1 * r1 * Math.pow(r2, 3) * t + l1 * Math.pow(r1, 3) * r2 * t +
                    l2 * r1 * Math.pow(r2, 3) * t + 3 * l2 * Math.pow(r1, 3) * r2 * t -
                    12 * a_1 * Math.pow(r1, 2) * Math.pow(r2, 2) * t + 3 * l1 * Math.pow(r1, 2) * Math.pow(r2, 2) * t +
                    2 * Math.pow(l1, 2) * r1 * Math.pow(r2, 2) * t + 2 * Math.pow(l1, 2) * Math.pow(r1, 2) * r2 * t +
                    3 * l2 * Math.pow(r1, 2) * Math.pow(r2, 2) * t + 2 * Math.pow(l2, 2) * r1 * Math.pow(r2, 2) * t +
                    2 * Math.pow(l2, 2) * Math.pow(r1, 2) * r2 * t + 8 * a_1 * l1 * Math.pow(r1, 2) * sinhExp -
                    8 * a_1 * l2 * Math.pow(r1, 2) * sinhExp - 4 * Math.pow(l1, 2) * r1 * r2 * sinhExp -
                    4 * Math.pow(l2, 2) * r1 * r2 * sinhExp - 4 * a_1 * l1 * r1 * Math.pow(r2, 2) * t -
                    8 * a_1 * l1 * Math.pow(r1, 2) * r2 * t + 4 * a_1 * l2 * r1 * Math.pow(r2, 2) * t +
                    8 * a_1 * l2 * Math.pow(r1, 2) * r2 * t - 4 * l1 * l2 * r1 * Math.pow(r2, 2) * t -
                    4 * l1 * l2 * Math.pow(r1, 2) * r2 * t + 8 * a_1 * l1 * r1 * r2 * sinhExp -
                    8 * a_1 * l2 * r1 * r2 * sinhExp + 8 * l1 * l2 * r1 * r2 * sinhExp;

            double denominator2 = 4 * (r1 + r2) * (Math.pow(l2, 2) * Math.pow(r1, 2) * t -
                    2 * Math.pow(l2, 2) * r1 * sinhExp - l1 * l2 * Math.pow(r1, 2) * t +
                    Math.pow(l2, 2) * r1 * r2 * t + 2 * l1 * l2 * r1 * sinhExp -
                    l1 * l2 * r1 * r2 * t);

            q[1][i] = numerator2 / denominator2;
        }

        for (int i = 0; i < m - 1; i++) {
            Matrix diag = new Matrix(2, 2);
            diag.set(0, 0, q[0][i]);
            diag.set(1, 1, q[1][i]);
            mmmpp[2 + i] = diag.elementMult(mmpp[1]);
        }

        Matrix diag = new Matrix(2, 2);
        double q0Sum = 0.0;
        double q1Sum = 0.0;
        for (int idx = 0; idx < q[0].length; idx++) {
            q0Sum += q[0][idx];
            q1Sum += q[1][idx];
        }
        diag.set(0, 0, 1 - q0Sum);
        diag.set(1, 1, 1 - q1Sum);
        mmmpp[2 + m - 1] = diag.elementMult(mmpp[1]);

        return mmmpp;
    }
}
