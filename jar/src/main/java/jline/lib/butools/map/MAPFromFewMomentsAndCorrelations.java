/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 *
 * Reference:
 * G Horvath, "Matching marginal moments and lag autocorrelations
 * with MAPs," ValueTools 2013, Torino, Italy (2013).
 */
package jline.lib.butools.map;

import jline.lang.processes.APH;
import jline.lib.butools.APHFrom3Moments;
import jline.util.matrix.Matrix;

public final class MAPFromFewMomentsAndCorrelations {
    private MAPFromFewMomentsAndCorrelations() {}

    public static Matrix[] mapFromFewMomentsAndCorrelations(double[] moms) {
        return mapFromFewMomentsAndCorrelations(moms, 0.0, null);
    }

    public static Matrix[] mapFromFewMomentsAndCorrelations(double[] moms, double corr1) {
        return mapFromFewMomentsAndCorrelations(moms, corr1, null);
    }

    /**
     * Creates a Markovian arrival process that has the given 2 or 3 marginal
     * moments and lag-1 autocorrelation.
     *
     * @return Array {D0, D1}
     */
    public static Matrix[] mapFromFewMomentsAndCorrelations(double[] moms, double corr1, Double r) {
        double m1 = moms[0];
        double c2 = moms[1] / moms[0] / moms[0] - 1.0;
        Double l3 = (moms.length > 2) ? Double.valueOf(moms[2] * moms[0] / moms[1] / moms[1] - 1.0) : null;

        Matrix alpha1;
        Matrix A1;
        Matrix alpha2;
        Matrix A2;
        double p1;
        double p2;

        if (corr1 >= 0) {
            Double rVal = r;
            double p1Calc;
            double p2Calc;

            if (rVal != null && (rVal <= 0 || rVal >= 1)) {
                throw new IllegalArgumentException("Parameter r is out of range");
            }

            if (rVal == null) {
                rVal = 2.0 * corr1 / (1.0 + corr1);
                p1Calc = (1.0 - (1.0 + corr1) / 2.0) / (1.0 + c2);
                p2Calc = (1.0 - (1.0 + corr1) / 2.0) * c2 / (1.0 + c2);
            } else {
                p1Calc = (1.0 - corr1 / rVal) / (1.0 + c2);
                p2Calc = (1.0 - corr1 / rVal) * c2 / (1.0 + c2);
            }
            p1 = p1Calc;
            p2 = p2Calc;

            double m11 = m1 * (1.0 - Math.sqrt(rVal));
            double m12 = m1 * (1.0 + c2 * Math.sqrt(rVal));

            if (l3 == null) {
                double cv21 = (Math.sqrt(c2) * (1.0 + c2) * (1.0 + Math.sqrt(rVal)))
                        / (1.0 - Math.sqrt(c2) * (-1.0 + Math.sqrt(rVal)) + c2 * Math.sqrt(rVal));
                double cv22 = -(c2 * (1.0 + c2) * (-1.0 + rVal))
                        / ((1.0 + c2 * Math.sqrt(rVal))
                                * (1.0 - Math.sqrt(c2) * (-1.0 + Math.sqrt(rVal)) + c2 * Math.sqrt(rVal)));
                double m21 = (cv21 + 1.0) * m11 * m11;
                double m22 = (cv22 + 1.0) * m12 * m12;
                Matrix[] aph1Result = aphFrom2Moments(new double[] { m11, m21 });
                alpha1 = aph1Result[0];
                A1 = aph1Result[1];
                Matrix[] aph2Result = aphFrom2Moments(new double[] { m12, m22 });
                alpha2 = aph2Result[0];
                A2 = aph2Result[1];
            } else {
                double cv21 = (c2 + Math.sqrt(rVal)) / (1.0 - Math.sqrt(rVal));
                double cv22 = c2 * (1.0 - Math.sqrt(rVal)) / (1.0 + c2 * Math.sqrt(rVal));
                double l31 = ((1.0 + c2) * l3)
                        / (c2 * (1.0 - Math.sqrt(rVal))
                                + Math.sqrt((1.0 + c2 * Math.sqrt(rVal)) * c2 * (1.0 - Math.sqrt(rVal))));
                double l32 = ((1.0 + c2) * l3)
                        / ((1.0 + c2 * Math.sqrt(rVal))
                                + Math.sqrt((1.0 + c2 * Math.sqrt(rVal)) * c2 * (1.0 - Math.sqrt(rVal))));
                double m21 = (cv21 + 1.0) * m11 * m11;
                double m22 = (cv22 + 1.0) * m12 * m12;
                double m31 = (l31 + 1.0) * m21 * m21 / m11;
                double m32 = (l32 + 1.0) * m22 * m22 / m12;
                APH aph1Result = APHFrom3Moments.APHFrom3Moments(new double[] { m11, m21, m31 });
                alpha1 = aph1Result.getInitProb();
                A1 = (Matrix) aph1Result.getParam(3).getValue();
                APH aph2Result = APHFrom3Moments.APHFrom3Moments(new double[] { m12, m22, m32 });
                alpha2 = aph2Result.getInitProb();
                A2 = (Matrix) aph2Result.getParam(3).getValue();
            }
        } else {
            Double rVal = r;

            if (c2 >= 1) {
                if (rVal != null && (rVal <= 0 || rVal >= 1.0 / c2)) {
                    throw new IllegalArgumentException("Parameter r is out of range");
                }
                if (rVal == null) {
                    rVal = -(2.0 * corr1) / (1.0 - c2 * corr1);
                    p1 = 0.5 * (1.0 + (1.0 - c2 * corr1) / 2.0);
                } else {
                    p1 = 0.5 * (1.0 - corr1 / rVal);
                }
                p2 = p1;
            } else {
                if (rVal != null && (rVal <= 0 || rVal >= 1)) {
                    throw new IllegalArgumentException("Parameter r is out of range");
                }
                if (rVal == null) {
                    rVal = -(2.0 * corr1) / (1.0 - corr1);
                    p1 = 0.5 * (1.0 + (1.0 - corr1) / 2.0);
                } else {
                    p1 = 0.5 * (1.0 - corr1 / rVal);
                }
                p2 = p1;
            }

            double m11 = m1 * (1.0 - Math.sqrt(c2 * rVal));
            double m12 = m1 * (1.0 + Math.sqrt(c2 * rVal));

            if (l3 == null) {
                double cv21 = c2 * (1.0 - rVal) / (1.0 - Math.sqrt(c2 * rVal));
                double cv22 = c2 * (1.0 - rVal) / (1.0 + Math.sqrt(c2 * rVal));
                double m21 = (cv21 + 1.0) * m11 * m11;
                double m22 = (cv22 + 1.0) * m12 * m12;
                Matrix[] aph1Result = aphFrom2Moments(new double[] { m11, m21 });
                alpha1 = aph1Result[0];
                A1 = aph1Result[1];
                Matrix[] aph2Result = aphFrom2Moments(new double[] { m12, m22 });
                alpha2 = aph2Result[0];
                A2 = aph2Result[1];
            } else {
                double cv21 = (c2 + Math.sqrt(c2 * rVal)) / (1.0 - Math.sqrt(c2 * rVal));
                double cv22 = (c2 - Math.sqrt(c2 * rVal)) / (1.0 + Math.sqrt(c2 * rVal));
                double l31 = 2.0 * l3 / (1.0 - Math.sqrt(c2 * rVal) + Math.sqrt(1.0 - c2 * rVal));
                double l32 = 2.0 * l3 / (1.0 + Math.sqrt(c2 * rVal) + Math.sqrt(1.0 - c2 * rVal));
                double m21 = (cv21 + 1.0) * m11 * m11;
                double m22 = (cv22 + 1.0) * m12 * m12;
                double m31 = (l31 + 1.0) * m21 * m21 / m11;
                double m32 = (l32 + 1.0) * m22 * m22 / m12;
                APH aph1Result = APHFrom3Moments.APHFrom3Moments(new double[] { m11, m21, m31 });
                alpha1 = aph1Result.getInitProb();
                A1 = (Matrix) aph1Result.getParam(3).getValue();
                APH aph2Result = APHFrom3Moments.APHFrom3Moments(new double[] { m12, m22, m32 });
                alpha2 = aph2Result.getInitProb();
                A2 = (Matrix) aph2Result.getParam(3).getValue();
            }
        }

        int N1 = A1.getNumRows();
        int N2 = A2.getNumRows();

        // Build D0 block diagonal
        Matrix D0 = Matrix.zeros(N1 + N2, N1 + N2);
        for (int i = 0; i < N1; i++) {
            for (int j = 0; j < N1; j++) {
                D0.set(i, j, A1.get(i, j));
            }
        }
        for (int i = 0; i < N2; i++) {
            for (int j = 0; j < N2; j++) {
                D0.set(N1 + i, N1 + j, A2.get(i, j));
            }
        }

        Matrix D1 = Matrix.zeros(N1 + N2, N1 + N2);

        double[] a1exit = new double[N1];
        for (int i = 0; i < N1; i++) {
            double sum = 0.0;
            for (int j = 0; j < N1; j++) {
                sum += A1.get(i, j);
            }
            a1exit[i] = -sum;
        }

        double[] a2exit = new double[N2];
        for (int i = 0; i < N2; i++) {
            double sum = 0.0;
            for (int j = 0; j < N2; j++) {
                sum += A2.get(i, j);
            }
            a2exit[i] = -sum;
        }

        for (int i = 0; i < N1; i++) {
            for (int j = 0; j < N1; j++) {
                D1.set(i, j, a1exit[i] * alpha1.get(0, j) * (1.0 - p1));
            }
        }
        for (int i = 0; i < N1; i++) {
            for (int j = 0; j < N2; j++) {
                D1.set(i, N1 + j, a1exit[i] * alpha2.get(0, j) * p1);
            }
        }
        for (int i = 0; i < N2; i++) {
            for (int j = 0; j < N1; j++) {
                D1.set(N1 + i, j, a2exit[i] * alpha1.get(0, j) * p2);
            }
        }
        for (int i = 0; i < N2; i++) {
            for (int j = 0; j < N2; j++) {
                D1.set(N1 + i, N1 + j, a2exit[i] * alpha2.get(0, j) * (1.0 - p2));
            }
        }

        return new Matrix[] { D0, D1 };
    }

    private static Matrix[] aphFrom2Moments(double[] moms) {
        double m1 = moms[0];
        double m2 = moms[1];
        double cv2 = m2 / m1 / m1 - 1.0;
        double lambd = 1.0 / m1;
        int N = Math.max((int) Math.ceil(1.0 / cv2), 2);
        double p = 1.0 / (cv2 + 1.0 + (cv2 - 1.0) / (N - 1));

        Matrix A = new Matrix(N, N);
        double diagVal = -lambd * p * N;
        for (int i = 0; i < N; i++) {
            A.set(i, i, diagVal);
        }
        for (int i = 0; i < N - 1; i++) {
            A.set(i, i + 1, -A.get(i, i));
        }
        A.set(N - 1, N - 1, -lambd * N);

        Matrix alpha = Matrix.zeros(1, N);
        alpha.set(0, 0, p);
        alpha.set(0, N - 1, 1.0 - p);

        return new Matrix[] { alpha, A };
    }
}
