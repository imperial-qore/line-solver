/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dph;

import jline.lib.butools.fitting.PHFromTrace;
import jline.lib.butools.ph.PH3Representation;
import jline.util.matrix.Matrix;
import org.apache.commons.math3.complex.Complex;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

public final class CanonicalFromDPH3 {
    private CanonicalFromDPH3() {}

    public static class DPH3Representation {
        public final Matrix beta;
        public final Matrix B;

        public DPH3Representation(Matrix beta, Matrix B) {
            this.beta = beta;
            this.B = B;
        }
    }

    public static DPH3Representation canonicalFromDPH3(Matrix alpha, Matrix A) {
        return canonicalFromDPH3(alpha, A, 1e-14);
    }

    public static DPH3Representation canonicalFromDPH3(Matrix alpha, Matrix A, double prec) {
        if (A.getNumRows() != 3 || A.getNumCols() != 3) {
            throw new IllegalArgumentException("CanonicalFromDPH3: Dimension must be 3!");
        }
        if (!checkMGRepresentation(alpha, A, prec)) {
            throw new IllegalArgumentException("CanonicalFromDPH3: Input isn't a valid MG distribution!");
        }

        List<Complex> evList = A.eig();
        evList.sort(new Comparator<Complex>() {
            public int compare(Complex a, Complex b) {
                double ma = Math.abs(a.getReal() * a.getReal() + a.getImaginary() * a.getImaginary());
                double mb = Math.abs(b.getReal() * b.getReal() + b.getImaginary() * b.getImaginary());
                return Double.compare(mb, ma);
            }
        });
        final List<Complex> lambda = evList;

        double a0 = -lambda.get(0).getReal() * lambda.get(1).getReal() * lambda.get(2).getReal();
        double a1 = lambda.get(0).getReal() * lambda.get(1).getReal()
                + lambda.get(0).getReal() * lambda.get(2).getReal()
                + lambda.get(1).getReal() * lambda.get(2).getReal();
        double a2 = -lambda.get(0).getReal() - lambda.get(1).getReal() - lambda.get(2).getReal();

        int N = A.getNumRows();
        Matrix e = new Matrix(N, 1);
        for (int i = 0; i < N; i++) e.set(i, 0, 1.0);

        Matrix alphaOut;
        Matrix Aout;

        double L0 = lambda.get(0).getReal(), L1 = lambda.get(1).getReal(), L2 = lambda.get(2).getReal();

        if (L0 > 0 && L1 >= 0 && L2 >= 0) {
            Matrix I = Matrix.eye(3);
            PH3Representation ph3Rep = canonicalFromPH3Adapter(alpha, A.sub(I), prec);
            alphaOut = ph3Rep.alpha;
            Aout = ph3Rep.A.add(I);
        } else if (L0 > 0 && L1 >= 0 && L2 < 0) {
            double x1 = L0;
            double x2 = L1 + L2;
            double x3 = L1 * L2 / (L1 + L2 - 1);

            Aout = new Matrix(3, 3);
            Aout.set(0, 0, x1);
            Aout.set(0, 1, 1 - x1);
            Aout.set(0, 2, 0.0);
            Aout.set(1, 0, 0.0);
            Aout.set(1, 1, x2);
            Aout.set(1, 2, 1 - x2);
            Aout.set(2, 0, 0.0);
            Aout.set(2, 1, x3);
            Aout.set(2, 2, 0.0);

            Matrix Ae = A.mult(e);
            Matrix eMinusAe = e.sub(Ae);
            Matrix b3 = eMinusAe.scale(1.0 / (1 - x3));
            Matrix b2 = A.mult(b3).scale(1.0 / (1 - x2));
            Matrix b1 = e.sub(b2).sub(b3);

            Matrix B = new Matrix(3, 3);
            for (int i = 0; i < 3; i++) {
                B.set(i, 0, b1.get(i, 0));
                B.set(i, 1, b2.get(i, 0));
                B.set(i, 2, b3.get(i, 0));
            }
            alphaOut = alpha.mult(B);
        } else if (L0 > 0 && L1 < 0 && L2 >= 0) {
            double x1 = -a2;
            double x2 = (a0 - a1 * a2) / (a2 * (1 + a2));
            double x3 = a0 * (1 + a2) / (a0 - a2 - a1 * a2 - a2 * a2);

            Aout = new Matrix(3, 3);
            Aout.set(0, 0, x1); Aout.set(0, 1, 1 - x1); Aout.set(0, 2, 0.0);
            Aout.set(1, 0, x2); Aout.set(1, 1, 0.0); Aout.set(1, 2, 1 - x2);
            Aout.set(2, 0, 0.0); Aout.set(2, 1, x3); Aout.set(2, 2, 0.0);

            Matrix Ae = A.mult(e);
            Matrix eMinusAe = e.sub(Ae);
            Matrix b3 = eMinusAe.scale(1.0 / (1 - x3));
            Matrix b2 = A.mult(b3).scale(1.0 / (1 - x2));
            Matrix b1 = e.sub(b2).sub(b3);

            double alphab1 = alpha.mult(b1).get(0, 0);
            if (alphab1 >= 0) {
                Matrix B = new Matrix(3, 3);
                for (int i = 0; i < 3; i++) {
                    B.set(i, 0, b1.get(i, 0));
                    B.set(i, 1, b2.get(i, 0));
                    B.set(i, 2, b3.get(i, 0));
                }
                alphaOut = alpha.mult(B);
            } else {
                boolean foundValid = false;
                double x33 = 0.0;
                double validX1 = 0.0, validX2 = 0.0, validX3 = 0.0;
                Matrix validB = new Matrix(3, 3);
                double validA1 = 0.0;
                while (x33 <= 1 && !foundValid) {
                    double[] sortEigs = new double[]{L0, L1, L2};
                    FirstInitElemResult result = firstInitElem(x33, sortEigs, alpha, A);
                    validA1 = result.a1;
                    validX1 = result.m1;
                    validX2 = result.m2;
                    validX3 = result.m3;
                    validB = result.B;
                    if (validA1 >= 0 && validX1 >= 0 && validX2 >= 0 && validX3 >= 0 && validX3 + x33 < 1) {
                        foundValid = true;
                        break;
                    }
                    x33 += 0.01;
                }
                if (foundValid) {
                    Aout.set(0, 0, validX1); Aout.set(0, 1, 1 - validX1); Aout.set(0, 2, 0.0);
                    Aout.set(1, 0, validX2); Aout.set(1, 1, 0.0); Aout.set(1, 2, 1 - validX2);
                    Aout.set(2, 0, 0.0); Aout.set(2, 1, validX3); Aout.set(2, 2, x33);
                    alphaOut = alpha.mult(validB);
                } else {
                    double px1 = L2;
                    double px2 = L0 + L1;
                    double px3 = L0 * L1 / (L0 + L1 - 1);
                    Aout.set(0, 0, px1); Aout.set(0, 1, 0.0); Aout.set(0, 2, 0.0);
                    Aout.set(1, 0, 0.0); Aout.set(1, 1, px2); Aout.set(1, 2, 1 - px2);
                    Aout.set(2, 0, 0.0); Aout.set(2, 1, px3); Aout.set(2, 2, 0.0);

                    Matrix Ae2 = A.mult(e);
                    Matrix eMinusAe2 = e.sub(Ae2);
                    double p1 = alpha.mult(eMinusAe2).get(0, 0);
                    double p2 = alpha.mult(A).mult(eMinusAe2).get(0, 0);

                    double l1v = L0;
                    double l2v = L1;
                    double l3v = L2;

                    double d1 = (1 - l1v) * ((1 - l2v) * (1 - l3v) + (-1 + l2v + l3v) * p1 - p2) / ((l1v - l2v) * (l1v - l3v));
                    double d2 = (l2v - 1) * ((1 - l1v) * (1 - l3v) + (-1 + l1v + l3v) * p1 - p2) / ((l1v - l2v) * (l2v - l3v));
                    double d3 = (l3v - 1) * ((1 - l1v) * (1 - l2v) + (-1 + l1v + l2v) * p1 - p2) / ((l2v - l3v) * (l3v - l1v));

                    alphaOut = new Matrix(1, 3);
                    alphaOut.set(0, 0, d3 / (1 - l3v));
                    alphaOut.set(0, 1, (d1 * l1v + d2 * l2v) / ((1 - l1v) * (1 - l2v)));
                    alphaOut.set(0, 2, (d1 + d2) * (1 - l1v - l2v) / ((1 - l1v) * (1 - l2v)));

                    if (alphaOut.elementMin() < 0 || Aout.elementMin() < 0) {
                        throw new IllegalArgumentException("CanonicalFromDPH3: Unhandled PNP case!");
                    }
                }
            }
        } else if (L0 > 0 && L1 < 0 && L2 < 0) {
            double absLambda2 = Math.abs(L1 * L1 + lambda.get(1).getImaginary() * lambda.get(1).getImaginary());
            boolean allReal = true;
            for (Complex ev : lambda) {
                if (Math.abs(ev.getImaginary()) >= prec) { allReal = false; break; }
            }
            if (allReal || absLambda2 <= 2 * L0 * (-L1)) {
                double x1 = -a2;
                double x2 = -a1 / (1 + a2);
                double x3 = -a0 / (1 + a1 + a2);

                Aout = new Matrix(3, 3);
                Aout.set(0, 0, x1); Aout.set(0, 1, 1 - x1); Aout.set(0, 2, 0.0);
                Aout.set(1, 0, x2); Aout.set(1, 1, 0.0); Aout.set(1, 2, 1 - x2);
                Aout.set(2, 0, x3); Aout.set(2, 1, 0.0); Aout.set(2, 2, 0.0);

                Matrix Ae = A.mult(e);
                Matrix eMinusAe = e.sub(Ae);
                Matrix b3 = eMinusAe.scale(1.0 / (1 - x3));
                Matrix b2 = A.mult(b3).scale(1.0 / (1 - x2));
                Matrix b1 = e.sub(b2).sub(b3);

                Matrix B = new Matrix(3, 3);
                for (int i = 0; i < 3; i++) {
                    B.set(i, 0, b1.get(i, 0));
                    B.set(i, 1, b2.get(i, 0));
                    B.set(i, 2, b3.get(i, 0));
                }
                alphaOut = alpha.mult(B);
            } else {
                Matrix I = Matrix.eye(3);
                PH3Representation ph3Rep = canonicalFromPH3Adapter(alpha, A.sub(I), prec);
                alphaOut = ph3Rep.alpha;
                Aout = ph3Rep.A.add(I);
            }
        } else {
            throw new IllegalArgumentException("CanonicalFromDPH3: Unhandled eigenvalue configuration!");
        }
        return new DPH3Representation(alphaOut, Aout);
    }

    public static DPH3Representation canonicalFromDPH3(double[] alpha, Matrix A) {
        return canonicalFromDPH3(new Matrix(alpha), A, 1e-14);
    }

    public static DPH3Representation canonicalFromDPH3(double[] alpha, Matrix A, double prec) {
        return canonicalFromDPH3(new Matrix(alpha), A, prec);
    }

    private static class FirstInitElemResult {
        final double a1;
        final double m1;
        final double m2;
        final double m3;
        final Matrix B;
        FirstInitElemResult(double a1, double m1, double m2, double m3, Matrix B) {
            this.a1 = a1; this.m1 = m1; this.m2 = m2; this.m3 = m3; this.B = B;
        }
    }

    private static FirstInitElemResult firstInitElem(double m33, double[] sortEigs, Matrix alpha, Matrix A) {
        double l1 = sortEigs[0], l2 = sortEigs[1], l3 = sortEigs[2];
        double m1 = -m33 + l1 + l2 + l3;

        double m2Num = -((l2 - l3) * (l1 * l1 - l1 * l2 - l1 * l3 + l2 * l3) *
                (m33 * m33 * m33 - 2 * m33 * m33 * l1 + m33 * l1 * l1 - 2 * m33 * m33 * l2 + 3 * m33 * l1 * l2 - l1 * l1 * l2 +
                        m33 * l2 * l2 - l1 * l2 * l2 - 2 * m33 * m33 * l3 + 3 * m33 * l1 * l3 - l1 * l1 * l3 + 3 * m33 * l2 * l3 -
                        2 * l1 * l2 * l3 - l2 * l2 * l3 + m33 * l3 * l3 - l1 * l3 * l3 - l2 * l3 * l3));
        double m2Den = 2 * m33 * l1 * l1 * l2 + 2 * m33 * m33 * l1 * l1 * l2 - l1 * l1 * l1 * l2 - 3 * m33 * l1 * l1 * l1 * l2 +
                l1 * l1 * l1 * l1 * l2 - 2 * m33 * l1 * l2 * l2 - 2 * m33 * m33 * l1 * l2 * l2 + l1 * l1 * l1 * l2 * l2 +
                l1 * l2 * l2 * l2 + 3 * m33 * l1 * l2 * l2 * l2 - l1 * l1 * l2 * l2 * l2 - l1 * l2 * l2 * l2 * l2 -
                2 * m33 * l1 * l1 * l3 - 2 * m33 * m33 * l1 * l1 * l3 + l1 * l1 * l1 * l3 + 3 * m33 * l1 * l1 * l1 * l3 -
                l1 * l1 * l1 * l1 * l3 + 2 * m33 * l2 * l2 * l3 + 2 * m33 * m33 * l2 * l2 * l3 - l2 * l2 * l2 * l3 -
                3 * m33 * l2 * l2 * l2 * l3 + l2 * l2 * l2 * l2 * l3 + 2 * m33 * l1 * l3 * l3 + 2 * m33 * m33 * l1 * l3 * l3 -
                l1 * l1 * l1 * l3 * l3 - 2 * m33 * l2 * l3 * l3 - 2 * m33 * m33 * l2 * l3 * l3 + l2 * l2 * l2 * l3 * l3 -
                l1 * l3 * l3 * l3 - 3 * m33 * l1 * l3 * l3 * l3 + l1 * l1 * l3 * l3 * l3 + l2 * l3 * l3 * l3 +
                3 * m33 * l2 * l3 * l3 * l3 - l2 * l2 * l3 * l3 * l3 + l1 * l3 * l3 * l3 * l3 - l2 * l3 * l3 * l3 * l3;
        double m2 = (Math.abs(m2Den) > 1e-14) ? m2Num / m2Den : 0.0;

        double m3Num = (l2 - l3) * (l1 * l1 - l1 * l2 - l1 * l3 + l2 * l3) *
                (m33 * m33 * m33 + m33 * m33 * m33 * m33 - m33 * m33 * l1 - 2 * m33 * m33 * m33 * l1 + m33 * m33 * l1 * l1 -
                        m33 * m33 * l2 - 2 * m33 * m33 * m33 * l2 + m33 * l1 * l2 + 3 * m33 * m33 * l1 * l2 - m33 * l1 * l1 * l2 +
                        m33 * m33 * l2 * l2 - m33 * l1 * l2 * l2 - m33 * m33 * l3 - 2 * m33 * m33 * m33 * l3 + m33 * l1 * l3 +
                        3 * m33 * m33 * l1 * l3 - m33 * l1 * l1 * l3 + m33 * l2 * l3 + 3 * m33 * m33 * l2 * l3 - l1 * l2 * l3 -
                        4 * m33 * l1 * l2 * l3 + l1 * l1 * l2 * l3 - m33 * l2 * l2 * l3 + l1 * l2 * l2 * l3 + m33 * m33 * l3 * l3 -
                        m33 * l1 * l3 * l3 - m33 * l2 * l3 * l3 + l1 * l2 * l3 * l3);
        double m3Den = -2 * m33 * l1 * l1 * l2 - 2 * m33 * m33 * l1 * l1 * l2 - m33 * m33 * m33 * l1 * l1 * l2 +
                l1 * l1 * l1 * l2 + 3 * m33 * l1 * l1 * l1 * l2 + 2 * m33 * m33 * l1 * l1 * l1 * l2 - l1 * l1 * l1 * l1 * l2 -
                m33 * l1 * l1 * l1 * l1 * l2 + 2 * m33 * l1 * l2 * l2 + 2 * m33 * m33 * l1 * l2 * l2 + m33 * m33 * m33 * l1 * l2 * l2 -
                l1 * l1 * l1 * l2 * l2 - 2 * m33 * l1 * l1 * l1 * l2 * l2 + l1 * l1 * l1 * l1 * l2 * l2 - l1 * l2 * l2 * l2 -
                3 * m33 * l1 * l2 * l2 * l2 - 2 * m33 * m33 * l1 * l2 * l2 * l2 + l1 * l1 * l2 * l2 * l2 +
                2 * m33 * l1 * l1 * l2 * l2 * l2 + l1 * l2 * l2 * l2 * l2 + m33 * l1 * l2 * l2 * l2 * l2 - l1 * l1 * l2 * l2 * l2 * l2 +
                2 * m33 * l1 * l1 * l3 + 2 * m33 * m33 * l1 * l1 * l3 + m33 * m33 * m33 * l1 * l1 * l3 - l1 * l1 * l1 * l3 -
                3 * m33 * l1 * l1 * l1 * l3 - 2 * m33 * m33 * l1 * l1 * l1 * l3 + l1 * l1 * l1 * l1 * l3 + m33 * l1 * l1 * l1 * l1 * l3 -
                2 * m33 * l2 * l2 * l3 - 2 * m33 * m33 * l2 * l2 * l3 - m33 * m33 * m33 * l2 * l2 * l3 + l2 * l2 * l2 * l3 +
                3 * m33 * l2 * l2 * l2 * l3 + 2 * m33 * m33 * l2 * l2 * l2 * l3 - l2 * l2 * l2 * l2 * l3 - m33 * l2 * l2 * l2 * l2 * l3 -
                2 * m33 * l1 * l3 * l3 - 2 * m33 * m33 * l1 * l3 * l3 - m33 * m33 * m33 * l1 * l3 * l3 + l1 * l1 * l1 * l3 * l3 +
                2 * m33 * l1 * l1 * l1 * l3 * l3 - l1 * l1 * l1 * l1 * l3 * l3 + 2 * m33 * l2 * l3 * l3 + 2 * m33 * m33 * l2 * l3 * l3 +
                m33 * m33 * m33 * l2 * l3 * l3 - l2 * l2 * l2 * l3 * l3 - 2 * m33 * l2 * l2 * l2 * l3 * l3 + l2 * l2 * l2 * l2 * l3 * l3 +
                l1 * l3 * l3 * l3 + 3 * m33 * l1 * l3 * l3 * l3 + 2 * m33 * m33 * l1 * l3 * l3 * l3 - l1 * l1 * l3 * l3 * l3 -
                2 * m33 * l1 * l1 * l3 * l3 * l3 - l2 * l3 * l3 * l3 - 3 * m33 * l2 * l3 * l3 * l3 - 2 * m33 * m33 * l2 * l3 * l3 * l3 +
                l2 * l2 * l3 * l3 * l3 + 2 * m33 * l2 * l2 * l3 * l3 * l3 - l1 * l3 * l3 * l3 * l3 - m33 * l1 * l3 * l3 * l3 * l3 +
                l1 * l1 * l3 * l3 * l3 * l3 + l2 * l3 * l3 * l3 * l3 + m33 * l2 * l3 * l3 * l3 * l3 - l2 * l2 * l3 * l3 * l3 * l3;
        double m3 = (Math.abs(m3Den) > 1e-14) ? m3Num / m3Den : 0.0;

        int N = A.getNumRows();
        Matrix e = new Matrix(N, 1);
        for (int i = 0; i < N; i++) e.set(i, 0, 1.0);
        Matrix I = Matrix.eye(N);

        Matrix eMinusAe = I.sub(A).mult(e);
        Matrix b3 = eMinusAe.scale(1.0 / (1 - m3 - m33));
        Matrix m33I = I.scale(m33);
        Matrix b2 = Matrix.negative(A.sub(m33I)).mult(b3).scale(1.0 / (1 - m2));
        Matrix b1 = A.mult(b2).sub(b3.scale(m3)).scale(1.0 / (1 - m1));

        Matrix B = new Matrix(3, 3);
        for (int i = 0; i < 3; i++) {
            B.set(i, 0, b1.get(i, 0));
            B.set(i, 1, b2.get(i, 0));
            B.set(i, 2, b3.get(i, 0));
        }
        double a1 = alpha.mult(b1).get(0, 0);
        return new FirstInitElemResult(a1, m1, m2, m3, B);
    }

    private static boolean checkMGRepresentation(Matrix alpha, Matrix A, double prec) {
        // Mirror utility from butools - delegated implementation
        return jline.lib.butools.dph.CheckMGRepresentation.checkMGRepresentation(alpha, A, prec);
    }

    private static PH3Representation canonicalFromPH3Adapter(Matrix alpha, Matrix A, double prec) {
        return jline.lib.butools.ph.CanonicalFromPH3.canonicalFromPH3(alpha, A, prec);
    }
}
