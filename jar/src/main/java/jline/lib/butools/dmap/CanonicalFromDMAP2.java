/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.apache.commons.math3.complex.Complex;

import jline.lib.butools.map.CanonicalFromMAP2;
import jline.lib.butools.mc.CRPSolve;
import jline.util.Pair;
import jline.util.matrix.Matrix;

public final class CanonicalFromDMAP2 {
    private CanonicalFromDMAP2() {}

    /**
     * Returns the canonical form of an order-2 discrete Markovian arrival process.
     */
    public static Pair<Matrix, Matrix> canonicalFromDMAP2(Matrix D0, Matrix D1, double prec) {
        if (!CheckDRAPRepresentation.checkDRAPRepresentation(D0, D1, prec)) {
            throw new IllegalArgumentException("CanonicalFromDMAP2: Input isn't a valid DRAP representation!");
        }

        if (D0.getNumRows() != 2) {
            throw new IllegalArgumentException("CanonicalFromDMAP2: Size is not 2!");
        }

        List<Complex> eigenvalues = new ArrayList<Complex>(D0.eig());
        Collections.sort(eigenvalues, new Comparator<Complex>() {
            @Override
            public int compare(Complex a, Complex b) {
                return Double.compare(b.getReal(), a.getReal());
            }
        });
        double s1 = eigenvalues.get(0).getReal();
        double s2 = eigenvalues.get(1).getReal();

        if (s2 >= 0) {
            Matrix I = Matrix.eye(2);
            Pair<Matrix, Matrix> sub = CanonicalFromMAP2.canonicalFromMAP2(D0.sub(I), D1, prec);
            return new Pair<Matrix, Matrix>(sub.getLeft().add(I), sub.getRight());
        }

        Matrix I = Matrix.eye(2);
        Matrix av = CRPSolve.drpSolve(I.sub(D0).inv().mult(D1));

        Matrix P = I.sub(D0).inv().mult(D1);
        List<Complex> gammaEig = new ArrayList<Complex>(P.eig());
        Collections.sort(gammaEig, new Comparator<Complex>() {
            @Override
            public int compare(Complex a, Complex b) {
                return Double.compare(Math.abs(b.getReal()), Math.abs(a.getReal()));
            }
        });
        double gamma = gammaEig.get(1).getReal();

        Matrix rowSum = new Matrix(2, 1);
        rowSum.set(0, 0, D0.get(0, 0) + D0.get(0, 1));
        rowSum.set(1, 0, D0.get(1, 0) + D0.get(1, 1));

        Matrix w1 = new Matrix(2, 1);
        w1.set(0, 0, (rowSum.get(0, 0) - s2) / (s1 - s2));
        w1.set(1, 0, (rowSum.get(1, 0) - s2) / (s1 - s2));

        Matrix w2 = new Matrix(2, 1);
        w2.set(0, 0, 1.0 - w1.get(0, 0));
        w2.set(1, 0, 1.0 - w1.get(1, 0));

        Matrix W = new Matrix(2, 2);
        W.set(0, 0, w1.get(0, 0));
        W.set(0, 1, w2.get(0, 0));
        W.set(1, 0, w1.get(1, 0));
        W.set(1, 1, w2.get(1, 0));

        Matrix avW = av.mult(W);
        Matrix A = avW.scale(1.0 - s1);
        double a1 = A.get(0, 0);

        Matrix G0;
        Matrix G1;

        if (gamma >= 0) {
            double innerExpr1 = -1.0 + s1 * s1 * (-2.0 + gamma) + gamma + s2 * (1.0 + a1 - a1 * gamma) + s1 * (3.0 - a1 - s2 - 2.0 * gamma + a1 * gamma);
            double innerExpr2 = -s1 * s1 * s1 * (-1.0 + gamma) + a1 * (-1.0 + s2) * s2 * (-1.0 + gamma) + s1 * s1 * (-2.0 + a1 + s2 + 2.0 * gamma - a1 * gamma) + s1 * (1.0 - a1 - s2 - gamma + a1 * gamma);
            double sqrtExpr = (-1.0 + s1 + s2) * (-1.0 + s1 + s2) * (innerExpr1 * innerExpr1 - 4.0 * (-1.0 + s1) * innerExpr2);

            double a = -(1.0 / (2.0 * (-1.0 + s1) * (-1.0 + s1 + s2) * (-1.0 + s1 + s2))) * (1.0 - 4.0 * s1 + a1 * s1 + 5.0 * s1 * s1 - a1 * s1 * s1 - 2.0 * s1 * s1 * s1 - 2.0 * s2 - a1 * s2 + 5.0 * s1 * s2 - 3.0 * s1 * s1 * s2 + s2 * s2 + a1 * s2 * s2 - s1 * s2 * s2 - gamma + 3.0 * s1 * gamma - a1 * s1 * gamma - 3.0 * s1 * s1 * gamma + a1 * s1 * s1 * gamma + s1 * s1 * s1 * gamma + s2 * gamma + a1 * s2 * gamma - 2.0 * s1 * s2 * gamma + s1 * s1 * s2 * gamma - a1 * s2 * s2 * gamma + Math.sqrt(sqrtExpr));
            double b = 1.0 + (a * (-1.0 + s1 + s2 - s1 * s2) * gamma) / ((a - 1.0) * (-s1 * s2 + a * (-1.0 + s1 + s2)));

            G0 = new Matrix(2, 2);
            G0.set(0, 0, s1 + s2);
            G0.set(0, 1, a * (1.0 - s1 - s2));
            G0.set(1, 0, s1 * s2 / (a * (s1 + s2 - 1.0)));
            G0.set(1, 1, 0.0);

            G1 = new Matrix(2, 2);
            G1.set(0, 0, (1.0 - a) * (1.0 - s1 - s2));
            G1.set(0, 1, 0.0);
            G1.set(1, 0, b * (1.0 + s1 * s2 / (a * (1.0 - s1 - s2))));
            G1.set(1, 1, (1.0 - b) * (1.0 + s1 * s2 / (a * (1.0 - s1 - s2))));
        } else {
            double innerNum = a1 * s1 - a1 * s1 * s1 + s2 - a1 * s2 - 3.0 * s1 * s2 + 2.0 * s1 * s1 * s2 - s2 * s2 + a1 * s2 * s2 + s1 * s2 * s2 + s1 * gamma - a1 * s1 * gamma - 2.0 * s1 * s1 * gamma + a1 * s1 * s1 * gamma + s1 * s1 * s1 * gamma + a1 * s2 * gamma - a1 * s2 * s2 * gamma;
            double sqrtArg = -4.0 * (-1.0 + s1) * s1 * s2 * (-1.0 + s1 + s2) * (a1 * (s1 - s2) * (-1.0 + gamma) + (-1.0 + s1) * (s2 + (-1.0 + s1) * gamma));
            double innerExpr = a1 * (-s1 + s1 * s1 + s2 - s2 * s2) * (-1.0 + gamma) + (-1.0 + s1) * ((-1.0 + 2.0 * s1) * s2 + s2 * s2 + (-1.0 + s1) * s1 * gamma);
            double sqrtFull = sqrtArg + innerExpr * innerExpr;
            double denom = 2.0 * (-1.0 + s1 + s2) * (a1 * (s1 - s2) * (-1.0 + gamma) + (-1.0 + s1) * (s2 + (-1.0 + s1) * gamma));

            double a = (innerNum + Math.sqrt(sqrtFull)) / denom;
            double b = -((a * (1.0 - s1) * (1.0 - s2) * gamma) / ((a - 1.0) * (-a + a * s1 + a * s2 - s1 * s2)));

            G0 = new Matrix(2, 2);
            G0.set(0, 0, s1 + s2);
            G0.set(0, 1, a * (1.0 - s1 - s2));
            G0.set(1, 0, s1 * s2 / (a * (s1 + s2 - 1.0)));
            G0.set(1, 1, 0.0);

            G1 = new Matrix(2, 2);
            G1.set(0, 0, 0.0);
            G1.set(0, 1, (1.0 - a) * (1.0 - s1 - s2));
            G1.set(1, 0, b * (1.0 - s1 * s2 / (a * (s1 + s2 - 1.0))));
            G1.set(1, 1, (1.0 - b) * (1.0 - s1 * s2 / (a * (s1 + s2 - 1.0))));
        }

        return new Pair<Matrix, Matrix>(G0, G1);
    }

    public static Pair<Matrix, Matrix> canonicalFromDMAP2(Matrix D0, Matrix D1) {
        return canonicalFromDMAP2(D0, D1, 1e-14);
    }
}
