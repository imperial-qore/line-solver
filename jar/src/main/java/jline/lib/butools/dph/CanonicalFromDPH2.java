/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dph;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Objects;

import org.apache.commons.math3.complex.Complex;

import jline.util.matrix.Matrix;

public final class CanonicalFromDPH2 {
    private CanonicalFromDPH2() {}

    /**
     * Result class for CanonicalFromDPH2 containing both beta and B.
     */
    public static final class DPH2Representation {
        public final Matrix beta;
        public final Matrix B;

        public DPH2Representation(Matrix beta, Matrix B) {
            this.beta = beta;
            this.B = B;
        }

        public Matrix component1() { return beta; }
        public Matrix component2() { return B; }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof DPH2Representation)) return false;
            DPH2Representation that = (DPH2Representation) o;
            return Objects.equals(beta, that.beta) && Objects.equals(B, that.B);
        }

        @Override
        public int hashCode() {
            return Objects.hash(beta, B);
        }
    }

    public static DPH2Representation canonicalFromDPH2(Matrix alpha, Matrix A) {
        return canonicalFromDPH2(alpha, A, 1e-14);
    }

    public static DPH2Representation canonicalFromDPH2(Matrix alpha, Matrix A, double prec) {
        if (A.getNumRows() != 2 || A.getNumCols() != 2) {
            throw new IllegalArgumentException("CanonicalFromDPH2: Dimension must be 2!");
        }

        if (!CheckMGRepresentation.checkMGRepresentation(alpha, A, prec)) {
            throw new IllegalArgumentException("CanonicalFromDPH2: Input isn't a valid MG distribution!");
        }

        // Get eigenvalues sorted by absolute value (descending)
        List<Complex> eigenvalues = A.eig();
        List<Complex> sorted = new ArrayList<Complex>(eigenvalues);
        Collections.sort(sorted, new Comparator<Complex>() {
            @Override
            public int compare(Complex a, Complex b) {
                double absA = Math.abs(a.getReal() * a.getReal() + a.getImaginary() * a.getImaginary());
                double absB = Math.abs(b.getReal() * b.getReal() + b.getImaginary() * b.getImaginary());
                return Double.compare(absB, absA);
            }
        });
        double[] lambda = new double[sorted.size()];
        for (int i = 0; i < sorted.size(); i++) {
            lambda[i] = sorted.get(i).getReal();
        }

        // e = [1; 1]
        // p1 = alpha * (e - A*e)
        int N = A.getNumRows();
        Matrix e = new Matrix(N, 1);
        for (int i = 0; i < N; i++) {
            e.set(i, 0, 1.0);
        }
        Matrix Ae = A.mult(e);
        Matrix eMinusAe = e.sub(Ae);
        double p1 = alpha.mult(eMinusAe).get(0, 0);

        Matrix beta;
        Matrix B;

        if (lambda[0] > 0 && lambda[1] > 0 && Math.abs(lambda[0] - lambda[1]) > prec) {
            // Case: Two distinct positive eigenvalues
            double d1 = (1 - lambda[0]) * (1 - p1 - lambda[1]) / (lambda[0] - lambda[1]);
            double d2 = p1 - d1;

            beta = new Matrix(1, 2);
            beta.set(0, 0, d1 * (lambda[0] - lambda[1]) / ((1 - lambda[0]) * (1 - lambda[1])));
            beta.set(0, 1, (d1 + d2) / (1 - lambda[1]));

            B = new Matrix(2, 2);
            B.set(0, 0, lambda[0]);
            B.set(0, 1, 1 - lambda[0]);
            B.set(1, 0, 0.0);
            B.set(1, 1, lambda[1]);
        } else if (lambda[0] > 0 && Math.abs(lambda[0] - lambda[1]) <= prec) {
            // Case: Two equal positive eigenvalues
            double d2 = p1;
            double d1 = (1 - lambda[0]) * (1 - d2 - lambda[0]) / lambda[0];

            beta = new Matrix(1, 2);
            beta.set(0, 0, d1 * lambda[0] / ((1 - lambda[0]) * (1 - lambda[0])));
            beta.set(0, 1, d2 / (1 - lambda[0]));

            B = new Matrix(2, 2);
            B.set(0, 0, lambda[0]);
            B.set(0, 1, 1 - lambda[0]);
            B.set(1, 0, 0.0);
            B.set(1, 1, lambda[0]);
        } else if (lambda[0] > 0) {
            // Case: One positive, one non-positive eigenvalue
            double d1 = (1 - lambda[0]) * (1 - p1 - lambda[1]) / (lambda[0] - lambda[1]);
            double d2 = p1 - d1;

            beta = new Matrix(1, 2);
            beta.set(0, 0, (d1 * lambda[0] + d2 * lambda[1]) / ((1 - lambda[0]) * (1 - lambda[1])));
            beta.set(0, 1, (d1 + d2) * (1 - lambda[0] - lambda[1]) / ((1 - lambda[0]) * (1 - lambda[1])));

            B = new Matrix(2, 2);
            B.set(0, 0, lambda[0] + lambda[1]);
            B.set(0, 1, 1 - lambda[0] - lambda[1]);
            B.set(1, 0, lambda[0] * lambda[1] / (lambda[0] + lambda[1] - 1));
            B.set(1, 1, 0.0);
        } else {
            throw new IllegalArgumentException("CanonicalFromDPH2: Cannot convert to canonical form!");
        }

        return new DPH2Representation(beta, B);
    }

    /**
     * Overload for double[] alpha.
     */
    public static DPH2Representation canonicalFromDPH2(double[] alpha, Matrix A, double prec) {
        return canonicalFromDPH2(new Matrix(alpha), A, prec);
    }

    public static DPH2Representation canonicalFromDPH2(double[] alpha, Matrix A) {
        return canonicalFromDPH2(new Matrix(alpha), A, 1e-14);
    }
}
