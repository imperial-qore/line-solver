/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 *
 * Reference:
 * M. Telek and A. Heindl, "Moment bounds for acyclic discrete and continuous
 * phase-type distributions of second order," in Proc. of UK Performance
 * Evaluation Workshop, UKPEW, 2002
 */
package jline.lib.butools.ph;

import java.util.Objects;

import jline.lib.butools.APH2ndMomentLowerBound;
import jline.lib.butools.APH3rdMomentLowerBound;
import jline.lib.butools.APH3rdMomentUpperBound;
import jline.util.matrix.Matrix;

public final class PH2From3Moments {
    private PH2From3Moments() {}

    /**
     * Result class for PH2From3Moments.
     */
    public static final class PH2Representation {
        public final Matrix alpha;
        public final Matrix A;

        public PH2Representation(Matrix alpha, Matrix A) {
            this.alpha = alpha;
            this.A = A;
        }

        public Matrix component1() { return alpha; }
        public Matrix component2() { return A; }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof PH2Representation)) return false;
            PH2Representation that = (PH2Representation) o;
            return Objects.equals(alpha, that.alpha) && Objects.equals(A, that.A);
        }

        @Override
        public int hashCode() {
            return Objects.hash(alpha, A);
        }

        @Override
        public String toString() {
            return "PH2Representation(alpha=" + alpha + ", A=" + A + ")";
        }
    }

    public static PH2Representation ph2From3Moments(double[] moms) {
        return ph2From3Moments(moms, 1e-14);
    }

    /**
     * Returns a PH(2) which has the same 3 moments as given.
     */
    public static PH2Representation ph2From3Moments(double[] moms, double prec) {
        double m1 = moms[0];
        double m2 = moms[1];
        double m3 = moms[2];

        // Check moment bounds
        double m2l = APH2ndMomentLowerBound.APH2ndMomentLowerBound(m1, 2);
        double m3l = APH3rdMomentLowerBound.APH3rdMomentLowerBound(m1, m2, 2);
        double m3u = APH3rdMomentUpperBound.APH3rdMomentUpperBound(m1, m2, 2);

        if (m2 < m2l) {
            throw new IllegalArgumentException("The given second moment is not feasible!");
        }
        if (m3 < m3l) {
            throw new IllegalArgumentException("The given third moment is not feasible (too small)!");
        }
        if (m3 > m3u) {
            throw new IllegalArgumentException("The given third moment is not feasible (too large)!");
        }

        // Check if we have an exponential distribution
        if (Math.abs(m2 / m1 / m1 - 2.0) < prec) {
            Matrix alpha = new Matrix(1, 1);
            alpha.set(0, 0, 1.0);
            Matrix A = new Matrix(1, 1);
            A.set(0, 0, -1.0 / m1);
            return new PH2Representation(alpha, A);
        }

        // Calculate parameters
        double b = 3.0 * m1 * m2 - m3;
        double c = 3.0 * m2 * m2 - 2.0 * m1 * m3;
        double e = -2.0 * m1 * m1 + m2;
        double a = b * b + 6.0 * c * e;
        if (a < 0) {
            a = 0.0;
        }
        a = Math.sqrt(a);

        double lambda1;
        double lambda2;
        double p;

        if (c > 0) {
            lambda1 = (b - a) / c;
            lambda2 = (b + a) / c;
            p = (-b - 6.0 * m1 * e + a) / (b + a);
        } else if (c < 0) {
            lambda1 = (b + a) / c;
            lambda2 = (b - a) / c;
            p = (b + 6.0 * m1 * e + a) / (-b + a);
        } else {
            lambda1 = 0.0;
            lambda2 = 1.0 / m1;
            p = 0.0;
        }

        Matrix alpha = new Matrix(1, 2);
        alpha.set(0, 0, p);
        alpha.set(0, 1, 1.0 - p);

        Matrix A = new Matrix(2, 2);
        A.set(0, 0, -lambda1);
        A.set(0, 1, lambda1);
        A.set(1, 0, 0.0);
        A.set(1, 1, -lambda2);

        return new PH2Representation(alpha, A);
    }
}
