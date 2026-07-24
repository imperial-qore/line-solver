/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 *
 * Reference:
 * A. van de Liefvoort. The moment problem for continuous distributions.
 * Technical report, University of Missouri, WP-CM-1990-02, Kansas City, 1990.
 */
package jline.lib.butools.dph;

import java.util.Objects;

import jline.lib.butools.FactorialMomsFromMoms;
import jline.lib.butools.ReducedMomsFromMoms;
import jline.lib.butools.ph.MEFromMoments;
import jline.lib.butools.ph.MERepresentation;
import jline.util.matrix.Matrix;

public final class MGFromMoments {
    private MGFromMoments() {}

    /**
     * Result class for MGFromMoments containing both alpha and A.
     */
    public static final class MGRepresentation {
        public final Matrix alpha;
        public final Matrix A;

        public MGRepresentation(Matrix alpha, Matrix A) {
            this.alpha = alpha;
            this.A = A;
        }

        public Matrix component1() {
            return alpha;
        }

        public Matrix component2() {
            return A;
        }

        public Matrix getAlpha() {
            return alpha;
        }

        public Matrix getA() {
            return A;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof MGRepresentation)) return false;
            MGRepresentation that = (MGRepresentation) o;
            return Objects.equals(alpha, that.alpha) && Objects.equals(A, that.A);
        }

        @Override
        public int hashCode() {
            return Objects.hash(alpha, A);
        }

        @Override
        public String toString() {
            return "MGRepresentation(alpha=" + alpha + ", A=" + A + ")";
        }
    }

    /**
     * Creates a matrix-geometric distribution that has the same moments as given.
     *
     * @param moms The list of moments. The order of the resulting matrix-geometric distribution
     *        is determined based on the number of moments given. To obtain a matrix-geometric
     *        distribution of order M, 2*M-1 moments are required.
     * @return The MGRepresentation containing alpha (initial vector) and A (matrix parameter).
     */
    public static MGRepresentation mgFromMoments(double[] moms) {
        // Convert raw moments to factorial moments, then to reduced moments
        Matrix fmoms = FactorialMomsFromMoms.factorialMomsFromMoms(new Matrix(moms));
        Matrix rfmoms = ReducedMomsFromMoms.reducedMomsFromMoms(fmoms);

        // Prepend 1 to rfmoms
        double[] rfmomsArr = new double[moms.length + 1];
        rfmomsArr[0] = 1.0;
        for (int i = 0; i < moms.length; i++) {
            rfmomsArr[i + 1] = rfmoms.get(i);
        }

        // Compute vlist
        double[] vlist = new double[moms.length];
        double[] tmpVec = new double[moms.length + 1];
        double k = 1.0;
        tmpVec[0] = rfmomsArr[0];

        for (int i = 1; i <= moms.length; i++) {
            double sign = (i % 2 == 0) ? 1.0 : -1.0;
            tmpVec[i] = sign * rfmomsArr[i];
            k *= i;
            double sum = 0.0;
            for (int j = 0; j <= i; j++) {
                sum += tmpVec[j];
            }
            vlist[i - 1] = k * sum;
        }

        // Use MEFromMoments to get alpha and C
        MERepresentation meRep = MEFromMoments.meFromMoments(vlist);
        Matrix alpha = meRep.alpha;
        Matrix C = meRep.A;

        // A = inv(C) * inv(inv(C) + I)
        Matrix iC = C.inv();
        Matrix I = Matrix.eye(C.getNumRows());
        Matrix A = iC.mult(iC.add(I).inv());

        return new MGRepresentation(alpha, A);
    }
}
