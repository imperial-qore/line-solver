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

import org.apache.commons.math3.complex.Complex;

import jline.lib.butools.reptrans.SimilarityMatrix;
import jline.util.matrix.Matrix;

public final class AcyclicDPHFromMG {
    private AcyclicDPHFromMG() {}

    /**
     * Transforms a matrix-geometric representation to an acyclic DPH representation
     * of the same size, if possible.
     *
     * @param alpha Initial vector of the distribution
     * @param A Matrix parameter of the distribution
     * @param prec Vector and matrix entries smaller than the precision are considered to be zeros.
     * @return The MGRepresentation containing beta and B
     */
    public static MGFromMoments.MGRepresentation acyclicDPHFromMG(Matrix alpha, Matrix A, double prec) {
        if (!CheckMGRepresentation.checkMGRepresentation(alpha, A, prec)) {
            throw new IllegalArgumentException("AcyclicDPHFromMG: Input isn't a valid MG distribution!");
        }

        int N = A.getNumRows();

        // Get eigenvalues
        List<Complex> eigenvalues = A.eig();

        // Check for complex eigenvalues
        for (Complex ev : eigenvalues) {
            if (Math.abs(ev.getImaginary()) > prec) {
                throw new IllegalArgumentException("AcyclicDPHFromMG: The input matrix has complex eigenvalue!");
            }
        }

        // Sort eigenvalues by real part (descending)
        List<Double> lambda = new ArrayList<Double>();
        for (Complex ev : eigenvalues) {
            lambda.add(ev.getReal());
        }
        Collections.sort(lambda, new Comparator<Double>() {
            @Override
            public int compare(Double a, Double b) {
                return Double.compare(b, a);
            }
        });

        // Create target acyclic matrix: mx = diag(lambda) + diag(1-lambda(1:end-1), 1)
        Matrix mx = Matrix.zeros(N, N);
        for (int i = 0; i < N; i++) {
            mx.set(i, i, lambda.get(i));
        }
        for (int i = 0; i < N - 1; i++) {
            mx.set(i, i + 1, 1.0 - lambda.get(i));
        }

        // Find similarity transformation matrix
        Matrix T = SimilarityMatrix.similarityMatrix(A, mx);

        // beta = alpha * T
        Matrix beta = alpha.mult(T);
        Matrix B = mx;

        // Verify the result is a valid DPH representation
        if (!CheckDPHRepresentation.checkDPHRepresentation(beta, B, prec)) {
            throw new IllegalArgumentException("AcyclicDPHFromMG: No acyclic representation found!");
        }

        return new MGFromMoments.MGRepresentation(beta, B);
    }

    public static MGFromMoments.MGRepresentation acyclicDPHFromMG(Matrix alpha, Matrix A) {
        return acyclicDPHFromMG(alpha, A, 1e-14);
    }

    public static MGFromMoments.MGRepresentation acyclicDPHFromMG(double[] alpha, Matrix A, double prec) {
        return acyclicDPHFromMG(new Matrix(alpha), A, prec);
    }

    public static MGFromMoments.MGRepresentation acyclicDPHFromMG(double[] alpha, Matrix A) {
        return acyclicDPHFromMG(new Matrix(alpha), A, 1e-14);
    }
}
