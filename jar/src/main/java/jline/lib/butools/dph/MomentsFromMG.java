/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dph;

import jline.lib.butools.MomsFromFactorialMoms;
import jline.util.matrix.Matrix;

public final class MomentsFromMG {
    private MomentsFromMG() {}

    /**
     * Returns the first K moments of a matrix-geometric distribution.
     *
     * @param alpha The initial vector of the matrix-geometric distribution.
     * @param A The matrix parameter of the matrix-geometric distribution.
     * @param K Number of moments to compute. If K=0, 2*M-1 moments are computed.
     * @return The vector of moments.
     */
    public static double[] momentsFromMG(Matrix alpha, Matrix A, int K) {
        int N = A.getNumRows();

        Matrix rowAlpha = alpha;
        if (alpha.getNumRows() == N && alpha.getNumCols() == 1) {
            rowAlpha = alpha.transpose();
        }

        int numMoments = (K == 0) ? 2 * N - 1 : K;

        // iA = inv(I - A)
        Matrix I = Matrix.eye(N);
        Matrix iA = I.sub(A).inv();

        // Compute factorial moments: fmoms(i) = i! * alpha * iA^i * A^(i-1)
        double[] fmoms = new double[numMoments];
        double factorial = 1.0;
        Matrix iAPower = iA.copy();
        Matrix APower = Matrix.eye(N);

        for (int i = 1; i <= numMoments; i++) {
            factorial *= i;
            Matrix term = rowAlpha.mult(iAPower).mult(APower);
            fmoms[i - 1] = factorial * term.elementSum();

            iAPower = iAPower.mult(iA);
            APower = APower.mult(A);
        }

        // Convert factorial moments to raw moments
        Matrix fmomsMatrix = new Matrix(fmoms);
        Matrix moms = MomsFromFactorialMoms.MomsFromFactorialMoms(fmomsMatrix);

        double[] result = new double[numMoments];
        for (int i = 0; i < numMoments; i++) {
            result[i] = moms.get(i);
        }
        return result;
    }

    public static double[] momentsFromMG(Matrix alpha, Matrix A) {
        return momentsFromMG(alpha, A, 0);
    }

    public static double[] momentsFromMG(double[] alpha, Matrix A, int K) {
        return momentsFromMG(new Matrix(alpha), A, K);
    }

    public static double[] momentsFromMG(double[] alpha, Matrix A) {
        return momentsFromMG(new Matrix(alpha), A, 0);
    }
}
