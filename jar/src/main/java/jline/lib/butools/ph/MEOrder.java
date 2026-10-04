/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.ph;

import jline.util.matrix.Matrix;

public final class MEOrder {
    private MEOrder() {}

    /**
     * Returns the order of the ME distribution (which is not necessarily equal to
     * the size of the representation).
     *
     * @param alpha The initial vector of the matrix-exponential distribution.
     * @param A     The matrix parameter of the matrix-exponential distribution.
     * @param kind  Determines which order is computed: "obs", "cont", "obscont", "moment".
     * @param prec  Precision for rank/determinant computations.
     * @return The order of the ME distribution.
     */
    public static int meOrder(Matrix alpha, Matrix A, String kind, double prec) {
        int N = A.getNumRows();
        if ("moment".equals(kind)) {
            double[] moms = MomentsFromME.momentsFromME(alpha, A, 2 * N - 1);
            return MEOrderFromMoments.meOrderFromMoments(moms, prec);
        } else if ("obs".equals(kind)) {
            // Observability matrix: row n = alpha * A^n, for n = 0, 1, ..., N-1
            Matrix re = new Matrix(N, N);
            Matrix alphaAn = alpha.copy();
            for (int n = 0; n < N; n++) {
                for (int j = 0; j < N; j++) {
                    re.set(n, j, alphaAn.get(0, j));
                }
                alphaAn = alphaAn.mult(A);
            }
            return re.rank();
        } else if ("cont".equals(kind)) {
            // Controllability matrix: row n = sum((A^T)^n, axis=0)
            Matrix AT = A.transpose();
            Matrix re = new Matrix(N, N);
            Matrix ATn = Matrix.eye(N);
            for (int n = 0; n < N; n++) {
                for (int j = 0; j < N; j++) {
                    double colSum = 0.0;
                    for (int i = 0; i < N; i++) {
                        colSum += ATn.get(i, j);
                    }
                    re.set(n, j, colSum);
                }
                ATn = ATn.mult(AT);
            }
            return re.rank();
        } else if ("obscont".equals(kind)) {
            int obsOrder = meOrder(alpha, A, "obs", prec);
            int contOrder = meOrder(alpha, A, "cont", prec);
            return Math.min(obsOrder, contOrder);
        } else {
            throw new IllegalArgumentException("MEOrder: Invalid 'kind' parameter '" + kind + "'!");
        }
    }

    public static int meOrder(Matrix alpha, Matrix A, String kind) {
        return meOrder(alpha, A, kind, 1e-10);
    }

    public static int meOrder(Matrix alpha, Matrix A) {
        return meOrder(alpha, A, "moment", 1e-10);
    }
}
