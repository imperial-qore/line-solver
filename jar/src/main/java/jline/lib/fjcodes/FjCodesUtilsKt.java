/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 *
 * Java translation of FJ_codes Kotlin top-level utility functions.
 * Based on FJ_codes MATLAB toolkit:
 * Z. Qiu, J.F. Perez, P. Harrison, "Beyond the Mean in Fork-Join Queues:
 * Efficient Approximation for Response-Time Tails", IFIP Performance 2015.
 */
package jline.lib.fjcodes;

import jline.util.matrix.Matrix;

/**
 * FJ_codes utility helpers translated from FJUtils.kt's top-level functions.
 * The class name preserves the Kotlin-generated facade name so existing
 * call sites (e.g. {@code FjCodesUtilsKt.setSubMatrix(...)}) compile unchanged.
 */
public final class FjCodesUtilsKt {

    private FjCodesUtilsKt() {}

    /**
     * Compute Kronecker sum of two matrices: A xor B = A x I + I x B.
     */
    public static Matrix kronsum(Matrix A, Matrix B) {
        Matrix IA = Matrix.eye(A.getNumRows());
        Matrix IB = Matrix.eye(B.getNumRows());
        return A.kron(IB).add(1.0, IA.kron(B));
    }

    /**
     * Find the 1-based index of the row in {@code matrix} that matches
     * {@code row} elementwise (within 1e-10), or -1 if no row matches.
     */
    public static int vectmatch(double[] row, Matrix matrix) {
        int m = matrix.getNumRows();
        int n = matrix.getNumCols();
        for (int outer = 0; outer < m; outer++) {
            boolean matches = true;
            for (int inner = 0; inner < n; inner++) {
                if (Math.abs(matrix.get(outer, inner) - row[inner]) > 1e-10) {
                    matches = false;
                    break;
                }
            }
            if (matches) {
                return outer + 1;
            }
        }
        return -1;
    }

    /**
     * Build combinatorial index patterns: enumerate all ways to distribute
     * {@code cr} items into {@code m} bins. Each row of the returned matrix
     * is one distribution pattern.
     */
    public static Matrix build_index(int m, int cr) {
        int totalDim = binomialCoeff(cr + m - 1, cr);
        Matrix indexes = new Matrix(totalDim, m);
        indexes.set(0, 0, cr);

        for (int row = 1; row < totalDim; row++) {
            int k = -1;
            for (int col = 0; col < m; col++) {
                if (indexes.get(row - 1, col) > 0.0) {
                    k = col;
                    break;
                }
            }

            if (k >= 0 && k < m - 1) {
                for (int col = 0; col < m; col++) {
                    indexes.set(row, col, indexes.get(row - 1, col));
                }
                double currentVal = indexes.get(row, k + 1);
                indexes.set(row, k + 1, currentVal + 1.0);
                double kVal = indexes.get(row, k);
                indexes.set(row, 0, kVal - 1.0);
                for (int col = 1; col <= k; col++) {
                    indexes.set(row, col, 0.0);
                }
            }
        }
        return indexes;
    }

    /**
     * Extract row {@code row} of {@code matrix} as a {@code double[]}.
     */
    public static double[] getRowAsArray(Matrix matrix, int row) {
        Matrix rowMatrix = matrix.getRow(row);
        int len = rowMatrix.length();
        double[] result = new double[len];
        for (int i = 0; i < len; i++) {
            result[i] = rowMatrix.get(i);
        }
        return result;
    }

    /**
     * Copy {@code source} into {@code target} starting at the given
     * (0-based) row/column position.
     */
    public static void setSubMatrix(Matrix target, int startRow, int startCol, Matrix source) {
        for (int i = 0; i < source.getNumRows(); i++) {
            for (int j = 0; j < source.getNumCols(); j++) {
                target.set(startRow + i, startCol + j, source.get(i, j));
            }
        }
    }

    /**
     * Binomial coefficient C(n, k) = n! / (k! * (n - k)!).
     */
    private static int binomialCoeff(int n, int k) {
        if (k > n) return 0;
        if (k == 0 || k == n) return 1;
        long result = 1L;
        int kMin = Math.min(k, n - k);
        for (int i = 0; i < kMin; i++) {
            result = result * (long) (n - i) / (long) (i + 1);
        }
        return (int) result;
    }
}
