/**
 * @file M3PP(2,m) interleaved MMAP construction
 *
 * @since LINE 3.0
 */
package jline.api.mam.m3pp;

import java.util.List;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class M3pp2m_interleave {
    private M3pp2m_interleave() {}

    /**
     * Computes the interleaved MMAP obtained by multiple M3PP(2,m).
     */
    public static MatrixCell m3pp2m_interleave(List<MatrixCell> m3pps) {
        if (m3pps.isEmpty()) {
            throw new IllegalArgumentException("Cannot interleave empty list of M3PPs");
        }

        if (m3pps.size() == 1) {
            return m3pps.get(0);
        }

        int L = m3pps.size();

        double[][] r = new double[2][L];

        // Compute r[0][i] for i = L down to 1
        r[0][L - 1] = m3pps.get(L - 1).get(0).get(0, 1);
        for (int i = L - 2; i >= 0; i--) {
            r[0][i] = m3pps.get(i).get(0).get(0, 1);
            for (int j = i + 1; j < L; j++) {
                r[0][i] -= r[0][j];
            }
        }

        // Compute r[1][i] for i = 1 to L
        r[1][0] = m3pps.get(0).get(0).get(1, 0);
        for (int i = 1; i < L; i++) {
            r[1][i] = m3pps.get(i).get(0).get(1, 0);
            for (int j = 0; j < i; j++) {
                r[1][i] -= r[1][j];
            }
        }

        // Compute total number of class matrices M
        int M = 0;
        for (int i = 0; i < L; i++) {
            M += m3pps.get(i).size() - 2;
        }

        // Create interleaved MMAP with n = 2 + (L-1) states
        int n = 2 + (L - 1);
        MatrixCell interleavedMmap = new MatrixCell(2 + M);

        // Initialize D0 matrix
        Matrix D0 = new Matrix(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (j > i) {
                    D0.set(i, j, r[0][j - 1]);
                } else if (j < i) {
                    D0.set(i, j, r[1][j]);
                }
            }
        }
        interleavedMmap.set(0, D0);

        // Create class matrices D1c
        int classIndex = 0;
        for (int i = 0; i < L; i++) {
            MatrixCell currentM3pp = m3pps.get(i);
            int numClasses = currentM3pp.size() - 2;

            for (int j = 0; j < numClasses; j++) {
                Matrix Dic = new Matrix(n, n);

                for (int h = 0; h < n; h++) {
                    if (h <= i) {
                        Dic.set(h, h, currentM3pp.get(2 + j).get(0, 0));
                    } else {
                        Dic.set(h, h, currentM3pp.get(2 + j).get(1, 1));
                    }
                }

                interleavedMmap.set(2 + classIndex, Dic);
                classIndex++;
            }
        }

        // Compute total D1 matrix
        Matrix D1 = new Matrix(n, n);
        for (int i = 0; i < M; i++) {
            Matrix classMatrix = interleavedMmap.get(2 + i);
            for (int row = 0; row < n; row++) {
                for (int col = 0; col < n; col++) {
                    D1.set(row, col, D1.get(row, col) + classMatrix.get(row, col));
                }
            }
        }
        interleavedMmap.set(1, D1);

        // Set diagonal elements of D0 to ensure stochastic property
        for (int h = 0; h < n; h++) {
            double rowSum = 0.0;
            for (int j = 0; j < n; j++) {
                if (h != j) {
                    rowSum += D0.get(h, j);
                }
            }
            rowSum += D1.get(h, h);
            D0.set(h, h, -rowSum);
        }

        return interleavedMmap;
    }
}
