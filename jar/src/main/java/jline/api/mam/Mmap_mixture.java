/**
 * @file Marked Markovian Arrival Process mixture modeling
 *
 * Creates probabilistic mixtures of MMAP processes with specified weights.
 * Essential for modeling heterogeneous multiclass traffic patterns and aggregation.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import java.util.Map;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Mmap_mixture {
    private Mmap_mixture() {}

    /**
     * Creates a mixture of MMAPs using the given weights (alpha) and MAPs.
     */
    public static MatrixCell mmap_mixture(Matrix alpha, Map<Integer, MatrixCell> MAPs) {
        MatrixCell Dk = new MatrixCell();
        int I = MAPs.size();

        // Initialize all matrices
        for (int j = 0; j < I + 2; j++) {
            Dk.set(j, new Matrix(0, 0));
        }

        // Replace empty MAPs with exponential
        for (int i = 0; i < I; i++) {
            if (MAPs.get(i).isEmpty()) {
                MAPs.put(i, Map_exponential.map_exponential(1e6));
            }
        }

        // Build D0 as block diagonal matrix
        for (int i = 0; i < I; i++) {
            if (i == 0) {
                Dk.set(0, MAPs.get(i).get(0).copy());
            } else {
                Dk.set(0, Dk.get(0).createBlockDiagonal(MAPs.get(i).get(0)));
            }
        }

        // Build arrival matrices
        for (int i = 0; i < I; i++) {
            MatrixCell mapI = MAPs.get(i);
            Matrix D0i = mapI.get(0);
            Matrix D1i_base = mapI.get(1);
            int numStatesI = D0i.getNumRows();
            Matrix e = Matrix.ones(numStatesI, 1);

            // Build D1i for MAP i
            Matrix D1i = new Matrix(numStatesI, 0);

            for (int j = 0; j < I; j++) {
                MatrixCell mapJ = MAPs.get(j);
                Matrix pieJ = Map_pie.map_pie(mapJ);
                int numStatesJ = mapJ.get(0).getNumRows();

                // alpha(j) * D1i_base * e * pie(MAP_j)
                Matrix term = D1i_base.mult(e).mult(pieJ).scale(alpha.get(j));
                D1i = Matrix.concatColumns(D1i, term, null);
            }

            // Add D1i to the total arrival matrix D1 (index 1)
            if (i == 0) {
                Dk.set(1, D1i.copy());
            } else {
                Dk.set(1, Matrix.concatRows(Dk.get(1), D1i, null));
            }

            // Add to class-specific matrices D{2+j}
            for (int j = 0; j < I; j++) {
                if (i == j) {
                    if (i == 0) {
                        Dk.set(2 + j, D1i.copy());
                    } else {
                        Dk.set(2 + j, Matrix.concatRows(Dk.get(2 + j), D1i, null));
                    }
                } else {
                    Matrix zeroMatrix = new Matrix(D1i.getNumRows(), D1i.getNumCols());
                    zeroMatrix.zero();
                    if (i == 0) {
                        Dk.set(2 + j, zeroMatrix);
                    } else {
                        Dk.set(2 + j, Matrix.concatRows(Dk.get(2 + j), zeroMatrix, null));
                    }
                }
            }
        }

        return Mmap_normalize.mmap_normalize(Dk);
    }
}
