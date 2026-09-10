/**
 * @file Inter-departure MAP of a two-resource closed subnetwork
 *
 * @since LINE 3.0
 */
package jline.api.fes;

import java.util.List;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * Builds the MAP (T0,T1) of the inter-departure times of a closed subnetwork made of one
 * MAP station and one MAP flow-equivalent server.
 *
 * Implements the block bidiagonal construction of Casale, Mi, Cherkasova and Smirni,
 * IEEE Trans. Soft. Eng. 37(5), 2011, Section 5.2.2. Level k is the population of the
 * flow-equivalent server, so the station holds n-k jobs and both processes may be load
 * dependent. Marked transitions are the completions of the station, which are the
 * departures fed to the rest of the model.
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
public final class Fes_map_interdeparture {
    private Fes_map_interdeparture() {}

    /**
     * Builds the inter-departure MAP of a subnetwork holding n jobs.
     *
     * @param maps per-level service processes of the station, index j-1 holding j jobs
     * @param fes  per-level processes of the flow-equivalent server, index k-1 holding k jobs
     * @param n    number of jobs circulating in the subnetwork
     * @return the pair (T0,T1) of the inter-departure MAP
     */
    public static MatrixCell fes_map_interdeparture(List<MatrixCell> maps, List<MatrixCell> fes, int n) {
        if (n < 1) {
            throw new IllegalArgumentException("The subnetwork population n must be at least 1.");
        }
        List<MatrixCell> mapsLev = Fes_map_levels.fes_map_levels(maps, n);
        List<MatrixCell> fesLev = Fes_map_levels.fes_map_levels(fes, n);

        int ms = mapsLev.get(0).get(0).getNumRows();
        int mf = fesLev.get(0).get(0).getNumRows();
        int blk = ms * mf;
        int dim = (n + 1) * blk;

        Matrix Ims = Matrix.eye(ms);
        Matrix Imf = Matrix.eye(mf);
        Matrix T0 = new Matrix(dim, dim);
        Matrix T1 = new Matrix(dim, dim);

        for (int k = 0; k <= n; k++) {
            int off = k * blk;
            int j = n - k;
            Matrix diag = null;
            if (j > 0) {
                diag = mapsLev.get(j - 1).get(0).kron(Imf);
            }
            if (k > 0) {
                Matrix fesDiag = Ims.kron(fesLev.get(k - 1).get(0));
                diag = (diag == null) ? fesDiag : diag.add(1.0, fesDiag);
            }
            addBlock(T0, diag, off, off);

            if (k > 0) {
                addBlock(T0, Ims.kron(fesLev.get(k - 1).get(1)), off, off - blk);
            }
            if (j > 0) {
                addBlock(T1, mapsLev.get(j - 1).get(1).kron(Imf), off, off + blk);
            }
        }

        return new MatrixCell(T0, T1);
    }

    private static void addBlock(Matrix target, Matrix block, int rowOff, int colOff) {
        for (int r = 0; r < block.getNumRows(); r++) {
            for (int c = 0; c < block.getNumCols(); c++) {
                double v = block.get(r, c);
                if (v != 0) {
                    target.set(rowOff + r, colOff + c, target.get(rowOff + r, colOff + c) + v);
                }
            }
        }
    }
}
