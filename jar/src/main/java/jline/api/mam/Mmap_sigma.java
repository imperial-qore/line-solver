package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Mmap_sigma {
    private Mmap_sigma() {}

    /**
     * Computes one-step class transition probabilities for a Marked Markovian Arrival Process (MMAP).
     *
     * The function computes the class transition probabilities p_{i,j} = P(C_k = j | C_{k-1} = i),
     * which represent the probability of transitioning from class i to class j in one step.
     *
     * @param MMAP the MMAP represented as a MatrixCell where MMAP[0] = D0, MMAP[1] = aggregate D1,
     *             and MMAP[2+i] = D_{i+1}
     * @return matrix of class transition probabilities
     */
    public static Matrix mmap_sigma(MatrixCell MMAP) {
        int C = MMAP.size() - 2;  // Number of classes
        Matrix sigma = Matrix.zeros(C, C);

        // Get the stationary probability vector of the underlying MAP
        Matrix alpha = Map_pie.map_pie(MMAP.get(0), MMAP.get(1));

        // Compute (-D0)^{-1} once for efficiency
        Matrix invNegD0 = MMAP.get(0).scale(-1.0).inv();

        for (int i = 0; i < C; i++) {
            // Compute alpha * (-D0)^{-1} * D_{i+1}
            Matrix start = alpha.mult(invNegD0).mult(MMAP.get(2 + i));

            for (int j = 0; j < C; j++) {
                // Compute start * (-D0)^{-1} * D_{j+1} * 1
                Matrix result = start.mult(invNegD0).mult(MMAP.get(2 + j));
                sigma.set(i, j, result.elementSum());
            }
        }

        return sigma;
    }

    /**
     * Computes one-step class transition probabilities for an MMAP given as Matrix[].
     */
    public static Matrix mmap_sigma(Matrix[] mmap) {
        // Convert to MatrixCell and delegate to main implementation
        MatrixCell mmapCell = new MatrixCell(mmap.length);
        for (int i = 0; i < mmap.length; i++) {
            mmapCell.set(i, mmap[i]);
        }
        return mmap_sigma(mmapCell);
    }
}
