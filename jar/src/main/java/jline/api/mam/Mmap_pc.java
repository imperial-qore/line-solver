package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Mmap_pc {
    private Mmap_pc() {}

    /**
     * Computes the proportion of counts (PC) for each type in a Markovian Arrival Process with marked arrivals (MMAP).
     *
     * <p>This method calculates the proportion of arrivals attributed to each type in the MMAP. It does so by computing the
     * stationary distribution of the underlying Markov chain and using it to weight the arrivals from each type.
     *
     * @param MMAP the MMAP from which to compute the proportions
     * @return a matrix where each element represents the proportion of counts for a type
     */
    public static Matrix mmap_pc(MatrixCell MMAP) {
        int m = MMAP.size() - 2;

        Matrix neg_D0 = MMAP.get(0).copy();
        neg_D0.scaleEq(-1.0);
        Matrix PC = new Matrix(m, 1, m);
        for (int i = 0; i < m; i++) {
            double v = Map_pie.map_pie(MMAP.get(0), MMAP.get(1)).mult(neg_D0.inv().mult(MMAP.get(2 + i))).elementSum();
            PC.set(i, v);
        }
        return PC;
    }
}
