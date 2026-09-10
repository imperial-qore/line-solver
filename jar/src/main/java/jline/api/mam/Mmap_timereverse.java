package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Mmap_timereverse {
    private Mmap_timereverse() {}

    /**
     * Computes the time-reversed version of a Markovian Arrival Process with marked arrivals (MMAP).
     *
     * This method takes an MMAP and returns its time-reversed version. The time-reversed MMAP is computed by transposing
     * the transition matrices and adjusting them using the stationary distribution of the underlying Markov chain.
     * The resulting matrices represent the same process viewed in reverse time.
     *
     * @param mmap the MatrixCell containing the transition matrices of the original MMAP
     * @return a MatrixCell representing the time-reversed MMAP
     */
    public static MatrixCell mmap_timereverse(MatrixCell mmap) {
        int K = mmap.size();
        Matrix piq = Map_piq.map_piq(mmap.get(0), mmap.get(1));
        Matrix D = Matrix.diag(piq.toArray1D());
        Matrix iD = D.inv();
        Matrix[] DK = new Matrix[K];
        for (int k = 0; k < K; k++) {
            DK[k] = new Matrix(iD.mult(mmap.get(k).transpose()).mult(D));
        }
        return new MatrixCell(DK);
    }
}
