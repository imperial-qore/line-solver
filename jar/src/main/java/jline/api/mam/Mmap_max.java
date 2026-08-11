/**
 * @file Marked Markovian Arrival Process maximum (fork-join synchronization)
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import java.util.List;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * MMAP of the maximum of two independent marked arrival flows, as used to
 * synchronize sibling branches at a fork-join Join node.
 *
 * <p>Port of matlab/lib/m3a/m3a/mmap/mmap_max.m. The construction augments the
 * Kronecker state space of the two flows with a synchronization queue of length
 * k, so the resulting MMAP has order na*nb*(1+2k): one "balanced" level in which
 * neither flow leads, plus 2k levels tracking which flow is ahead and by how
 * many arrivals (up to k).</p>
 */
public final class Mmap_max {
    private Mmap_max() {}

    /**
     * Computes the MMAP representing the maximum of two MMAPs, synchronizing
     * through a queue of length k.
     *
     * @param MMAPa first flow, {D0, D1, D1_c1, ...}
     * @param MMAPb second flow, {D0, D1, D1_c1, ...}
     * @param k     length of the synchronization queue (k >= 1)
     * @return the synchronized MMAP, with the same number of marks as MMAPa
     */
    public static MatrixCell mmap_max(MatrixCell MMAPa, MatrixCell MMAPb, int k) {
        if (k < 1) {
            throw new IllegalArgumentException("mmap_max: synchronization queue length k must be >= 1");
        }
        Matrix D0a = MMAPa.get(0);
        Matrix D0b = MMAPb.get(0);
        Matrix D1a = MMAPa.get(1);
        Matrix D1b = MMAPb.get(1);

        int na = D0a.getNumRows();
        int nb = D0b.getNumRows();

        Matrix Ia = Matrix.eye(na);
        Matrix Ib = Matrix.eye(nb);

        Matrix A0B0 = D0a.krons(D0b);
        Matrix A1IB = D1a.kron(Ib);
        Matrix IAB1 = Ia.kron(D1b);
        Matrix IAB0 = Ia.kron(D0b);
        Matrix A0IB = D0a.kron(Ib);

        int iRows = A0B0.getNumRows();
        int iCols = A0B0.getNumCols();
        int nBlocks = 1 + 2 * k;
        int dim = iRows * nBlocks;

        Matrix M0 = new Matrix(dim, dim);
        Matrix M1 = new Matrix(dim, dim);

        // MATLAB: M0(1:iRows, 1:iCols*3) = [A0B0 A1IB IAB1];
        setBlock(M0, 0, 0, iRows, iCols, A0B0);
        setBlock(M0, 0, 1, iRows, iCols, A1IB);
        setBlock(M0, 0, 2, iRows, iCols, IAB1);

        // MATLAB: for i=2:1+(k-1)*2 -> diagonal A0B0 on blocks 1..2k-2
        for (int i = 2; i <= 1 + (k - 1) * 2; i++) {
            setBlock(M0, i - 1, i - 1, iRows, iCols, A0B0);
        }

        // MATLAB: the two terminal blocks 2k-1 and 2k
        setBlock(M0, 2 * k - 1, 2 * k - 1, iRows, iCols, IAB0);
        setBlock(M0, 2 * k, 2 * k, iRows, iCols, A0IB);

        // MATLAB: for i=2:k -> lead-increasing transitions without arrival marks
        for (int i = 2; i <= k; i++) {
            int r = 1 + 2 * (i - 2);
            int c = 3 + 2 * (i - 2);
            setBlock(M0, r, c, iRows, iCols, A1IB);
            setBlock(M0, r + 1, c + 1, iRows, iCols, IAB1);
        }

        // MATLAB: M1(iRows+1:iRows*3, 1:iCols) = [IAB1; A1IB];
        setBlock(M1, 1, 0, iRows, iCols, IAB1);
        setBlock(M1, 2, 0, iRows, iCols, A1IB);

        for (int i = 2; i <= k; i++) {
            int r = 1 + 2 * (i - 1);
            int c = 1 + 2 * (i - 2);
            setBlock(M1, r, c, iRows, iCols, IAB1);
            setBlock(M1, r + 1, c + 1, iRows, iCols, A1IB);
        }

        MatrixCell MMAP = new MatrixCell(MMAPa.size());
        MMAP.set(0, M0);
        MMAP.set(1, M1);

        // Per-mark matrices: same block structure as M1, driven by the
        // per-class arrival matrices of each flow (MATLAB: for cls = 3:end)
        for (int cls = 2; cls < MMAPa.size(); cls++) {
            Matrix Mc = new Matrix(dim, dim);
            Matrix D1ca = MMAPa.get(cls);
            Matrix D1cb = MMAPb.get(cls);
            Matrix IAB1c = Ia.kron(D1cb);
            Matrix A1IBc = D1ca.kron(Ib);

            setBlock(Mc, 1, 0, iRows, iCols, IAB1c);
            setBlock(Mc, 2, 0, iRows, iCols, A1IBc);

            for (int i = 2; i <= k; i++) {
                int r = 1 + 2 * (i - 1);
                int c = 1 + 2 * (i - 2);
                setBlock(Mc, r, c, iRows, iCols, IAB1c);
                setBlock(Mc, r + 1, c + 1, iRows, iCols, A1IBc);
            }
            MMAP.set(cls, Mc);
        }

        return MMAP;
    }

    /**
     * Iteratively synchronizes a list of MMAPs through a queue of length k.
     */
    public static MatrixCell mmap_max_multiple(List<MatrixCell> mmaps, int k) {
        if (mmaps.isEmpty()) {
            throw new IllegalArgumentException("Must provide at least one MMAP");
        }
        MatrixCell result = mmaps.get(0);
        for (int i = 1; i < mmaps.size(); i++) {
            result = mmap_max(result, mmaps.get(i), k);
        }
        return result;
    }

    /**
     * Copies src into the (blockRow, blockCol) block of dest.
     */
    private static void setBlock(Matrix dest, int blockRow, int blockCol,
                                 int iRows, int iCols, Matrix src) {
        int r0 = blockRow * iRows;
        int c0 = blockCol * iCols;
        for (int i = 0; i < iRows; i++) {
            for (int j = 0; j < iCols; j++) {
                double v = src.get(i, j);
                if (v != 0.0) {
                    dest.set(r0 + i, c0 + j, v);
                }
            }
        }
    }
}
