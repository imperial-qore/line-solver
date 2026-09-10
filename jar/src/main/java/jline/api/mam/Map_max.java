/**
 * @file Markovian Arrival Process maximum operation
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Map_max {
    private Map_max() {}

    /**
     * Computes the MAP that represents the maximum of two independent MAPs.
     *
     * The phase space is ordered as [(i,j) pairs, B-only phases, A-only phases]:
     * in the first block both A and B are still running, in the second block A
     * has already completed and B is awaited, in the third block B has completed
     * and A is awaited. An arrival is recorded when the second of the two
     * completes, i.e. only out of the last two blocks.
     *
     * @param A the first MAP, as {D0,D1}
     * @param B the second MAP, as {D0,D1}
     * @return a MAP whose inter-arrival times are distributed as max(X,Y)
     */
    public static MatrixCell map_max(MatrixCell A, MatrixCell B) {
        Matrix A0 = A.get(0);
        Matrix B0 = B.get(0);
        int na = A0.getNumRows();
        int nb = B0.getNumRows();
        int np = na * nb;
        int n = np + nb + na;

        // completion-rate vectors of each process out of each of its phases
        Matrix a = A0.scale(-1.0).mult(Matrix.ones(na, 1));
        Matrix b = B0.scale(-1.0).mult(Matrix.ones(nb, 1));

        Matrix M0 = new Matrix(n, n);

        // pair block: both still running; state (i,j) has index i*nb+j
        Matrix pairs = A0.krons(B0);
        // A completes out of phase i while B stays in phase j -> B-only state j
        Matrix toBonly = a.kron(Matrix.eye(nb));
        // B completes out of phase j while A stays in phase i -> A-only state i
        Matrix toAonly = Matrix.eye(na).kron(b);
        for (int i = 0; i < np; i++) {
            for (int j = 0; j < np; j++) M0.set(i, j, pairs.get(i, j));
            for (int j = 0; j < nb; j++) M0.set(i, np + j, toBonly.get(i, j));
            for (int j = 0; j < na; j++) M0.set(i, np + nb + j, toAonly.get(i, j));
        }

        // B-only block: B evolves on its own until it completes
        for (int i = 0; i < nb; i++) {
            for (int j = 0; j < nb; j++) M0.set(np + i, np + j, B0.get(i, j));
        }

        // A-only block: A evolves on its own until it completes
        for (int i = 0; i < na; i++) {
            for (int j = 0; j < na; j++) M0.set(np + nb + i, np + nb + j, A0.get(i, j));
        }

        // both processes restart from their embedded equilibrium at each arrival
        Matrix pieKron = Map_pie.map_pie(A).kron(Map_pie.map_pie(B));
        Matrix pie = new Matrix(1, n);
        for (int j = 0; j < np; j++) pie.set(0, j, pieKron.get(0, j));

        // rate at which the second of the two processes completes; it is zero in
        // the pair block, where only the first of the two can still complete
        Matrix d = new Matrix(n, 1);
        for (int i = 0; i < nb; i++) d.set(np + i, 0, b.get(i, 0));
        for (int i = 0; i < na; i++) d.set(np + nb + i, 0, a.get(i, 0));

        MatrixCell result = new MatrixCell();
        result.set(0, M0);
        result.set(1, d.mult(pie));
        return result;
    }
}
