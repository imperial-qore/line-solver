/**
 * @file Quasi-Birth-Death process BMAP/BMAP/1 queue analysis
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import java.util.ArrayList;
import java.util.List;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Qbd_bmapbmap1 {
    private Qbd_bmapbmap1() {}

    /**
     * Result containing QBD matrices for BMAP/BMAP/1 queue.
     */
    public static final class QbdBmapResult {
        public final Matrix A0;
        public final Matrix A_1;
        public final List<Matrix> A1;
        public final Matrix B0;
        public final List<Matrix> B1;

        public QbdBmapResult(Matrix A0, Matrix A_1, List<Matrix> A1, Matrix B0, List<Matrix> B1) {
            this.A0 = A0;
            this.A_1 = A_1;
            this.A1 = A1;
            this.B0 = B0;
            this.B1 = B1;
        }
    }

    /**
     * Set up QBD matrices for BMAP/BMAP/1 queue analysis.
     */
    public static QbdBmapResult qbd_bmapbmap1(MatrixCell MAPa, Matrix pbatcha, MatrixCell MAPs) {
        int na = MAPa.get(0).getNumRows();
        int ns = MAPs.get(0).getNumRows();
        int maxbatch = pbatcha.length();

        Matrix IA = Matrix.eye(na);
        Matrix IS = Matrix.eye(ns);

        // Upward transition blocks for each batch size
        List<Matrix> A1 = new ArrayList<Matrix>();
        for (int b = 0; b < maxbatch; b++) {
            Matrix A1_b = MAPa.get(1).scale(pbatcha.get(b)).kron(IS);
            A1.add(A1_b);
        }

        Matrix A0 = MAPa.get(0).kron(IS).add(IA.kron(MAPs.get(0)));
        Matrix A_1 = IA.kron(MAPs.get(1));
        Matrix B0 = MAPa.get(0).kron(IS).add(Matrix.eye(na * ns));

        List<Matrix> B1 = new ArrayList<Matrix>();
        for (int b = 0; b < maxbatch; b++) {
            Matrix B1_b = MAPa.get(1).scale(pbatcha.get(b)).kron(IS);
            B1.add(B1_b);
        }

        return new QbdBmapResult(A0, A_1, A1, B0, B1);
    }
}
