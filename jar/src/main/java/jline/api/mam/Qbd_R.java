package jline.api.mam;

import jline.util.matrix.Matrix;

/**
 * QBD R-matrix computation algorithms.
 *
 * Provides methods for computing the R-matrix in Quasi-Birth-Death (QBD) processes using
 * matrix analytic methods. The R-matrix is fundamental to analyzing QBD processes and
 * represents the conditional probability that the level increases before decreasing.
 *
 * @since LINE 3.0
 */
public final class Qbd_R {
    private Qbd_R() {}

    /**
     * Compute R matrix using successive substitutions method for QBD processes.
     *
     * @param B Backward matrix
     * @param L Local matrix
     * @param F Forward matrix
     * @param iterMax Maximum number of iterations (default: 100000)
     * @return R matrix
     */
    public static Matrix qbd_R(Matrix B, Matrix L, Matrix F, int iterMax) {
        Matrix LInv = L.inv();
        Matrix Fil = F.mult(LInv);
        Matrix BiL = B.mult(LInv);

        Matrix R = Fil.scale(-1.0);
        Matrix Rprime = Fil.scale(-1.0).sub(R.mult(R).mult(BiL));

        for (int iter = 1; iter <= iterMax; iter++) {
            R = Rprime.copy();
            Rprime = Fil.scale(-1.0).sub(R.mult(R).mult(BiL));

            Matrix diff = R.sub(Rprime);
            if (diff.norm() <= 1e-12) {
                break;
            }
        }

        return Rprime;
    }

    public static Matrix qbd_R(Matrix B, Matrix L, Matrix F) {
        return qbd_R(B, L, F, 100000);
    }
}
