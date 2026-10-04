/**
 * @file Quasi-Birth-Death process R-matrix logarithmic reduction
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;

public final class Qbd_R_logred {
    private Qbd_R_logred() {}

    /**
     * Compute R matrix using logarithmic reduction method for QBD processes.
     *
     * @param B Backward matrix
     * @param L Local matrix
     * @param F Forward matrix
     * @param iterMax Maximum number of iterations
     * @return R matrix
     */
    public static Matrix qbd_R_logred(Matrix B, Matrix L, Matrix F, int iterMax) {
        int r = L.getNumRows();
        Matrix LInv = L.inv();

        Matrix iLF = LInv.mult(F).scale(-1.0);
        Matrix iLB = LInv.mult(B).scale(-1.0);

        Matrix T = iLF.copy();
        Matrix S = iLB.copy();

        for (int iter = 1; iter <= iterMax; iter++) {
            Matrix D = iLF.mult(iLB).add(iLB.mult(iLF));
            Matrix eyeMinusD = Matrix.eye(r).sub(D);
            Matrix eyeMinusDInv = eyeMinusD.inv();

            iLF = eyeMinusDInv.mult(iLF.mult(iLF));
            iLB = eyeMinusDInv.mult(iLB.mult(iLB));

            S = S.add(T.mult(iLB));
            T.multEq(iLF);

            Matrix ones = Matrix.ones(r, 1);
            Matrix convergenceTest = ones.sub(S.mult(ones));
            if (convergenceTest.norm() <= 1e-12) {
                break;
            }
        }

        Matrix U = L.add(F.mult(S));
        return F.scale(-1.0).mult(U.inv());
    }

    public static Matrix qbd_R_logred(Matrix B, Matrix L, Matrix F) {
        return qbd_R_logred(B, L, F, 100000);
    }
}
