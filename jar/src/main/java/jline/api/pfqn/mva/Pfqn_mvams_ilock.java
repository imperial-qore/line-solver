/**
 * MVA entry point for models carrying the interlocked-flow correction.
 *
 * The interlock of Franks (1999), Ch. 4, Eq. (4.7) is defined only for closed single-server
 * models, so that is the one shape accepted here; anything else is refused rather than served
 * without the correction. Models with no interlock go to {@link Pfqn_mvams}.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva;

import jline.io.Ret;
import jline.util.Utils;
import jline.util.matrix.Matrix;

public final class Pfqn_mvams_ilock {
    private Pfqn_mvams_ilock() {}

    /**
     * As {@link Pfqn_mvams#pfqn_mvams}, but carrying the interlock matrix IL through to
     * {@link Pfqn_mva_ilock}. IL is required.
     */
    public static Ret.pfqnMVA pfqn_mvams_ilock(Matrix lambda, Matrix L, Matrix N, Matrix Z, Matrix mi, Matrix S, Matrix IL) {
        if (IL == null || IL.isEmpty()) {
            throw new RuntimeException("pfqn_mvams_ilock: an interlock matrix is required; use pfqn_mvams for the standard arrival theorem");
        }
        Matrix Zlocal = (Z == null || Z.isEmpty()) ? new Matrix(1, L.getNumCols()) : Z.copy();
        if (Zlocal.getNumRows() > 1) {
            Zlocal = Zlocal.sumCols();
        }
        Matrix miLocal = (mi == null) ? null : mi.copy();
        Matrix Slocal = (S == null) ? null : S.copy();
        int M = L.getNumRows();
        boolean NhasInf = false;
        for (int i = 0; i < N.getNumRows(); i++) {
            for (int j = 0; j < N.getNumCols(); j++) {
                if (Utils.isInf(N.get(i, j))) {
                    NhasInf = true;
                }
            }
        }
        if (Slocal == null) {
            Slocal = Matrix.ones(M, 1);
        }
        if (miLocal == null) {
            miLocal = Matrix.ones(M, 1);
        }
        boolean hasMultiServer = false;
        for (int i = 0; i < Slocal.getNumRows() && !hasMultiServer; i++) {
            for (int j = 0; j < Slocal.getNumCols(); j++) {
                double num = Slocal.get(i, j);
                if (Double.isFinite(num) && num > 1.0) {
                    hasMultiServer = true;
                    break;
                }
            }
        }
        if (NhasInf || hasMultiServer) {
            throw new RuntimeException("pfqn_mvams_ilock: the interlock correction is available in exact MVA for closed single-server models only; use an AMVA method for this model.");
        }
        if (lambda != null && !lambda.isEmpty() && lambda.any()) {
            throw new RuntimeException("pfqn_mvams_ilock: the interlock correction is available in exact MVA for closed single-server models only; use an AMVA method for this model.");
        }
        Ret.pfqnMVA retMVA = Pfqn_mva_ilock.pfqn_mva_ilock(L, N, Zlocal, miLocal, IL);
        return new Ret.pfqnMVA(retMVA.X, retMVA.Q, retMVA.U, retMVA.R, retMVA.lGN);
    }
}
