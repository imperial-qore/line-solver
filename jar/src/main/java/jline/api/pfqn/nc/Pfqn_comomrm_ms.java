/**
 * @file CoMoM normalizing constant method for multiserver repairman models
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import jline.GlobalConstants;
import jline.api.pfqn.ld.Pfqn_mu_ms;
import jline.io.Ret;
import jline.util.matrix.Matrix;
import org.apache.commons.math3.util.FastMath;

public final class Pfqn_comomrm_ms {
    private Pfqn_comomrm_ms() {}

    /**
     * Four-argument form. The fourth positional argument is the replication
     * factor m, matching MATLAB pfqn_comomrm_ms(L,N,Z,m); it read as S here,
     * so the same call transcribed from MATLAB silently solved a different
     * model. S defaults to 1, as in MATLAB.
     */
    public static Ret.pfqnComomrmMs pfqn_comomrm_ms(Matrix L, Matrix N, Matrix Z, int m) {
        return pfqn_comomrm_ms(L, N, Z, m, 1);
    }

    /**
     * Compute the normalizing constant of a multiserver repairman model using CoMoM.
     */
    public static Ret.pfqnComomrmMs pfqn_comomrm_ms(Matrix L, Matrix N, Matrix Z, int m, int S) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        if (M != 1) {
            throw new RuntimeException("pfqn_comomrm_ms: The solver accepts at most a single queueing station.");
        }

        double atol = GlobalConstants.FineTol;
        Matrix lambda = new Matrix(1, R);
        lambda.fill(0.0);
        Ret.pfqnNcSanitize ret = Pfqn_nc_sanitize.pfqn_nc_sanitize(lambda, L, N, Z, atol);
        Matrix L_new = ret.L;
        Matrix N_new = ret.N;
        Matrix Z_new = ret.Z;
        double lG0 = ret.lGremaind;
        // see _kb/03-api-layer.md for rationale
        R = L_new.getNumCols();

        int Nt = (int) N_new.elementSum();

        Matrix mu;
        if (m > 1) {
            mu = Pfqn_mu_ms.pfqn_mu_ms(Nt, m, S);
        } else {
            mu = new Matrix(1, Nt);
            for (int i = 0; i < Nt; i++) {
                mu.set(i, FastMath.min((double) S, (double) (i + 1)));
            }
        }

        Matrix h = new Matrix(Nt + 1, 1);
        h.set(Nt, 1.0);
        Matrix scale = new Matrix(Nt, 1);
        int nt = 0;

        for (int r = 0; r < R; r++) {
            Matrix Tr = Matrix.eye(Nt + 1).scale(Z_new.get(r));
            for (int i = 0; i < Nt; i++) {
                Tr.set(i, i + 1, L_new.get(r) * (Nt - i) / mu.get(Nt - i - 1));
            }
            for (int nr = 1; nr <= (int) N_new.get(r); nr++) {
                Matrix hT = Tr.copy().scale(1.0 / (double) nr);
                h = hT.mult(h);
                scale.set(nt, FastMath.abs(h.sort().elementSum()));
                h.absEq();
                h.scaleEq(1.0 / scale.get(nt));
                nt++;
            }
        }

        double lG = lG0 + scale.log().elementSum();
        double G = FastMath.exp(lG);
        Matrix prob = h.reverse().scale(1.0 / G);
        prob.divideEq(prob.elementSum());

        return new Ret.pfqnComomrmMs(G, lG, prob);
    }
}
