/**
 * @file Knessl-Tier expansion with the Stirling-remainder correction (BKT)
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.FastMath;

import jline.GlobalConstants;
import jline.io.Ret;
import jline.util.matrix.Matrix;

public final class Pfqn_bkt {
    private Pfqn_bkt() {}

    /**
     * The Stirling remainder s(N) = log(N!) - (N log N - N + log(2 pi N)/2), exactly, for
     * N >= 1. s(1) is the constant pfqn_ble adds per station direction, and
     * s(N) = 1/(12 N) + O(N^-2).
     */
    public static double stirlingRemainder(double n) {
        return Gamma.logGamma(n + 1.0) - (n + 0.5) * FastMath.log(n) + n
                - FastMath.log(2 * FastMath.PI) / 2;
    }

    /**
     * Knessl-Tier expansion corrected for the Stirling remainder that steepest descent
     * drops in each class direction.
     *
     * pfqn_kt extracts N from the generating function of G by steepest descent. On the
     * demand-free integral the exact coefficient is [u^N] exp(Z u) = Z^N/N!, whereas the
     * expansion returns N log Z - (N log N - N + log(2 pi N)/2), Stirling's approximation
     * of log(N!) in place of log(N!). So KT lies ABOVE the exact value by s(N) per Laplaced
     * class direction, and BKT subtracts sum_r s(N_r). The remainder is evaluated exactly
     * from logGamma: truncating it at 1/(12 N) loses an order of magnitude (on the 1562
     * models of Cas17 sec5.3.1 the median |error| is 0.083 nats for KT, 1.9e-4 for the
     * truncation and 1.4e-5 for the exact remainder). With a think time BKT is the SAME
     * estimator as pfqn_ble, to the accuracy of the two saddle-point solvers.
     *
     * Only the classes pfqn_kt actually Laplaces are corrected: a class with no jobs is
     * dropped by its recursion, and a self-looping class (one nonzero demand and no think
     * time) has its coefficient extracted exactly, so neither carries a remainder. The
     * predicate below is pfqn_kt's own. See _kb/03-api-layer.md.
     */
    public static Ret.pfqnNc pfqn_bkt(Matrix L, Matrix N, Matrix Z) {
        Matrix Zc = (Z == null || Z.isEmpty()) ? new Matrix(1, L.getNumCols()) : Z;
        Ret.pfqnNc base = Pfqn_kt.pfqn_kt(L, N, Zc);
        int M = L.getNumRows();
        int R = L.getNumCols();
        int nKeep = 0;
        for (int r = 0; r < R; r++) {
            if (N.get(0, r) > GlobalConstants.Zero) nKeep++;
        }
        double corr = 0.0;
        for (int r = 0; r < R; r++) {
            double n = N.get(0, r);
            if (n <= GlobalConstants.Zero) continue;          // dropped by pfqn_kt's recursion
            if (nKeep > 1) {                                   // pfqn_kt folds self-loops only with Rorig > 1
                int nnz = 0;
                for (int k = 0; k < M; k++) {
                    if (L.get(k, r) > GlobalConstants.Zero) nnz++;
                }
                if (nnz == 1 && Zc.get(0, r) == 0.0) continue; // extracted exactly, no remainder
            }
            corr += stirlingRemainder(n);
        }
        double lGn = base.lG - corr;
        return new Ret.pfqnNc(Double.valueOf(FastMath.exp(lGn)), Double.valueOf(lGn));
    }

    public static Ret.pfqnNc pfqn_bkt(Matrix L, Matrix N) {
        return pfqn_bkt(L, N, new Matrix(1, L.getNumCols()));
    }
}
