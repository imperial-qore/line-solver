/**
 * @file Logistic expansion with the eps->0 bias correction (BLE)
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import org.apache.commons.math3.util.FastMath;

import jline.GlobalConstants;
import jline.io.Ret;
import jline.util.matrix.Matrix;

public final class Pfqn_ble {
    private Pfqn_ble() {}

    /** The eps->0 bias in nats per Laplaced direction. */
    public static final double BLE_CORRECTION = 1.0 - FastMath.log(2 * FastMath.PI) / 2;

    /**
     * Logistic expansion with an additive (1-log(2*pi)/2) bias correction per
     * Gaussian direction the branch actually takes.
     *
     * Cas17 Theorem 4.1 holds for eps >= eps_N > 0; the K(1+eps*N) self-looping
     * populations are what make the integrand concentrate. Evaluated at eps->0, as
     * pfqn_le does, the curvature at the saddle tends to 1 rather than growing with
     * N, so Laplace's method has no asymptotic regime there and carries an O(1)
     * relative bias of e/sqrt(2*pi) PER LAPLACED DIRECTION. The count is the exponent
     * on sqrt(2*pi) in the branch taken: M-1 with Z=0, where the radial integral is
     * exact as Gamma(N+M), and M with Z>0, where the radius is Laplaced too. Measured
     * over the 1562 models of the Cas17 dataset (Zenodo 546873, sec5.3.1, sigma=100),
     * the Z>0 deficit is M to within 0.01 units. The published expansion is NOT in
     * error and the correction is empirical; see _kb/03-api-layer.md.
     */
    public static Ret.pfqnNc pfqn_ble(Matrix L, Matrix N, Matrix Z) {
        Ret.pfqnNc base = Pfqn_le.pfqn_le(L, N, Z);
        if (L.isEmpty() || N.isEmpty() || N.elementSum() == 0.0
                || L.elementSum() < GlobalConstants.CoarseTol) {
            // Degenerate branch: the delay term is exact, there is no Laplace step to correct.
            return base;
        }
        boolean noDelay = Z.isEmpty() || Z.elementSum() < GlobalConstants.Zero;
        int nGauss = noDelay ? L.getNumRows() - 1 : L.getNumRows();
        double lGn = base.lG + nGauss * BLE_CORRECTION;
        return new Ret.pfqnNc(Double.valueOf(FastMath.exp(lGn)), Double.valueOf(lGn));
    }

    public static Ret.pfqnNc pfqn_ble(Matrix L, Matrix N) {
        return pfqn_ble(L, N, new Matrix(1, L.getNumCols()));
    }
}
