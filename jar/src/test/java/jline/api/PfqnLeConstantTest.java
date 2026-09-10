package jline.api;

import jline.api.pfqn.nc.Pfqn_ca;
import jline.api.pfqn.nc.Pfqn_le;
import jline.api.pfqn.nc.Pfqn_ble;
import jline.io.Ret;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * The logistic expansion evaluated at eps->0 under-reports log G by exactly
 * (1 - log(2 pi)/2) nats PER LAPLACED DIRECTION: M-1 with Z = 0, where the radial
 * integral is exact as Gamma(N+M), and M with Z > 0, where the radius is Laplaced
 * too. The deficit is a constant: independent of the number of classes, of the
 * demands and of rank(L), and it does NOT vanish as the population grows, so it
 * cannot be dismissed as the method's O(1/N) error. pfqn_le keeps Cas17 eq. (34)
 * as published, pfqn_ble adds the correction.
 */
public class PfqnLeConstantTest {

    private static final double C = 1.0 - Math.log(2 * Math.PI) / 2;

    private static Matrix demands(int M, int R, long seed) {
        Matrix L = new Matrix(M, R);
        long x = seed;
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                x = (x * 6364136223846793005L + 1442695040888963407L);
                double u = ((x >>> 11) / (double) (1L << 53));
                L.set(i, r, 0.5 + u);
            }
        }
        return L;
    }

    private static Matrix pop(int R, double n) {
        Matrix N = new Matrix(1, R);
        for (int r = 0; r < R; r++) N.set(0, r, n);
        return N;
    }

    @Test
    public void logisticExpansionDeficitIsTheConstant() {
        double worst = 0.0;
        for (int M = 2; M <= 8; M++) {
            for (int R = 1; R <= 2; R++) {
                Matrix L = demands(M, R, 12345L + M);
                Matrix N = pop(R, 200.0);
                Matrix Z = new Matrix(1, R);
                double lLe = Pfqn_le.pfqn_le(L, N, Z).lG;
                double lCa = Pfqn_ca.pfqn_ca(L, N, Z).lG;
                worst = Math.max(worst, Math.abs((lCa - lLe) - (M - 1) * C));
            }
        }
        // The raw gap peaks at (M-1)*0.0811 = 0.568 nats at M=8 and grows linearly in M.
        // Net of the constant, what is left is the method's genuine O(1/N) error, 0.055
        // nats on this sweep, which unlike the constant does shrink with the population.
        assertTrue(worst < 0.1, "worst |(lG_ca - lG_le) - (M-1)*c| was " + worst + " nats");
    }

    @Test
    public void bleMatchesConvolutionAtLargePopulation() {
        double worst = 0.0;
        for (int M = 2; M <= 8; M++) {
            for (int R = 1; R <= 2; R++) {
                Matrix L = demands(M, R, 12345L + M);
                Matrix N = pop(R, 200.0);
                Matrix Z = new Matrix(1, R);
                Ret.pfqnNc ble = Pfqn_ble.pfqn_ble(L, N, Z);
                double lCa = Pfqn_ca.pfqn_ca(L, N, Z).lG;
                worst = Math.max(worst, Math.abs(ble.lG - lCa));
            }
        }
        assertTrue(worst < 0.1, "worst |lG_ble - lG_ca| was " + worst + " nats");
    }

    /**
     * With a think time the expansion Laplaces the radius as well, so the deficit is M
     * units, not M-1. Measured over the 1562 models of the Cas17 dataset (Zenodo
     * 546873, sec5.3.1, sigma = 100) the count is M to within 0.01 units, and using
     * M-1 there leaves a residual of exactly one unit on every model.
     */
    @Test
    public void theDelayBranchLosesOneUnitPerStationIncludingTheRadius() {
        double worst = 0.0;
        for (int M = 2; M <= 6; M++) {
            Matrix L = demands(M, 2, 999L + M);
            Matrix N = pop(2, 40.0);
            Matrix Z = new Matrix(1, 2);
            Z.set(0, 0, 100.0);
            Z.set(0, 1, 100.0);
            double lLe = Pfqn_le.pfqn_le(L, N, Z).lG;
            worst = Math.max(worst, Math.abs((Pfqn_ca.pfqn_ca(L, N, Z).lG - lLe) - M * C));
            assertEquals(M * C, Pfqn_ble.pfqn_ble(L, N, Z).lG - lLe, 1e-12);
        }
        assertTrue(worst < 0.1, "worst |(lG_ca - lG_le) - M*c| was " + worst + " nats");
    }

    /** An all-zero think time is the Z = 0 branch, so the count stays M-1. */
    @Test
    public void anAllZeroThinkTimeKeepsTheZeroBranchCount() {
        Matrix L = demands(4, 2, 4242L);
        Matrix N = pop(2, 20.0);
        Matrix Z = new Matrix(1, 2);
        assertEquals(3 * C, Pfqn_ble.pfqn_ble(L, N, Z).lG - Pfqn_le.pfqn_le(L, N, Z).lG, 1e-12);
    }
}
