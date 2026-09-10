package jline.api;

import jline.api.pfqn.nc.Pfqn_aghq;
import jline.api.pfqn.nc.Pfqn_ca;
import jline.api.pfqn.nc.Pfqn_le;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * The adaptive Gauss-Hermite rule over the simplex factor of the McKenna-Mitra integral.
 *
 * At Z = 0 the integrand is homogeneous, so the radius integrates exactly to gamma(N+M)
 * and every bit of the error of Pfqn_le sits in the closure of the simplex integral that
 * remains. Pfqn_aghq reads the mode and curvature of that closure as a quadrature rule
 * whose q = 1 member IS Pfqn_le. With Z &gt; 0 it integrates the radius rather than closing
 * it. The reference values are MATLAB pfqn_aghq on the same models, and agree with the
 * python and cpp ports to 1e-12.
 */
public class PfqnAghqTest {

    private static Matrix demands3x2() {
        Matrix L = new Matrix(3, 2);
        L.set(0, 0, 1.0); L.set(0, 1, 0.5);
        L.set(1, 0, 0.7); L.set(1, 1, 1.2);
        L.set(2, 0, 0.3); L.set(2, 1, 0.9);
        return L;
    }

    private static Matrix row(double... v) {
        Matrix m = new Matrix(1, v.length);
        for (int i = 0; i < v.length; i++) {
            m.set(0, i, v[i]);
        }
        return m;
    }

    @Test
    public void aghqAtOneNodeIsTheLogisticExpansion() {
        // A single node sits at the mode with weight sqrt(2 pi), so the rule collapses to
        // Cas17 eq. (34), up to the tolerance of the fixed point the two share.
        Matrix L = demands3x2();
        Matrix N = row(2.0, 3.0);
        Matrix Z0 = row(0.0, 0.0);
        assertEquals(Pfqn_le.pfqn_le(L, N, Z0).lG, Pfqn_aghq.pfqn_aghq(L, N, Z0, 1).lG, 1e-9);
    }

    @Test
    public void aghqConvergesInQ() {
        Matrix L = demands3x2();
        Matrix N = row(2.0, 3.0);
        Matrix Z0 = row(0.0, 0.0);
        double exact = Pfqn_ca.pfqn_ca(L, N, Z0).lG;
        double e2 = Math.abs(Pfqn_aghq.pfqn_aghq(L, N, Z0, 2).lG - exact);
        double e7 = Math.abs(Pfqn_aghq.pfqn_aghq(L, N, Z0, 7).lG - exact);
        assertTrue(e7 <= e2, "the rule must not get worse with more nodes: " + e2 + " -> " + e7);
        assertTrue(e7 < 0.02, "q=7 should be within 0.02 nats, was " + e7);
    }

    @Test
    public void theRuleIsExactAtASingleStation() {
        // M = 1 leaves a point mass on the simplex, so there is nothing left to close.
        Matrix L = new Matrix(1, 2);
        L.set(0, 0, 2.0);
        L.set(0, 1, 3.0);
        Matrix N = row(4.0, 1.0);
        double exact = Pfqn_ca.pfqn_ca(L, N, row(0.0, 0.0)).lG;
        assertEquals(exact, Pfqn_aghq.pfqn_aghq(L, N).lG, 1e-10);
    }

    @Test
    public void anAllZeroThinkTimeIsTheZeroBranch() {
        // Pfqn_nc always passes Z.sumCols(), so this is the shape a delay-free model
        // arrives in. Pfqn_le used to branch on isEmpty alone, which sent it down the
        // Z > 0 branch while Pfqn_ble still corrected it as if Z were zero.
        Matrix L = demands3x2();
        Matrix N = row(2.0, 3.0);
        Matrix Z0 = row(0.0, 0.0);
        assertEquals(Pfqn_le.pfqn_le(L, N).lG, Pfqn_le.pfqn_le(L, N, Z0).lG, 1e-12);
        assertEquals(Pfqn_aghq.pfqn_aghq(L, N).lG, Pfqn_aghq.pfqn_aghq(L, N, Z0).lG, 1e-12);
    }
}
