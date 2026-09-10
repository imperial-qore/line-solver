package jline.api;

import jline.api.pfqn.nc.Pfqn_bkt;
import jline.api.pfqn.nc.Pfqn_le;
import jline.api.pfqn.nc.Pfqn_lekt;
import jline.api.pfqn.nc.Pfqn_ble;
import jline.util.matrix.Matrix;
import org.apache.commons.math3.special.Gamma;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * pfqn_lekt computes the estimator pfqn_ble and pfqn_bkt share, on the cheaper side.
 * With a think time the two corrected expansions are one estimator (dual saddle points,
 * Sylvester on the Hessians), so the value must not depend on the route; without one the
 * LE side carries M kappa - r(N+M) in place of ble's (M-1) kappa, which lands it on the
 * KT value.
 */
public class PfqnLektTest {

    private static final double KAPPA = 1.0 - Math.log(2 * Math.PI) / 2;

    private static double r(double a) {
        return Gamma.logGamma(a) - (a - 0.5) * Math.log(a) + a - 0.5 * Math.log(2 * Math.PI);
    }

    private static Matrix row(double... v) {
        Matrix m = new Matrix(1, v.length);
        for (int i = 0; i < v.length; i++) m.set(0, i, v[i]);
        return m;
    }

    /** L = [1 .5; .7 1.2; .3 .9], 3 x 2: the KT side. */
    private static Matrix modelA() {
        Matrix L = new Matrix(3, 2);
        L.set(0, 0, 1.0); L.set(0, 1, 0.5);
        L.set(1, 0, 0.7); L.set(1, 1, 1.2);
        L.set(2, 0, 0.3); L.set(2, 1, 0.9);
        return L;
    }

    /** Its transpose, 2 x 3: the LE side. */
    private static Matrix wide() {
        Matrix L = new Matrix(2, 3);
        L.set(0, 0, 1.0); L.set(0, 1, 0.7); L.set(0, 2, 0.3);
        L.set(1, 0, 0.5); L.set(1, 1, 1.2); L.set(1, 2, 0.9);
        return L;
    }

    @Test
    public void theRouteFollowsTheSmallerDimensionAndSelfLoops() {
        assertEquals("kt", Pfqn_lekt.route(modelA(), row(2, 3), row(1, 0.5)));
        assertEquals("le", Pfqn_lekt.route(wide(), row(2, 3, 4), row(1, 0.5, 0.2)));
        assertEquals("le", Pfqn_lekt.route(wide(), row(2, 3, 4), new Matrix(1, 3)));
        Matrix Ls = new Matrix(2, 3);
        Ls.set(0, 0, 1.0); Ls.set(0, 2, 0.4); Ls.set(1, 0, 0.7); Ls.set(1, 1, 1.2);
        assertEquals("kt", Pfqn_lekt.route(Ls, row(3, 2, 2), new Matrix(1, 3)));
        assertEquals("le", Pfqn_lekt.route(Ls, row(3, 2, 2), row(0, 1, 1)));  // both single-station classes think
    }

    @Test
    public void theKtSideIsBkt() {
        Matrix N = row(2, 3), Z = row(1, 0.5);
        assertEquals(Pfqn_bkt.pfqn_bkt(modelA(), N, Z).lG, Pfqn_lekt.pfqn_lekt(modelA(), N, Z).lG, 0.0);
        assertEquals(Pfqn_bkt.pfqn_bkt(modelA(), N).lG, Pfqn_lekt.pfqn_lekt(modelA(), N).lG, 0.0);
    }

    @Test
    public void theLeSideWithAThinkTimeIsBleAndLandsOnBkt() {
        Matrix N = row(2, 3, 4), Z = row(1, 0.5, 0.2);
        double v = Pfqn_lekt.pfqn_lekt(wide(), N, Z).lG;
        assertEquals(Pfqn_ble.pfqn_ble(wide(), N, Z).lG, v, 0.0);
        double d = Math.abs(v - Pfqn_bkt.pfqn_bkt(wide(), N, Z).lG);
        assertTrue(d < 1e-5, "|lekt(le) - bkt| = " + d);   // solver floors, not the identity
    }

    @Test
    public void matchesMatlab() {
        // MATLAB pfqn_lekt on the transposed modelA (2 x 3, the LE side) and on a self-looping model
        assertEquals(7.64923112714508, Pfqn_lekt.pfqn_lekt(wide(), row(2, 3, 4), row(1, 0.5, 0.2)).lG, 1e-8);
        assertEquals(7.09110182626673, Pfqn_lekt.pfqn_lekt(wide(), row(2, 3, 4)).lG, 1e-8);
        Matrix Ls = new Matrix(2, 3);
        Ls.set(0, 0, 1.0); Ls.set(0, 2, 0.4); Ls.set(1, 0, 0.7); Ls.set(1, 1, 1.2);
        assertEquals(2.08635455377189, Pfqn_lekt.pfqn_lekt(Ls, row(3, 2, 2)).lG, 1e-8);
    }

    @Test
    public void theLeSideWithoutAThinkTimeCarriesTheCommonConstant() {
        Matrix N = row(2, 3, 4), Z0 = new Matrix(1, 3);
        double eta = 9.0 + 2;
        double v = Pfqn_lekt.pfqn_lekt(wide(), N, Z0).lG;
        assertEquals(Pfqn_le.pfqn_le(wide(), N, Z0).lG + 2 * KAPPA - r(eta), v, 1e-12);
        assertEquals(Pfqn_ble.pfqn_ble(wide(), N, Z0).lG + KAPPA - r(eta), v, 1e-12);
        double d = Math.abs(v - Pfqn_bkt.pfqn_bkt(wide(), N, Z0).lG);
        assertTrue(d < 1e-5, "|lekt(le) - bkt| at Z=0 = " + d);  // bkt - LE = M kappa - r(eta)
    }
}
