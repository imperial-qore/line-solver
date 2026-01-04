package jline.lib;

import jline.lang.processes.MAP;
import jline.lib.kpctoolbox.MAPCatalog;
import jline.lib.lti.euler;
import jline.lib.lti.laguerre;
import jline.lib.lti.talbot;
import jline.lib.perm.BethePermanent;
import jline.util.matrix.Matrix;
import org.apache.commons.math3.complex.Complex;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Coverage tests for library classes with 0% coverage.
 */
public class LibCoverageTest {

    private static final double TOL = 1e-6;

    // ========== lib.lti.euler ==========

    @Test
    void testEuler() {
        // Test getalpha
        ArrayList<Complex> alpha = euler.INSTANCE.getalpha(5);
        assertNotNull(alpha);
        assertEquals(5, alpha.size());

        // Check that alpha values are complex
        assertNotNull(alpha.get(0));
        assertTrue(alpha.get(0).getReal() > 0);  // Should be positive real part

        // Test geteta
        double[] eta = euler.INSTANCE.geteta(5);
        assertNotNull(eta);
        assertEquals(5, eta.length);
        assertEquals(0.5, eta[0], TOL);  // First element should be 0.5

        // Test getomega
        ArrayList<Complex> omega = euler.INSTANCE.getomega(5);
        assertNotNull(omega);
        assertEquals(5, omega.size());
    }

    // ========== lib.lti.talbot ==========

    @Test
    void testTalbot() {
        // Test getalpha
        ArrayList<Complex> alpha = talbot.INSTANCE.getalpha(5);
        assertNotNull(alpha);
        assertEquals(5, alpha.size());

        // First element has 0 imaginary part
        assertEquals(0.0, alpha.get(0).getImaginary(), TOL);

        // Test getomega
        ArrayList<Complex> omega = talbot.INSTANCE.getomega(5, alpha);
        assertNotNull(omega);
        assertEquals(5, omega.size());
    }

    // ========== lib.lti.laguerre ==========

    @Test
    void testLaguerre() {
        // Test getLaguerreCoefficients
        double[] coeffs = laguerre.INSTANCE.getLaguerreCoefficients(5);
        assertNotNull(coeffs);
        assertEquals(6, coeffs.length);  // n+1 coefficients

        // First coefficient should be 1 (divided by 0! = 1)
        assertEquals(1.0, coeffs[0], TOL);

        // Test getLaguerreRoots
        double[] roots = laguerre.INSTANCE.getLaguerreRoots(coeffs);
        assertNotNull(roots);
        assertEquals(5, roots.length);  // n roots
    }

    // ========== lib.kpctoolbox.MAPCatalog ==========

    @Test
    void testMAPCatalog() {
        // Test various MAP factory methods
        MAP mama1 = MAPCatalog.mama1();
        assertNotNull(mama1);
        assertEquals(2, mama1.getNumberOfPhases());
        assertTrue(mama1.getMean() > 0);

        MAP mama3 = MAPCatalog.mama3();
        assertNotNull(mama3);
        assertEquals(2, mama3.getNumberOfPhases());

        // Test MMPP2 models
        MAP mmpp2Lrd = MAPCatalog.mmpp2_lrd();
        assertNotNull(mmpp2Lrd);

        MAP mmpp2Srd = MAPCatalog.mmpp2_srd();
        assertNotNull(mmpp2Srd);

        // Test hyperexponential MAP
        MAP hyper8 = MAPCatalog.hyper8();
        assertNotNull(hyper8);
        assertEquals(8, hyper8.getNumberOfPhases());

        // Test hypothetical model
        MAP hypo = MAPCatalog.hypo065_095();
        assertNotNull(hypo);

        // Test example models
        MAP example1 = MAPCatalog.map_example1();
        assertNotNull(example1);

        MAP poisson = MAPCatalog.poisson_example();
        assertNotNull(poisson);
        assertEquals(1, poisson.getNumberOfPhases());  // Poisson is 1-state
    }


    // ========== lib.lti.cme (complex, test simpler functionality) ==========
    // Note: cme requires UnivariateFunction and complex high-precision arithmetic.
    // Testing simpler helper methods.

    @Test
    void testCmeHelpers() {
        // Test deepcopy helper
        ArrayList<Double> original = new ArrayList<>();
        original.add(1.0);
        original.add(2.0);
        original.add(3.0);

        ArrayList<Double> copy = jline.lib.lti.cme.INSTANCE.deepcopy(original);
        assertNotNull(copy);
        assertEquals(3, copy.size());
        assertEquals(1.0, copy.get(0), TOL);
        assertEquals(2.0, copy.get(1), TOL);
        assertEquals(3.0, copy.get(2), TOL);

        // Verify it's a deep copy
        original.set(0, 99.0);
        assertEquals(1.0, copy.get(0), TOL);  // copy unchanged

        // Test getnormalrandom
        ArrayList<Double> randoms = jline.lib.lti.cme.INSTANCE.getnormalrandom(10);
        assertNotNull(randoms);
        assertEquals(10, randoms.size());

        // Test binaryCheck
        boolean[] binary = jline.lib.lti.cme.INSTANCE.binaryCheck(5, 4);  // 5 = 101 in binary
        assertNotNull(binary);
        assertEquals(4, binary.length);
    }
}
