package jline.lib;

import jline.lib.butools.mc.CheckGeneratorKt;
import jline.lib.butools.mc.CheckProbMatrixKt;
import jline.lib.butools.ph.PH2From3MomentsKt;
import jline.lib.butools.ph.PH2Representation;
import jline.lib.lti.customromberg;
import jline.lib.lti.function_wrapper;
import jline.lib.perm.AdaPartSampler;
import jline.lib.perm.BethePermanent;
import jline.lib.perm.HeuristicPermanent;
import jline.lib.perm.HuberLawSampler;
import jline.lib.perm.NaivePermanent;
import jline.lib.perm.RyzerPermanent;
import jline.util.matrix.Matrix;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apfloat.Apcomplex;
import org.apfloat.Apfloat;
import org.junit.jupiter.api.Test;

import java.math.BigDecimal;
import java.util.function.UnaryOperator;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Coverage tests for utility functions with 0% coverage.
 */
public class UtilityCoverageTest {

    private static final double TOL = 1e-6;

    // ========== CheckGenerator ==========

    @Test
    void testCheckGenerator() {
        // Create a valid 2x2 generator matrix
        // Row sums must be 0, diagonal negative, off-diagonal non-negative
        Matrix validGen = new Matrix(2, 2);
        validGen.set(0, 0, -1.0);
        validGen.set(0, 1, 1.0);
        validGen.set(1, 0, 0.5);
        validGen.set(1, 1, -0.5);

        assertTrue(CheckGeneratorKt.checkGenerator(validGen, false, 1e-14));

        // Create an invalid generator (positive diagonal)
        Matrix invalidGen = new Matrix(2, 2);
        invalidGen.set(0, 0, 1.0);  // Invalid: diagonal should be negative
        invalidGen.set(0, 1, -1.0);
        invalidGen.set(1, 0, 0.5);
        invalidGen.set(1, 1, -0.5);

        assertFalse(CheckGeneratorKt.checkGenerator(invalidGen, false, 1e-14));

        // Create an invalid generator (non-zero row sum)
        Matrix invalidRowSum = new Matrix(2, 2);
        invalidRowSum.set(0, 0, -1.0);
        invalidRowSum.set(0, 1, 2.0);  // Row sum = 1.0 (should be 0)
        invalidRowSum.set(1, 0, 0.5);
        invalidRowSum.set(1, 1, -0.5);

        assertFalse(CheckGeneratorKt.checkGenerator(invalidRowSum, false, 1e-14));

        // Non-square matrix should fail
        Matrix nonSquare = new Matrix(2, 3);
        assertFalse(CheckGeneratorKt.checkGenerator(nonSquare, false, 1e-14));
    }

    // ========== CheckProbMatrix ==========

    @Test
    void testCheckProbMatrix() {
        // Create a valid stochastic probability matrix (row sums = 1)
        Matrix validProb = new Matrix(2, 2);
        validProb.set(0, 0, 0.3);
        validProb.set(0, 1, 0.7);
        validProb.set(1, 0, 0.4);
        validProb.set(1, 1, 0.6);

        assertTrue(CheckProbMatrixKt.checkProbMatrix(validProb, false, 1e-14));

        // Create an invalid probability matrix (negative element)
        Matrix invalidProb = new Matrix(2, 2);
        invalidProb.set(0, 0, -0.1);  // Invalid: negative
        invalidProb.set(0, 1, 1.1);
        invalidProb.set(1, 0, 0.4);
        invalidProb.set(1, 1, 0.6);

        assertFalse(CheckProbMatrixKt.checkProbMatrix(invalidProb, false, 1e-14));

        // Create an invalid probability matrix (row sum != 1)
        Matrix invalidRowSum = new Matrix(2, 2);
        invalidRowSum.set(0, 0, 0.3);
        invalidRowSum.set(0, 1, 0.5);  // Row sum = 0.8 (should be 1)
        invalidRowSum.set(1, 0, 0.4);
        invalidRowSum.set(1, 1, 0.6);

        assertFalse(CheckProbMatrixKt.checkProbMatrix(invalidRowSum, false, 1e-14));

        // Non-square matrix should fail
        Matrix nonSquare = new Matrix(2, 3);
        assertFalse(CheckProbMatrixKt.checkProbMatrix(nonSquare, false, 1e-14));
    }

    // ========== customromberg ==========

    @Test
    void testCustomRomberg() {
        // Test numerical integration: integral of exp(-x) from 1 to 2
        // Analytical result: -exp(-2) - (-exp(-1)) = exp(-1) - exp(-2) ≈ 0.2325
        UnivariateFunction f = x -> Math.exp(-x);

        BigDecimal result = customromberg.INSTANCE.starter(f, 1.0, 2.0);
        assertNotNull(result);

        double expected = Math.exp(-1.0) - Math.exp(-2.0);
        assertEquals(expected, result.doubleValue(), 1e-4);  // Romberg numerical tolerance

        // Test hn method
        BigDecimal hn = customromberg.INSTANCE.hn(0.0, 1.0, 4);
        assertEquals(0.25, hn.doubleValue(), TOL);
    }

    // ========== function_wrapper ==========

    @Test
    void testFunctionWrapper() {
        // Test add function
        UnaryOperator<Apcomplex> f1 = x -> x;  // identity
        UnaryOperator<Apcomplex> f2 = x -> new Apcomplex(new Apfloat(1.0));  // constant 1

        UnaryOperator<Apcomplex> sum = function_wrapper.INSTANCE.add(f1, f2);
        Apcomplex input = new Apcomplex(new Apfloat(2.0));
        Apcomplex result = sum.apply(input);
        // f1(2) + f2(2) = 2 + 1 = 3
        assertEquals(3.0, result.real().doubleValue(), TOL);

        // Test multiply function
        UnaryOperator<Apcomplex> product = function_wrapper.INSTANCE.multiply(f1, f2);
        result = product.apply(input);
        // f1(2) * f2(2) = 2 * 1 = 2
        assertEquals(2.0, result.real().doubleValue(), TOL);
    }

    // ========== NaivePermanent ==========

    @Test
    void testNaivePermanent() {
        // Test with a simple 2x2 matrix
        // perm([[1, 2], [3, 4]]) = 1*4 + 2*3 = 10
        Matrix m = new Matrix(2, 2);
        m.set(0, 0, 1.0);
        m.set(0, 1, 2.0);
        m.set(1, 0, 3.0);
        m.set(1, 1, 4.0);

        NaivePermanent np = new NaivePermanent(m, true);  // solve immediately
        assertEquals(10.0, np.getValue(), TOL);

        // Test with solve=false then call solve()
        NaivePermanent np2 = new NaivePermanent(m, false);
        np2.solve();
        assertEquals(10.0, np2.getValue(), TOL);
    }

    // ========== RyzerPermanent ==========

    @Test
    void testRyzerPermanent() {
        // Test with the same 2x2 matrix
        Matrix m = new Matrix(2, 2);
        m.set(0, 0, 1.0);
        m.set(0, 1, 2.0);
        m.set(1, 0, 3.0);
        m.set(1, 1, 4.0);

        // Test graycode mode (default)
        RyzerPermanent rp1 = new RyzerPermanent(m, "graycode", true);
        assertEquals(10.0, rp1.getValue(), TOL);

        // Test naive mode
        RyzerPermanent rp2 = new RyzerPermanent(m, "naive", true);
        assertEquals(10.0, rp2.getValue(), TOL);

        // Test without immediate solve
        RyzerPermanent rp3 = new RyzerPermanent(m, "graycode", false);
        rp3.solve();
        assertEquals(10.0, rp3.getValue(), TOL);
    }

    // ========== HeuristicPermanent ==========

    @Test
    void testHeuristicPermanent() {
        // Test with a positive 2x2 matrix
        Matrix m = new Matrix(2, 2);
        m.set(0, 0, 1.0);
        m.set(0, 1, 2.0);
        m.set(1, 0, 3.0);
        m.set(1, 1, 4.0);

        // Heuristic permanent is an approximation
        HeuristicPermanent hp = new HeuristicPermanent(m, 1e-10, 1000, true);
        double result = hp.getValue();
        assertTrue(Double.isFinite(result));
        // The approximation should be in the ballpark
        assertTrue(result > 5 && result < 15);

        // Test with constructor that only takes matrix
        HeuristicPermanent hp2 = new HeuristicPermanent(m);
        hp2.solve();
        assertTrue(Double.isFinite(hp2.getValue()));

        // Test with solve=true shorthand
        HeuristicPermanent hp3 = new HeuristicPermanent(m, true);
        assertTrue(Double.isFinite(hp3.getValue()));

        // Test that negative matrix throws
        Matrix negative = new Matrix(2, 2);
        negative.set(0, 0, -1.0);  // negative element
        negative.set(0, 1, 2.0);
        negative.set(1, 0, 3.0);
        negative.set(1, 1, 4.0);

        assertThrows(IllegalArgumentException.class, () -> new HeuristicPermanent(negative, true));
    }

    // ========== PH2From3Moments ==========

    @Test
    void testPH2From3Moments() {
        // Test with moments that correspond to a valid PH(2) distribution
        // For an exponential with rate 1: m1=1, m2=2, m3=6
        double[] expMoments = {1.0, 2.0, 6.0};
        PH2Representation expResult = PH2From3MomentsKt.ph2From3Moments(expMoments, 1e-14);

        assertNotNull(expResult);
        assertNotNull(expResult.getAlpha());
        assertNotNull(expResult.getA());

        // For an exponential, we expect a 1x1 matrix
        assertEquals(1, expResult.getAlpha().getNumCols());
        assertEquals(1, expResult.getA().getNumRows());

        // Test with moments that require a 2-phase distribution
        // Using moments from a hyper-exponential
        double[] hyperMoments = {1.0, 3.0, 15.0};
        PH2Representation hyperResult = PH2From3MomentsKt.ph2From3Moments(hyperMoments, 1e-14);

        assertNotNull(hyperResult);
        assertEquals(2, hyperResult.getAlpha().getNumCols());  // 2-phase
        assertEquals(2, hyperResult.getA().getNumRows());

        // Test that infeasible moments throw
        double[] invalidMoments = {1.0, 0.5, 1.0};  // m2 too small
        assertThrows(IllegalArgumentException.class,
                () -> PH2From3MomentsKt.ph2From3Moments(invalidMoments, 1e-14));
    }

    // ========== BethePermanent ==========

    @Test
    void testBethePermanent() {
        // Create a simple 3x3 matrix
        Matrix m = new Matrix(3, 3);
        m.set(0, 0, 1.0);
        m.set(0, 1, 2.0);
        m.set(0, 2, 3.0);
        m.set(1, 0, 4.0);
        m.set(1, 1, 5.0);
        m.set(1, 2, 6.0);
        m.set(2, 0, 7.0);
        m.set(2, 1, 8.0);
        m.set(2, 2, 9.0);

        // Create BethePermanent solver
        BethePermanent bp = new BethePermanent(m, 0.001, 1000, false);
        assertNotNull(bp);

        // Solve
        bp.solve();
        double result = bp.getValue();

        // The Bethe permanent approximates the permanent
        // For this matrix, permanent = 1*5*9 + 2*6*7 + 3*4*8 - 3*5*7 - 2*4*9 - 1*6*8
        // = 45 + 84 + 96 - 105 - 72 - 48 = 0
        // But Bethe permanent is an approximation, so just check it returns a value
        assertTrue(Double.isFinite(result));
    }

    // ========== HuberLawSampler ==========

    @Test
    void testHuberLawSampler_sampleMode() {
        // Create a simple 2x2 positive matrix
        Matrix m = new Matrix(2, 2);
        m.set(0, 0, 1.0);
        m.set(0, 1, 2.0);
        m.set(1, 0, 3.0);
        m.set(1, 1, 4.0);

        // Use sample mode with limited samples for fast test
        HuberLawSampler sampler = new HuberLawSampler(
                m,
                0.1,        // delta
                0.000001,   // alpha2
                0.1,        // epsilon
                "sample",   // mode
                50,         // numberOfSamples
                30000,      // maximumTime
                true        // solve
        );

        double result = sampler.getValue();
        assertTrue(Double.isFinite(result));
        assertTrue(result > 0);

        // Actual permanent is 1*4 + 2*3 = 10
        // Sampler should be in ballpark (may not be exact due to sampling)
        assertTrue(result > 1 && result < 100);

        // Check tracking variables are populated
        assertFalse(sampler.getSampleAccepted().isEmpty());
        assertFalse(sampler.getSampleTime().isEmpty());
        assertFalse(sampler.getPermStep().isEmpty());
    }

    @Test
    void testHuberLawSampler_timeMode() {
        // Create a simple 2x2 positive matrix
        Matrix m = new Matrix(2, 2);
        m.set(0, 0, 1.0);
        m.set(0, 1, 1.0);
        m.set(1, 0, 1.0);
        m.set(1, 1, 1.0);

        // Use time mode with short time limit
        HuberLawSampler sampler = new HuberLawSampler(
                m,
                0.1,        // delta
                0.000001,   // alpha2
                0.1,        // epsilon
                "time",     // mode
                1000,       // numberOfSamples
                100,        // maximumTime (100ms)
                true        // solve
        );

        double result = sampler.getValue();
        assertTrue(Double.isFinite(result));

        // Permanent of all-ones 2x2 matrix is 2
        // Result should be reasonable approximation
        assertTrue(result >= 0);
    }

    @Test
    void testHuberLawSampler_noAutoSolve() {
        // Test constructor without automatic solve
        Matrix m = new Matrix(2, 2);
        m.set(0, 0, 1.0);
        m.set(0, 1, 2.0);
        m.set(1, 0, 3.0);
        m.set(1, 1, 4.0);

        HuberLawSampler sampler = new HuberLawSampler(
                m,
                0.1,
                0.000001,
                0.1,
                "sample",
                20,
                30000,
                false  // don't solve immediately
        );

        // Value should be uninitialized (0.0)
        assertEquals(0.0, sampler.getValue());

        // Now solve manually
        sampler.solve();
        assertTrue(sampler.getValue() > 0);
    }

    @Test
    void testHuberLawSampler_3x3Matrix() {
        // Test with a 3x3 matrix
        Matrix m = new Matrix(3, 3);
        m.set(0, 0, 1.0);
        m.set(0, 1, 2.0);
        m.set(0, 2, 3.0);
        m.set(1, 0, 4.0);
        m.set(1, 1, 5.0);
        m.set(1, 2, 6.0);
        m.set(2, 0, 7.0);
        m.set(2, 1, 8.0);
        m.set(2, 2, 9.0);

        HuberLawSampler sampler = new HuberLawSampler(
                m,
                0.2,
                0.00001,
                0.2,
                "sample",
                30,
                30000,
                true
        );

        double result = sampler.getValue();
        assertTrue(Double.isFinite(result));
    }

    // ========== AdaPartSampler ==========

    @Test
    void testAdaPartSampler_classicMode() {
        // Create a simple 2x2 positive matrix
        Matrix m = new Matrix(2, 2);
        m.set(0, 0, 1.0);
        m.set(0, 1, 2.0);
        m.set(1, 0, 3.0);
        m.set(1, 1, 4.0);

        // Use classic mode with limited accepted samples
        AdaPartSampler sampler = new AdaPartSampler(
                m,
                10,        // maximumAcceptedSamples
                30000,     // maximumTime
                100,       // maximumSamples
                "classic", // mode
                true       // solve
        );

        double result = sampler.getValue();
        assertTrue(Double.isFinite(result));
        assertTrue(result > 0);

        // Check tracking variables are populated
        assertFalse(sampler.getSampleAccepted().isEmpty());
        assertFalse(sampler.getSampleTime().isEmpty());
        assertFalse(sampler.getPermStep().isEmpty());
    }

    @Test
    void testAdaPartSampler_sampleMode() {
        // Create a simple 2x2 matrix
        Matrix m = new Matrix(2, 2);
        m.set(0, 0, 1.0);
        m.set(0, 1, 1.0);
        m.set(1, 0, 1.0);
        m.set(1, 1, 1.0);

        // Use sample mode
        AdaPartSampler sampler = new AdaPartSampler(
                m,
                100,       // maximumAcceptedSamples
                30000,     // maximumTime
                20,        // maximumSamples (limited)
                "sample",  // mode
                true       // solve
        );

        double result = sampler.getValue();
        assertTrue(Double.isFinite(result));
        assertTrue(result > 0);

        // All-ones 2x2 matrix has permanent = 2
        // Result should be reasonable
        assertTrue(result > 0 && result < 10);
    }

    @Test
    void testAdaPartSampler_timeMode() {
        // Test time mode with short time limit
        Matrix m = new Matrix(2, 2);
        m.set(0, 0, 2.0);
        m.set(0, 1, 1.0);
        m.set(1, 0, 1.0);
        m.set(1, 1, 2.0);

        AdaPartSampler sampler = new AdaPartSampler(
                m,
                100,      // maximumAcceptedSamples
                100,      // maximumTime (100ms)
                1000,     // maximumSamples
                "time",   // mode
                true      // solve
        );

        double result = sampler.getValue();
        assertTrue(Double.isFinite(result));
        // Permanent = 2*2 + 1*1 = 5
        assertTrue(result >= 0);
    }

    @Test
    void testAdaPartSampler_noAutoSolve() {
        // Test constructor without automatic solve
        Matrix m = new Matrix(2, 2);
        m.set(0, 0, 1.0);
        m.set(0, 1, 2.0);
        m.set(1, 0, 3.0);
        m.set(1, 1, 4.0);

        AdaPartSampler sampler = new AdaPartSampler(
                m,
                10,
                30000,
                50,
                "sample",
                false  // don't solve immediately
        );

        // Value should be uninitialized (0.0)
        assertEquals(0.0, sampler.getValue());

        // Now solve manually
        sampler.solve();
        assertTrue(sampler.getValue() > 0);
    }

    @Test
    void testAdaPartSampler_3x3Matrix() {
        // Test with a 3x3 matrix
        Matrix m = new Matrix(3, 3);
        m.set(0, 0, 1.0);
        m.set(0, 1, 2.0);
        m.set(0, 2, 1.0);
        m.set(1, 0, 2.0);
        m.set(1, 1, 1.0);
        m.set(1, 2, 2.0);
        m.set(2, 0, 1.0);
        m.set(2, 1, 2.0);
        m.set(2, 2, 1.0);

        AdaPartSampler sampler = new AdaPartSampler(
                m,
                5,
                30000,
                30,
                "classic",
                true
        );

        double result = sampler.getValue();
        assertTrue(Double.isFinite(result));
        assertTrue(result > 0);
    }
}
