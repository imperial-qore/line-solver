package jline.lib.m3a;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.BeforeEach;
import static org.junit.jupiter.api.Assertions.*;
import static jline.TestTools.FINE_TOL;

/**
 * Tests for M3A library functions.
 */
public class M3ATest {

    private MatrixCell sampleMMAP;

    @BeforeEach
    public void setUp() {
        // Create a sample 2x2 MMAP for testing
        sampleMMAP = new MatrixCell();
        
        // D0 matrix (generator matrix) - properly constructed to ensure invertibility
        Matrix D0 = new Matrix(2, 2);
        D0.set(0, 0, -3.0);  // Larger negative diagonal to ensure stability
        D0.set(0, 1, 0.8);
        D0.set(1, 0, 0.5);
        D0.set(1, 1, -2.0);
        sampleMMAP.set(0, D0);
        
        // D1 matrix (sum of marking matrices) - ensure D0 + D1 has zero row sums
        Matrix D1 = new Matrix(2, 2);
        D1.set(0, 0, 1.2);   // 3.0 - 0.8 = 2.2, so we need 1.2 + 1.0 = 2.2
        D1.set(0, 1, 1.0);
        D1.set(1, 0, 0.8);   // 2.0 - 0.5 = 1.5, so we need 0.8 + 0.7 = 1.5
        D1.set(1, 1, 0.7);
        sampleMMAP.set(1, D1);
          
        // D11 matrix (marking matrix for class 1)
        Matrix D11 = new Matrix(2, 2);
        D11.set(0, 0, 0.8);
        D11.set(0, 1, 0.3);
        D11.set(1, 0, 0.4);
        D11.set(1, 1, 0.2);
        sampleMMAP.set(2, D11);
        
        // D12 matrix (marking matrix for class 2)
        Matrix D12 = new Matrix(2, 2);
        D12.set(0, 0, 0.4);
        D12.set(0, 1, 0.7);
        D12.set(1, 0, 0.4);
        D12.set(1, 1, 0.5);
        sampleMMAP.set(3, D12);
    }

    @Test
    public void testMomentComputation() {
        // Test moment computation for different orders
        Double[] moments1 = M3aUtils.Companion.computeMoments(sampleMMAP, 1);
        Double[] moments2 = M3aUtils.Companion.computeMoments(sampleMMAP, 2);
        Double[] moments3 = M3aUtils.Companion.computeMoments(sampleMMAP, 3);
        
        // Verify that moments are computed
        assertEquals(1, moments1.length, "Should compute 1 moment");
        assertEquals(2, moments2.length, "Should compute 2 moments");
        assertEquals(3, moments3.length, "Should compute 3 moments");
        
        // First moment should be positive (mean inter-arrival time)
        assertTrue(moments1[0] > 0, "First moment should be positive");
        assertTrue(moments2[0] > 0, "First moment should be positive");
        assertTrue(moments3[0] > 0, "First moment should be positive");
        
        // Higher moments should be positive
        assertTrue(moments2[1] > 0, "Second moment should be positive");
        assertTrue(moments3[1] > 0, "Second moment should be positive");
        assertTrue(moments3[2] > 0, "Third moment should be positive");
    }

    @Test
    public void testAutocorrelationComputation() {
        // Test autocorrelation computation
        Double[] acf = M3aUtils.Companion.computeAutocorrelation(sampleMMAP, 5);
        
        // Verify autocorrelation function properties
        assertEquals(5, acf.length, "Should compute 5 autocorrelation lags");
        
        // Autocorrelation should be bounded
        for (int i = 0; i < acf.length; i++) {
            assertTrue(acf[i] >= -1.0, "Autocorrelation should be >= -1");
            assertTrue(acf[i] <= 1.0, "Autocorrelation should be <= 1");
        }
        
        // For a well-behaved MMAP, autocorrelation should decay
        // (This is a statistical property, not guaranteed mathematically)
        assertTrue(Math.abs(acf[0]) >= Math.abs(acf[acf.length - 1]) * 0.1,
            "Autocorrelation should generally decay with lag");
    }

    @Test
    public void testCoefficientOfVariation() {
        // Test coefficient of variation computation
        double cv = M3aUtils.Companion.computeCoeffVar(sampleMMAP);
        
        // Coefficient of variation should be positive
        assertTrue(cv > 0, "Coefficient of variation should be positive");
        
        // For a reasonable MMAP, CV should be finite
        assertTrue(cv < Double.POSITIVE_INFINITY, "Coefficient of variation should be finite");
        assertTrue(!Double.isNaN(cv), "Coefficient of variation should not be NaN");
    }

    @Test
    public void testIDCComputation() {
        // Test Index of Dispersion for Counts computation
        // Note: time window parameter is ignored in revised implementation (uses asymptotic IDC)
        double idc1 = M3aUtils.Companion.computeIDC(sampleMMAP, 1.0);
        double idc2 = M3aUtils.Companion.computeIDC(sampleMMAP, 2.0);
        
        // IDC should be positive
        assertTrue(idc1 > 0, "IDC should be positive");
        assertTrue(idc2 > 0, "IDC should be positive");
        
        // Since using asymptotic IDC, both values should be the same
        assertEquals(idc1, idc2, FINE_TOL, "Asymptotic IDC should be independent of time window");
        
        // IDC should be finite
        assertTrue(idc1 < Double.POSITIVE_INFINITY, "IDC should be finite");
        assertTrue(!Double.isNaN(idc1), "IDC should not be NaN");
    }

    @Test
    public void testSpectralGapComputation() {
        // Test spectral gap computation
        double spectralGap = M3aUtils.Companion.computeSpectralGap(sampleMMAP);
        
        // Spectral gap should be non-negative
        assertTrue(spectralGap >= 0, "Spectral gap should be non-negative");
        
        // For a stable MMAP, spectral gap should be finite
        assertTrue(spectralGap < Double.POSITIVE_INFINITY, "Spectral gap should be finite");
        assertTrue(!Double.isNaN(spectralGap), "Spectral gap should not be NaN");
    }

    @Test
    public void testKLDivergence() {
        // Test KL divergence computation between two MMAPs
        MatrixCell mmap2 = new MatrixCell();
        
        // Create a slightly different MMAP
        Matrix D0_2 = new Matrix(2, 2);
        D0_2.set(0, 0, -1.8);
        D0_2.set(0, 1, 0.4);
        D0_2.set(1, 0, 0.25);
        D0_2.set(1, 1, -1.3);
        mmap2.set(0, D0_2);
        
        Matrix D1_2 = new Matrix(2, 2);
        D1_2.set(0, 0, 0.7);
        D1_2.set(0, 1, 0.7);
        D1_2.set(1, 0, 0.55);
        D1_2.set(1, 1, 0.5);
        mmap2.set(1, D1_2);
        
        double kl = M3aUtils.Companion.computeKLDivergence(sampleMMAP, mmap2, 1000, 23000);
        
        // KL divergence should be non-negative
        assertTrue(kl >= 0, "KL divergence should be non-negative");
        
        // Should be finite
        assertTrue(kl < Double.POSITIVE_INFINITY, "KL divergence should be finite");
        assertTrue(!Double.isNaN(kl), "KL divergence should not be NaN");
        
        // KL divergence of MMAP with itself should be close to 0
        double selfKL = M3aUtils.Companion.computeKLDivergence(sampleMMAP, sampleMMAP, 1000, 23000);
        assertTrue(selfKL < 0.1, "Self KL divergence should be close to 0");
    }

    @Test
    public void testParameterOptimization() {
        // Test COBYLA optimization functionality
        Double[] initialParams = {1.0, 0.5};
        
        // Simple quadratic objective function: (x-2)^2 + (y-1)^2
        kotlin.jvm.functions.Function1<Double[], Double> objective = params -> {
            double x = params[0];
            double y = params[1];
            return (x - 2.0) * (x - 2.0) + (y - 1.0) * (y - 1.0);
        };
        
        // Constraint: x >= 0, y >= 0
        @SuppressWarnings("unchecked")
        kotlin.jvm.functions.Function1<Double[], Double>[] constraints = new kotlin.jvm.functions.Function1[2];
        constraints[0] = params -> params[0]; // x >= 0
        constraints[1] = params -> params[1]; // y >= 0
        
        Double[] result = M3aUtils.Companion.optimizeParameters(initialParams, objective, constraints, 1e-6);
        
        // Result should be close to optimal (2, 1)
        assertEquals(2, result.length, "Should return 2 parameters");
        assertTrue(result[0] >= 0, "First parameter should be non-negative");
        assertTrue(result[1] >= 0, "Second parameter should be non-negative");
        
        // Should be reasonably close to optimum
        assertTrue(Math.abs(result[0] - 2.0) < 1.0, "First parameter should be close to 2");
        assertTrue(Math.abs(result[1] - 1.0) < 1.0, "Second parameter should be close to 1");
    }

    @Test
    public void testMmapValidation() {
        // Test valid MMAP
        assertTrue(M3aUtils.Companion.validateMMAP(sampleMMAP), "Sample MMAP should be valid");
        
        // Test invalid MMAP - empty
        MatrixCell emptyMMAP = new MatrixCell();
        assertFalse(M3aUtils.Companion.validateMMAP(emptyMMAP), "Empty MMAP should be invalid");
        
        // Test invalid MMAP - only D0
        MatrixCell incompleteMAP = new MatrixCell();
        incompleteMAP.set(0, sampleMMAP.get(0));
        assertFalse(M3aUtils.Companion.validateMMAP(incompleteMAP), "Incomplete MMAP should be invalid");
        
        // Test invalid MMAP - mismatched dimensions
        MatrixCell mismatchedMMAP = new MatrixCell();
        mismatchedMMAP.set(0, sampleMMAP.get(0));
        Matrix wrongSizeD1 = new Matrix(3, 3); // Different size
        wrongSizeD1.set(0, 0, 1.0);
        mismatchedMMAP.set(1, wrongSizeD1);
        assertFalse(M3aUtils.Companion.validateMMAP(mismatchedMMAP), "Mismatched dimensions should be invalid");
        
        // Test invalid MMAP - positive diagonal in D0
        MatrixCell invalidDiagMMAP = new MatrixCell();
        Matrix invalidD0 = new Matrix(2, 2);
        invalidD0.set(0, 0, 1.0); // Positive diagonal (invalid)
        invalidD0.set(0, 1, 0.5);
        invalidD0.set(1, 0, 0.3);
        invalidD0.set(1, 1, -1.5);
        invalidDiagMMAP.set(0, invalidD0);
        invalidDiagMMAP.set(1, sampleMMAP.get(1));
        assertFalse(M3aUtils.Companion.validateMMAP(invalidDiagMMAP), "Positive diagonal in D0 should be invalid");
        
        // Test invalid MMAP - negative off-diagonal in D0
        MatrixCell negativeOffDiagMMAP = new MatrixCell();
        Matrix negativeD0 = new Matrix(2, 2);
        negativeD0.set(0, 0, -2.0);
        negativeD0.set(0, 1, -0.5); // Negative off-diagonal (invalid)
        negativeD0.set(1, 0, 0.3);
        negativeD0.set(1, 1, -1.5);
        negativeOffDiagMMAP.set(0, negativeD0);
        negativeOffDiagMMAP.set(1, sampleMMAP.get(1));
        assertFalse(M3aUtils.Companion.validateMMAP(negativeOffDiagMMAP), "Negative off-diagonal in D0 should be invalid");
        
        // Test invalid MMAP - negative elements in D1
        MatrixCell negativeD1MMAP = new MatrixCell();
        negativeD1MMAP.set(0, sampleMMAP.get(0));
        Matrix negativeD1 = new Matrix(2, 2);
        negativeD1.set(0, 0, -0.5); // Negative element (invalid)
        negativeD1.set(0, 1, 0.7);
        negativeD1.set(1, 0, 0.6);
        negativeD1.set(1, 1, 0.6);
        negativeD1MMAP.set(1, negativeD1);
        assertFalse(M3aUtils.Companion.validateMMAP(negativeD1MMAP), "Negative elements in D1 should be invalid");
    }

    @Test
    public void testEdgeCases() {
        // Test with 1x1 MMAP
        MatrixCell oneDimMMAP = new MatrixCell();
        Matrix d0_1x1 = new Matrix(1, 1);
        d0_1x1.set(0, 0, -1.0);
        oneDimMMAP.set(0, d0_1x1);
        
        Matrix d1_1x1 = new Matrix(1, 1);
        d1_1x1.set(0, 0, 1.0);
        oneDimMMAP.set(1, d1_1x1);
        
        // Should handle 1x1 case gracefully
        assertTrue(M3aUtils.Companion.validateMMAP(oneDimMMAP), "1x1 MMAP should be valid");
        
        Double[] moments = M3aUtils.Companion.computeMoments(oneDimMMAP, 2);
        assertTrue(moments[0] > 0, "1x1 MMAP should have positive first moment");
        
        double cv = M3aUtils.Companion.computeCoeffVar(oneDimMMAP);
        assertTrue(cv > 0, "1x1 MMAP should have positive coefficient of variation");
        
        double spectralGap = M3aUtils.Companion.computeSpectralGap(oneDimMMAP);
        assertTrue(spectralGap >= 0, "1x1 MMAP should have non-negative spectral gap");
    }
}
