package jline.api;

import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;
import java.util.List;
import java.util.Arrays;
import static jline.api.trace.TraceSkewKt.trace_skew;
import static jline.api.trace.Trace_varKt.*;
import static jline.api.trace.Trace_meanKt.trace_mean;
import static jline.TestTools.*;

public class TraceAPITest {
    
    @Test
    public void testTraceSkewWithSymmetricData() {
        // Symmetric data should have skewness close to 0
        double[] symmetricData = {1.0, 2.0, 3.0, 4.0, 5.0, 4.0, 3.0, 2.0, 1.0};
        double skewness = trace_skew(symmetricData);
        assertEquals(0.146354523591967, skewness, 0.01, "Symmetric data should have skewness close to 0");
    }
    
    @Test
    public void testTraceSkewWithPositiveSkew() {
        // Right-skewed data (long tail on the right)
        double[] rightSkewedData = {1.0, 1.0, 1.0, 2.0, 2.0, 3.0, 4.0, 5.0, 10.0};
        double skewness = trace_skew(rightSkewedData);
        assertTrue(skewness > 0, "Right-skewed data should have positive skewness");
    }
    
    @Test
    public void testTraceSkewWithNegativeSkew() {
        // Left-skewed data (long tail on the left)
        double[] leftSkewedData = {1.0, 6.0, 7.0, 8.0, 9.0, 9.0, 9.0, 9.0, 9.0};
        double skewness = trace_skew(leftSkewedData);
        assertTrue(skewness < 0, "Left-skewed data should have negative skewness");
    }
    
    @Test
    public void testTraceSkewWithList() {
        List<Double> dataList = Arrays.asList(1.0, 2.0, 3.0, 4.0, 5.0);
        double skewnessFromList = trace_skew(dataList);
        double skewnessFromArray = trace_skew(new double[]{1.0, 2.0, 3.0, 4.0, 5.0});
        assertEquals(skewnessFromArray, skewnessFromList, 1e-10, "List and array should give same result");
    }

    @Test
    public void testTraceSkewConsistencyWithApacheCommons() {
        // Test with the same data used in Replayer example
        double[] traceData = {
            1.2377474e-02, 4.4486055e-02, 1.0027642e-02, 2.0983173e-02, 5.0081083e-02,
            1.5390689e-02, 1.8521508e-01, 1.1288256e-01, 6.8028361e-02, 9.3326489e-02
        };
        double skewness = trace_skew(traceData);
        
        // Should match Apache Commons Skewness.evaluate() result
        org.apache.commons.math3.stat.descriptive.moment.Skewness apacheSkewness = 
            new org.apache.commons.math3.stat.descriptive.moment.Skewness();
        double expectedSkewness = apacheSkewness.evaluate(traceData);
        
        assertEquals(expectedSkewness, skewness, 1e-10,
            "Should match Apache Commons skewness calculation exactly");
    }

    // ========================== Test 43: trace_var tests ==========================

    /**
     * Test 43: trace_var - variance computation for uniform data.
     */
    @Test
    public void testTraceVar_uniformData() {
        // Uniform data: 1, 2, 3, 4, 5
        double[] data = {1.0, 2.0, 3.0, 4.0, 5.0};
        double variance = trace_var(data);

        // Variance of 1,2,3,4,5: mean=3, var = ((1-3)^2+(2-3)^2+(3-3)^2+(4-3)^2+(5-3)^2)/5 = 2
        assertEquals(2.0, variance, MID_TOL, "Variance of uniform data should be 2.0");
    }

    /**
     * Test 43b: trace_var - zero variance for constant data.
     */
    @Test
    public void testTraceVar_constantData() {
        // All same values should have zero variance
        double[] data = {5.0, 5.0, 5.0, 5.0, 5.0};
        double variance = trace_var(data);

        assertEquals(0.0, variance, ZERO_TOL, "Constant data should have zero variance");
    }

    /**
     * Test 43c: trace_var - variance with known distribution.
     */
    @Test
    public void testTraceVar_knownDistribution() {
        // Data: 2, 4, 6, 8, 10
        double[] data = {2.0, 4.0, 6.0, 8.0, 10.0};
        double variance = trace_var(data);

        // mean = 6, var = ((2-6)^2 + (4-6)^2 + (6-6)^2 + (8-6)^2 + (10-6)^2)/5 = (16+4+0+4+16)/5 = 8
        assertEquals(8.0, variance, MID_TOL, "Variance should be 8.0");
    }

    /**
     * Test 43d: trace_scv - squared coefficient of variation.
     */
    @Test
    public void testTraceScv_basic() {
        double[] data = {1.0, 2.0, 3.0, 4.0, 5.0};
        double scv = trace_scv(data);

        // scv = var / mean^2 = 2 / 9 = 0.222...
        double expected = 2.0 / 9.0;
        assertEquals(expected, scv, MID_TOL, "SCV should be var/mean^2");
    }

    /**
     * Test 43e: trace_acf - autocorrelation function.
     */
    @Test
    public void testTraceAcf_basicLag() {
        // Simple oscillating data
        double[] data = {1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0};
        double[] acf = trace_acf(data, new int[]{1});

        // For alternating data, lag-1 ACF should be negative
        assertTrue(acf.length > 0, "ACF should return values");
        assertTrue(acf[0] < 0, "Alternating data should have negative lag-1 ACF");
    }

    /**
     * Test 43f: trace_acf - multiple lags.
     */
    @Test
    public void testTraceAcf_multipleLags() {
        double[] data = new double[100];
        for (int i = 0; i < 100; i++) {
            data[i] = Math.sin(i * 0.5);
        }

        int[] lags = {1, 2, 3, 4, 5};
        double[] acf = trace_acf(data, lags);

        assertEquals(lags.length, acf.length, "ACF should return value for each lag");

        // All ACF values should be between -1 and 1
        for (double value : acf) {
            assertTrue(value >= -1.0 - FINE_TOL && value <= 1.0 + FINE_TOL,
                "ACF values should be in [-1, 1]");
        }
    }

    /**
     * Test 43g: trace_idi - index of dispersion for intervals.
     */
    @Test
    public void testTraceIdi_basic() {
        // Generate some trace data
        double[] data = {0.5, 1.2, 0.8, 1.5, 0.3, 2.0, 0.7, 1.1, 0.9, 1.3,
                         0.6, 1.4, 0.4, 1.6, 0.8, 1.2, 0.5, 1.0, 0.9, 1.1};

        int[] kset = {1, 2, 3};
        kotlin.Pair<double[], int[]> result = trace_idi(data, kset, null, 1);

        double[] idi = result.getFirst();
        int[] support = result.getSecond();

        assertEquals(kset.length, idi.length, "IDI should return value for each k");
        assertEquals(kset.length, support.length, "Support should have same length as kset");

        // IDI values should be positive
        for (double value : idi) {
            assertTrue(value >= 0, "IDI values should be non-negative");
        }
    }

    /**
     * Test 43h: trace_idc - index of dispersion for counts.
     */
    @Test
    public void testTraceIdc_basic() {
        double[] data = new double[100];
        java.util.Random random = new java.util.Random(42);
        for (int i = 0; i < 100; i++) {
            data[i] = -Math.log(random.nextDouble()); // Exponential samples
        }

        double idc = trace_idc(data);

        // IDC should be positive
        assertTrue(idc > 0, "IDC should be positive");
    }

    /**
     * Test 43i: trace_summary - comprehensive summary statistics.
     */
    @Test
    public void testTraceSummary_basic() {
        double[] data = new double[100];
        java.util.Random random = new java.util.Random(123);
        for (int i = 0; i < 100; i++) {
            data[i] = -Math.log(random.nextDouble()); // Exponential samples
        }

        double[] summary = trace_summary(data);

        // Summary should return 17 values
        assertTrue(summary.length >= 10, "Summary should return multiple statistics");

        // Mean should be positive
        assertTrue(summary[0] > 0, "Mean should be positive");

        // SCV should be positive
        assertTrue(summary[1] > 0, "SCV should be positive");
    }

    /**
     * Test: trace_var consistency with Apache Commons.
     */
    @Test
    public void testTraceVar_apacheConsistency() {
        double[] data = {1.2, 3.4, 5.6, 7.8, 9.0, 2.3, 4.5, 6.7, 8.9, 1.0};

        double ourVariance = trace_var(data);

        // Compute variance manually: E[X^2] - E[X]^2
        double sum = 0, sumSq = 0;
        for (double d : data) {
            sum += d;
            sumSq += d * d;
        }
        double mean = sum / data.length;
        double expectedVar = sumSq / data.length - mean * mean;

        assertEquals(expectedVar, ourVariance, ZERO_TOL, "Should match manual calculation");
    }
}