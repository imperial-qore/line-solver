/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.util.GMMUtils;
import jline.util.Maths;
import jline.util.RandomManager;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.util.Random;

import static org.junit.jupiter.api.Assertions.*;

public class GMMTest {
    
    private static final double TOLERANCE = 1e-6;
    
    @BeforeAll
    public static void setUp() {
        // Ensure MATLAB-compatible random number generation with fixed seed
        Maths.setRandomNumbersMatlab(true);
        RandomManager.setMasterSeed(23000);
    }
    
    @Test
    public void testGMMCreation() {
        double[] weights = {0.3, 0.7};
        double[] means = {0.0, 5.0};
        double[] sigmas = {1.0, 2.0};
        
        GMM gmm = new GMM(weights, means, sigmas);
        
        assertEquals(2, gmm.getNumberOfComponents());
        
        // Test mean calculation
        double expectedMean = 0.3 * 0.0 + 0.7 * 5.0; // 3.5
        assertEquals(expectedMean, gmm.getMean(), TOLERANCE);
    }
    
    @Test
    public void testGMMVariance() {
        double[] weights = {0.3, 0.7};
        double[] means = {0.0, 5.0};
        double[] sigmas = {1.0, 2.0};
        
        GMM gmm = new GMM(weights, means, sigmas);
        
        // Calculate expected variance
        double mean = gmm.getMean();
        double e2 = 0.3 * (1.0 + 0.0) + 0.7 * (4.0 + 25.0); // E[X^2]
        double expectedVar = e2 - mean * mean;
        
        assertEquals(expectedVar, gmm.getVar(), TOLERANCE);
    }
    
    @Test
    public void testGMMConvolve() {
        double[] weights1 = {1.0};
        double[] means1 = {2.0};
        double[] sigmas1 = {1.0};
        GMM gmm1 = new GMM(weights1, means1, sigmas1);
        
        double[] weights2 = {1.0};
        double[] means2 = {3.0};
        double[] sigmas2 = {2.0};
        GMM gmm2 = new GMM(weights2, means2, sigmas2);
        
        GMM result = gmm1.convolve(gmm2);
        
        // Mean of sum should be sum of means
        assertEquals(5.0, result.getMean(), TOLERANCE);
        
        // Variance of sum should be sum of variances
        assertEquals(5.0, result.getVar(), TOLERANCE); // 1^2 + 2^2 = 5
    }
    
    @Test
    public void testGMMMixture() {
        double[] weights1 = {1.0};
        double[] means1 = {0.0};
        double[] sigmas1 = {1.0};
        GMM gmm1 = new GMM(weights1, means1, sigmas1);
        
        double[] weights2 = {1.0};
        double[] means2 = {10.0};
        double[] sigmas2 = {1.0};
        GMM gmm2 = new GMM(weights2, means2, sigmas2);
        
        GMM mixture = GMM.mixture(gmm1, gmm2, 0.5, 0.5);
        
        assertEquals(2, mixture.getNumberOfComponents());
        assertEquals(5.0, mixture.getMean(), TOLERANCE); // 0.5 * 0 + 0.5 * 10
    }
    
    @Test
    public void testGMMMerge() {
        double[] weights = {0.3, 0.3, 0.4};
        double[] means = {0.0, 1.0, 5.0};
        double[] sigmas = {1.0, 1.0, 2.0};
        
        GMM gmm = new GMM(weights, means, sigmas);
        GMM merged = gmm.merge(0, 1);
        
        assertEquals(2, merged.getNumberOfComponents());
        
        // The merged component should have weight 0.6
        double[] mergedWeights = merged.getWeights();
        assertEquals(0.6, mergedWeights[0], TOLERANCE);
        assertEquals(0.4, mergedWeights[1], TOLERANCE);
    }
    
    @Test
    public void testGMMAutoMerge() {
        double[] weights = {0.2, 0.2, 0.3, 0.3};
        double[] means = {0.0, 1.0, 5.0, 6.0};
        double[] sigmas = {1.0, 1.0, 1.0, 1.0};
        
        GMM gmm = new GMM(weights, means, sigmas);
        GMM merged = gmm.autoMerge(2);
        
        assertEquals(2, merged.getNumberOfComponents());
    }
    
    @Test
    public void testMatrixConversion() {
        Matrix sgmm = new Matrix(2, 3);
        sgmm.set(0, 0, 0.3);  // weight
        sgmm.set(0, 1, 0.0);  // mean
        sgmm.set(0, 2, 1.0);  // variance
        sgmm.set(1, 0, 0.7);  // weight
        sgmm.set(1, 1, 5.0);  // mean
        sgmm.set(1, 2, 4.0);  // variance
        
        GMM gmm = GMM.fromMatrix(sgmm);
        
        assertEquals(2, gmm.getNumberOfComponents());
        assertEquals(3.5, gmm.getMean(), TOLERANCE);
        
        // Convert back to matrix
        Matrix sgmm2 = gmm.toMatrix();
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 3; j++) {
                assertEquals(sgmm.get(i, j), sgmm2.get(i, j), TOLERANCE);
            }
        }
    }
    
    @Test
    public void testGMMUtils() {
        double[] weights = {0.5, 0.5};
        double[] means = {-2.0, 2.0};
        double[] sigmas = {1.0, 1.0};
        
        GMM gmm = new GMM(weights, means, sigmas);
        
        // Test median (should be close to 0 for symmetric mixture)
        double median = GMMUtils.gmmMedian(gmm);
        assertEquals(0.0, median, 0.1);
        
        // Test that min < mean < max
        double min = GMMUtils.gmmMin(gmm);
        double max = GMMUtils.gmmMax(gmm);
        assertTrue(min < gmm.getMean());
        assertTrue(gmm.getMean() < max);
    }
    
    @Test
    public void testSampling() {
        double[] weights = {0.3, 0.7};
        double[] means = {0.0, 10.0};
        double[] sigmas = {1.0, 2.0};
        
        GMM gmm = new GMM(weights, means, sigmas);
        
        // Use fixed seed for reproducible test results
        Random random = new Random(42);
        int n = 10000;
        double[] samples = gmm.sample(n, random);
        
        // Check that sample mean is close to theoretical mean
        double sampleMean = 0;
        for (double s : samples) {
            sampleMean += s;
        }
        sampleMean /= n;
        
        assertEquals(gmm.getMean(), sampleMean, 0.1);
    }
}