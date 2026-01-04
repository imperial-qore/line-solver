/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.BeforeEach;

import java.util.Random;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Unit tests for MultivariateNormal distribution.
 */
class MultivariateNormalTest {

    private static final double TOL = 1e-6;
    private static final double CONVERGENCE_TOL = 0.1; // 10% tolerance for sampling
    private Random random;

    @BeforeEach
    void setUp() {
        random = new Random(42);
    }

    // =================== CONSTRUCTION TESTS ===================

    @Test
    void testConstruction2D() {
        double[] mu = {1.0, 2.0};
        Matrix Sigma = new Matrix(new double[][]{{1.0, 0.5}, {0.5, 1.0}});

        MultivariateNormal mvn = new MultivariateNormal(mu, Sigma);

        assertEquals(2, mvn.getDimension());
        assertEquals(1.0, mvn.getMeanVector().get(0, 0), TOL);
        assertEquals(2.0, mvn.getMeanVector().get(1, 0), TOL);
    }

    @Test
    void testConstruction3D() {
        double[] mu = {1.0, 2.0, 3.0};
        Matrix Sigma = new Matrix(new double[][]{
            {1.0, 0.5, 0.2},
            {0.5, 2.0, 0.3},
            {0.2, 0.3, 1.5}
        });

        MultivariateNormal mvn = new MultivariateNormal(mu, Sigma);

        assertEquals(3, mvn.getDimension());
        assertEquals(1.0, mvn.getMean(), TOL);  // First component
        assertEquals(1.0, mvn.getVar(), TOL);  // Variance of first component
    }

    @Test
    void testConstructionMatrixInput() {
        Matrix mu = new Matrix(2, 1);
        mu.set(0, 0, 1.0);
        mu.set(1, 0, 2.0);

        Matrix Sigma = new Matrix(new double[][]{{1.0, 0.5}, {0.5, 1.0}});

        MultivariateNormal mvn = new MultivariateNormal(mu, Sigma);

        assertEquals(2, mvn.getDimension());
    }

    @Test
    void testConstructionInvalidNonPD() {
        double[] mu = {0, 0};
        Matrix Sigma = new Matrix(new double[][]{{1.0, 2.0}, {2.0, 1.0}});  // Not PD

        assertThrows(IllegalArgumentException.class, () -> {
            new MultivariateNormal(mu, Sigma);
        });
    }

    @Test
    void testConstructionDimensionMismatch() {
        double[] mu = {0, 0};
        Matrix Sigma = new Matrix(new double[][]{{1.0, 0.5}, {0.5, 1.0}, {0.1, 0.1}});  // 3x2

        assertThrows(IllegalArgumentException.class, () -> {
            new MultivariateNormal(mu, Sigma);
        });
    }

    @Test
    void testConstructionEmptyMean() {
        double[] mu = {};
        Matrix Sigma = new Matrix(0, 0);

        assertThrows(IllegalArgumentException.class, () -> {
            new MultivariateNormal(mu, Sigma);
        });
    }

    // =================== ACCESSOR TESTS ===================

    @Test
    void testGetDimension() {
        double[] mu = {1, 2, 3};
        Matrix Sigma = Matrix.eye(3);

        MultivariateNormal mvn = new MultivariateNormal(mu, Sigma);

        assertEquals(3, mvn.getDimension());
    }

    @Test
    void testGetMeanVector() {
        double[] mu = {1.5, 2.5};
        Matrix Sigma = Matrix.eye(2);

        MultivariateNormal mvn = new MultivariateNormal(mu, Sigma);
        Matrix mean = mvn.getMeanVector();

        assertEquals(1.5, mean.get(0, 0), TOL);
        assertEquals(2.5, mean.get(1, 0), TOL);
    }

    @Test
    void testGetCovariance() {
        double[] mu = {0, 0};
        Matrix Sigma = new Matrix(new double[][]{{2.0, 0.5}, {0.5, 1.5}});

        MultivariateNormal mvn = new MultivariateNormal(mu, Sigma);
        Matrix cov = mvn.getCovariance();

        assertEquals(2.0, cov.get(0, 0), TOL);
        assertEquals(0.5, cov.get(0, 1), TOL);
        assertEquals(0.5, cov.get(1, 0), TOL);
        assertEquals(1.5, cov.get(1, 1), TOL);
    }

    @Test
    void testGetCorrelation() {
        double[] mu = {0, 0};
        Matrix Sigma = new Matrix(new double[][]{{1.0, 0.5}, {0.5, 2.0}});

        MultivariateNormal mvn = new MultivariateNormal(mu, Sigma);
        Matrix corr = mvn.getCorrelation();

        assertEquals(1.0, corr.get(0, 0), TOL);  // Diagonal elements are 1
        assertEquals(1.0, corr.get(1, 1), TOL);
        assertEquals(0.5 / sqrt(1.0 * 2.0), corr.get(0, 1), TOL);  // Off-diagonal correlation
        assertEquals(0.5 / sqrt(1.0 * 2.0), corr.get(1, 0), TOL);
    }

    // =================== STATISTICAL PROPERTIES TESTS ===================

    @Test
    void testGetMean() {
        double[] mu = {3.5, 4.5};
        Matrix Sigma = Matrix.eye(2);

        MultivariateNormal mvn = new MultivariateNormal(mu, Sigma);

        assertEquals(3.5, mvn.getMean(), TOL);  // First component
    }

    @Test
    void testGetVar() {
        double[] mu = {0, 0};
        Matrix Sigma = new Matrix(new double[][]{{2.0, 0.5}, {0.5, 1.5}});

        MultivariateNormal mvn = new MultivariateNormal(mu, Sigma);

        assertEquals(2.0, mvn.getVar(), TOL);  // Variance of first component
    }

    @Test
    void testGetSCV() {
        double[] mu = {1, 2};
        Matrix Sigma = Matrix.eye(2);

        MultivariateNormal mvn = new MultivariateNormal(mu, Sigma);

        assertTrue(Double.isNaN(mvn.getSCV()));
    }

    @Test
    void testGetSkewness() {
        double[] mu = {0, 0};
        Matrix Sigma = Matrix.eye(2);

        MultivariateNormal mvn = new MultivariateNormal(mu, Sigma);

        assertEquals(0.0, mvn.getSkewness(), TOL);
    }

    // =================== SAMPLING TESTS ===================

    @Test
    void testSamplingDimensions() {
        double[] mu = {0, 0};
        Matrix Sigma = Matrix.eye(2);
        MultivariateNormal mvn = new MultivariateNormal(mu, Sigma);

        Matrix samples = mvn.sampleMatrix(100, random);

        assertEquals(100, samples.getNumRows());
        assertEquals(2, samples.getNumCols());
    }

    @Test
    void testSamplingMeanConvergence() {
        double[] mu = {1.0, 2.0};
        Matrix Sigma = new Matrix(new double[][]{{1.0, 0.3}, {0.3, 2.0}});
        MultivariateNormal mvn = new MultivariateNormal(mu, Sigma);

        Matrix samples = mvn.sampleMatrix(10000, random);

        // Compute empirical mean
        double[] empMean = new double[2];
        for (int i = 0; i < 10000; i++) {
            empMean[0] += samples.get(i, 0);
            empMean[1] += samples.get(i, 1);
        }
        empMean[0] /= 10000;
        empMean[1] /= 10000;

        // Check convergence (10% tolerance)
        assertEquals(mu[0], empMean[0], CONVERGENCE_TOL * mu[0]);
        assertEquals(mu[1], empMean[1], CONVERGENCE_TOL * mu[1]);
    }

    @Test
    void testSamplingCovarianceConvergence() {
        double[] mu = {0, 0};
        Matrix Sigma = new Matrix(new double[][]{{1.0, 0.5}, {0.5, 2.0}});
        MultivariateNormal mvn = new MultivariateNormal(mu, Sigma);

        Matrix samples = mvn.sampleMatrix(10000, random);

        // Compute empirical covariance
        double[] means = new double[2];
        for (int i = 0; i < 10000; i++) {
            means[0] += samples.get(i, 0);
            means[1] += samples.get(i, 1);
        }
        means[0] /= 10000;
        means[1] /= 10000;

        double[][] empCov = new double[2][2];
        for (int i = 0; i < 10000; i++) {
            double x0 = samples.get(i, 0) - means[0];
            double x1 = samples.get(i, 1) - means[1];
            empCov[0][0] += x0 * x0;
            empCov[0][1] += x0 * x1;
            empCov[1][0] += x0 * x1;
            empCov[1][1] += x1 * x1;
        }
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                empCov[i][j] /= 10000;
            }
        }

        // Check covariance convergence
        assertEquals(1.0, empCov[0][0], 0.15);
        assertEquals(0.5, empCov[0][1], 0.15);
        assertEquals(2.0, empCov[1][1], 0.15);
    }

    // =================== PDF EVALUATION TESTS ===================

    @Test
    void testPDFEvaluationStandard() {
        double[] mu = {0, 0};
        Matrix Sigma = Matrix.eye(2);
        MultivariateNormal mvn = new MultivariateNormal(mu, Sigma);

        // At mean, PDF should be 1/(2π)
        double pdfAtMean = mvn.evalPDF(new double[]{0, 0});
        double expected = 1.0 / (2 * Math.PI);

        assertEquals(expected, pdfAtMean, TOL);
    }

    @Test
    void testPDFEvaluation1D() {
        double[] mu = {0};
        Matrix Sigma = new Matrix(new double[][]{{1.0}});
        MultivariateNormal mvn = new MultivariateNormal(mu, Sigma);

        // Compare with univariate normal: f(0) = 1/sqrt(2π)
        double pdfAt0 = mvn.evalPDF(new double[]{0});
        double expected = 1.0 / Math.sqrt(2 * Math.PI);

        assertEquals(expected, pdfAt0, TOL);
    }

    @Test
    void testPDFEvaluationNonZeroMean() {
        double[] mu = {1, 2};
        Matrix Sigma = Matrix.eye(2);
        MultivariateNormal mvn = new MultivariateNormal(mu, Sigma);

        // PDF at mean should be 1/(2π)
        double pdfAtMean = mvn.evalPDF(new double[]{1, 2});
        double expected = 1.0 / (2 * Math.PI);

        assertEquals(expected, pdfAtMean, TOL);
    }

    @Test
    void testPDFEvaluationMatrix() {
        double[] mu = {0, 0};
        Matrix Sigma = Matrix.eye(2);
        MultivariateNormal mvn = new MultivariateNormal(mu, Sigma);

        Matrix x = new Matrix(2, 1);
        x.set(0, 0, 0.0);
        x.set(1, 0, 0.0);

        double pdfAtMean = mvn.evalPDF(x);
        double expected = 1.0 / (2 * Math.PI);

        assertEquals(expected, pdfAtMean, TOL);
    }

    @Test
    void testPDFEvaluationDimensionMismatch() {
        double[] mu = {0, 0};
        Matrix Sigma = Matrix.eye(2);
        MultivariateNormal mvn = new MultivariateNormal(mu, Sigma);

        assertThrows(IllegalArgumentException.class, () -> {
            mvn.evalPDF(new double[]{0, 0, 0});  // Wrong dimension
        });
    }

    // =================== CDF/LST TESTS ===================

    @Test
    void testEvalCDFThrows() {
        double[] mu = {0, 0};
        Matrix Sigma = Matrix.eye(2);
        MultivariateNormal mvn = new MultivariateNormal(mu, Sigma);

        assertThrows(UnsupportedOperationException.class, () -> {
            mvn.evalCDF(0.5);
        });
    }

    @Test
    void testEvalLSTThrows() {
        double[] mu = {0, 0};
        Matrix Sigma = Matrix.eye(2);
        MultivariateNormal mvn = new MultivariateNormal(mu, Sigma);

        assertThrows(UnsupportedOperationException.class, () -> {
            mvn.evalLST(0.5);
        });
    }

    // =================== MARGINAL DISTRIBUTION TESTS ===================

    @Test
    void testGetMarginal2DFrom3D() {
        double[] mu = {1, 2, 3};
        Matrix Sigma = new Matrix(new double[][]{
            {1.0, 0.5, 0.2},
            {0.5, 2.0, 0.3},
            {0.2, 0.3, 1.5}
        });
        MultivariateNormal mvn = new MultivariateNormal(mu, Sigma);

        MultivariateNormal marginal = mvn.getMarginal(new int[]{0, 2});

        assertEquals(2, marginal.getDimension());
        assertEquals(1.0, marginal.getMeanVector().get(0, 0), TOL);
        assertEquals(3.0, marginal.getMeanVector().get(1, 0), TOL);
        assertEquals(1.0, marginal.getCovariance().get(0, 0), TOL);
        assertEquals(0.2, marginal.getCovariance().get(0, 1), TOL);
        assertEquals(1.5, marginal.getCovariance().get(1, 1), TOL);
    }

    @Test
    void testGetMarginalUniv() {
        double[] mu = {1, 2, 3};
        Matrix Sigma = new Matrix(new double[][]{
            {1.0, 0.5, 0.2},
            {0.5, 2.0, 0.3},
            {0.2, 0.3, 1.5}
        });
        MultivariateNormal mvn = new MultivariateNormal(mu, Sigma);

        Normal norm = mvn.getMarginalUniv(1);

        assertEquals(2.0, norm.getMean(), TOL);
        assertEquals(Math.sqrt(2.0), norm.getStd(), TOL);
    }

    @Test
    void testGetMarginalUnivFirst() {
        double[] mu = {1.5, 2.5};
        Matrix Sigma = new Matrix(new double[][]{{3.0, 1.0}, {1.0, 2.0}});
        MultivariateNormal mvn = new MultivariateNormal(mu, Sigma);

        Normal norm = mvn.getMarginalUniv(0);

        assertEquals(1.5, norm.getMean(), TOL);
        assertEquals(Math.sqrt(3.0), norm.getStd(), TOL);
    }

    @Test
    void testGetMarginalUnivIndexOutOfBounds() {
        double[] mu = {1, 2};
        Matrix Sigma = Matrix.eye(2);
        MultivariateNormal mvn = new MultivariateNormal(mu, Sigma);

        assertThrows(IllegalArgumentException.class, () -> {
            mvn.getMarginalUniv(5);
        });
    }

    // =================== PRIOR INTEGRATION TESTS ===================

    @Test
    void testPriorIntegrationTwoAlternatives() {
        MultivariateNormal mvn1 = new MultivariateNormal(
            new double[]{0, 0},
            Matrix.eye(2)
        );
        MultivariateNormal mvn2 = new MultivariateNormal(
            new double[]{2, 2},
            new Matrix(new double[][]{{0.5, 0.2}, {0.2, 0.5}})
        );

        java.util.List<Distribution> dists = new java.util.ArrayList<Distribution>();
        dists.add(mvn1);
        dists.add(mvn2);
        double[] probs = {0.6, 0.4};

        Prior prior = new Prior(dists, probs);

        assertEquals(2, prior.getNumAlternatives());
        // Prior.getMean() = sum of first component means weighted by probabilities
        double expectedMean = 0.6 * 0.0 + 0.4 * 2.0;
        assertEquals(expectedMean, prior.getMean(), TOL);
    }

    // =================== UTILITY TESTS ===================

    @Test
    void testGetProcess() {
        double[] mu = {1, 2};
        Matrix Sigma = new Matrix(new double[][]{{1.0, 0.5}, {0.5, 2.0}});
        MultivariateNormal mvn = new MultivariateNormal(mu, Sigma);

        jline.util.matrix.MatrixCell cell = mvn.getProcess();

        // Should contain mu and Sigma
        assertNotNull(cell);
    }

    @Test
    void testToString() {
        double[] mu = {1, 2, 3};
        Matrix Sigma = Matrix.eye(3);
        MultivariateNormal mvn = new MultivariateNormal(mu, Sigma);

        String str = mvn.toString();
        assertTrue(str.contains("MultivariateNormal"));
        assertTrue(str.contains("d=3"));
    }

    // =================== HELPER METHODS ===================

    private double sqrt(double x) {
        return Math.sqrt(x);
    }
}
