/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.basic;

import jline.lang.processes.MultivariateNormal;
import jline.lang.processes.Normal;
import jline.lang.processes.Prior;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.List;

/**
 * Example demonstrating MultivariateNormal distribution usage.
 *
 * This example shows:
 * 1. Creating a multivariate normal distribution
 * 2. Accessing statistical properties
 * 3. Sampling from the distribution
 * 4. Evaluating the PDF
 * 5. Extracting marginal distributions
 * 6. Using MultivariateNormal within a Prior for mixture models
 */
public class MultivariateNormalExample {

    public static void main(String[] args) {
        System.out.println("=== MultivariateNormal Distribution Example ===\n");

        // Example 1: Basic 2D normal distribution
        example1_Basic2DDistribution();

        // Example 2: Sampling and statistical properties
        example2_SamplingAndStatistics();

        // Example 3: PDF evaluation
        example3_PDFEvaluation();

        // Example 4: Marginal distributions
        example4_MarginalDistributions();

        // Example 5: Mixture model using Prior
        example5_MixtureModelWithPrior();
    }

    /**
     * Example 1: Create and inspect a basic 2D multivariate normal distribution.
     */
    static void example1_Basic2DDistribution() {
        System.out.println("Example 1: Basic 2D Distribution");
        System.out.println("---------------------------------");

        // Create mean vector
        double[] mu = {1.0, 2.0};

        // Create covariance matrix
        Matrix Sigma = new Matrix(new double[][]{
            {1.0, 0.5},
            {0.5, 1.0}
        });

        // Create multivariate normal
        MultivariateNormal mvn = new MultivariateNormal(mu, Sigma);

        // Inspect properties
        System.out.printf("Dimensionality: %d\n", mvn.getDimension());
        System.out.printf("Mean: [%.2f, %.2f]\n", mvn.getMeanVector().get(0, 0), mvn.getMeanVector().get(1, 0));
        System.out.printf("Variance (dim 1): %.4f\n", mvn.getVar());
        System.out.printf("Skewness: %.4f (always 0 for MVN)\n\n", mvn.getSkewness());
    }

    /**
     * Example 2: Generate samples and verify convergence of sample statistics.
     */
    static void example2_SamplingAndStatistics() {
        System.out.println("Example 2: Sampling and Statistics");
        System.out.println("-----------------------------------");

        double[] mu = {0.0, 0.0};
        Matrix Sigma = new Matrix(new double[][]{{2.0, 0.8}, {0.8, 1.0}});
        MultivariateNormal mvn = new MultivariateNormal(mu, Sigma);

        // Generate large sample
        int n = 10000;
        Matrix samples = mvn.sampleMatrix(n, new java.util.Random(42));

        // Compute empirical statistics
        double[] empMean = new double[2];
        double[][] empCov = new double[2][2];

        for (int i = 0; i < n; i++) {
            empMean[0] += samples.get(i, 0);
            empMean[1] += samples.get(i, 1);
        }
        empMean[0] /= n;
        empMean[1] /= n;

        for (int i = 0; i < n; i++) {
            double x0 = samples.get(i, 0) - empMean[0];
            double x1 = samples.get(i, 1) - empMean[1];
            empCov[0][0] += x0 * x0;
            empCov[0][1] += x0 * x1;
            empCov[1][0] += x0 * x1;
            empCov[1][1] += x1 * x1;
        }
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                empCov[i][j] /= n;
            }
        }

        System.out.printf("True mean: [%.4f, %.4f]\n", mu[0], mu[1]);
        System.out.printf("Sample mean: [%.4f, %.4f]\n\n", empMean[0], empMean[1]);

        System.out.println("True covariance:");
        System.out.printf("  [%.4f  %.4f]\n", Sigma.get(0, 0), Sigma.get(0, 1));
        System.out.printf("  [%.4f  %.4f]\n\n", Sigma.get(1, 0), Sigma.get(1, 1));

        System.out.println("Sample covariance:");
        System.out.printf("  [%.4f  %.4f]\n", empCov[0][0], empCov[0][1]);
        System.out.printf("  [%.4f  %.4f]\n\n", empCov[1][0], empCov[1][1]);
    }

    /**
     * Example 3: Evaluate the probability density function.
     */
    static void example3_PDFEvaluation() {
        System.out.println("Example 3: PDF Evaluation");
        System.out.println("------------------------");

        // Standard 2D normal
        MultivariateNormal mvn = new MultivariateNormal(
            new double[]{0.0, 0.0},
            new Matrix(new double[][]{{1.0, 0.0}, {0.0, 1.0}})
        );

        // At mean: f(0,0) = 1/(2π) ≈ 0.1592
        double pdf_at_mean = mvn.evalPDF(new double[]{0.0, 0.0});
        System.out.printf("PDF at mean (0,0): %.6f (expected: %.6f)\n", pdf_at_mean, 1.0 / (2 * Math.PI));

        // Away from mean
        double pdf_away = mvn.evalPDF(new double[]{1.0, 1.0});
        System.out.printf("PDF at (1,1): %.6f\n", pdf_away);

        // Far from mean
        double pdf_far = mvn.evalPDF(new double[]{3.0, 3.0});
        System.out.printf("PDF at (3,3): %.6f\n\n", pdf_far);
    }

    /**
     * Example 4: Extract and work with marginal distributions.
     */
    static void example4_MarginalDistributions() {
        System.out.println("Example 4: Marginal Distributions");
        System.out.println("--------------------------------");

        // Create 3D distribution
        double[] mu = {1.0, 2.0, 3.0};
        Matrix Sigma = new Matrix(new double[][]{
            {1.0, 0.5, 0.2},
            {0.5, 2.0, 0.3},
            {0.2, 0.3, 1.5}
        });
        MultivariateNormal mvn3d = new MultivariateNormal(mu, Sigma);

        System.out.printf("3D Distribution: dimension = %d\n\n", mvn3d.getDimension());

        // Extract 2D marginal (dimensions 0 and 2)
        MultivariateNormal mvn2d = mvn3d.getMarginal(new int[]{0, 2});
        System.out.printf("2D Marginal (dims 0,2): dimension = %d\n", mvn2d.getDimension());
        System.out.printf("Marginal mean: [%.2f, %.2f]\n\n",
            mvn2d.getMeanVector().get(0, 0), mvn2d.getMeanVector().get(1, 0));

        // Extract univariate marginals
        Normal norm0 = mvn3d.getMarginalUniv(0);
        Normal norm1 = mvn3d.getMarginalUniv(1);
        Normal norm2 = mvn3d.getMarginalUniv(2);

        System.out.println("Univariate marginals:");
        System.out.printf("  Dim 0: N(%.2f, %.4f)\n", norm0.getMean(), norm0.getVar());
        System.out.printf("  Dim 1: N(%.2f, %.4f)\n", norm1.getMean(), norm1.getVar());
        System.out.printf("  Dim 2: N(%.2f, %.4f)\n\n", norm2.getMean(), norm2.getVar());
    }

    /**
     * Example 5: Use MultivariateNormal within a Prior for mixture models.
     */
    static void example5_MixtureModelWithPrior() {
        System.out.println("Example 5: Mixture Model with Prior");
        System.out.println("------------------------------------");

        // Create two alternative MVN distributions
        MultivariateNormal mvn1 = new MultivariateNormal(
            new double[]{0.0, 0.0},
            new Matrix(new double[][]{{1.0, 0.3}, {0.3, 1.0}})
        );

        MultivariateNormal mvn2 = new MultivariateNormal(
            new double[]{2.0, 2.0},
            new Matrix(new double[][]{{0.5, 0.1}, {0.1, 0.5}})
        );

        // Create Prior with 60% probability for mvn1, 40% for mvn2
        List<jline.lang.processes.Distribution> dists = new ArrayList<jline.lang.processes.Distribution>();
        dists.add(mvn1);
        dists.add(mvn2);
        double[] probs = {0.6, 0.4};

        Prior prior = new Prior(dists, probs);

        System.out.printf("Prior mixture model:\n");
        System.out.printf("  Alternative 1 (60%%): MVN(μ=[0,0], Σ=[1,0.3;0.3,1])\n");
        System.out.printf("  Alternative 2 (40%%): MVN(μ=[2,2], Σ=[0.5,0.1;0.1,0.5])\n\n");

        // Properties of the mixture
        System.out.printf("Mixture properties:\n");
        System.out.printf("  Number of alternatives: %d\n", prior.getNumAlternatives());
        System.out.printf("  Mean (first component): %.4f\n", prior.getMean());
        System.out.printf("  (Expected: 0.6*0 + 0.4*2 = 0.8)\n\n");

        // Generate mixture samples
        System.out.println("Generating 10 samples from the mixture:");
        double[] samples = prior.sample(10);
        for (int i = 0; i < 10; i++) {
            System.out.printf("  Sample %d: %.4f\n", i + 1, samples[i]);
        }
    }
}
