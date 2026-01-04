/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import static jline.GlobalConstants.Inf;
import static jline.GlobalConstants.NegInf;

import jline.GlobalConstants;
import jline.util.Pair;
import jline.util.RandomManager;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.io.Serializable;
import java.util.Random;

import static org.apache.commons.math3.util.FastMath.*;
import org.apache.commons.math3.linear.CholeskyDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularMatrixException;

/**
 * A multivariate normal (Gaussian) distribution.
 *
 * Represents a d-dimensional normal distribution with mean vector mu and
 * covariance matrix Sigma. The distribution can be used standalone or within
 * a Prior for mixture models.
 */
public class MultivariateNormal extends ContinuousDistribution implements Serializable {

    private int dimension;
    private Matrix L; // Cholesky decomposition (cached)
    private double detSigma; // Determinant of Sigma (cached)
    private Matrix invSigma; // Inverse of Sigma (cached)

    /**
     * Creates a multivariate normal distribution.
     *
     * @param mu    Mean vector (d-dimensional)
     * @param Sigma Covariance matrix (d x d, must be positive definite)
     * @throws IllegalArgumentException if Sigma is not positive definite or dimensions don't match
     */
    public MultivariateNormal(double[] mu, Matrix Sigma) {
        super("MultivariateNormal", 2, new Pair<Double, Double>(NegInf, Inf));
        // Convert array to column vector Matrix
        Matrix muMatrix = new Matrix(mu.length, 1);
        for (int i = 0; i < mu.length; i++) {
            muMatrix.set(i, 0, mu[i]);
        }
        initializeFields(muMatrix, Sigma);
    }

    /**
     * Creates a multivariate normal distribution.
     *
     * @param mu    Mean vector (d x 1 Matrix)
     * @param Sigma Covariance matrix (d x d, must be positive definite)
     * @throws IllegalArgumentException if Sigma is not positive definite or dimensions don't match
     */
    public MultivariateNormal(Matrix mu, Matrix Sigma) {
        super("MultivariateNormal", 2, new Pair<Double, Double>(NegInf, Inf));
        initializeFields(mu, Sigma);
    }

    /**
     * Initializes the fields for the distribution.
     */
    private void initializeFields(Matrix mu, Matrix Sigma) {
        // Validate inputs
        if (mu.getNumRows() == 0) {
            throw new IllegalArgumentException("Mean vector cannot be empty");
        }
        if (mu.getNumCols() > 1 && mu.getNumRows() == 1) {
            // Row vector - transpose to column vector
            mu = mu.transpose();
        }
        if (mu.getNumCols() != 1) {
            throw new IllegalArgumentException("Mean must be a column vector (d x 1)");
        }

        int d = mu.getNumRows();
        if (Sigma.getNumRows() != d || Sigma.getNumCols() != d) {
            throw new IllegalArgumentException(
                String.format("Covariance matrix dimensions (%d x %d) don't match mean length (%d)",
                    Sigma.getNumRows(), Sigma.getNumCols(), d)
            );
        }

        // Validate positive definiteness and compute Cholesky decomposition
        validateAndDecompose(Sigma);

        this.dimension = d;
        this.setParam(1, "mu", mu);
        this.setParam(2, "Sigma", Sigma);

        // Compute mean of first component for getMean() compatibility
        this.mean = mu.get(0, 0);
        this.immediate = false;
    }

    /**
     * Validates that Sigma is positive definite and computes Cholesky decomposition.
     */
    private void validateAndDecompose(Matrix Sigma) {
        try {
            RealMatrix sigmaMatrix = convertToRealMatrix(Sigma);
            CholeskyDecomposition cholesky = new CholeskyDecomposition(sigmaMatrix);

            // Store Cholesky factor L
            RealMatrix L_real = cholesky.getL();
            this.L = convertFromRealMatrix(L_real);

            // Compute and cache determinant
            this.detSigma = cholesky.getDeterminant();

            // Compute and cache inverse
            try {
                RealMatrix invSigma_real = MatrixUtils.inverse(sigmaMatrix);
                this.invSigma = convertFromRealMatrix(invSigma_real);
            } catch (SingularMatrixException e) {
                throw new IllegalArgumentException("Covariance matrix is singular");
            }
        } catch (SingularMatrixException e) {
            throw new IllegalArgumentException("Covariance matrix is not positive definite: " + e.getMessage());
        }
    }

    /**
     * Converts LINE Matrix to Apache Commons RealMatrix.
     */
    private RealMatrix convertToRealMatrix(Matrix m) {
        double[][] data = new double[m.getNumRows()][m.getNumCols()];
        for (int i = 0; i < m.getNumRows(); i++) {
            for (int j = 0; j < m.getNumCols(); j++) {
                data[i][j] = m.get(i, j);
            }
        }
        return MatrixUtils.createRealMatrix(data);
    }

    /**
     * Converts Apache Commons RealMatrix to LINE Matrix.
     */
    private Matrix convertFromRealMatrix(RealMatrix m) {
        double[][] data = m.getData();
        return new Matrix(data);
    }

    // =================== ACCESSORS ===================

    /**
     * Gets the dimensionality of this distribution.
     *
     * @return the number of dimensions d
     */
    public int getDimension() {
        return dimension;
    }

    /**
     * Gets the mean vector.
     *
     * @return d x 1 Matrix containing the mean
     */
    public Matrix getMeanVector() {
        return (Matrix) this.getParam(1).getValue();
    }

    /**
     * Gets the covariance matrix.
     *
     * @return d x d Matrix containing the covariance
     */
    public Matrix getCovariance() {
        return (Matrix) this.getParam(2).getValue();
    }

    /**
     * Gets the correlation matrix from the covariance matrix.
     *
     * @return d x d Matrix containing the correlation coefficients
     */
    public Matrix getCorrelation() {
        Matrix Sigma = getCovariance();
        Matrix R = Sigma.copy();

        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                double std_i = sqrt(Sigma.get(i, i));
                double std_j = sqrt(Sigma.get(j, j));
                if (std_i > GlobalConstants.FineTol && std_j > GlobalConstants.FineTol) {
                    R.set(i, j, Sigma.get(i, j) / (std_i * std_j));
                } else {
                    R.set(i, j, (i == j) ? 1.0 : 0.0);
                }
            }
        }

        return R;
    }

    // =================== DISTRIBUTION METHODS ===================

    /**
     * Gets the mean of the first component (for Prior compatibility).
     * For multivariate use, prefer getMeanVector().
     *
     * @return the mean of the first dimension
     */
    @Override
    public double getMean() {
        return getMeanVector().get(0, 0);
    }

    /**
     * Gets the variance of the first component.
     *
     * @return the variance of the first dimension
     */
    @Override
    public double getVar() {
        return getCovariance().get(0, 0);
    }

    /**
     * Gets the squared coefficient of variation (not meaningful for multivariate).
     *
     * @return Double.NaN (not defined for multivariate distributions)
     */
    @Override
    public double getSCV() {
        return Double.NaN;
    }

    /**
     * Gets the skewness.
     *
     * @return 0.0 (multivariate normal is symmetric)
     */
    @Override
    public double getSkewness() {
        return 0.0;
    }

    /**
     * Generates random samples from the distribution.
     * Uses the Cholesky decomposition: X = mu + L*Z where Z ~ N(0, I).
     *
     * @param n      the number of samples to generate
     * @param random the random number generator to use
     * @return flattened array of n*d samples (each dimension is a column)
     */
    @Override
    public double[] sample(int n, Random random) {
        Matrix samples = sampleMatrix(n, random);
        double[] result = new double[n * dimension];

        int idx = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < dimension; j++) {
                result[idx++] = samples.get(i, j);
            }
        }

        return result;
    }

    /**
     * Generates random samples from the distribution.
     * Returns samples as an n x d matrix.
     *
     * @param n      the number of samples to generate
     * @param random the random number generator to use
     * @return n x d Matrix of samples
     */
    public Matrix sampleMatrix(int n, Random random) {
        Matrix mu = getMeanVector();
        Matrix samples = new Matrix(n, dimension);

        for (int i = 0; i < n; i++) {
            // Generate standard normal vector Z (d x 1)
            Matrix Z = new Matrix(dimension, 1);
            for (int j = 0; j < dimension; j++) {
                Z.set(j, 0, random.nextGaussian());
            }

            // Compute X = mu + L*Z
            Matrix sample = mu.add(L.mult(Z));

            // Store as row in samples matrix
            for (int j = 0; j < dimension; j++) {
                samples.set(i, j, sample.get(j, 0));
            }
        }

        return samples;
    }

    /**
     * Evaluates the probability density function.
     *
     * @param x the point at which to evaluate (d x 1 Matrix)
     * @return the PDF value at x
     */
    public double evalPDF(Matrix x) {
        // Validate input dimension
        if (x.getNumRows() != dimension || x.getNumCols() != 1) {
            throw new IllegalArgumentException(
                String.format("Point dimension (%d x %d) doesn't match distribution dimension (%d x 1)",
                    x.getNumRows(), x.getNumCols(), dimension)
            );
        }

        // Compute x - mu
        Matrix diff = x.add(-1.0, getMeanVector());

        // Compute quadratic form: (x-mu)' * Sigma^-1 * (x-mu)
        Matrix temp = invSigma.mult(diff);
        double quadForm = diff.transpose().mult(temp).get(0, 0);

        // Compute normalization constant: 1 / sqrt((2π)^d * |Sigma|)
        double normConst = 1.0 / sqrt(pow(2 * PI, dimension) * abs(detSigma));

        // Return PDF value
        return normConst * exp(-0.5 * quadForm);
    }

    /**
     * Evaluates the probability density function.
     *
     * @param x the point at which to evaluate (d-dimensional array)
     * @return the PDF value at x
     */
    public double evalPDF(double[] x) {
        if (x.length != dimension) {
            throw new IllegalArgumentException(
                String.format("Point dimension (%d) doesn't match distribution dimension (%d)",
                    x.length, dimension)
            );
        }

        Matrix xMatrix = new Matrix(dimension, 1);
        for (int i = 0; i < dimension; i++) {
            xMatrix.set(i, 0, x[i]);
        }

        return evalPDF(xMatrix);
    }

    /**
     * CDF is not well-defined for multivariate distributions.
     *
     * @param t ignored
     * @return never returns
     * @throws UnsupportedOperationException always
     */
    @Override
    public double evalCDF(double t) {
        throw new UnsupportedOperationException(
            "CDF is not defined for multivariate distributions");
    }

    /**
     * LST is not well-defined for multivariate distributions.
     *
     * @param s ignored
     * @return never returns
     * @throws UnsupportedOperationException always
     */
    @Override
    public double evalLST(double s) {
        throw new UnsupportedOperationException(
            "Laplace-Stieltjes Transform (LST) is not defined for multivariate distributions");
    }

    // =================== MARGINAL DISTRIBUTIONS ===================

    /**
     * Extracts a marginal distribution for a subset of dimensions.
     *
     * @param indices array of dimension indices to keep (0-based)
     * @return MultivariateNormal distribution for the marginal
     */
    public MultivariateNormal getMarginal(int[] indices) {
        Matrix mu = getMeanVector();
        Matrix Sigma = getCovariance();

        // Extract marginal mean
        Matrix mu_marg = new Matrix(indices.length, 1);
        for (int i = 0; i < indices.length; i++) {
            mu_marg.set(i, 0, mu.get(indices[i], 0));
        }

        // Extract marginal covariance
        Matrix Sigma_marg = new Matrix(indices.length, indices.length);
        for (int i = 0; i < indices.length; i++) {
            for (int j = 0; j < indices.length; j++) {
                Sigma_marg.set(i, j, Sigma.get(indices[i], indices[j]));
            }
        }

        return new MultivariateNormal(mu_marg, Sigma_marg);
    }

    /**
     * Extracts a univariate marginal distribution.
     *
     * @param index dimension index to extract (0-based)
     * @return Normal distribution for that dimension
     */
    public Normal getMarginalUniv(int index) {
        if (index < 0 || index >= dimension) {
            throw new IllegalArgumentException(
                String.format("Index %d out of range [0, %d)", index, dimension));
        }

        Matrix mu = getMeanVector();
        Matrix Sigma = getCovariance();

        double mean = mu.get(index, 0);
        double variance = Sigma.get(index, index);
        double std = sqrt(max(GlobalConstants.FineTol, variance));

        return new Normal(mean, std);
    }

    // =================== SERIALIZATION ===================

    /**
     * Gets the process representation with actual distribution parameters.
     *
     * @return MatrixCell containing mu and Sigma
     */
    @Override
    public MatrixCell getProcess() {
        MatrixCell representation = new MatrixCell();
        representation.set(0, getMeanVector());
        representation.set(1, getCovariance());
        return representation;
    }

    @Override
    public String toString() {
        return String.format("jline.MultivariateNormal(d=%d)", dimension);
    }
}
