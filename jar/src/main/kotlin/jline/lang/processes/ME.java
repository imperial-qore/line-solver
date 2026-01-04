/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import static jline.lib.butools.ph.CheckMERepresentationKt.checkMERepresentation;

/**
 * A Matrix Exponential (ME) distribution.
 *
 * ME distributions are characterized by an initial vector alpha and a matrix parameter A.
 * They generalize Phase-Type (PH) distributions by allowing alpha to have entries outside [0,1]
 * and A to have arbitrary structure (not necessarily a valid sub-generator).
 *
 * Representation:
 * - alpha: initial vector (may have negative entries or sum != 1)
 * - A: matrix parameter (must have all eigenvalues with negative real parts)
 * - Dominant eigenvalue of A must be negative and real
 *
 * The distribution function is: F(t) = 1 - alpha * exp(A*t) * e (under certain conditions)
 * Moments: m_k = k! * alpha * (-A)^(-k) * e
 */
public class ME extends Markovian {

    /**
     * Creates a Matrix Exponential distribution with specified initial vector and matrix parameter.
     *
     * @param alpha the initial vector (row vector as Matrix)
     * @param A the matrix parameter (must be square with negative real eigenvalues)
     * @throws IllegalArgumentException if the representation is invalid
     */
    public ME(Matrix alpha, Matrix A) {
        super("ME", 2);

        // Ensure alpha is a row vector (1 x n)
        // If passed as column vector (n x 1), transpose it
        Matrix alphaRow = alpha;
        if (alpha.getNumRows() > 1 && alpha.getNumCols() == 1) {
            alphaRow = alpha.transpose();
        }

        // Strict validation using BuTools
        if (!checkMERepresentation(alphaRow, A, 1e-14)) {
            throw new IllegalArgumentException("Invalid ME representation: " +
                    "Check that A is square, alpha and A have compatible dimensions, " +
                    "all eigenvalues of A have negative real parts, and " +
                    "the dominant eigenvalue is real.");
        }

        nPhases = alphaRow.getNumElements();

        // Store parameters (always store as row vector)
        this.setParam(1, "alpha", alphaRow);
        this.setParam(2, "A", A);

        // Build MatrixCell process representation for compatibility with map_* functions
        // Process format: {D0=A, D1=-A*e*alpha'}
        MatrixCell rep = new MatrixCell();
        rep.set(0, A);  // D0 = A

        // D1 = -A * e * alpha' = outer product of (-A*e) and alpha
        // where e is column vector of ones
        Matrix ones = Matrix.ones(nPhases, 1);
        Matrix Ae = A.mult(ones).scale(-1);  // -A * e (column vector, n x 1)
        Matrix D1 = Ae.mult(alphaRow);  // outer product: (n x 1) * (1 x n) = (n x n)
        rep.set(1, D1);

        this.setProcess(rep);
    }

    /**
     * Gets the initial vector alpha.
     *
     * @return the initial vector as a Matrix
     */
    public Matrix getAlpha() {
        return (Matrix) this.getParam(1).getValue();
    }

    /**
     * Gets the matrix parameter A.
     *
     * @return the matrix parameter A
     */
    public Matrix getA() {
        return (Matrix) this.getParam(2).getValue();
    }

    @Override
    public long getNumberOfPhases() {
        return nPhases;
    }

    @Override
    public MatrixCell getProcess() {
        return this.process;
    }

    /**
     * Creates an ME distribution by fitting the given moments.
     * Uses BuTools MEFromMoments algorithm.
     *
     * @param moments array of moments (requires 2*M-1 moments for order M ME distribution)
     * @return an ME distribution matching the given moments
     * @throws IllegalArgumentException if moments are invalid or fitting fails
     */
    public static ME fitMoments(double[] moments) {
        // Use the MEFromMoments algorithm to fit the distribution
        jline.lib.butools.ph.MERepresentation rep = jline.lib.butools.ph.MEFromMomentsKt.meFromMoments(moments);
        return new ME(rep.getAlpha(), rep.getA());
    }

    /**
     * Creates an ME distribution from an exponential distribution.
     * This is a convenience method showing that Exp is a special case of ME.
     *
     * @param rate the rate parameter (lambda)
     * @return an ME distribution equivalent to Exp(rate)
     */
    public static ME fromExp(double rate) {
        Matrix alpha = new Matrix(new double[]{1.0});
        Matrix A = new Matrix(new double[][]{{-rate}});
        return new ME(alpha, A);
    }

    /**
     * Creates an ME distribution from an Erlang distribution.
     * This is a convenience method showing that Erlang is a special case of ME.
     *
     * @param k number of phases
     * @param rate rate parameter for each phase
     * @return an ME distribution equivalent to Erlang(k, rate)
     */
    public static ME fromErlang(int k, double rate) {
        Matrix alpha = new Matrix(1, k);
        alpha.set(0, 0, 1.0);  // alpha = [1, 0, 0, ..., 0]

        Matrix A = new Matrix(k, k);
        for (int i = 0; i < k; i++) {
            A.set(i, i, -rate);  // diagonal
            if (i < k - 1) {
                A.set(i, i + 1, rate);  // super-diagonal
            }
        }

        return new ME(alpha, A);
    }

    /**
     * Creates an ME distribution from a HyperExponential distribution.
     * This is a convenience method showing that HyperExp is a special case of ME.
     *
     * @param p array of probabilities for each branch
     * @param rates array of rates for each branch
     * @return an ME distribution equivalent to HyperExp(p, rates)
     */
    public static ME fromHyperExp(double[] p, double[] rates) {
        if (p.length != rates.length) {
            throw new IllegalArgumentException("p and rates must have the same length");
        }

        int k = p.length;
        // Create alpha as row vector (1 x k), not column vector
        Matrix alpha = new Matrix(1, k);
        for (int i = 0; i < k; i++) {
            alpha.set(0, i, p[i]);
        }

        Matrix A = new Matrix(k, k);
        for (int i = 0; i < k; i++) {
            A.set(i, i, -rates[i]);  // diagonal matrix of rates
        }

        return new ME(alpha, A);
    }
}
