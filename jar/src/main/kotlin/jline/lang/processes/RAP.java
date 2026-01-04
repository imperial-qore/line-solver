/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import static jline.lib.butools.ph.CheckRAPRepresentationKt.checkRAPRepresentation;

/**
 * A Rational Arrival Process (RAP) distribution.
 *
 * RAP is a generalization of the Markovian Arrival Process (MAP) where the matrices
 * H0 and H1 represent hidden and visible transitions respectively, but with relaxed
 * constraints compared to MAP.
 *
 * Representation:
 * - H0: matrix of hidden transition rates (transitions without arrivals)
 * - H1: matrix of visible transition rates (transitions with arrivals)
 * - H0 + H1 must form a valid infinitesimal generator (row sums = 0)
 * - All eigenvalues of H0 must have negative real parts
 * - Dominant eigenvalue of H0 must be negative and real
 *
 * The marginal distribution of inter-arrival times is a Matrix Exponential (ME) distribution.
 */
public class RAP extends Markovian {

    /**
     * Creates a Rational Arrival Process with specified H0 and H1 matrices.
     *
     * @param H0 matrix of hidden transition rates (square matrix)
     * @param H1 matrix of visible transition rates (square matrix, same size as H0)
     * @throws IllegalArgumentException if the representation is invalid
     */
    public RAP(Matrix H0, Matrix H1) {
        super("RAP", 2);

        // Strict validation using BuTools
        if (!checkRAPRepresentation(H0, H1, 1e-14)) {
            throw new IllegalArgumentException("Invalid RAP representation: " +
                    "Check that H0 and H1 are square matrices of the same size, " +
                    "H0 + H1 forms a valid infinitesimal generator (row sums = 0), " +
                    "all eigenvalues of H0 have negative real parts, and " +
                    "the dominant eigenvalue of H0 is real.");
        }

        nPhases = H0.getNumRows();

        // Store parameters
        this.setParam(1, "H0", H0);
        this.setParam(2, "H1", H1);

        // Build MatrixCell process representation
        // RAP process format is same as MAP: {D0=H0, D1=H1}
        MatrixCell rep = new MatrixCell();
        rep.set(0, H0);  // D0 = H0
        rep.set(1, H1);  // D1 = H1

        this.setProcess(rep);
    }

    /**
     * Gets the H0 matrix (hidden transition rates).
     *
     * @return the H0 matrix
     */
    public Matrix getH0() {
        return (Matrix) this.getParam(1).getValue();
    }

    /**
     * Gets the H1 matrix (visible transition rates).
     *
     * @return the H1 matrix
     */
    public Matrix getH1() {
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
     * Creates a RAP by fitting the given moments and correlations.
     * Uses BuTools RAPFromMomentsAndCorrelations algorithm if available.
     *
     * @param moments array of moments
     * @param correlations array of lag-k correlations
     * @return a RAP matching the given moments and correlations
     * @throws UnsupportedOperationException if fitting algorithm is not yet implemented
     */
    public static RAP fitMomentsAndCorrelations(double[] moments, double[] correlations) {
        // This will be implemented after RAPFromMomentsAndCorrelations is available
        throw new UnsupportedOperationException(
                "RAP.fitMomentsAndCorrelations() requires RAPFromMomentsAndCorrelations from BuTools. " +
                "Use the RAP(H0, H1) constructor directly for now.");
    }

    /**
     * Creates a RAP by fitting the given moments.
     * Uses BuTools RAPFromMoments algorithm if available.
     *
     * @param moments array of moments
     * @return a RAP matching the given moments
     * @throws UnsupportedOperationException if fitting algorithm is not yet implemented
     */
    public static RAP fitMoments(double[] moments) {
        // This will be implemented after RAPFromMoments is available
        throw new UnsupportedOperationException(
                "RAP.fitMoments() requires RAPFromMoments from BuTools. " +
                "Use the RAP(H0, H1) constructor directly for now.");
    }

    /**
     * Creates a RAP from a Markovian Arrival Process (MAP).
     * This shows that MAP is a special case of RAP.
     *
     * @param map the MAP to convert to RAP
     * @return a RAP equivalent to the given MAP
     */
    public static RAP fromMAP(MAP map) {
        Matrix D0 = map.D(0);
        Matrix D1 = map.D(1);
        return new RAP(D0, D1);
    }

    /**
     * Creates a RAP from an exponential renewal process (Poisson process).
     *
     * @param rate the arrival rate (lambda)
     * @return a RAP representing a Poisson process
     */
    public static RAP fromPoisson(double rate) {
        Matrix H0 = new Matrix(new double[][]{{-rate}});
        Matrix H1 = new Matrix(new double[][]{{rate}});
        return new RAP(H0, H1);
    }

    /**
     * Creates a RAP from an Erlang renewal process.
     *
     * @param k number of phases
     * @param rate rate parameter for each phase
     * @return a RAP representing an Erlang renewal process
     */
    public static RAP fromErlang(int k, double rate) {
        Matrix H0 = new Matrix(k, k);
        Matrix H1 = new Matrix(k, k);

        for (int i = 0; i < k; i++) {
            H0.set(i, i, -rate);  // diagonal
            if (i < k - 1) {
                H0.set(i, i + 1, rate);  // transitions to next phase (no arrival)
            }
        }
        // Last phase transition with arrival back to first phase
        H1.set(k - 1, 0, rate);

        return new RAP(H0, H1);
    }

    /**
     * Creates a RAP from a Hyper-exponential distribution.
     *
     * @param p probability of using the first phase
     * @param lambda1 rate parameter for the first phase
     * @param lambda2 rate parameter for the second phase
     * @return a RAP representing a Hyper-exponential renewal process
     */
    public static RAP fromHyperExp(double p, double lambda1, double lambda2) {
        HyperExp hyperExp = new HyperExp(p, lambda1, lambda2);
        Matrix D0 = hyperExp.D(0);
        Matrix D1 = hyperExp.D(1);
        MAP map = new MAP(D0, D1);
        return fromMAP(map);
    }

    /**
     * Creates a RAP by fitting the given mean and squared coefficient of variation (SCV).
     * Uses Erlang fitting for SCV < 1 and HyperExp fitting for SCV >= 1.
     *
     * @param mean target mean inter-arrival time
     * @param scv target squared coefficient of variation
     * @return a RAP matching the given mean and SCV
     */
    public static RAP fitMeanAndSCV(double mean, double scv) {
        if (scv < 1.0) {
            // For SCV < 1, use Erlang distribution
            // Erlang-k has SCV = 1/k, so k = 1/scv
            int k = (int) Math.round(1.0 / scv);
            if (k < 1) {
                k = 1;
            }
            // Erlang-k mean = k/rate, so rate = k/mean
            double rate = k / mean;
            return fromErlang(k, rate);
        } else {
            // For SCV >= 1, use HyperExp distribution
            HyperExp hyperExp = HyperExp.fitMeanAndSCV(mean, scv);
            Matrix D0 = hyperExp.D(0);
            Matrix D1 = hyperExp.D(1);
            MAP map = new MAP(D0, D1);
            return fromMAP(map);
        }
    }
}
