/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.mva;

import jline.util.matrix.Matrix;

/**
 * Result class for probability calculations in solvers.
 *
 * <p>Supports both scalar probabilities (single value) and probability distributions
 * (vectors/matrices) for marginal and aggregated state probabilities.</p>
 *
 * <p>For marginal probabilities (getProbMarg), the probability vector contains
 * P(n jobs of class r) for n=0,1,...,N where N is the population.</p>
 *
 * <p>For aggregated probabilities (getProbAggr), the probability is a scalar
 * representing P(n1 jobs of class 1, n2 jobs of class 2, ...).</p>
 */
public class ProbabilityResult {
    /**
     * The scalar probability value (for getProbAggr-style results)
     */
    public final double probability;

    /**
     * The logarithm of the scalar probability value
     */
    public final double logProbability;

    /**
     * The probability distribution matrix (for getProbMarg-style results).
     * For marginal probabilities, this is a row vector where element k = P(k jobs).
     */
    private Matrix probabilityVector;

    /**
     * Storage for multi-station, multi-class marginal probabilities.
     * Indexed as [station][class].
     */
    private Matrix[][] marginals;

    /**
     * Number of stations for marginal storage
     */
    private int nstations;

    /**
     * Number of classes for marginal storage
     */
    private int nclasses;

    /**
     * Constructor for scalar probability result.
     *
     * @param probability the probability value
     * @param logProbability the logarithm of the probability value
     */
    public ProbabilityResult(double probability, double logProbability) {
        this.probability = probability;
        this.logProbability = logProbability;
        this.probabilityVector = null;
        this.marginals = null;
    }

    /**
     * Constructor for probability distribution result.
     *
     * @param probabilityVector the probability distribution as a row vector
     */
    public ProbabilityResult(Matrix probabilityVector) {
        this.probabilityVector = probabilityVector;
        // Compute scalar probability as sum (should be ~1.0 for valid distribution)
        double sum = 0;
        double maxLogP = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < probabilityVector.length(); i++) {
            double p = probabilityVector.get(0, i);
            sum += p;
            if (p > 0) {
                double logP = Math.log(p);
                if (logP > maxLogP) {
                    maxLogP = logP;
                }
            }
        }
        this.probability = sum;
        this.logProbability = maxLogP;
        this.marginals = null;
    }

    /**
     * Constructor for multi-station, multi-class marginal storage.
     *
     * @param nstations number of stations
     * @param nclasses number of job classes
     */
    public ProbabilityResult(int nstations, int nclasses) {
        this.nstations = nstations;
        this.nclasses = nclasses;
        this.marginals = new Matrix[nstations][nclasses];
        this.probability = 0;
        this.logProbability = Double.NEGATIVE_INFINITY;
        this.probabilityVector = null;
    }

    /**
     * Set the marginal probability distribution for a specific station and class.
     *
     * @param station station index
     * @param jobclass job class index
     * @param probVector probability distribution vector
     */
    public void setMarginal(int station, int jobclass, Matrix probVector) {
        if (marginals == null) {
            throw new IllegalStateException("ProbabilityResult not initialized for marginal storage");
        }
        if (station >= nstations || jobclass >= nclasses) {
            throw new IndexOutOfBoundsException("Station or class index out of bounds");
        }
        marginals[station][jobclass] = probVector;
        // Also set as primary probability vector for single-query access
        this.probabilityVector = probVector;
    }

    /**
     * Get the marginal probability distribution for a specific station and class.
     *
     * @param station station index
     * @param jobclass job class index
     * @return probability distribution vector, or null if not set
     */
    public Matrix getMarginal(int station, int jobclass) {
        if (marginals == null) {
            return probabilityVector;  // Single-station case
        }
        if (station >= nstations || jobclass >= nclasses) {
            throw new IndexOutOfBoundsException("Station or class index out of bounds");
        }
        return marginals[station][jobclass];
    }

    /**
     * Get the probability distribution vector.
     *
     * @return the probability vector, or null if this is a scalar result
     */
    public Matrix getProbabilityVector() {
        return probabilityVector;
    }

    /**
     * Check if this result contains a probability distribution (vs scalar).
     *
     * @return true if this result contains a probability vector
     */
    public boolean hasDistribution() {
        return probabilityVector != null;
    }

    /**
     * Get probability at a specific state index.
     *
     * @param stateIndex the state index (e.g., number of jobs)
     * @return the probability at that state
     */
    public double getProbability(int stateIndex) {
        if (probabilityVector == null) {
            throw new IllegalStateException("No probability distribution available");
        }
        if (stateIndex >= probabilityVector.length()) {
            return 0.0;  // States beyond computed range have zero probability
        }
        return probabilityVector.get(0, stateIndex);
    }

    /**
     * Get the length of the probability distribution.
     *
     * @return number of states in the distribution, or 1 for scalar results
     */
    public int length() {
        if (probabilityVector == null) {
            return 1;
        }
        return (int) probabilityVector.length();
    }
}