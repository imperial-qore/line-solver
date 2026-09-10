/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.util.matrix.Matrix;
import jline.api.mc.Dtmc_stochcomp;

import static jline.api.mc.Dtmc_hitting_time.dtmc_hitting_time;
import static jline.api.mc.Dtmc_isfeasible.dtmc_isfeasible;
import static jline.api.mc.Dtmc_makestochastic.dtmc_makestochastic;
import static jline.api.mc.Dtmc_rand.dtmc_rand;
import static jline.api.mc.Dtmc_solve.dtmc_solve;
import static jline.api.mc.Dtmc_stochcomp.dtmc_stochcomp;
import static jline.api.mc.Dtmc_stochcomp.dtmc_stochcomp_full;
import static jline.api.mc.Dtmc_uniformization.dtmc_uniformization;
import static jline.api.mc.Dtmc_timereverse.dtmc_timereverse;
import static jline.api.mc.Dtmc_transient.dtmc_transient;

/**
 * A class for a discrete time Markov chain
 */
public class MarkovChain extends Process {
    protected Matrix transMat;
    protected Matrix stateSpace;
    protected boolean isfinite;

    /**
     * Creates a DTMC with the specified transition matrix
     * @param transMat the transition matrix
     */
    public MarkovChain(Matrix transMat) {
        this(transMat, true);
    }

    /**
     * Creates a DTMC with the specified transition matrix and finite flag
     * @param transMat the transition matrix
     * @param isFinite whether the DTMC is finite
     */
    public MarkovChain(Matrix transMat, boolean isFinite) {
        super("DTMC", 1);
        this.transMat = dtmc_makestochastic(transMat);
        this.stateSpace = null;
        this.isfinite = isFinite;
    }

    /**
     * Convert to CTMC by subtracting identity matrix
     * @return the equivalent CTMC
     */
    public MarkovProcess toCTMC() {
        Matrix Q = transMat.sub(Matrix.eye(transMat.getNumRows()));
        MarkovProcess ctmc = new MarkovProcess(Q);
        if (stateSpace != null) {
            ctmc.setStateSpace(stateSpace);
        }
        return ctmc;
    }

    /**
     * Get time-reversed DTMC
     * @return time-reversed DTMC
     */
    public MarkovChain toTimeReversed() {
        return new MarkovChain(dtmc_timereverse(transMat));
    }

    /**
     * Get the transition matrix
     * @return the transition matrix
     */
    public Matrix getTransMat() {
        return transMat;
    }

    /**
     * Solve the DTMC for steady-state probabilities. Twin of MarkovProcess.solve.
     * @return the steady-state probability vector
     */
    public Matrix solve() {
        return dtmc_solve(transMat);
    }

    /**
     * Distribution at each step 0,...,steps.
     *
     * @param pi0   initial distribution, uniform when null
     * @param steps number of steps
     * @return matrix with one row per step
     */
    public Matrix transientProb(Matrix pi0, int steps) {
        return dtmc_transient(transMat, pi0, steps);
    }

    /**
     * Mean number of steps to reach any target state.
     *
     * @param targetStates 0-based indices of the target states
     * @return column vector of mean hitting times
     */
    public Matrix hittingTime(int[] targetStates) {
        return dtmc_hitting_time(transMat, targetStates);
    }

    /**
     * Stochastic complement of a subset of states.
     *
     * @param I 0-based indices of the states to retain
     * @return the transition matrix of the complement on those states
     */
    public Matrix stochComp(java.util.List<Integer> I) {
        return dtmc_stochcomp(transMat, I);
    }

    /**
     * Stochastic complement of a subset of states together with the four blocks
     * of the transition matrix partitioned by I and its complement. Twin of the
     * MATLAB [S,P11,P12,P21,P22] = dtmc.stochCompFull(I) and of the Python
     * stochCompFull.
     *
     * @param I 0-based indices of the states to retain
     * @return the complement and the blocks it was built from
     */
    public Dtmc_stochcomp.DtmcStochCompResult stochCompFull(java.util.List<Integer> I) {
        return dtmc_stochcomp_full(transMat, I);
    }

    /**
     * Distribution at time t of the DTMC seen through uniformization. The chain
     * is read as the randomized image of a CTMC, so t is continuous here,
     * unlike the step count taken by transientProb.
     *
     * @param pi0 initial distribution, uniform when null
     * @param t   time point
     * @return the distribution at time t
     */
    public Matrix transientUnif(Matrix pi0, double t) {
        int n = transMat.getNumRows();
        Matrix start = new Matrix(1, n);
        for (int i = 0; i < n; i++) {
            start.set(0, i, (pi0 == null || pi0.length() != n) ? 1.0 / n : pi0.get(i));
        }
        return dtmc_uniformization(start, transMat, t).pi;
    }

    /**
     * @return true when the transition matrix is stochastic
     */
    public boolean isFeasible() {
        return dtmc_isfeasible(transMat) > 0;
    }

    /**
     * Set the state space
     * @param stateSpace the state space matrix
     */
    public void setStateSpace(Matrix stateSpace) {
        this.stateSpace = stateSpace;
    }

    /**
     * Get the state space
     * @return the state space matrix
     */
    public Matrix getStateSpace() {
        return stateSpace;
    }

    /**
     * Check if the DTMC is finite
     * @return true if finite
     */
    public boolean isFinite() {
        return isfinite;
    }

    /**
     * Sample from the DTMC
     * @return sampled value
     */
    @Override
    public Matrix sample() {
        return sample(1);
    }

    /**
     * Sample n steps from the DTMC starting from initial distribution
     * @param n number of steps to sample
     * @return matrix containing sampled states (column vector)
     */
    public Matrix sample(int n) {
        return sample(n, null, new java.util.Random());
    }

    /**
     * Sample n steps from the DTMC
     * @param n number of steps to sample
     * @param pi0 initial state distribution (if null, uses uniform distribution)
     * @param random random number generator
     * @return matrix containing sampled states (column vector)
     */
    public Matrix sample(int n, Matrix pi0, java.util.Random random) {
        int numStates = transMat.getNumRows();
        
        // Set initial distribution if not provided
        if (pi0 == null) {
            pi0 = new Matrix(numStates, 1);
            for (int i = 0; i < numStates; i++) {
                pi0.set(i, 0, 1.0 / numStates);
            }
        }
        
        // Sample initial state
        double r = random.nextDouble();
        double cumulative = 0.0;
        int currentState = 0;
        for (int i = 0; i < numStates; i++) {
            cumulative += pi0.get(i, 0);
            if (r < cumulative) {
                currentState = i;
                break;
            }
        }
        
        Matrix result = new Matrix(n, 1);
        
        // Sample n steps
        for (int step = 0; step < n; step++) {
            result.set(step, 0, currentState);
            
            // Sample next state based on current state's transition probabilities
            if (step < n - 1) { // Don't sample next state for the last step
                r = random.nextDouble();
                cumulative = 0.0;
                for (int j = 0; j < numStates; j++) {
                    cumulative += transMat.get(currentState, j);
                    if (r < cumulative) {
                        currentState = j;
                        break;
                    }
                }
            }
        }
        
        return result;
    }

    /**
     * Create a random DTMC
     * @param nStates number of states
     * @return random DTMC
     */
    public static MarkovChain rand(int nStates) {
        return new MarkovChain(dtmc_rand(nStates));
    }

    /**
     * Create DTMC from sample system aggregation
     * @param samples matrix where each row is a sample trajectory, each column is a time step
     * @return DTMC constructed from samples by estimating transition probabilities
     */
    public static MarkovChain fromSampleSysAggr(Matrix samples) {
        return fromSampleSysAggr(samples, -1);
    }

    /**
     * Create DTMC from sample system aggregation
     * @param samples matrix where each row is a sample trajectory, each column is a time step
     * @param numStates number of states (if -1, infer from data)
     * @return DTMC constructed from samples by estimating transition probabilities
     */
    public static MarkovChain fromSampleSysAggr(Matrix samples, int numStates) {
        if (samples == null || samples.getNumRows() == 0 || samples.getNumCols() < 2) {
            throw new IllegalArgumentException("Invalid sample data: need at least 1 trajectory with 2+ time steps");
        }

        // Determine number of states if not provided
        if (numStates <= 0) {
            double maxState = samples.elementMax();
            double minState = samples.elementMin();
            numStates = (int) (maxState - minState) + 1;
        }

        // Count transitions
        Matrix transitionCounts = new Matrix(numStates, numStates);
        
        for (int trajectory = 0; trajectory < samples.getNumRows(); trajectory++) {
            for (int t = 0; t < samples.getNumCols() - 1; t++) {
                int fromState = (int) samples.get(trajectory, t);
                int toState = (int) samples.get(trajectory, t + 1);
                
                // Validate state indices
                if (fromState >= 0 && fromState < numStates && toState >= 0 && toState < numStates) {
                    transitionCounts.set(fromState, toState, 
                        transitionCounts.get(fromState, toState) + 1.0);
                }
            }
        }

        // Convert counts to probabilities
        Matrix transitionMatrix = new Matrix(numStates, numStates);
        for (int i = 0; i < numStates; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < numStates; j++) {
                rowSum += transitionCounts.get(i, j);
            }
            
            if (rowSum > 0) {
                for (int j = 0; j < numStates; j++) {
                    transitionMatrix.set(i, j, transitionCounts.get(i, j) / rowSum);
                }
            } else {
                // If no transitions observed from state i, make it absorbing
                transitionMatrix.set(i, i, 1.0);
            }
        }

        return new MarkovChain(transitionMatrix);
    }

    /**
     * Create DTMC from sample system aggregation (legacy interface)
     * @param sa the sample aggregation object
     * @return DTMC constructed from samples
     */
    public static MarkovChain fromSampleSysAggr(Object sa) {
        if (sa instanceof Matrix) {
            return fromSampleSysAggr((Matrix) sa);
        } else {
            throw new UnsupportedOperationException("fromSampleSysAggr only supports Matrix input currently");
        }
    }
}