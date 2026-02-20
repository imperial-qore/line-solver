/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.util.matrix.Matrix;
import jline.util.RandomManager;
import java.util.Random;

import static jline.api.mc.Ctmc_makeinfgenKt.ctmc_makeinfgen;
import static jline.api.mc.Ctmc_randKt.ctmc_rand;
import static jline.api.mc.Ctmc_solveKt.ctmc_solve;
import static jline.api.mc.Ctmc_timereverseKt.ctmc_timereverse;

/**
 * A class for a continuous time Markov chain
 */
public class MarkovProcess extends Process {
    protected Matrix infGen;
    protected Matrix stateSpace;
    protected boolean isfinite;

    /**
     * Creates a CTMC with the specified infinitesimal generator
     * @param infGen the infinitesimal generator matrix
     */
    public MarkovProcess(Matrix infGen) {
        this(infGen, true, null);
    }

    /**
     * Creates a CTMC with the specified infinitesimal generator and finite flag
     * @param infGen the infinitesimal generator matrix
     * @param isFinite whether the CTMC is finite
     */
    public MarkovProcess(Matrix infGen, boolean isFinite) {
        this(infGen, isFinite, null);
    }

    /**
     * Creates a CTMC with the specified infinitesimal generator, finite flag, and state space
     * @param infGen the infinitesimal generator matrix
     * @param isFinite whether the CTMC is finite
     * @param stateSpace the state space representation
     */
    public MarkovProcess(Matrix infGen, boolean isFinite, Matrix stateSpace) {
        super("CTMC", 1);
        this.infGen = ctmc_makeinfgen(infGen);
        this.isfinite = isFinite;
        this.stateSpace = stateSpace;
    }

    /**
     * Convert to DTMC using uniformization
     * @return the equivalent DTMC
     */
    public MarkovChain toDTMC() {
        double q = infGen.elementMax() + RandomManager.nextDouble();
        return toDTMC(q);
    }

    public MarkovChain toDTMC(double q) {
        Matrix P = infGen.scale(1.0/q).add(Matrix.eye(infGen.getNumRows()));
        MarkovChain dtmc = new MarkovChain(P);
        if (stateSpace != null) {
            dtmc.setStateSpace(stateSpace);
        }
        return dtmc;
    }

    /**
     * Get time-reversed CTMC
     * @return time-reversed CTMC
     */
    public MarkovProcess toTimeReversed() {
        return new MarkovProcess(ctmc_timereverse(infGen));
    }

    /**
     * Set the state space
     * @param stateSpace the state space matrix
     */
    public void setStateSpace(Matrix stateSpace) {
        this.stateSpace = stateSpace;
    }

    /**
     * Get the infinitesimal generator matrix
     * @return the generator matrix
     */
    public Matrix getGenerator() {
        return infGen;
    }

    /**
     * Get probability of a specific state using steady-state analysis
     * @param state the state index (0-based)
     * @return the probability of the specified state
     */
    public double getProbState(int state) {
        Matrix steadyState = solve();
        if (state < 0 || state >= steadyState.getNumRows()) {
            throw new IndexOutOfBoundsException("State index " + state + " is out of bounds");
        }
        return steadyState.get(state, 0);
    }

    /**
     * Get probability of a specific state using Cramer's rule
     * @param state the state vector (if state space is defined) or single state index
     * @return the probability matrix containing the state probability
     */
    public Matrix getProbState(Matrix state) {
        if (stateSpace != null && state.getNumRows() == stateSpace.getNumCols()) {
            // Find matching state in state space
            for (int i = 0; i < stateSpace.getNumRows(); i++) {
                boolean match = true;
                for (int j = 0; j < stateSpace.getNumCols(); j++) {
                    if (Math.abs(stateSpace.get(i, j) - state.get(j, 0)) > 1e-10) {
                        match = false;
                        break;
                    }
                }
                if (match) {
                    Matrix result = new Matrix(1, 1);
                    result.set(0, 0, getProbState(i));
                    return result;
                }
            }
            // State not found in state space
            Matrix result = new Matrix(1, 1);
            result.set(0, 0, 0.0);
            return result;
        } else {
            // Assume state is a single index
            int stateIndex = (int) state.get(0, 0);
            Matrix result = new Matrix(1, 1);
            result.set(0, 0, getProbState(stateIndex));
            return result;
        }
    }

    /**
     * Solve the CTMC for steady-state probabilities
     * @return the steady-state probability vector
     */
    public Matrix solve() {
        return ctmc_solve(infGen);
    }

    /**
     * Sample from the CTMC
     * @return sampled value
     */
    @Override
    public Matrix sample() {
        return sample(1);
    }

    /**
     * Sample n state transitions from the CTMC starting from steady-state
     * @param n number of state transitions to sample
     * @return matrix containing sampled states (column vector)
     */
    public Matrix sample(int n) {
        Matrix pi0 = solve(); // Start from steady-state
        return sample(n, pi0, RandomManager.getThreadRandomAsRandom());
    }

    /**
     * Sample n state transitions from the CTMC
     * @param n number of state transitions to sample
     * @param pi0 initial state distribution (if null, uses steady-state)
     * @param random random number generator
     * @return matrix containing sampled states (column vector)
     */
    public Matrix sample(int n, Matrix pi0, Random random) {
        int numStates = infGen.getNumRows();
        
        // Set initial distribution if not provided
        if (pi0 == null) {
            pi0 = solve(); // Use steady-state distribution
        }
        
        // Sample initial state from pi0
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
        
        // Precompute cumulative transition probabilities for each state
        Matrix F = new Matrix(numStates, numStates);
        for (int i = 0; i < numStates; i++) {
            double rowSum = 0.0;
            // Calculate sum of off-diagonal elements (transition rates)
            for (int j = 0; j < numStates; j++) {
                if (i != j) {
                    rowSum += infGen.get(i, j);
                }
            }
            
            if (rowSum > 0) {
                double cumulativeProb = 0.0;
                for (int j = 0; j < numStates; j++) {
                    if (i != j) {
                        cumulativeProb += infGen.get(i, j) / rowSum;
                    }
                    F.set(i, j, cumulativeProb);
                }
            }
        }
        
        Matrix result = new Matrix(n, 1);
        
        // Sample n state transitions
        for (int step = 0; step < n; step++) {
            result.set(step, 0, currentState);
            
            // Sample next state based on current state's transition rates
            if (step < n - 1) { // Don't sample next state for the last step
                r = random.nextDouble();
                for (int j = 0; j < numStates; j++) {
                    if (r < F.get(currentState, j)) {
                        currentState = j;
                        break;
                    }
                }
            }
        }
        
        return result;
    }

    /**
     * Get the state space
     * @return the state space matrix
     */
    public Matrix getStateSpace() {
        return stateSpace;
    }

    /**
     * Check if the CTMC is finite
     * @return true if finite
     */
    public boolean isFinite() {
        return isfinite;
    }

    /**
     * Create a random CTMC
     * @param nStates number of states
     * @return random CTMC
     */
    public static MarkovProcess rand(int nStates) {
        return new MarkovProcess(ctmc_rand(nStates));
    }

    /**
     * Create CTMC from sample system aggregation
     * @param samples matrix where each row is a sample trajectory, each column is a time step
     * @param sojournTimes matrix of sojourn times corresponding to each state in samples
     * @return CTMC constructed from samples by estimating transition rates
     */
    public static MarkovProcess fromSampleSysAggr(Matrix samples, Matrix sojournTimes) {
        return fromSampleSysAggr(samples, sojournTimes, -1);
    }

    /**
     * Create CTMC from sample system aggregation
     * @param samples matrix where each row is a sample trajectory, each column is a time step
     * @param sojournTimes matrix of sojourn times corresponding to each state in samples
     * @param numStates number of states (if -1, infer from data)
     * @return CTMC constructed from samples by estimating transition rates
     */
    public static MarkovProcess fromSampleSysAggr(Matrix samples, Matrix sojournTimes, int numStates) {
        if (samples == null || samples.getNumRows() == 0 || samples.getNumCols() < 2) {
            throw new IllegalArgumentException("Invalid sample data: need at least 1 trajectory with 2+ time steps");
        }
        if (sojournTimes == null || sojournTimes.getNumRows() != samples.getNumRows() || 
            sojournTimes.getNumCols() != samples.getNumCols()) {
            throw new IllegalArgumentException("Sojourn times matrix must have same dimensions as samples matrix");
        }

        // Determine number of states if not provided
        if (numStates <= 0) {
            double maxState = samples.elementMax();
            double minState = samples.elementMin();
            numStates = (int) (maxState - minState) + 1;
        }

        // Count transitions and accumulate sojourn times
        Matrix transitionCounts = new Matrix(numStates, numStates);
        Matrix totalSojournTimes = new Matrix(numStates, 1);
        
        for (int trajectory = 0; trajectory < samples.getNumRows(); trajectory++) {
            for (int t = 0; t < samples.getNumCols() - 1; t++) {
                int fromState = (int) samples.get(trajectory, t);
                int toState = (int) samples.get(trajectory, t + 1);
                double sojournTime = sojournTimes.get(trajectory, t);
                
                // Validate state indices and positive sojourn time
                if (fromState >= 0 && fromState < numStates && toState >= 0 && toState < numStates && sojournTime > 0) {
                    // Count transition
                    if (fromState != toState) { // Only count actual state changes
                        transitionCounts.set(fromState, toState, 
                            transitionCounts.get(fromState, toState) + 1.0);
                    }
                    
                    // Accumulate sojourn time in from state
                    totalSojournTimes.set(fromState, 0, 
                        totalSojournTimes.get(fromState, 0) + sojournTime);
                }
            }
        }

        // Convert counts and times to rates
        Matrix infGen = new Matrix(numStates, numStates);
        for (int i = 0; i < numStates; i++) {
            double totalTime = totalSojournTimes.get(i, 0);
            
            if (totalTime > 0) {
                // Calculate transition rates
                double totalRate = 0.0;
                for (int j = 0; j < numStates; j++) {
                    if (i != j) {
                        double rate = transitionCounts.get(i, j) / totalTime;
                        infGen.set(i, j, rate);
                        totalRate += rate;
                    }
                }
                
                // Set diagonal element (negative sum of off-diagonal rates)
                infGen.set(i, i, -totalRate);
            }
            // If no sojourn time observed, leave rates as zero (absorbing state)
        }

        return new MarkovProcess(infGen);
    }

    /**
     * Create CTMC from sample system aggregation (assuming unit sojourn times)
     * @param samples matrix where each row is a sample trajectory, each column is a time step
     * @return CTMC constructed from samples with unit time steps
     */
    public static MarkovProcess fromSampleSysAggr(Matrix samples) {
        // Create unit sojourn times matrix
        Matrix unitTimes = new Matrix(samples.getNumRows(), samples.getNumCols());
        for (int i = 0; i < samples.getNumRows(); i++) {
            for (int j = 0; j < samples.getNumCols(); j++) {
                unitTimes.set(i, j, 1.0);
            }
        }
        return fromSampleSysAggr(samples, unitTimes);
    }

    /**
     * Create CTMC from sample system aggregation (legacy interface)
     * @param sa the sample aggregation object
     * @return CTMC constructed from samples
     */
    public static MarkovProcess fromSampleSysAggr(Object sa) {
        if (sa instanceof Matrix) {
            return fromSampleSysAggr((Matrix) sa);
        } else {
            throw new UnsupportedOperationException("fromSampleSysAggr only supports Matrix input currently");
        }
    }
}