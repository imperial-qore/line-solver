/*
 * Copyright (c) 2012-2025, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.fluid.handlers;

import jline.GlobalConstants;
import jline.lang.NetworkStruct;
import jline.solvers.SolverOptions;
import jline.solvers.ctmc.SolverCTMC;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import static jline.api.mc.Ctmc_stochcompKt.ctmc_stochcomp;
import static jline.io.InputOutputKt.line_warning;
import jline.VerboseLevel;

/**
 * Eliminates immediate transitions from Fluid ODE system using stochastic complementation
 * to reduce stiffness and improve performance
 */
public class ImmediateElimination {

    /**
     * Result of immediate elimination containing reduced system and state mapping
     */
    public static class EliminationResult {
        public Matrix allJumpsReduced;
        public Matrix rateBaseReduced;
        public Matrix eventIdxReduced;
        public int[] stateMap;

        public EliminationResult(Matrix allJumpsReduced, Matrix rateBaseReduced,
                                Matrix eventIdxReduced, int[] stateMap) {
            this.allJumpsReduced = allJumpsReduced;
            this.rateBaseReduced = rateBaseReduced;
            this.eventIdxReduced = eventIdxReduced;
            this.stateMap = stateMap;
        }
    }

    /**
     * Eliminate immediate transitions from ODE system
     *
     * @param allJumps  [n_states x n_transitions] matrix of state changes
     * @param rateBase  [n_transitions x 1] vector of base rates
     * @param eventIdx  [n_transitions x 1] vector of state indices for rate gating
     * @param sn        Model structure
     * @param options   Solver options
     * @return EliminationResult containing reduced system
     */
    public static EliminationResult eliminateImmediate(
            Matrix allJumps,
            Matrix rateBase,
            Matrix eventIdx,
            NetworkStruct sn,
            SolverOptions options) {

        // Get immediate detection threshold
        double immTol = GlobalConstants.Immediate / 10.0; // Default: 1e7
        if (options.config != null && options.config.containsKey("immediate_tol")) {
            immTol = (Double) options.config.get("immediate_tol");
        }

        // Identify immediate transitions
        List<Integer> immIdx = new ArrayList<Integer>();
        for (int i = 0; i < rateBase.getNumRows(); i++) {
            if (rateBase.get(i, 0) >= immTol) {
                immIdx.add(i);
            }
        }

        // If no immediate transitions, return unchanged
        if (immIdx.isEmpty()) {
            int[] identityMap = new int[(int) allJumps.getNumRows()];
            for (int i = 0; i < identityMap.length; i++) {
                identityMap[i] = i;
            }
            return new EliminationResult(allJumps, rateBase, eventIdx, identityMap);
        }

        // Report detected immediate transitions
        if (options.verbose == VerboseLevel.STD || options.verbose == VerboseLevel.DEBUG) {
            System.out.println(String.format("Eliminating %d immediate transitions (out of %d total)",
                    immIdx.size(), rateBase.getNumRows()));
        }

        try {
            // Apply stochastic complement elimination
            return eliminateViaStochcomp(allJumps, rateBase, eventIdx, immIdx, sn, options);
        } catch (Exception e) {
            // Elimination failed - fall back to original system
            line_warning("ImmediateElimination",
                    String.format("Immediate transition elimination failed: %s. Using original system.",
                            e.getMessage()));
            int[] identityMap = new int[(int) allJumps.getNumRows()];
            for (int i = 0; i < identityMap.length; i++) {
                identityMap[i] = i;
            }
            return new EliminationResult(allJumps, rateBase, eventIdx, identityMap);
        }
    }

    private static EliminationResult eliminateViaStochcomp(
            Matrix allJumps,
            Matrix rateBase,
            Matrix eventIdx,
            List<Integer> immIdx,
            NetworkStruct sn,
            SolverOptions options) {

        int nStates = (int) allJumps.getNumRows();
        int nTransitions = (int) rateBase.getNumRows();

        // Step 1: Construct infinitesimal generator matrix from jumps and rates
        Matrix W = new Matrix(nStates, nStates);

        for (int t = 0; t < nTransitions; t++) {
            int srcState = (int) eventIdx.get(t, 0);

            // Find destination state (where jump > 0)
            for (int s = 0; s < nStates; s++) {
                double jump = allJumps.get(s, t);
                if (jump > 0) {
                    // Add transition rate from src to dst
                    double currentRate = W.get(srcState, s);
                    W.set(srcState, s, currentRate + rateBase.get(t, 0) * jump);
                }
            }
        }

        // Make W a proper generator (row sums = 0)
        for (int i = 0; i < nStates; i++) {
            double rowSum = 0;
            for (int j = 0; j < nStates; j++) {
                if (i != j) {
                    rowSum += W.get(i, j);
                }
            }
            W.set(i, i, -rowSum);
        }

        // Step 2: Identify immediate vs timed states
        Set<Integer> immStatesSet = new HashSet<Integer>();
        for (int idx : immIdx) {
            immStatesSet.add((int) eventIdx.get(idx, 0));
        }

        List<Integer> timedStatesList = new ArrayList<Integer>();
        for (int i = 0; i < nStates; i++) {
            if (!immStatesSet.contains(i)) {
                timedStatesList.add(i);
            }
        }

        // Check that we have some timed states remaining
        if (timedStatesList.isEmpty()) {
            throw new RuntimeException("All states have immediate transitions - cannot eliminate");
        }

        // Convert to array for ctmc_stochcomp
        int[] timedStates = new int[timedStatesList.size()];
        for (int i = 0; i < timedStatesList.size(); i++) {
            timedStates[i] = timedStatesList.get(i);
        }

        // Step 3: Apply stochastic complement
        // Convert int[] to MutableList<Double?>
        java.util.ArrayList<Double> timedStatesList2 = new java.util.ArrayList<>();
        for (int state : timedStates) {
            timedStatesList2.add((double) state);
        }
        SolverCTMC.StochCompResult result = ctmc_stochcomp(W, timedStatesList2);
        Matrix WReduced = result.S;

        // Step 4: Convert reduced generator back to jump/rate representation
        JumpRateResult jumpRateResult = generatorToJumps(WReduced);

        // Step 5: Expand to original state space dimension
        // This ensures compatibility with ode_rates_closing which uses original q_indices
        Matrix allJumpsReduced = new Matrix(nStates, jumpRateResult.allJumps.getNumCols());
        for (int t = 0; t < jumpRateResult.allJumps.getNumCols(); t++) {
            for (int s = 0; s < timedStates.length; s++) {
                double jump = jumpRateResult.allJumps.get(s, t);
                if (jump != 0) {
                    allJumpsReduced.set(timedStates[s], t, jump);
                }
            }
        }

        // Map eventIdx from reduced to original state space
        Matrix eventIdxReduced = new Matrix(jumpRateResult.eventIdx.getNumRows(), 1);
        for (int i = 0; i < jumpRateResult.eventIdx.getNumRows(); i++) {
            int localIdx = (int) jumpRateResult.eventIdx.get(i, 0);
            eventIdxReduced.set(i, 0, timedStates[localIdx]);
        }

        // Report reduction
        if (options.verbose == VerboseLevel.STD || options.verbose == VerboseLevel.DEBUG) {
            double reductionPct = 100.0 * (1.0 - (double) timedStates.length / nStates);
            System.out.println(String.format("State space reduced from %d to %d states (%.1f%% reduction)",
                    nStates, timedStates.length, reductionPct));
        }

        return new EliminationResult(allJumpsReduced, jumpRateResult.rateBase,
                                     eventIdxReduced, timedStates);
    }

    /**
     * Result of generator to jumps conversion
     */
    private static class JumpRateResult {
        Matrix allJumps;
        Matrix rateBase;
        Matrix eventIdx;

        JumpRateResult(Matrix allJumps, Matrix rateBase, Matrix eventIdx) {
            this.allJumps = allJumps;
            this.rateBase = rateBase;
            this.eventIdx = eventIdx;
        }
    }

    /**
     * Convert infinitesimal generator matrix to jump/rate representation
     *
     * @param W [n_states x n_states] infinitesimal generator matrix
     * @return JumpRateResult containing jump matrix, rate vector, and event indices
     */
    private static JumpRateResult generatorToJumps(Matrix W) {
        int nStates = (int) W.getNumRows();

        // Find all non-zero off-diagonal elements
        List<Integer> srcStates = new ArrayList<Integer>();
        List<Integer> dstStates = new ArrayList<Integer>();
        List<Double> rates = new ArrayList<Double>();

        for (int i = 0; i < nStates; i++) {
            for (int j = 0; j < nStates; j++) {
                if (i != j && W.get(i, j) > 0) {
                    srcStates.add(i);
                    dstStates.add(j);
                    rates.add(W.get(i, j));
                }
            }
        }

        int nTransitions = srcStates.size();

        // Pre-allocate outputs
        Matrix allJumps = new Matrix(nStates, nTransitions);
        Matrix rateBase = new Matrix(nTransitions, 1);
        Matrix eventIdx = new Matrix(nTransitions, 1);

        // Construct jump vectors
        for (int t = 0; t < nTransitions; t++) {
            int src = srcStates.get(t);
            int dst = dstStates.get(t);

            // Jump: decrease at source, increase at destination
            allJumps.set(src, t, -1);
            allJumps.set(dst, t, +1);

            rateBase.set(t, 0, rates.get(t));
            eventIdx.set(t, 0, src);
        }

        return new JumpRateResult(allJumps, rateBase, eventIdx);
    }
}
