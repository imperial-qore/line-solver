/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.jmt;

import jline.lang.Metric;
import jline.solvers.SolverResult;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.List;

public class JMTResult extends SolverResult {

    // The following field is used only by the JMT solver.
    public List<Metric> metrics = new ArrayList<>();
    public double logNormConstAggr;

    /**
     * Cache node throughputs [nnodes x nclasses].
     * Stores the throughput values for cache nodes directly from simulation results.
     * Non-cache nodes will have 0 values. This is used to override the computed
     * values in getAvgNode() since cache throughputs from simulation are more accurate.
     */
    public Matrix cacheTN;

    /**
     * Cache node arrival rates [nnodes x nclasses].
     * Stores the arrival rate values for cache nodes directly from simulation results.
     */
    public Matrix cacheAN;

    /**
     * List of cache node indices in the network.
     */
    public List<Integer> cacheNodeIndices = new ArrayList<>();
    
    /**
     * Result class for transient probability analysis of aggregated states.
     */
    public static class TransientProbabilityResult {
        public Matrix Pi_t;      // Time-probability matrix [time, probability]
        public Matrix SSnode_a;  // Aggregated state space for the node
        
        public TransientProbabilityResult() {
            this.Pi_t = new Matrix(0, 0);
            this.SSnode_a = new Matrix(0, 0);
        }
        
        public TransientProbabilityResult(Matrix Pi_t, Matrix SSnode_a) {
            this.Pi_t = Pi_t;
            this.SSnode_a = SSnode_a;
        }
    }
}
