/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

// Copyright (c) 2012-2026, Imperial College London
// All rights reserved.

package jline.solvers;

import jline.util.matrix.Matrix;

/**
 * Result container for layered queueing network solver analysis.
 * 
 * <p>LayeredSolverResult extends the base SolverResult to include metrics specific
 * to layered queueing networks (LQNs). In addition to standard performance metrics,
 * it provides processor utilization and service time measurements essential for
 * analyzing software system models with nested service requests.</p>
 * 
 * <p>LQN-specific metrics:
 * <ul>
 *   <li>Processor utilization for hardware resource modeling</li>
 *   <li>Service time breakdowns across software layers</li>
 *   <li>End-to-end response time analysis</li>
 *   <li>Resource bottleneck identification</li>
 * </ul>
 * </p>
 * 
 * @see jline.solvers.ln.SolverLN
 * @see jline.solvers.lqns.SolverLQNS
 * @see SolverResult
 * @since 1.0
 */
public class LayeredSolverResult extends SolverResult {

    /** Mean processor utilization [processors x 1] */
    public Matrix PN;

    /** Mean service time [tasks x 1] */
    public Matrix SN;

    // Raw average metrics for detailed tables (aligned with MATLAB RawAvg structure)

    /** Raw utilization by node [nodes x 1] */
    public Matrix rawUtilization;

    /** Phase 1 utilization by node [nodes x 1] */
    public Matrix rawPhase1Utilization;

    /** Phase 2 utilization by node [nodes x 1] */
    public Matrix rawPhase2Utilization;

    /** Phase 1 service time by node [nodes x 1] */
    public Matrix rawPhase1ServiceTime;

    /** Phase 2 service time by node [nodes x 1] */
    public Matrix rawPhase2ServiceTime;

    /** Throughput by node [nodes x 1] */
    public Matrix rawThroughput;

    /** Processor waiting time by node [nodes x 1] */
    public Matrix rawProcWaiting;

    /** Processor utilization by node [nodes x 1] */
    public Matrix rawProcUtilization;

    /** Waiting time by call edge [calls x 1] */
    public Matrix rawEdgesWaiting;

    /**
     * Creates a deep copy of this layered solver result.
     * 
     * @return a new LayeredSolverResult with copied data
     */
    public LayeredSolverResult deepCopy() {

        LayeredSolverResult copy = new LayeredSolverResult();

        copy.method = this.method;

        copy.QN = this.QN.copy();
        copy.UN = this.UN.copy();
        copy.RN = this.RN.copy();
        copy.WN = this.WN.copy();
        copy.TN = this.TN.copy();
        copy.AN = this.AN.copy();
        copy.PN = this.PN.copy();
        copy.SN = this.SN.copy();
        copy.runtime = this.runtime;

        // Copy raw metrics if present
        if (this.rawUtilization != null) {
            copy.rawUtilization = this.rawUtilization.copy();
        }
        if (this.rawPhase1Utilization != null) {
            copy.rawPhase1Utilization = this.rawPhase1Utilization.copy();
        }
        if (this.rawPhase2Utilization != null) {
            copy.rawPhase2Utilization = this.rawPhase2Utilization.copy();
        }
        if (this.rawPhase1ServiceTime != null) {
            copy.rawPhase1ServiceTime = this.rawPhase1ServiceTime.copy();
        }
        if (this.rawPhase2ServiceTime != null) {
            copy.rawPhase2ServiceTime = this.rawPhase2ServiceTime.copy();
        }
        if (this.rawThroughput != null) {
            copy.rawThroughput = this.rawThroughput.copy();
        }
        if (this.rawProcWaiting != null) {
            copy.rawProcWaiting = this.rawProcWaiting.copy();
        }
        if (this.rawProcUtilization != null) {
            copy.rawProcUtilization = this.rawProcUtilization.copy();
        }
        if (this.rawEdgesWaiting != null) {
            copy.rawEdgesWaiting = this.rawEdgesWaiting.copy();
        }

        return copy;
    }
}
