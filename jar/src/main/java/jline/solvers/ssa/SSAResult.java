/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.ssa;

import jline.lang.NetworkStruct;
import jline.lang.nodes.Station;
import jline.solvers.SolverResult;
import jline.util.matrix.Matrix;

import java.util.Map;

/**
 * Result container for Stochastic Simulation Algorithm (SSA) solver computations.
 *
 * <p>SSAResult extends the base SolverResult to include SSA-specific simulation
 * results and state-space information. SSA uses stochastic simulation to generate
 * sample paths and estimate performance metrics with statistical confidence intervals.</p>
 * 
 * <p>SSA analysis results include:
 * <ul>
 *   <li>Transient system state evolution over time</li>
 *   <li>State space representations for each station</li>
 *   <li>Synchronization matrices for multi-class systems</li>
 *   <li>Statistical estimates from simulation runs</li>
 * </ul>
 * </p>
 * 
 * @see jline.solvers.ssa.SolverSSA
 * @see SolverResult
 * @since 1.0
 */
public class SSAResult extends SolverResult {

    /** Transient system state matrices indexed by time points */
    public Map<Integer, Matrix> tranSysState;
    
    /** Transient synchronization matrix for multi-class coordination */
    public Matrix tranSync;
    
    /** Network structure used for the analysis */
    public NetworkStruct sn;
    
    /** State space matrices for each station in the network */
    public Map<Station, Matrix> space;

    /** Confidence interval half-widths for queue lengths [stations x classes] */
    public Matrix QNCI;

    /** Confidence interval half-widths for utilizations [stations x classes] */
    public Matrix UNCI;

    /** Confidence interval half-widths for response times [stations x classes] */
    public Matrix RNCI;

    /** Confidence interval half-widths for throughputs [stations x classes] */
    public Matrix TNCI;

    /** Confidence interval half-widths for arrival rates [stations x classes] */
    public Matrix ANCI;

    /** Confidence interval half-widths for residence times [stations x classes] */
    public Matrix WNCI;

    /**
     * Constructs an empty SSAResult to allow field population from different sources.
     */
    public SSAResult() {

    }

    /**
     * Constructs an SSAResult with the specified performance metrics and state information.
     * 
     * @param QN queue lengths [stations x classes]
     * @param UN utilizations [stations x classes]  
     * @param RN response times [stations x classes]
     * @param TN throughputs [stations x classes]
     * @param CN visit counts [stations x classes]
     * @param XN arrival rates [stations x classes]
     * @param tranSysState transient system states indexed by time
     * @param tranSync transient synchronization matrix
     * @param sn network structure information
     */
    public SSAResult(Matrix QN, Matrix UN, Matrix RN, Matrix TN, Matrix CN, Matrix XN,
                     Map<Integer, Matrix> tranSysState, Matrix tranSync, NetworkStruct sn) {
        this.QN = QN;
        this.UN = UN;
        this.RN = RN;
        this.TN = TN;
        this.CN = CN;
        this.XN = XN;
        this.tranSysState = tranSysState;
        this.tranSync = tranSync;
        this.sn = sn;
    }


}
