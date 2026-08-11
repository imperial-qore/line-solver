/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers;

/**
 * Container for organized performance metric handles computed by queueing network solvers.
 * 
 * <p>SolverAvgHandles groups all the average performance metric handles into a single
 * container for convenient access and management. Each handle corresponds to a different
 * type of performance measure and provides structured access to metrics by station
 * and job class.</p>
 * 
 * <p>Standard performance metrics included:
 * <ul>
 *   <li><strong>Q</strong> - Queue lengths (mean number of jobs)</li>
 *   <li><strong>U</strong> - Utilization (fraction of time busy)</li>
 *   <li><strong>R</strong> - Response times (time from arrival to departure)</li>
 *   <li><strong>W</strong> - Residence times (time spent at specific stations)</li>
 *   <li><strong>T</strong> - Throughputs (job completion rates)</li>
 *   <li><strong>A</strong> - Arrival rates (job arrival rates)</li>
 * </ul>
 * </p>
 * 
 * @see AvgHandle
 * @see SolverTranHandles
 * @since 1.0
 */
public class SolverAvgHandles {

    /** Queue length handles [stations x classes] */
    public AvgHandle Q;
    
    /** Utilization handles [stations x classes] */
    public AvgHandle U;
    
    /** Response time handles [stations x classes] */
    public AvgHandle R;
    
    /** Residence time handles [stations x classes] */
    public AvgHandle W;
    
    /** Throughput handles [stations x classes] */
    public AvgHandle T;
    
    /** Arrival rate handles [stations x classes] */
    public AvgHandle A;

    /** Tardiness handles [stations x classes] */
    public AvgHandle Tard;

    /** System tardiness handles [1 x classes] */
    public AvgHandle SysTard;

    /** Transient queue length handles */
    public AvgHandle Qt;

    /** Transient utilization handles */
    public AvgHandle Ut;

    /** Transient throughput handles */
    public AvgHandle Tt;

    /**
     * Constructs a new SolverAvgHandles container with the specified metric handles.
     *
     * @param Q queue length handles
     * @param U utilization handles
     * @param R response time handles
     * @param W residence time handles
     * @param T throughput handles
     * @param A arrival rate handles
     */
    public SolverAvgHandles(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        this.Q = Q;
        this.U = U;
        this.R = R;
        this.T = T;
        this.W = W;
        this.A = A;
        this.Tard = null;
        this.SysTard = null;
    }

    /**
     * Constructs a new SolverAvgHandles container with all metric handles including tardiness.
     *
     * @param Q queue length handles
     * @param U utilization handles
     * @param R response time handles
     * @param W residence time handles
     * @param T throughput handles
     * @param A arrival rate handles
     * @param Tard tardiness handles
     * @param SysTard system tardiness handles
     */
    public SolverAvgHandles(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, AvgHandle Tard, AvgHandle SysTard) {
        this.Q = Q;
        this.U = U;
        this.R = R;
        this.T = T;
        this.W = W;
        this.A = A;
        this.Tard = Tard;
        this.SysTard = SysTard;
    }

    public AvgHandle getAvgArvRHandles() {
        return A;
    }

    public AvgHandle getAvgQLenHandles() {
        return Q;
    }

    public AvgHandle getAvgResidTHandles() {
        return W;
    }

    public AvgHandle getAvgRespTHandles() {
        return R;
    }

    public AvgHandle getAvgTputHandles() {
        return T;
    }

    public AvgHandle getAvgUtilHandles() {
        return U;
    }

}
