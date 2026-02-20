/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */


package jline.solvers;

import jline.util.matrix.Matrix;

/**
 * Container for storing performance metrics computed by queueing network solvers.
 * <p>
 * This class holds both steady-state and transient performance results including
 * queue lengths, utilizations, response times, throughputs, and arrival rates.
 * Results are organized as matrices with dimensions [stations x classes] for
 * station-level metrics and [1 x classes] for system-level metrics.
 * <p>
 * For solver-specific results, create separate classes (e.g., SolverFluidResult).
 *
 * @see NetworkSolver
 * @see Matrix
 */
public class SolverResult {

    /**
     * The solution method used to compute these results
     */
    public String method;

    /**
     * The name of the solver that computed these results
     */
    public String solver;

    /**
     * Number of iterations performed by iterative solvers
     */
    public int iter;

    /**
     * Mean queue lengths [stations x classes]
     */
    public Matrix QN;

    /**
     * Mean server utilizations [stations x classes]
     */
    public Matrix UN;

    /**
     * Mean response times [stations x classes]
     */
    public Matrix RN;

    /**
     * Mean throughputs [stations x classes]
     */
    public Matrix TN;

    /**
     * Mean arrival rates [stations x classes]
     */
    public Matrix AN;

    /**
     * Mean residence times [stations x classes]
     */
    public Matrix WN;

    /**
     * Mean tardiness [stations x classes]
     */
    public Matrix TardN;

    /**
     * Mean system tardiness [1 x classes]
     */
    public Matrix SysTardN;

    /**
     * Mean system response times [1 x chains]
     */
    public Matrix CN;

    /**
     * Mean system throughputs [1 x chains]
     */
    public Matrix XN;

    /**
     * FCR mean queue lengths (jobs waiting at ingress) [regions x classes]
     */
    public Matrix QNfcr;

    /**
     * FCR mean utilizations (capacity utilization) [regions x classes]
     */
    public Matrix UNfcr;

    /**
     * FCR mean response times (wait time at ingress) [regions x classes]
     */
    public Matrix RNfcr;

    /**
     * FCR mean throughputs (entry rate) [regions x classes]
     */
    public Matrix TNfcr;

    /**
     * FCR mean arrival rates (attempted entries) [regions x classes]
     */
    public Matrix ANfcr;

    /**
     * FCR mean residence times [regions x classes]
     */
    public Matrix WNfcr;

    /**
     * Transient queue lengths [time_points][stations x classes].
     * Time steps are stored separately in matrix 't' for efficiency.
     */
    public Matrix[][] QNt;

    /**
     * Transient response times [time_points][stations x classes]
     */
    public Matrix[][] RNt;

    /**
     * Transient utilizations [time_points][stations x classes]
     */
    public Matrix[][] UNt;

    /**
     * Transient throughputs [time_points][stations x classes]
     */
    public Matrix[][] TNt;

    /**
     * Transient system throughputs [time_points][1 x chains]
     */
    public Matrix[][] XNt;

    /**
     * Transient system response times [time_points][1 x chains]
     */
    public Matrix[][] CNt;

    /**
     * Time points for transient analysis
     */
    public Matrix t;

    /**
     * Transient probability distributions over time
     */
    public Matrix pi_t;

    /**
     * Steady-state probability distribution
     */
    public Matrix SS;

    /**
     * Computation time in seconds
     */
    public double runtime;

    /**
     * Creates a deep copy of this SolverResult instance.
     * All matrix fields are cloned to ensure independence.
     *
     * @return a deep copy of this result object
     */
    public SolverResult deepCopy() {

        SolverResult clone = new SolverResult();

        clone.method = this.method;

        clone.QN = this.QN == null ? null : this.QN.copy();
        clone.UN = this.UN == null ? null : this.UN.copy();
        clone.RN = this.RN == null ? null : this.RN.copy();
        clone.WN = this.WN == null ? null : this.WN.copy();
        clone.TN = this.TN == null ? null : this.TN.copy();
        clone.AN = this.AN == null ? null : this.AN.copy();
        clone.TardN = this.TardN == null ? null : this.TardN.copy();
        clone.SysTardN = this.SysTardN == null ? null : this.SysTardN.copy();
        clone.CN = this.CN == null ? null : this.CN.copy();
        clone.XN = this.XN == null ? null : this.XN.copy();

        // Copy FCR metrics
        clone.QNfcr = this.QNfcr == null ? null : this.QNfcr.copy();
        clone.UNfcr = this.UNfcr == null ? null : this.UNfcr.copy();
        clone.RNfcr = this.RNfcr == null ? null : this.RNfcr.copy();
        clone.TNfcr = this.TNfcr == null ? null : this.TNfcr.copy();
        clone.ANfcr = this.ANfcr == null ? null : this.ANfcr.copy();
        clone.WNfcr = this.WNfcr == null ? null : this.WNfcr.copy();

        if (this.QNt == null && this.UNt == null && this.TNt == null) {
            clone.QNt = null;
            clone.UNt = null;
            clone.TNt = null;
        } else {
            clone.QNt = new Matrix[this.QNt.length][this.QNt[0].length];
            clone.UNt = new Matrix[this.UNt.length][this.UNt[0].length];
            clone.TNt = new Matrix[this.TNt.length][this.TNt[0].length];
            for (int i = 0; i < this.QNt.length; i++) {
                for (int j = 0; j < this.QNt[0].length; j++) {
                    clone.QNt[i][j] = this.QNt[i][j].copy();
                    clone.UNt[i][j] = this.UNt[i][j].copy();
                    clone.TNt[i][j] = this.TNt[i][j].copy();
                }
            }
        }

        clone.t = this.t == null ? null : this.t.copy();
        clone.runtime = this.runtime;

        return clone;
    }

    /**
     * Resets all stored results to null and runtime to zero.
     */
    public void reset() {
        this.QN = null;
        this.UN = null;
        this.RN = null;
        this.TN = null;
        this.AN = null;
        this.WN = null;
        this.TardN = null;
        this.SysTardN = null;
        this.CN = null;
        this.XN = null;

        // Reset FCR metrics
        this.QNfcr = null;
        this.UNfcr = null;
        this.RNfcr = null;
        this.TNfcr = null;
        this.ANfcr = null;
        this.WNfcr = null;

        this.QNt = null;
        this.UNt = null;
        this.TNt = null;
        this.t = null;
        this.pi_t = null;
        this.SS = null;
        this.runtime = 0.0;
    }
}
