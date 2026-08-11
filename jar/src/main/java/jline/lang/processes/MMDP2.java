/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.io.Serializable;

/**
 * A 2-state Markov-Modulated Deterministic Process.
 *
 * Specialized MMDP with exactly 2 phases, using a convenient
 * parameterization analogous to MMPP2.
 *
 * Parameterization:
 * <ul>
 *   <li>r0, r1: Deterministic rates in states 0 and 1</li>
 *   <li>sigma0: Transition rate from state 0 to state 1</li>
 *   <li>sigma1: Transition rate from state 1 to state 0</li>
 * </ul>
 *
 * The generator matrix is:
 *   Q = [-sigma0, sigma0; sigma1, -sigma1]
 *
 * The rate matrix is:
 *   R = diag([r0, r1])
 */
public class MMDP2 extends MMDP implements Serializable {

    private double r0;
    private double r1;
    private double sigma0;
    private double sigma1;

    /**
     * Creates a 2-state Markov-Modulated Deterministic Process.
     *
     * @param r0     deterministic rate in state 0
     * @param r1     deterministic rate in state 1
     * @param sigma0 transition rate from state 0 to state 1
     * @param sigma1 transition rate from state 1 to state 0
     */
    public MMDP2(double r0, double r1, double sigma0, double sigma1) {
        super(buildQ(sigma0, sigma1), buildR(r0, r1));
        this.name = "MMDP2";
        this.r0 = r0;
        this.r1 = r1;
        this.sigma0 = sigma0;
        this.sigma1 = sigma1;

        // Override parameters with original scalar values for convenience
        this.setParam(1, "r0", r0);
        this.setParam(2, "r1", r1);
        this.setParam(3, "sigma0", sigma0);
        this.setParam(4, "sigma1", sigma1);
    }

    private static Matrix buildQ(double sigma0, double sigma1) {
        Matrix Q = new Matrix(2, 2);
        Q.set(0, 0, -sigma0);
        Q.set(0, 1, sigma0);
        Q.set(1, 0, sigma1);
        Q.set(1, 1, -sigma1);
        return Q;
    }

    private static Matrix buildR(double r0, double r1) {
        Matrix R = new Matrix(2, 2);
        R.set(0, 0, r0);
        R.set(1, 1, r1);
        return R;
    }

    /**
     * Returns the generator matrix Q (closed-form for 2 states).
     *
     * @return 2x2 generator matrix
     */
    @Override
    public Matrix Q() {
        return buildQ(sigma0, sigma1);
    }

    /**
     * Returns the rate matrix R (diagonal, closed-form for 2 states).
     *
     * @return 2x2 diagonal rate matrix
     */
    @Override
    public Matrix R() {
        return buildR(r0, r1);
    }

    /**
     * Returns the rate vector (diagonal of R).
     *
     * @return 2-vector of rates [r0; r1]
     */
    @Override
    public Matrix r() {
        Matrix r = new Matrix(2, 1);
        r.set(0, 0, r0);
        r.set(1, 0, r1);
        return r;
    }

    @Override
    public long getNumberOfPhases() {
        return 2;
    }

    /**
     * Computes the stationary mean rate (closed-form).
     *
     * For a 2-state MMDP, the mean rate has the closed form:
     *   E[r] = (r0*sigma1 + r1*sigma0) / (sigma0 + sigma1)
     *
     * @return stationary mean deterministic rate
     */
    @Override
    public double getMeanRate() {
        return (r0 * sigma1 + r1 * sigma0) / (sigma0 + sigma1);
    }

    /**
     * Computes the squared coefficient of variation (closed-form).
     *
     * For a 2-state MMDP, the SCV has a closed form based on
     * the variance of rates over the stationary distribution.
     *
     * @return squared coefficient of variation
     */
    @Override
    public double getSCV() {
        // Stationary probabilities
        double pi0 = sigma1 / (sigma0 + sigma1);
        double pi1 = sigma0 / (sigma0 + sigma1);

        // Mean and variance
        double meanRate = pi0 * r0 + pi1 * r1;
        double varRate = pi0 * r0 * r0 + pi1 * r1 * r1 - meanRate * meanRate;

        if (meanRate > 0) {
            return varRate / (meanRate * meanRate);
        }
        return Double.POSITIVE_INFINITY;
    }

    @Override
    public String toString() {
        return String.format("MMDP2(r0=%.6f, r1=%.6f, sigma0=%.6f, sigma1=%.6f)",
                r0, r1, sigma0, sigma1);
    }
}
