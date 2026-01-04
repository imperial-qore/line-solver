/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.api.mc.Ctmc_solveKt;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.io.Serializable;

/**
 * A Markov-Modulated Deterministic Process (MMDP) for fluid queue modeling.
 *
 * Models fluid flow with deterministic rates modulated by a background Markov chain.
 * Suitable for arrival and service processes in Markovian fluid queues analyzed
 * by the mfq method of SolverFLD.
 *
 * The (Q, R) parameterization follows BUTools conventions:
 * - Q: Generator matrix of the modulating CTMC (row sums = 0)
 * - R: Diagonal matrix of deterministic rates per state
 *
 * MMDP is the deterministic analogue of MMPP:
 * - MMPP: Poisson arrival rates modulated by a Markov chain
 * - MMDP: Deterministic rates modulated by a Markov chain
 */
public class MMDP extends Markovian implements Serializable {

    /**
     * Creates an empty MMDP instance.
     */
    public MMDP() {
        super("MMDP", 2);
    }

    /**
     * Creates an MMDP with specified generator Q and rate matrix R.
     *
     * @param Q n×n generator matrix (row sums = 0)
     * @param R n×n diagonal matrix of deterministic rates, or n-vector
     */
    public MMDP(Matrix Q, Matrix R) {
        super("MMDP", 2);

        // Convert vector to diagonal matrix if needed
        Matrix R_diag;
        if (R.getNumCols() == 1) {
            // Column vector - convert to diagonal matrix
            int n = R.getNumRows();
            R_diag = new Matrix(n, n);
            for (int i = 0; i < n; i++) {
                R_diag.set(i, i, R.get(i, 0));
            }
        } else if (R.getNumRows() == 1) {
            // Row vector - convert to diagonal matrix
            int n = R.getNumCols();
            R_diag = new Matrix(n, n);
            for (int i = 0; i < n; i++) {
                R_diag.set(i, i, R.get(0, i));
            }
        } else {
            R_diag = R;
        }

        this.nPhases = Q.getNumRows();
        this.setParam(1, "Q", Q);
        this.setParam(2, "R", R_diag);

        MatrixCell rep = new MatrixCell();
        rep.set(0, Q);
        rep.set(1, R_diag);
        this.setProcess(rep);
    }

    /**
     * Returns the generator matrix Q.
     *
     * @return n×n generator matrix of the modulating CTMC
     */
    public Matrix Q() {
        return (Matrix) this.getParam(1).getValue();
    }

    /**
     * Returns the rate matrix R (diagonal).
     *
     * @return n×n diagonal matrix of deterministic rates
     */
    public Matrix R() {
        return (Matrix) this.getParam(2).getValue();
    }

    /**
     * Returns the rate vector (diagonal of R).
     *
     * @return n-vector of deterministic rates per phase
     */
    public Matrix r() {
        Matrix R = R();
        int n = R.getNumRows();
        Matrix r = new Matrix(n, 1);
        for (int i = 0; i < n; i++) {
            r.set(i, 0, R.get(i, i));
        }
        return r;
    }

    @Override
    public long getNumberOfPhases() {
        return Q().getNumRows();
    }

    /**
     * Computes the stationary mean rate.
     *
     * For MMDP: E[r] = π * diag(R), where π is the stationary distribution
     * of the modulating CTMC with generator Q.
     *
     * @return stationary mean deterministic rate
     */
    public double getMeanRate() {
        Matrix Q = Q();
        Matrix r = r();
        Matrix pi = Ctmc_solveKt.ctmc_solve(Q);
        return pi.mult(r).get(0, 0);
    }

    @Override
    public double getMean() {
        double rate = getMeanRate();
        if (rate > 0) {
            return 1.0 / rate;
        }
        return Double.POSITIVE_INFINITY;
    }

    @Override
    public double getRate() {
        return getMeanRate();
    }

    /**
     * Computes the squared coefficient of variation of rates.
     *
     * Computes Var[r]/E[r]^2 where expectation is over the
     * stationary distribution of the modulating CTMC.
     *
     * @return squared coefficient of variation
     */
    @Override
    public double getSCV() {
        Matrix Q = Q();
        Matrix r = r();
        Matrix pi = Ctmc_solveKt.ctmc_solve(Q);

        double meanRate = pi.mult(r).get(0, 0);

        // Compute E[r^2]
        Matrix r2 = new Matrix(r.getNumRows(), 1);
        for (int i = 0; i < r.getNumRows(); i++) {
            r2.set(i, 0, r.get(i, 0) * r.get(i, 0));
        }
        double meanRateSq = pi.mult(r2).get(0, 0);

        double varRate = meanRateSq - meanRate * meanRate;

        if (meanRate > 0) {
            return varRate / (meanRate * meanRate);
        }
        return Double.POSITIVE_INFINITY;
    }

    @Override
    public MatrixCell getProcess() {
        MatrixCell res = new MatrixCell();
        res.set(0, Q());
        res.set(1, R());
        return res;
    }

    /**
     * Converts a MAP to MMDP (deterministic representation).
     *
     * Converts a Markovian Arrival Process to a Markov-Modulated
     * Deterministic Process by extracting the full generator and
     * using row sums of D1 as the deterministic rates.
     *
     * @param map MAP object to convert
     * @return MMDP representation of the MAP
     */
    public static MMDP fromMAP(MAP map) {
        Matrix D0 = map.D(0);
        Matrix D1 = map.D(1);
        Matrix Q = D0.add(1.0, D1);

        // Row sums of D1 as diagonal rate matrix
        int n = D1.getNumRows();
        Matrix R = new Matrix(n, n);
        for (int i = 0; i < n; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < D1.getNumCols(); j++) {
                rowSum += D1.get(i, j);
            }
            R.set(i, i, rowSum);
        }

        return new MMDP(Q, R);
    }

    /**
     * Creates an MMDP from MMPP2 parameters.
     *
     * Creates a 2-state MMDP using the same parameterization as MMPP2.
     *
     * @param r0 rate in state 0
     * @param r1 rate in state 1
     * @param sigma0 transition rate from state 0 to state 1
     * @param sigma1 transition rate from state 1 to state 0
     * @return 2-state MMDP
     */
    public static MMDP fromMMPP2(double r0, double r1, double sigma0, double sigma1) {
        Matrix Q = new Matrix(2, 2);
        Q.set(0, 0, -sigma0);
        Q.set(0, 1, sigma0);
        Q.set(1, 0, sigma1);
        Q.set(1, 1, -sigma1);

        Matrix R = new Matrix(2, 2);
        R.set(0, 0, r0);
        R.set(1, 1, r1);

        return new MMDP(Q, R);
    }

    /**
     * Checks if the given (Q, R) matrices define a valid MMDP.
     *
     * Requirements:
     * - Q must be a valid generator (square, row sums = 0, proper signs)
     * - R must be diagonal with non-negative entries
     * - Q and R must have compatible dimensions
     *
     * @param Q generator matrix
     * @param R rate matrix
     * @return true if valid MMDP, false otherwise
     */
    public static boolean isFeasible(Matrix Q, Matrix R) {
        double tol = 1e-10;

        int n = Q.getNumRows();
        int m = Q.getNumCols();

        // Q must be square
        if (n != m) {
            return false;
        }

        // Q must be a valid generator
        for (int i = 0; i < n; i++) {
            // Diagonal non-positive
            if (Q.get(i, i) > tol) {
                return false;
            }
            // Off-diagonal non-negative
            for (int j = 0; j < n; j++) {
                if (i != j && Q.get(i, j) < -tol) {
                    return false;
                }
            }
            // Row sums = 0
            double rowSum = 0.0;
            for (int j = 0; j < n; j++) {
                rowSum += Q.get(i, j);
            }
            if (Math.abs(rowSum) > tol) {
                return false;
            }
        }

        // R must be n×n
        if (R.getNumRows() != n || R.getNumCols() != n) {
            return false;
        }

        // R must be diagonal with non-negative entries
        for (int i = 0; i < n; i++) {
            if (R.get(i, i) < -tol) {
                return false;
            }
            for (int j = 0; j < n; j++) {
                if (i != j && Math.abs(R.get(i, j)) > tol) {
                    return false;
                }
            }
        }

        return true;
    }

    @Override
    public String toString() {
        Matrix Q = Q();
        double meanRate = getMeanRate();
        double scv = getSCV();
        return String.format("MMDP(%dx%d, meanRate=%.6f, scv=%.6f)",
                           Q.getNumRows(), Q.getNumCols(), meanRate, scv);
    }
}
