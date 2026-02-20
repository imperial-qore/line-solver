/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.util.Pair;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.apache.commons.math3.util.FastMath;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

/**
 * A Deterministic distribution taking a single constant value.
 */
public class Det extends Distribution implements Serializable {
    protected Map<Integer, Matrix> process; // <0, D0>, <1, D1>, <2, D2>

    /**
     * Creates a deterministic distribution with the specified constant value.
     * 
     * @param t the constant value for this deterministic distribution
     */
    public Det(double t) {
        super("Det", 1, new Pair<Double, Double>(t, t));
        this.setParam(1, "t", t);

        Matrix D0 = Matrix.singleton(-t);
        Matrix D1 = Matrix.singleton(t);
        Map<Integer, Matrix> res = new HashMap<Integer, Matrix>() {{
            put(1, D0);
            put(2, D1);
        }};
        setProcess(res);
    }

    /**
     * Creates a deterministic distribution with the specified mean value.
     * This is equivalent to new Det(mean).
     *
     * @param mean the mean value for the deterministic distribution
     * @return a new Det distribution with the specified mean
     */
    public static Det fitMean(double mean) {
        return new Det(mean);
    }

    /**
     * Evaluates the cumulative distribution function at the given point.
     * Returns 0 if t < value, 1 if t >= value.
     * 
     * @param t the point at which to evaluate the CDF
     * @return 0 or 1 depending on whether t has reached the deterministic value
     */
    public double evalCDF(double t) {
        // Evaluate the cumulative distribution function at t
        if (t < (double) this.getParam(1).getValue()) {
            return 0;
        } else {
            return 1;
        }
    }

    /**
     * Evaluates the Laplace-Stieltjes transform at the given point.
     * 
     * @param s the transform variable
     * @return e^(-s*t) where t is the deterministic value
     */
    public double evalLST(double s) {
        // Evaluate the Laplace-Stieltjes transform of the distribution function at t
        return FastMath.exp(-s * (double) this.getParam(1).getValue());
    }

    /**
     * Gets the mean of this deterministic distribution.
     * 
     * @return the constant value of the distribution
     */
    public double getMean() {
        return (double) this.getParam(1).getValue();
    }

    /**
     * Sets the mean (constant value) of this deterministic distribution.
     * 
     * @param t the new constant value
     */
    public void setMean(double t) {
        Matrix D0 = Matrix.singleton(-t);
        Matrix D1 = Matrix.singleton(t);
        Map<Integer, Matrix> res = new HashMap<Integer, Matrix>() {{
            put(1, D0);
            put(2, D1);
        }};
        setProcess(res);
    }

    /**
     * Gets the rate (inverse of the constant value).
     * 
     * @return 1/t where t is the deterministic value
     */
    public double getMu() {
        return 1.0 / (double) this.getParam(1).getValue();
    }

    /**
     * Gets the phase parameter.
     * 
     * @return always returns 1 for deterministic distribution
     */
    public double getPhi() {
        return 1;
    }

    /**
     * Gets the process representation with actual distribution parameters.
     * Returns MatrixCell{[t]} for deterministic value.
     *
     * @return MatrixCell containing the distribution parameter
     */
    public MatrixCell getRepresentation() {
        MatrixCell proc = new MatrixCell();
        proc.set(0, Matrix.singleton((Double) this.getParam(1).getValue()));  // t
        return proc;
    }

    /**
     * Gets the matrix representation of this deterministic process (for PH compatibility).
     *
     * @return map containing the process matrices
     */
    public Map<Integer, Matrix> getProcess() {
        return process;
    }

    /**
     * Sets the matrix representation of this deterministic process.
     *
     * @param proc the process matrices to set
     */
    public void setProcess(Map<Integer, Matrix> proc) {
        process = proc;
    }

    /**
     * Gets the rate of this deterministic distribution.
     * 
     * @return 1/t where t is the deterministic value
     */
    public double getRate() {
        return 1.0 / (double) this.getParam(1).getValue();
    }

    /**
     * Gets the squared coefficient of variation.
     * For deterministic distribution, SCV = 0 (no variance).
     * 
     * @return always returns 0
     */
    public double getSCV() {
        return 0;
    }

    /**
     * Gets the skewness of this deterministic distribution.
     * 
     * @return always returns 0 (no skewness)
     */
    @Override
    public double getSkewness() {
        return 0;
    }

    /**
     * Gets the variance of this deterministic distribution.
     * 
     * @return always returns 0 (no variance)
     */
    public double getVar() {
        return 0;
    }

    @Override
    public boolean isContinuous() {
        return true;
    }

    public boolean isDisabled() {
        return false;
    }

    @Override
    public boolean isDiscrete() {
        return true;
    }

    public boolean isImmediate() {
        return false;
    }

    @Override
    public double[] sample(int n) {
        return this.sample(n, null);
    }

    @Override
    public double[] sample(int n, Random random) {
        double[] ret_list = new double[(int) n];
        for (int i = 0; i < n; i++) {
            ret_list[i] = ((double) this.getParam(1).getValue());
        }
        return ret_list;
//        throw new RuntimeException("Not implemented");
    }

    // =================== KOTLIN-STYLE PROPERTY ALIASES ===================
    
    /**
     * Kotlin-style property alias for getMu()
     */
    public double mu() {
        return getMu();
    }
    
    /**
     * Kotlin-style property alias for getPhi()
     */
    public double phi() {
        return getPhi();
    }
    
    /**
     * Kotlin-style property alias for getProcess()
     */
    public Map<Integer, Matrix> process() {
        return getProcess();
    }
    
    /**
     * Kotlin-style property alias for getMean()
     */
    public double mean() {
        return getMean();
    }
    
    /**
     * Kotlin-style property alias for getRate()
     */
    public double rate() {
        return getRate();
    }
    
    /**
     * Kotlin-style property alias for getSCV()
     */
    public double scv() {
        return getSCV();
    }
    
    /**
     * Kotlin-style property alias for getSkewness()
     */
    public double skewness() {
        return getSkewness();
    }
    
    /**
     * Kotlin-style property alias for getVar()
     */
    public double var() {
        return getVar();
    }

}
