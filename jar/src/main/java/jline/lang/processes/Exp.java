/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.GlobalConstants;
import jline.util.RandomManager;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.io.Serializable;
import java.util.Random;

import static jline.api.mam.Map_sample.map_sample;
import static jline.util.Maths.max;
import static jline.util.Maths.min;
import static org.apache.commons.math3.util.FastMath.exp;

/**
 * An exponential distribution.
 */

public class Exp extends Markovian implements Serializable {
    /**
     * Creates an exponential distribution with the specified rate parameter.
     * 
     * @param lambda the rate parameter (must be positive)
     */
    public Exp(double lambda) {
        super("Exp", 1);
        nPhases = 1;
        this.setParam(1, "lambda", lambda);
        Matrix D0 = new Matrix(1, 1, 1);
        Matrix D1 = new Matrix(1, 1, 1);
        D0.set(0, 0, -lambda);
        D1.set(0, 0, lambda);
        MatrixCell rep = new MatrixCell();
        rep.set(0, D0);
        rep.set(1, D1);
        setProcess(rep);
    }

    /**
     * Creates an exponential distribution fitted to the specified mean.
     * 
     * @param MEAN the desired mean value
     * @return an exponential distribution with the specified mean
     */
    public static Exp fitMean(double MEAN) {
        return new Exp(min(GlobalConstants.Immediate, max(GlobalConstants.Zero, 1 / MEAN)));
    }

    /**
     * Creates an exponential distribution fitted to the specified rate.
     * 
     * @param RATE the desired rate value
     * @return an exponential distribution with the specified rate
     */
    public static Exp fitRate(double RATE) {
        return new Exp(min(GlobalConstants.Immediate, max(GlobalConstants.Zero, RATE)));
    }

    /**
     * Evaluates the cumulative distribution function at the given point.
     * 
     * @param t the point at which to evaluate the CDF
     * @return the CDF value F(t) = 1 - e^(-λt)
     */
    public double evalCDF(double t) {
        double lambda = (double) this.getParam(1).getValue();
        return 1 - exp(-lambda * t);
    }

    /**
     * Evaluates the Laplace-Stieltjes transform at the given point.
     * 
     * @param s the transform variable
     * @return the LST value λ/(λ+s)
     */
    public double evalLST(double s) {
        double lambda = (double) this.getParam(1).getValue();
        return (lambda / (lambda + s));
    }

    /**
     * Sets the mean by rescaling the process AND the rate parameter.
     *
     * Markovian.setMean rewrites the MAP representation only, but getRate,
     * getMean, evalCDF, evalLST and Network.refreshStruct all read the rate
     * parameter, so a rescaled Exp would otherwise keep reporting its old rate
     * and every estimator's model update would be silently lost.
     *
     * @param newMean the new mean value
     */
    @Override
    public void setMean(double newMean) {
        super.setMean(newMean);
        double lambda = min(GlobalConstants.Immediate,
                max(GlobalConstants.Zero, newMean == 0.0 ? GlobalConstants.Immediate : 1 / newMean));
        this.setParam(1, "lambda", lambda);
    }

    /**
     * Gets the mean of this exponential distribution.
     * 
     * @return the mean value (1/λ)
     */
    public double getMean() {
        return 1 / getRate();
    }

    /**
     * Gets the number of phases in this distribution.
     * Exponential distribution has a single phase.
     * 
     * @return always returns 1
     */
    public long getNumberOfPhases() {
        return 1;
    }

    /**
     * Gets the rate parameter of this exponential distribution.
     * 
     * @return the rate parameter λ
     */
    public double getRate() {
        return (double) this.getParam(1).getValue();
    }

    /**
     * Gets the squared coefficient of variation.
     * For exponential distribution, SCV = 1.
     * 
     * @return always returns 1
     */
    public double getSCV() {
        return 1;
    }

    /**
     * Gets the skewness of this exponential distribution.
     * 
     * @return always returns 2
     */
    public double getSkewness() {
        return 2;
    }

    /**
     * Gets the variance of this exponential distribution.
     * 
     * @return the variance (1/λ²)
     */
    public double getVar() {
        return 1 / (Math.pow(getRate(), 2));
    }

    @Override
    public double[] sample(int n) {
        return this.sample(n, RandomManager.getThreadRandomAsRandom());
    }

    @Override
    public double[] sample(int n, Random rand) {
        return map_sample(D(0), D(1), n, rand);
    }

    public String toString() {
        return String.format("jline.Exp(%f)", this.getRate());
    }

    /**
     * Updates the rate parameter of this exponential distribution.
     * 
     * @param rate the new rate parameter
     */
    public void updateRate(double rate) {
        this.setParam(1, "lambda", rate);
        this.mean = 1.0 / rate;
        this.immediate = 1.0 / rate < GlobalConstants.FineTol;
    }

    // =================== PROPERTY ALIASES ===================
    
    /**
     * Property alias for getNumberOfPhases
     */
    public long numberOfPhases() {
        return getNumberOfPhases();
    }
    
    /**
     * Property alias for getNumberOfPhases
     */
    public long numPhases() {
        return getNumberOfPhases();
    }
    
    /**
     * Property alias for getMean
     */
    public double mean() {
        return getMean();
    }
    
    /**
     * Property alias for getRate
     */
    public double rate() {
        return getRate();
    }
    
    /**
     * Property alias for getSCV
     */
    public double scv() {
        return getSCV();
    }
    
    /**
     * Property alias for getSkewness
     */
    public double skewness() {
        return getSkewness();
    }
    
    /**
     * Property alias for getVar
     */
    public double var() {
        return getVar();
    }
}
