/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.GlobalConstants;
import jline.util.NamedParam;
import jline.util.Pair;
import jline.util.RandomManager;
import jline.lang.Copyable;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import static jline.io.InputOutputKt.line_error;
import static jline.io.InputOutputKt.mfilename;

/**
 * An abstract class of a general distribution
 */

abstract public class Distribution implements Copyable {
    protected double mean;
    protected boolean immediate;

    protected String name;
    protected int numParam;
    protected Pair<Double, Double> support;
    protected List<NamedParam> params;

    /**
     * Creates a new distribution with the specified characteristics.
     *
     * @param name     the name of this distribution type
     * @param numParam the number of parameters required
     * @param support  the support range [min, max] for this distribution
     */
    public Distribution(String name, int numParam, Pair<Double, Double> support) {
        this.params = new ArrayList<NamedParam>();

        this.name = name;
        this.numParam = numParam;
        this.support = support;
        for (int i = 0; i < this.numParam; i++) {
            this.params.add(new NamedParam("NULL_PARAM", null));
        }
    }

    /**
     * Evaluates the cumulative distribution function (CDF) at the given point.
     * 
     * @param t the point at which to evaluate the CDF
     * @return the CDF value at point t
     */
    public abstract double evalCDF(double t);

    /**
     * Evaluates the probability of the distribution falling within the given interval.
     * Computes P(t0 <= X <= t1) = F(t1) - F(t0).
     * 
     * @param t0 the lower bound of the interval
     * @param t1 the upper bound of the interval
     * @return the probability of falling within [t0, t1]
     * @throws RuntimeException if t1 < t0
     */
    public double evalProbInterval(double t0, double t1) {
        double Ft1 = 0.0;
        double Ft0 = 0.0;
        if (t1 >= t0) {
            Ft1 = evalCDF(t1);
            Ft0 = evalCDF(t0);
        } else {
            line_error(mfilename(new Object() {
            }), "CDF interval incorrectly specified (t1<t0)");
        }
        return Ft1 - Ft0;
    }

    /**
     * Gets the mean (expected value) of this distribution.
     * 
     * @return the mean value
     */
    public abstract double getMean();

    /**
     * Gets the name of this distribution type.
     * 
     * @return the distribution name
     */
    public String getName() {
        return name;
    }

    /**
     * Gets the number of parameters for this distribution.
     * 
     * @param id parameter identifier (currently unused)
     * @return the number of parameters
     */
    public int getNumParams(int id) {
        return numParam;
    }

    /**
     * Gets a parameter by its ID.
     * 
     * @param id the parameter ID (1-based index)
     * @return the named parameter at the specified position
     */
    public NamedParam getParam(int id) {
        return this.params.get(id - 1);
    }

    /**
     * Gets the rate of this distribution (inverse of mean).
     * 
     * @return the rate value (1/mean)
     */
    public double getRate() {
        return 1.0 / getMean();
    }

    /**
     * Gets the squared coefficient of variation (SCV) of this distribution.
     * SCV = Var(X) / E[X]^2.
     * 
     * @return the squared coefficient of variation
     */
    public abstract double getSCV();

    /**
     * Gets the skewness of this distribution.
     * Skewness measures the asymmetry of the probability distribution.
     * 
     * @return the skewness value
     */
    public abstract double getSkewness();

    /**
     * Gets the support range of this distribution.
     * 
     * @return a pair containing [min, max] values where the distribution is defined
     */
    public Pair<Double, Double> getSupport() {
        return support;
    }

    /**
     * Gets the variance of this distribution.
     * Computed as SCV * mean^2.
     * 
     * @return the variance
     */
    public double getVar() {
        return getSCV() * getMean() * getMean();
    }

    /**
     * Checks if this is a continuous distribution.
     * 
     * @return true if this is a continuous distribution, false otherwise
     */
    public boolean isContinuous() {
        return this instanceof ContinuousDistribution;
    }

    /**
     * Checks if this is a disabled distribution.
     * 
     * @return true if this is a disabled distribution, false otherwise
     */
    public boolean isDisabled() {
        return this instanceof Disabled;
    }

    /**
     * Checks if this is a discrete distribution.
     * 
     * @return true if this is a discrete distribution, false otherwise
     */
    public boolean isDiscrete() {
        return this instanceof DiscreteDistribution;
    }

    /**
     * Checks if this distribution has immediate (zero) service time.
     * 
     * @return true if the distribution is immediate or has mean < Zero threshold
     */
    public boolean isImmediate() {
        return this.immediate || getMean() < GlobalConstants.Zero;
    }

    /**
     * Generates random samples from this distribution using default random generator.
     *
     * @param n the number of samples to generate
     * @return array of random samples
     */
    public double[] sample(int n) {
        return this.sample(n, RandomManager.getThreadRandomAsRandom());
    }

    /**
     * Generates random samples from this distribution using the specified random generator.
     *
     * @param n      the number of samples to generate
     * @param random the random number generator to use
     * @return array of random samples
     */
    public abstract double[] sample(int n, Random random);

    /**
     * Sets the number of parameters for this distribution.
     * 
     * @param num the number of parameters
     */
    public void setNumParams(int num) {
        numParam = num;
    }

    /**
     * Sets a parameter value for this distribution.
     * 
     * @param id the parameter ID (1-based index)
     * @param name the parameter name
     * @param value the parameter value
     */
    public void setParam(int id, String name, Object value) {
        if (id >= this.params.size()) {
            int shortfall = (id - this.params.size());
            for (int i = 0; i < shortfall; i++) {
                this.params.add(new NamedParam("NULL_PARAM", null));
            }
        }
        this.params.set(id - 1, new NamedParam(name, value));
    }

    /*
     * ===================================================================================
     * MISSING METHODS FROM MATLAB DISTRIBUTION IMPLEMENTATION - NOT YET MIGRATED
     * ===================================================================================
     *
     * Based on analysis of /matlab/src/lang/processes/Distribution.m
     */

    // =================== MISSING ABSTRACT METHODS ===================
    
    /**
     * Evaluate the Laplace-Stieltjes Transform at s
     * 
     * @param s the Laplace domain variable
     * @return the LST value at s
     */
    public abstract double evalLST(double s);

    // =================== ENHANCED PARAMETER MANAGEMENT ===================
    // Enhanced parameter structure with:
    // - Named parameters with structured storage
    // - Better parameter validation
    // - Support for complex parameter types

    // =================== COPYABLE INTERFACE SUPPORT ===================
    // public Distribution copy()          // Create deep copy of distribution
    // public Distribution copyElement()   // Protected method for deep copying

    // =================== MATLAB-SPECIFIC INITIALIZATION ===================
    // Enhanced constructor with:
    // - Structured parameter initialization with paramName/paramValue fields
    // - Proper immediate flag handling
    // - MATLAB cell array parameter storage


    /**
     * Check if this distribution is Markovian (has a matrix process representation)
     * 
     * @return true if this is a Markovian distribution, false otherwise
     */
    public boolean isMarkovian() {
        return this instanceof Markovian;
    }

    // =================== KOTLIN-STYLE PROPERTY ALIASES ===================
    
    /**
     * Kotlin-style property alias for getName()
     */
    public String name() {
        return getName();
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
    
    /**
     * Kotlin-style property alias for getSupport()
     */
    public Pair<Double, Double> support() {
        return getSupport();
    }
    
    /**
     * Kotlin-style property alias for getNumParams(int id)
     */
    public int numParams(int id) {
        return getNumParams(id);
    }
    
    /**
     * Kotlin-style property alias for getParam(int id)
     */
    public NamedParam param(int id) {
        return getParam(id);
    }

}
