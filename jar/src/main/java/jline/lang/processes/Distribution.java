/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.GlobalConstants;
import jline.util.NamedParam;
import org.apache.commons.math3.complex.Complex;
import jline.util.Pair;
import jline.util.RandomManager;
import jline.lang.Copyable;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import static jline.io.InputOutput.line_error;
import static jline.io.InputOutput.mfilename;

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
        this.immediate = false; // MATLAB Distribution sets it here, subclasses override
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
     * The registry name this distribution is marked under by
     * {@link jline.lang.Network#getUsedLangFeatures()}.
     *
     * Separate from {@link #getName()} because that one also resolves the
     * {@link jline.lang.constant.ProcessType} and drives the JSON wire type: a
     * subclass whose feature name is more specific than its process type -- a
     * Trace, or a two-phase Coxian -- must be able to say so to the feature gate
     * without moving to a different process type. Defaults to getName(), so a
     * distribution that needs no distinction inherits the old behaviour.
     *
     * @return the FeatureSet entry naming this distribution
     */
    public String getFeatureName() {
        return getName();
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
     * Gets the number of parameters of this distribution, twin of the MATLAB
     * {@code getNumParams()}, which reports {@code length(params)}.
     *
     * @return the number of parameters
     */
    public int getNumParams() {
        return this.params.size();
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
        // MATLAB setNumParams re-initializes the params cell, so the two stay
        // consistent: getNumParams there is length(params).
        this.params = new ArrayList<NamedParam>();
        for (int i = 0; i < num; i++) {
            this.params.add(new NamedParam("NULL_PARAM", null));
        }
    }

    /**
     * Sets a parameter value for this distribution.
     * 
     * @param id the parameter ID (1-based index)
     * @param name the parameter name
     * @param value the parameter value
     */
    public void setParam(int id, String name, Object value) {
        if (id < 1) {
            line_error(mfilename(new Object() {
            }), "Distribution parameter ids are 1-based, got " + id + ".");
        }
        if (name == null || name.trim().length() == 0) {
            line_error(mfilename(new Object() {
            }), "Distribution parameter " + id + " was set with an empty name.");
        }
        if (id >= this.params.size()) {
            int shortfall = (id - this.params.size());
            for (int i = 0; i < shortfall; i++) {
                this.params.add(new NamedParam("NULL_PARAM", null));
            }
            this.numParam = this.params.size(); // MATLAB getNumParams reads length(params)
        }
        this.params.set(id - 1, new NamedParam(name, value));
    }

    /**
     * Returns the parameter carrying a given name, twin of reading the MATLAB
     * {@code params} cell by its {@code paramName} field.
     *
     * @param name the parameter name
     * @return the parameter, or null when no parameter carries that name
     */
    public NamedParam getParam(String name) {
        for (NamedParam param : this.params) {
            if (param.getName() != null && param.getName().equals(name)) {
                return param;
            }
        }
        return null;
    }

    /**
     * Tests whether a parameter has been set at the given position, i.e. whether
     * it still holds the placeholder installed by the constructor.
     *
     * @param id the parameter id (1-based)
     * @return true when the parameter has been assigned
     */
    public boolean hasParam(int id) {
        if (id < 1 || id > this.params.size()) {
            return false;
        }
        NamedParam param = this.params.get(id - 1);
        return param != null && !"NULL_PARAM".equals(param.getName());
    }

    /**
     * Returns every parameter of this distribution, in declaration order. The
     * values are the stored objects, so a matrix-valued parameter (a PH
     * representation, a MAP block) is returned as it is held rather than
     * coerced to a scalar.
     *
     * @return the parameter list
     */
    public List<NamedParam> getParams() {
        return this.params;
    }

    /**
     * Evaluate the Laplace-Stieltjes Transform at s
     *
     * @param s the Laplace domain variable
     * @return the LST value at s
     */
    public abstract double evalLST(double s);

    /**
     * Evaluate the Laplace-Stieltjes Transform at a COMPLEX argument.
     *
     * A transform is evaluated off the real axis by everything that inverts it
     * or locates its roots: the Abate-Whitt Euler sum walks a vertical line, and
     * the matrix transform int exp(Ut) dF(t) is read off the spectrum of U,
     * which is complex in general. The real {@link #evalLST(double)} cannot
     * serve either.
     *
     * This default is the Riemann-Stieltjes sum over the law's own CDF with
     * complex exponentials: the weights are true probability increments, so it
     * is a proper measure for ANY law, including one with an atom or with no
     * density at all. Subclasses that own a closed form override it; the
     * phase-type families do so in {@link Markovian}.
     *
     * @param s the Laplace domain variable
     * @return the LST value at s
     */
    public Complex evalLST(Complex s) {
        // A REAL argument is what every real-line walk passes -- the G/M/1
        // sigma-root among them -- and there the exact evalLST(double) a
        // subclass already owns is authoritative. The Riemann-Stieltjes sum
        // below is a TRUNCATED, mass-renormalized quadrature, so reaching for it
        // ahead of a closed form silently replaced Pareto's and Replayer's exact
        // transforms with an approximation: the sigma-root landed 8.4e-4 off
        // MATLAB on Pareto, and on Replayer the distortion left NO sign change
        // to bracket, so sigma came back NaN and the solver fell through to the
        // G/G/1 KLB approximation without saying so. Delegating here reproduces
        // MATLAB to 12 significant figures on both.
        if (s.getImaginary() == 0.0) {
            try {
                return new Complex(evalLST(s.getReal()), 0.0);
            } catch (org.apache.commons.lang3.NotImplementedException e) {
                // This law has no closed form on the real line either, which is
                // exactly what the quadrature below exists for.
            }
        }
        final int nGrid = 2400;
        double hi = getMean() * 60.0;
        final double var = getVar();
        if (Double.isFinite(var) && var > 0.0) {
            hi = Math.max(hi, getMean() + 12.0 * Math.sqrt(var));
        }
        if (!Double.isFinite(hi) || hi <= 0.0) {
            return new Complex(1.0, 0.0);
        }
        final double step = hi / nGrid;
        Complex acc = new Complex(0.0, 0.0);
        double mass = 0.0;
        double prev = evalCDF(0.0);
        for (int i = 0; i < nGrid; i++) {
            final double right = (i + 1) * step;
            final double cur = evalCDF(right);
            final double w = cur - prev;
            prev = cur;
            if (w == 0.0) continue;
            mass += w;
            final double x = (i + 0.5) * step;
            final double mag = w * Math.exp(-s.getReal() * x);
            acc = acc.add(new Complex(mag * Math.cos(s.getImaginary() * x),
                    -mag * Math.sin(s.getImaginary() * x)));
        }
        // renormalize the truncated tail so the result is still a transform
        return mass > 0.0 ? acc.divide(mass) : new Complex(1.0, 0.0);
    }

    /**
     * Check if this distribution is Markovian (has a matrix process representation)
     * 
     * @return true if this is a Markovian distribution, false otherwise
     */
    public boolean isMarkovian() {
        return this instanceof Markovian;
    }

    // =================== PROPERTY ALIASES ===================
    
    /**
     * Property alias for getName
     */
    public String name() {
        return getName();
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
    
    /**
     * Property alias for getSupport
     */
    public Pair<Double, Double> support() {
        return getSupport();
    }
    
    /**
     * Property alias for getNumParams(int id)
     */
    public int numParams(int id) {
        return getNumParams(id);
    }
    
    /**
     * Property alias for getParam(int id)
     */
    public NamedParam param(int id) {
        return getParam(id);
    }

}
