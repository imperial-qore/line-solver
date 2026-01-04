/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import static jline.GlobalConstants.Inf;

import jline.util.Pair;
import jline.util.RandomManager;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.apache.commons.math3.distribution.LogNormalDistribution;
import org.apache.commons.math3.util.FastMath;

import java.io.Serializable;
import java.util.Random;

import static jline.io.InputOutputKt.line_error;
import static jline.io.InputOutputKt.mfilename;

/**
 * A Lognormal distribution.
 */

public class Lognormal extends ContinuousDistribution implements Serializable {
    public Lognormal(double mu, double sigma) {
        super("Lognormal", 2, new Pair<Double, Double>(0.0, Inf));
        if (sigma < 0) {
            line_error(mfilename(new Object() {
            }), "sigma parameter must be >= 0.0");
        }
        this.setParam(1, "mu", mu);
        this.setParam(2, "sigma", sigma);
    }

    public static Lognormal fitMeanAndSCV(double mean, double scv) {
        double c = FastMath.sqrt(scv);
        double mu = FastMath.log(mean / FastMath.sqrt(Math.pow(c, 2) + 1));
        double sigma = FastMath.sqrt(Math.log(Math.pow(c, 2) + 1));
        return new Lognormal(mu, sigma);
    }

    @Override
    public double evalCDF(double t) {
        double mu = (double) this.getParam(1).getValue();
        double sigma = (double) this.getParam(2).getValue();
        LogNormalDistribution logNormalDistribution = new LogNormalDistribution(mu, sigma);
        return logNormalDistribution.cumulativeProbability(t);
//        return 0.5 + 0.5 * UTIL.erf((Math.log(t) - mu) / (Math.sqrt(2) * sigma));
    }

    /**
     * Evaluates the probability density function (PDF) at the given point.
     *
     * @param t the point at which to evaluate the PDF
     * @return the PDF value at point t
     */
    public double evalPDF(double t) {
        if (t <= 0) {
            return 0.0;
        }
        double mu = (double) this.getParam(1).getValue();
        double sigma = (double) this.getParam(2).getValue();
        LogNormalDistribution logNormalDistribution = new LogNormalDistribution(mu, sigma);
        return logNormalDistribution.density(t);
    }

    @Override
    public double evalLST(double s) {
        // The LST of Lognormal distribution doesn't have a simple closed form.
        // We compute it numerically using the definition: E[e^(-sX)] = ∫₀^∞ e^(-sx) f(x) dx
        // where f(x) is the lognormal PDF
        
        double mu = (double) this.getParam(1).getValue();
        double sigma = (double) this.getParam(2).getValue();
        
        // Use numerical integration (adaptive method)
        // For practical purposes, we integrate up to a reasonable upper bound
        double upperBound = FastMath.exp(mu + 5 * sigma); // covers most of the distribution
        int n = 1000;
        double dx = upperBound / n;
        double sum = 0.0;
        
        for (int i = 1; i <= n; i++) {
            double x = i * dx;
            if (x > 0) {
                double logx = FastMath.log(x);
                double pdf = FastMath.exp(-FastMath.pow(logx - mu, 2) / (2 * sigma * sigma)) / 
                           (x * sigma * FastMath.sqrt(2 * Math.PI));
                sum += FastMath.exp(-s * x) * pdf;
            }
        }
        
        return sum * dx;
    }

    @Override
    public double getMean() {
        double mu = (double) this.getParam(1).getValue();
        double sigma = (double) this.getParam(2).getValue();
        return FastMath.exp(mu + FastMath.pow(sigma, 2) / 2.0);
    }

    @Override
    public double getRate() {
        return 1.0 / getMean();
    }

    @Override
    public double getSCV() {
        double mu = (double) this.getParam(1).getValue();
        double sigma = (double) this.getParam(2).getValue();
        double ex = FastMath.exp(mu + FastMath.pow(sigma, 2) / 2);
        double var = (Math.exp(Math.pow(sigma, 2)) - 1) * FastMath.exp(2 * mu + FastMath.pow(sigma, 2));
        return var / FastMath.pow(ex, 2);
    }

    @Override
    public double getSkewness() {
        double sigma = (double) this.getParam(2).getValue();
        return (Math.exp(Math.pow(sigma, 2)) + 2) * FastMath.sqrt(Math.exp(Math.pow(sigma, 2)) - 1);
    }

    @Override
    public double getVar() {
        double mu = (double) this.getParam(1).getValue();
        double sigma = (double) this.getParam(2).getValue();
        return (Math.exp(Math.pow(sigma, 2)) - 1) * FastMath.exp(2 * mu + FastMath.pow(sigma, 2));
    }

    @Override
    public double[] sample(int n, Random random) {
        double[] samples = new double[(int) n];
        double mu = (double) this.getParam(1).getValue();
        double sigma = (double) this.getParam(2).getValue();
        for (int i = 0; i < n; i++) {
            double z = random.nextGaussian(); // Standard Normal Distribution
            double value = FastMath.exp(mu + sigma * z); // Log-normal Distribution
            samples[i] = value;
        }
        return samples;
    }

    /**
     * Gets n samples from the distribution
     *
     * @param n - the number of samples
     * @return - n samples from the distribution
     */
    @Override
    public double[] sample(int n) {
        return this.sample(n, RandomManager.getThreadRandomAsRandom());
    }

    public MatrixCell getProcess() {
        // Returns {mu, sigma} for lognormal distribution
        MatrixCell representation = new MatrixCell();
        representation.set(0, Matrix.singleton((Double) this.getParam(1).getValue()));  // mu
        representation.set(1, Matrix.singleton((Double) this.getParam(2).getValue()));  // sigma
        return representation;
    }

    // =================== KOTLIN-STYLE PROPERTY ALIASES ===================
    
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
     * Kotlin-style property alias for getProcess()
     */
    public MatrixCell process() {
        return getProcess();
    }
}

