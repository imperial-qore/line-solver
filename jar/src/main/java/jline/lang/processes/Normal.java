/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import static jline.GlobalConstants.Inf;
import static jline.GlobalConstants.NegInf;

import jline.GlobalConstants;
import jline.util.Pair;
import jline.util.RandomManager;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.io.Serializable;
import java.util.Random;

import static jline.util.Maths.max;
import static org.apache.commons.math3.util.FastMath.*;
import static org.apache.commons.math3.special.Erf.erf;

/**
 * A normal (Gaussian) distribution.
 */
public class Normal extends ContinuousDistribution implements Serializable {
    
    public Normal(double mu, double sigma) {
        super("Normal", 2, new Pair<Double, Double>(NegInf, Inf));
        this.setParam(1, "mu", mu);
        this.setParam(2, "sigma", sigma);
        this.mean = mu;
        this.immediate = false;
    }
    
    public static Normal fitMean(double mean) {
        return new Normal(mean, 1.0);
    }
    
    public static Normal fitMeanAndStd(double mean, double std) {
        return new Normal(mean, max(GlobalConstants.FineTol, std));
    }
    
    public static Normal fitMeanAndVar(double mean, double var) {
        return new Normal(mean, max(GlobalConstants.FineTol, sqrt(var)));
    }
    
    @Override
    public double evalCDF(double x) {
        double mu = (double) this.getParam(1).getValue();
        double sigma = (double) this.getParam(2).getValue();
        return 0.5 * (1 + erf((x - mu) / (sigma * sqrt(2))));
    }
    
    @Override
    public double evalLST(double s) {
        double mu = (double) this.getParam(1).getValue();
        double sigma = (double) this.getParam(2).getValue();
        return exp(mu * s + 0.5 * sigma * sigma * s * s);
    }
    
    public double evalPDF(double x) {
        double mu = (double) this.getParam(1).getValue();
        double sigma = (double) this.getParam(2).getValue();
        double diff = x - mu;
        return (1.0 / (sigma * sqrt(2 * PI))) * exp(-0.5 * (diff * diff) / (sigma * sigma));
    }
    
    @Override
    public double getMean() {
        return (double) this.getParam(1).getValue();
    }
    
    public double getStd() {
        return (double) this.getParam(2).getValue();
    }
    
    @Override
    public double getVar() {
        double sigma = (double) this.getParam(2).getValue();
        return sigma * sigma;
    }
    
    @Override
    public double getSCV() {
        double mu = (double) this.getParam(1).getValue();
        double sigma = (double) this.getParam(2).getValue();
        if (abs(mu) < GlobalConstants.FineTol) {
            return Inf;
        }
        return (sigma * sigma) / (mu * mu);
    }
    
    @Override
    public double getSkewness() {
        return 0.0; // Normal distribution has zero skewness
    }
    
    @Override
    public double[] sample(int n) {
        return this.sample(n, RandomManager.getThreadRandomAsRandom());
    }
    
    @Override
    public double[] sample(int n, Random random) {
        double mu = (double) this.getParam(1).getValue();
        double sigma = (double) this.getParam(2).getValue();
        double[] samples = new double[n];
        
        for (int i = 0; i < n; i++) {
            samples[i] = mu + sigma * random.nextGaussian();
        }
        
        return samples;
    }
    
    @Override
    public MatrixCell getProcess() {
        // Returns {mu, sigma} for normal distribution
        MatrixCell representation = new MatrixCell();
        representation.set(0, Matrix.singleton((Double) this.getParam(1).getValue()));  // mu
        representation.set(1, Matrix.singleton((Double) this.getParam(2).getValue()));  // sigma
        return representation;
    }

    @Override
    public String toString() {
        return String.format("jline.Normal(%f,%f)", getMean(), getStd());
    }
}