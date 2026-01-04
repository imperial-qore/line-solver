/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.util.Pair;
import jline.util.RandomManager;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.apache.commons.math3.util.FastMath;

import java.io.Serializable;
import java.util.Random;

/**
 * A continuous Uniform distribution
 */
public class Uniform extends ContinuousDistribution implements Serializable {
    public Uniform(double minVal, double maxVal) {
        // Constructs a uniform distribution with specified minimum and maximum values
        super("Uniform", 2, new Pair<Double, Double>(minVal, maxVal));
        this.setParam(1, "min", minVal);
        this.setParam(2, "max", maxVal);
    }

    public double evalCDF(double t) {
        // Evaluate the cumulative distribution function at t
        double minVal = (double) this.getParam(1).getValue();
        double maxVal = (double) this.getParam(2).getValue();
        if (t < minVal) {
            return 0;
        } else if (t > maxVal) {
            return 1;
        } else {
            return (t - minVal) / (maxVal - minVal);
        }
    }

    /**
     * Evaluates the probability density function (PDF) at the given point.
     * Uniform PDF: 1/(max-min) for min <= x <= max, 0 otherwise
     *
     * @param t the point at which to evaluate the PDF
     * @return the PDF value at point t
     */
    public double evalPDF(double t) {
        double minVal = (double) this.getParam(1).getValue();
        double maxVal = (double) this.getParam(2).getValue();
        if (t >= minVal && t <= maxVal) {
            return 1.0 / (maxVal - minVal);
        }
        return 0.0;
    }

    public double evalLST(double s) {
        // Evaluate the Laplace-Stieltjes transform of the distribution function at t
        double minVal = (double) this.getParam(1).getValue();
        double maxVal = (double) this.getParam(2).getValue();
        return (Math.exp(-s * minVal) - FastMath.exp(-s * maxVal)) / (s * (maxVal - minVal));
    }

    public double getMean() {
        // Get distribution mean
        return ((double) this.getParam(1).getValue() + (double) this.getParam(2).getValue()) / 2.0;
    }

    public double getRate() {
        return 1.0 / getMean();
    }

    public double getSCV() {
        // Get distribution squared coefficient of variation (SCV = variance / mean^2)
        return getVar() / FastMath.pow(getMean(), 2);
    }

    public double getSkewness() {
        return 0;
    }

    public double getVar() {
        double minVal = (double) this.getParam(1).getValue();
        double maxVal = (double) this.getParam(2).getValue();
        return FastMath.pow(maxVal - minVal, 2) / 12.0;
    }

    @Override
    public double[] sample(int n, Random random) {
        double minVal = (double) this.getParam(1).getValue();
        double maxVal = (double) this.getParam(2).getValue();
        double[] samples = new double[(int) n];
        for (int i = 0; i < n; i++) {
            double randomValue = minVal + (maxVal - minVal) * random.nextDouble();
            samples[i] = randomValue;
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
        // Returns {min, max} for uniform distribution
        MatrixCell representation = new MatrixCell();
        representation.set(0, Matrix.singleton((Double) this.getParam(1).getValue()));  // min
        representation.set(1, Matrix.singleton((Double) this.getParam(2).getValue()));  // max
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

