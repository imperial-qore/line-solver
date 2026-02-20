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
 * A discrete distribution that samples uniformly among a set of elements.
 */

public class DiscreteUniform extends DiscreteDistribution implements Serializable {
    public DiscreteUniform(double minVal, double maxVal) {
        // Constructs a uniform distribution with specified minimum and maximum values
        super("DiscreteUniform", 2, new Pair<Double, Double>(minVal, maxVal));
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
            return 0;
        } else {
            return (Math.floor(t) - minVal + 1) / (maxVal - minVal + 1);
        }
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
        return (Math.pow(maxVal - minVal + 1, 2) - 1) / 12.0;
    }

    /**
     * Evaluates the Laplace-Stieltjes Transform at s.
     * For DiscreteUniform(a, b), LST(s) = (e^(-as) - e^(-(b+1)s)) / ((b-a+1)(1 - e^(-s)))
     *
     * @param s the Laplace domain variable
     * @return the LST value at s
     */
    public double evalLST(double s) {
        double a = (double) this.getParam(1).getValue();
        double b = (double) this.getParam(2).getValue();
        
        if (Math.abs(s) < 1e-10) {
            // Handle s ≈ 0 case to avoid division by zero
            return 1.0;
        }
        
        double e_neg_s = FastMath.exp(-s);
        double numerator = FastMath.exp(-a * s) - FastMath.exp(-(b + 1) * s);
        double denominator = (b - a + 1) * (1 - e_neg_s);
        
        return numerator / denominator;
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

    @Override
    public double[] sample(int n, Random random) {
        double minVal = (double) this.getParam(1).getValue();
        double maxVal = (double) this.getParam(2).getValue();
        double[] samples = new double[(int) n];
        for (int i = 0; i < n; i++) {
            double randomValue = FastMath.round(minVal + (maxVal - minVal) * random.nextDouble());
            samples[i] = randomValue;
        }
        return samples;
    }

    public MatrixCell getProcess() {
        MatrixCell representation = new MatrixCell();
        representation.set(0, Matrix.singleton(getMean()));
        representation.set(1, Matrix.singleton(getSCV()));
        return representation;
    }

}

