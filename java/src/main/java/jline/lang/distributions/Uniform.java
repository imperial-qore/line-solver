package jline.lang.distributions;

import jline.util.Matrix;
import jline.util.Pair;

import java.util.*;
import java.io.Serializable;

public class Uniform extends ContinuousDistribution implements Serializable {
    public Uniform(double minVal, double maxVal) {
        // Constructs a uniform distribution with specified minimum and maximum values
        super("Uniform", 2, new Pair<Double,Double>(minVal, maxVal));
        this.setParam(1, "min", minVal);
        this.setParam(2, "max", maxVal);
    }

    /**
     * Gets n samples from the distribution
     *
     * @param n - the number of samples
     * @return - n samples from the distribution
     */
    @Override
    public Matrix sample(long n) {
        return this.sample(n, new Random());
    }

    @Override
    public Matrix sample(long n, Random random) {
        double minVal = (double) this.getParam(1).getValue();
        double maxVal = (double) this.getParam(2).getValue();
        List<Double> samples = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            double randomValue = minVal + (maxVal - minVal) * random.nextDouble();
            samples.add(randomValue);
        }
        return new Matrix(samples);
    }

    public double getMean() {
        // Get distribution mean
        return ((double) this.getParam(1).getValue() + (double) this.getParam(2).getValue()) / 2.0;
    }

    public double getRate() {
        return 1 / getMean();
    }

    public double getSCV() {
        // Get distribution squared coefficient of variation (SCV = variance / mean^2)
        return getVar() / Math.pow(getMean(), 2);
    }

    public double getVar() {
        double minVal = (double) this.getParam(1).getValue();
        double maxVal = (double) this.getParam(2).getValue();
        return Math.pow(maxVal - minVal, 2) / 12.0;
    }

    public double getSkew() {
        return 0;
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
            return 1 / (maxVal - minVal);
        }
    }

    public double evalLST(double s) {
        // Evaluate the Laplace-Stieltjes transform of the distribution function at t
        double minVal = (double) this.getParam(1).getValue();
        double maxVal = (double) this.getParam(2).getValue();
        return (Math.exp(-s * minVal) - Math.exp(-s * maxVal)) / (s * (maxVal - minVal));
    }
}

