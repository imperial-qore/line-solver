/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import static jline.GlobalConstants.Inf;

import jline.util.Pair;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.apache.commons.math3.util.FastMath;

import java.io.Serializable;
import jline.util.RandomManager;
import java.util.Random;

/**
 * A Geometric distribution.
 */

public class Geometric extends DiscreteDistribution implements Serializable {

    // Geometric with support in 1,2,..., so the first success is counted.
    public Geometric(double probability) {
        super("Geometric", 1, new Pair<Double, Double>(0.0, Inf));
        this.setParam(1, "p", probability);
    }

    public static void main(String[] args) {
        Geometric g = new Geometric(0.5);
        double[] samples = g.sample(100000, RandomManager.getThreadRandomAsRandom());
        System.out.println(g.evalPMF(3));
    }

    public double evalCDF(double t) {
        double p = (double) this.getParam(1).getValue();
        return 1.0 - FastMath.pow(1.0 - p, t);
    }

    public double evalPMF(double t) {
        double p = (double) this.getParam(1).getValue();
        return FastMath.pow(1 - p, t - 1) * p;
    }

    /**
     * Evaluates the Laplace-Stieltjes Transform at s.
     * For Geometric(p), LST(s) = p*e^(-s) / (1 - (1-p)*e^(-s))
     *
     * @param s the Laplace domain variable
     * @return the LST value at s
     */
    public double evalLST(double s) {
        double p = (double) this.getParam(1).getValue();
        double e_neg_s = FastMath.exp(-s);
        return (p * e_neg_s) / (1 - (1 - p) * e_neg_s);
    }

    public double getMean() {
        double p = (double) this.getParam(1).getValue();
        return 1 / p;
    }

    public double getRate() {
        return (double) this.getParam(1).getValue();
    }

    public double getSCV() {
        return getVar() / Math.pow(getMean(), 2);
    }

    public double getSkewness() {
        double p = (double) this.getParam(1).getValue();
        return (2 - p) / Math.sqrt(1 - p);
    }

    public double getVar() {
        double p = (double) this.getParam(1).getValue();
        return (1 - p) / p / p;
    }

    /**
     * Gets n samples from the distribution
     *
     * @param n - the number of samples
     * @return - n samples from the distribution
     */
    @Override
    public double[] sample(int n, Random random) {
        double[] samples = new double[n];
        double p = (double) this.getParam(1).getValue();
        for (int i = 0; i < n; i++) {
            samples[i] = (int) FastMath.ceil(Math.log(1 - random.nextDouble()) / FastMath.log(1 - p));
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
