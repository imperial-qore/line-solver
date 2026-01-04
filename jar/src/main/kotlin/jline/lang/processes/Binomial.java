/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import static jline.GlobalConstants.Inf;

import jline.util.Maths;
import jline.util.Pair;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.apache.commons.math3.util.FastMath;

import java.io.Serializable;
import java.util.Random;

import static org.apache.commons.math3.special.Beta.regularizedBeta;

/**
 * A Binomial distribution
 */
public class Binomial extends DiscreteDistribution implements Serializable {

    public Binomial(int n, double prob) {
        super("Binomial", 2, new Pair<Double, Double>(0.0, Inf));
        this.setParam(1, "n", n); // number of trials
        this.setParam(2, "p", prob); // success probability
    }

    public double evalCDF(int t) {
        return evalCDF((double) t);
    }

    public double evalCDF(double t) {
        int n = (int) this.getParam(1).getValue();
        double p = (double) this.getParam(2).getValue();
        double ret = 0.0;
        if (t >= n) {
            ret = 1.0;
        } else if (t >= 0) {
            ret = 1.0 - regularizedBeta(p, t + 1.0, (n - t));
        }
        return ret;
    }

    public double evalPMF(int k) {
        int n = (int) this.getParam(1).getValue();
        double p = (double) this.getParam(2).getValue();
        return FastMath.exp(Maths.factln(n) - Maths.factln(k) - Maths.factln(n - k)) * FastMath.pow(p, k) * FastMath.pow(1 - p, n - k);
    }

    /**
     * Evaluates the Laplace-Stieltjes Transform at s.
     * For Binomial(n, p), LST(s) = (1 - p + p*e^(-s))^n
     *
     * @param s the Laplace domain variable
     * @return the LST value at s
     */
    public double evalLST(double s) {
        int n = (int) this.getParam(1).getValue();
        double p = (double) this.getParam(2).getValue();
        return FastMath.pow(1 - p + p * FastMath.exp(-s), n);
    }

    public double getMean() {
        return ((double) this.getParam(2).getValue()) * ((int) this.getParam(1).getValue());
    }

    public double getSCV() {
        return getVar() / getMean() / getMean();
    }

    public double getSkewness() {
        int n = (int) this.getParam(1).getValue();
        double p = (double) this.getParam(2).getValue();
        double q = 1 - p;

        return (q - p) / (Math.sqrt(n * p * q));
    }

    public double getVar() {
        int n = (int) this.getParam(1).getValue();
        double p = (double) this.getParam(2).getValue();
        return n * p * (1.0 - p);
    }

    @Override
    public double[] sample(int nsamples, Random random) {
        double[] ret = new double[(int) nsamples];
        int n = (int) this.getParam(1).getValue();
        double p = (double) this.getParam(2).getValue();
        for (int k = 0; k < nsamples; k++) {
            int acc = 0;
            for (int i = 0; i < n; i++) {
                if (random.nextDouble() <= p) {
                    acc++;
                }
            }
            ret[k] = acc;
        }
        return ret;
    }

    public MatrixCell getProcess() {
        MatrixCell representation = new MatrixCell();
        representation.set(0, Matrix.singleton(getMean()));
        representation.set(1, Matrix.singleton(getSCV()));
        return representation;
    }
}
