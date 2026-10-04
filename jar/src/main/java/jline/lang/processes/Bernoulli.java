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
 * A Bernoulli distribution for modeling binary random variables.
 * 
 * <p>The Bernoulli distribution models a single trial with two possible outcomes:
 * success (1) with probability p, or failure (0) with probability (1-p). It's
 * the fundamental building block for binomial processes and binary decision modeling.</p>
 * 
 * <p>Key characteristics:
 * <ul>
 *   <li>Support: {0, 1}</li>
 *   <li>Parameter: p ∈ [0,1] (success probability)</li>
 *   <li>Mean: p</li>
 *   <li>Variance: p(1-p)</li>
 *   <li>Special case of binomial distribution with n=1</li>
 * </ul>
 * </p>
 * 
 * <p>Common applications include modeling job completion states, server availability,
 * and routing decisions in queueing networks.</p>
 * 
 * @see Binomial
 * @see DiscreteDistribution
 * @since 1.0
 */
public class Bernoulli extends DiscreteDistribution implements Serializable {

    public Bernoulli(double prob) {
        // a Binomial with n=1
        super("Bernoulli", 2, new Pair<Double, Double>(0.0, Inf));
        this.setParam(1, "p", prob); // success probability
    }

    public double evalCDF(int t) {
        return evalCDF((double) t);
    }

    public double evalCDF(double t) {
        int n = 1;
        double p = (double) this.getParam(1).getValue();
        double ret = 0.0;
        if (t >= n) {
            ret = 1.0;
        } else if (t >= 0) {
            ret = 1.0 - regularizedBeta(p, t + 1.0, (n - t));
        }
        return ret;
    }

    public double evalPMF(int k) {
        int n = 1;
        double p = (double) this.getParam(1).getValue();
        return FastMath.exp(Maths.factln(n) - Maths.factln(k) - Maths.factln(n - k)) * FastMath.pow(p, k) * FastMath.pow(1 - p, n - k);
    }

    /**
     * Evaluates the Laplace-Stieltjes Transform at s.
     * For Bernoulli(p), LST(s) = (1 - p + p*e^(-s))
     *
     * @param s the Laplace domain variable
     * @return the LST value at s
     */
    public double evalLST(double s) {
        double p = (double) this.getParam(1).getValue();
        return 1 - p + p * FastMath.exp(-s);
    }

    public double getMean() {
        // Mean of Bernoulli(p) is simply p (n*p where n=1)
        return (double) this.getParam(1).getValue();
    }

    public double getSCV() {
        return getVar() / getMean() / getMean();
    }

    public double getSkewness() {
        int n = 1;
        double p = (double) this.getParam(1).getValue();
        double q = 1 - p;

        return (q - p) / (Math.sqrt(n * p * q));
    }

    public double getVar() {
        int n = 1;
        double p = (double) this.getParam(1).getValue();
        return n * p * (1.0 - p);
    }

    @Override
    public double[] sample(int nsamples, Random random) {
        double[] ret = new double[(int) nsamples];
        int n = 1;
        double p = (double) this.getParam(1).getValue();
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
