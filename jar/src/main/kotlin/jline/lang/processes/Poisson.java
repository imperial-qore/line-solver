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
import java.util.Random;

import static org.apache.commons.math3.util.CombinatoricsUtils.factorial;

/**
 * A Poisson discrete distribution
 */
public class Poisson extends DiscreteDistribution implements Serializable {
    public Poisson(double rate) {
        super("Poisson", 1, new Pair<Double, Double>(0.0, Inf));
        this.setParam(1, "lambda", rate);
    }

    public static void main(String[] args) {
        Poisson p = new Poisson(61.0);
        double[] samples = p.sample(1000);
        double mean = 0.0;
        for (int i = 0; i < samples.length; i++) {
            mean += samples[i];
        }
        mean /= samples.length;
        System.out.println("mean = " + mean);
    }

    /**
     * Evaluates the cumulative distribution function (CDF) at value t.
     * For a discrete Poisson distribution, this is P(X <= floor(t)).
     *
     * @param t the value to evaluate
     * @return the CDF value at t
     */
    public double evalCDF(double t) {
        if (t < 0) {
            return 0.0;
        }

        double lambda = getRate();
        int k = (int) Math.floor(t);

        // For Poisson distribution: CDF(k) = sum from i=0 to k of (lambda^i * exp(-lambda) / i!)
        double cdf = 0.0;
        double term = Math.exp(-lambda); // First term: lambda^0 * exp(-lambda) / 0!

        for (int i = 0; i <= k; i++) {
            cdf += term;
            if (i < k) {
                term *= lambda / (i + 1); // Next term: multiply by lambda/(i+1)
            }
        }

        return Math.min(1.0, cdf); // Ensure CDF doesn't exceed 1 due to numerical errors
    }

    public double evalPDF(int n) {
        double lambda = this.getRate();
        double num = FastMath.pow(lambda, n) * Math.exp(-lambda);
        return num / ((double) factorial(n));
    }

    /**
     * Evaluates the Laplace-Stieltjes Transform at s.
     * For Poisson(λ), LST(s) = exp(λ(e^(-s) - 1))
     *
     * @param s the Laplace domain variable
     * @return the LST value at s
     */
    public double evalLST(double s) {
        double lambda = this.getRate();
        return FastMath.exp(lambda * (FastMath.exp(-s) - 1));
    }

    public double getMean() {
        return this.getRate();
    }

    public double getRate() {
        return (double) this.getParam(1).getValue();
    }

    public double getSCV() {
        return 1.0 / this.getRate();
    }

    public double getSkewness() {
        return FastMath.sqrt(this.getRate());
    }

    public double getVar() {
        return this.getRate();
    }

    @Override
    public double[] sample(int nsamples, Random random) {
        double lambda = this.getRate();
        double[] samples = new double[nsamples];
        for (int i = 0; i < nsamples; i++) {
            double l = FastMath.exp(-lambda);
            double p = 1.0;
            int k = 0;
            do {
                k++;
                p *= random.nextDouble();
            } while (p > l);
            samples[i] = k - 1;
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
