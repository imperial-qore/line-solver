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

        double lambda = getMean();
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

    /**
     * Evaluates the probability mass function (PMF) at value k.
     * For Poisson(lambda), PMF(k) = lambda^k * exp(-lambda) / k!
     *
     * @param k the value to evaluate
     * @return the PMF value at k
     */
    public double evalPMF(double k) {
        if (k < 0 || k != Math.floor(k)) {
            return 0.0;
        }
        double lambda = this.getMean();
        int n = (int) k;
        // Computed in log-space for numerical stability at large k
        return FastMath.exp(n * FastMath.log(lambda) - lambda - Maths.factln(n));
    }

    public double evalPDF(int n) {
        return evalPMF((double) n);
    }

    /**
     * Evaluates the Laplace-Stieltjes Transform at s.
     * For Poisson(λ), LST(s) = exp(λ(e^(-s) - 1))
     *
     * @param s the Laplace domain variable
     * @return the LST value at s
     */
    public double evalLST(double s) {
        double lambda = this.getMean();
        return FastMath.exp(lambda * (FastMath.exp(-s) - 1));
    }

    public double getMean() {
        return (double) this.getParam(1).getValue();
    }

    // see _kb/01-model-classes.md (Java process-construction notes) for rationale

    public double getSCV() {
        return 1.0 / getMean();
    }

    public double getSkewness() {
        // Poisson skewness is lambda^(-1/2): E[(X-mu)^3] = lambda and
        // sigma^3 = lambda^(3/2). The former sqrt(lambda) inverted it.
        return 1.0 / FastMath.sqrt(getMean());
    }

    public double getVar() {
        return getMean();
    }

    @Override
    public double[] sample(int nsamples, Random random) {
        double lambda = this.getMean();
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
