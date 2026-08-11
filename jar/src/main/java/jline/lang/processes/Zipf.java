/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.util.Pair;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.apache.commons.math3.util.FastMath;

import java.util.ArrayList;
import java.util.List;
import jline.util.RandomManager;
import java.util.Random;

/**
 * A Zipf-like probability distribution
 */
public class Zipf extends DiscreteDistribution {
    /**
     * Construct a Zipf-like distribution
     *
     * @param s - the shape
     * @param n - number of items
     */
    public Zipf(double s, int n) {
        super("Zipf", 4, new Pair<Double, Double>(1.0, (double) n));
        Matrix p = new Matrix(1, n);
        double h = Zipf.genHarmonic(s, n);
        Matrix x = new Matrix(1, n);
        for (int i = 0; i < p.getNumCols(); i++) {
            p.set(i, 1 / Math.pow(i + 1, s) / h);
            x.set(i, i + 1);
        }
        setParam(1, "p", p);
        setParam(2, "x", x);
        setParam(3, "s", s);
        setParam(4, "n", n);
    }

    /**
     * Generate harmonic numbers to normalize a Zipf-like distribution on n items with shape s
     */
    public static double genHarmonic(double s, int n) {
        double Hnm = 0;
        for (int k = 1; k <= n; k++) {
            Hnm += 1.0 / Math.pow(k, s);
        }
        return Hnm;
    }

    public static void main(String[] args) {
        Zipf z = new Zipf(2.1, 10);
        System.out.println(z.getSkewness());
    }

    /**
     * Evaluates the cumulative distribution function at t
     *
     * @param t - the point where the cdf will be evaluated
     * @return - the cdf evaluated at t
     */
    @Override
    public double evalCDF(double t) {
        double s = (double) this.getParam(3).getValue();
        int n = (int) this.getParam(4).getValue();
        return genHarmonic(s, (int) t) / genHarmonic(s, n);
    }

    public Matrix evalPMF() {
        int n = (int) this.getParam(4).getValue();
        List<Double> t = new ArrayList<>();
        for (int i = 1; i <= n; i++) {
            t.add((double) i);
        }
        return evalPMF(t);
    }

    /**
     * Evaluates the probability mass function at t
     *
     * @param t - the point where the pmf will be evaluated
     * @return - the pfm evaluated at t
     */
    @Override
    public Matrix evalPMF(List<Double> t) {
        double s = (double) this.getParam(3).getValue();
        int n = (int) this.getParam(4).getValue();
        List<Double> retList = new ArrayList<>();
        double Hns = genHarmonic(s, n);
        for (double d : t) {
            retList.add(1 / Math.pow(d, s) / Hns);
        }
        return new Matrix(retList);
    }

    /**
     * Computes the distribution mean
     *
     * @return - the mean of the distribution
     */
    @Override
    public double getMean() {
        double s = (double) this.getParam(3).getValue();
        int n = (int) this.getParam(4).getValue();
        return genHarmonic(s - 1, n) / genHarmonic(s, n);
    }

    /**
     * Computes the squared coefficient of variation == variance/mean^2
     *
     * @return - the squared coefficient of variation
     */
    @Override
    public double getSCV() {
        double s = (double) this.getParam(3).getValue();
        int n = (int) this.getParam(4).getValue();
        double ex = getMean();
        double var = genHarmonic(s - 2, n) / genHarmonic(s, n) - ex * ex;
        return var / (ex * ex);
    }

    @Override
    public double getSkewness() {
        Matrix p = (Matrix) this.getParam(1).getValue();
        Matrix x = (Matrix) this.getParam(2).getValue();
        double e3 = 0.0;
        double e1 = getMean();
        double e2 = 0.0;
        for (int k = 0; k < x.length(); k++) {
            e2 += FastMath.pow(x.get(k), 2) * p.get(k);
        }

        for (int k = 0; k < x.length(); k++) {
            e3 += FastMath.pow(x.get(k), 3) * p.get(k);
        }

        return (e3 - 3 * e1 * e2 + 2 * e1 * e1 * e1) / FastMath.pow(getVar(), 1.5);
    }

    /**
     * Evaluates the Laplace-Stieltjes Transform at s.
     * For Zipf distribution, LST(s) = Σ(k=1 to n) P(X=k) * e^(-s*k)
     *
     * @param s the Laplace domain variable
     * @return the LST value at s
     */
    public double evalLST(double s) {
        double shape = (double) this.getParam(3).getValue();
        int n = (int) this.getParam(4).getValue();
        double Hns = genHarmonic(shape, n);
        
        double lst = 0.0;
        for (int k = 1; k <= n; k++) {
            double prob = 1.0 / Math.pow(k, shape) / Hns;
            lst += prob * FastMath.exp(-s * k);
        }
        
        return lst;
    }

    /**
     * Generates samples from the Zipf distribution using the specified random generator.
     * Uses the inverse transform sampling method.
     *
     * @param n      the number of samples to generate
     * @param random the random number generator to use
     * @return array of samples from the Zipf distribution
     */
    @Override
    public double[] sample(int n, Random random) {
        double[] samples = new double[n];
        double s = (double) this.getParam(3).getValue();
        int m = (int) this.getParam(4).getValue();

        // Precompute normalization constant (Hurwitz zeta function approximation)
        double norm = 0.0;
        for (int k = 1; k <= m; k++) {
            norm += 1.0 / Math.pow(k, s);
        }

        for (int i = 0; i < n; i++) {
            double u = random.nextDouble();
            double cdf = 0.0;

            // Find the value k such that CDF(k-1) < u <= CDF(k)
            for (int k = 1; k <= m; k++) {
                cdf += (1.0 / Math.pow(k, s)) / norm;
                if (u <= cdf) {
                    samples[i] = k;
                    break;
                }
            }
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
        MatrixCell representation = new MatrixCell();
        representation.set(0, Matrix.singleton(getMean()));
        representation.set(1, Matrix.singleton(getSCV()));
        return representation;
    }

}
