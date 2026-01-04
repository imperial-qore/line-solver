/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.GlobalConstants;
import jline.util.Pair;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import static jline.io.InputOutputKt.line_error;
import static jline.io.InputOutputKt.mfilename;

/**
 * Prior distribution representing parameter uncertainty over alternative distributions.
 *
 * Prior represents parameter uncertainty by specifying a discrete set of
 * alternative distributions with associated probabilities. When used with
 * setService or setArrival, it causes the Posterior solver to expand the
 * model into a family of networks, one for each alternative.
 *
 * This is NOT a mixture distribution - each alternative represents a
 * separate model realization with its associated prior probability.
 *
 * Example usage:
 * <pre>
 * List&lt;Distribution&gt; dists = Arrays.asList(new Exp(1.0), new Exp(2.0), new Erlang(1.5, 2));
 * double[] probs = {0.4, 0.35, 0.25};
 * Prior prior = new Prior(dists, probs);
 * queue.setService(jobclass, prior);
 * </pre>
 */
public class Prior extends Distribution implements Serializable {

    /** Cell array of alternative distributions */
    protected List<Distribution> distributions;

    /** Vector of prior probabilities (must sum to 1) */
    protected double[] probabilities;

    /**
     * Creates a Prior distribution with specified alternatives and probabilities.
     *
     * @param distributions list of alternative Distribution objects
     * @param probabilities array of probabilities (must sum to 1 and be non-negative)
     * @throws RuntimeException if validation fails
     */
    public Prior(List<Distribution> distributions, double[] probabilities) {
        super("Prior", 2, new Pair<Double, Double>(0.0, Double.POSITIVE_INFINITY));

        // Validate distributions
        if (distributions == null || distributions.isEmpty()) {
            line_error(mfilename(new Object(){}), "distributions cannot be null or empty");
        }

        // Validate probabilities
        if (probabilities == null || probabilities.length != distributions.size()) {
            line_error(mfilename(new Object(){}),
                "Number of distributions must match number of probabilities");
        }

        double sum = 0.0;
        for (double p : probabilities) {
            if (p < 0) {
                line_error(mfilename(new Object(){}), "Probabilities must be non-negative");
            }
            sum += p;
        }

        if (Math.abs(sum - 1.0) > GlobalConstants.CoarseTol) {
            line_error(mfilename(new Object(){}),
                String.format("Probabilities must sum to 1 (current sum: %f)", sum));
        }

        this.distributions = new ArrayList<Distribution>(distributions);
        this.probabilities = probabilities.clone();

        this.setParam(1, "distributions", this.distributions);
        this.setParam(2, "probabilities", this.probabilities);
    }

    /**
     * Convenience constructor for array-based input.
     *
     * @param distributions array of alternative Distribution objects
     * @param probabilities array of probabilities
     */
    public Prior(Distribution[] distributions, double[] probabilities) {
        this(arrayToList(distributions), probabilities);
    }

    private static List<Distribution> arrayToList(Distribution[] arr) {
        List<Distribution> list = new ArrayList<Distribution>();
        for (Distribution d : arr) {
            list.add(d);
        }
        return list;
    }

    /**
     * Returns the number of alternative distributions.
     *
     * @return number of alternatives
     */
    public int getNumAlternatives() {
        return distributions.size();
    }

    /**
     * Returns the distribution at the specified index.
     *
     * @param idx 0-based index
     * @return the distribution at index idx
     * @throws IndexOutOfBoundsException if idx is out of bounds
     */
    public Distribution getAlternative(int idx) {
        if (idx < 0 || idx >= distributions.size()) {
            line_error(mfilename(new Object(){}), "Index out of bounds");
        }
        return distributions.get(idx);
    }

    /**
     * Returns the probability of the alternative at the specified index.
     *
     * @param idx 0-based index
     * @return the probability at index idx
     * @throws IndexOutOfBoundsException if idx is out of bounds
     */
    public double getProbability(int idx) {
        if (idx < 0 || idx >= probabilities.length) {
            line_error(mfilename(new Object(){}), "Index out of bounds");
        }
        return probabilities[idx];
    }

    /**
     * Returns the array of all probabilities.
     *
     * @return copy of probability array
     */
    public double[] getProbabilities() {
        return probabilities.clone();
    }

    /**
     * Returns the list of all distributions.
     *
     * @return list of distributions (unmodifiable view recommended)
     */
    public List<Distribution> getDistributions() {
        return new ArrayList<Distribution>(distributions);
    }

    /**
     * Gets the prior-weighted mean (expected mean over alternatives).
     * E[X] = sum_i p_i * E[X_i]
     *
     * @return the weighted mean
     */
    @Override
    public double getMean() {
        double mean = 0.0;
        for (int i = 0; i < distributions.size(); i++) {
            mean += probabilities[i] * distributions.get(i).getMean();
        }
        return mean;
    }

    /**
     * Gets the prior-weighted SCV using law of total variance.
     * Var(X) = E[Var(X|D)] + Var(E[X|D])
     * SCV = Var(X) / E[X]^2
     *
     * @return the weighted SCV
     */
    @Override
    public double getSCV() {
        double E_mean = 0.0;      // E[E[X|D]]
        double E_var = 0.0;       // E[Var(X|D)]
        double E_mean_sq = 0.0;   // E[E[X|D]^2]

        for (int i = 0; i < distributions.size(); i++) {
            double m = distributions.get(i).getMean();
            double v = distributions.get(i).getSCV() * m * m;  // Var(X|D=i)
            E_mean += probabilities[i] * m;
            E_var += probabilities[i] * v;
            E_mean_sq += probabilities[i] * m * m;
        }

        // Total variance = E[Var(X|D)] + Var(E[X|D])
        // Var(E[X|D]) = E[E[X|D]^2] - E[E[X|D]]^2
        double total_var = E_var + (E_mean_sq - E_mean * E_mean);
        return total_var / (E_mean * E_mean);
    }

    /**
     * Gets the prior-weighted skewness.
     * Uses mixture formula with moment calculations.
     *
     * @return the weighted skewness
     */
    @Override
    public double getSkewness() {
        double mu = getMean();
        double sigma2 = getVar();
        double sigma = Math.sqrt(sigma2);

        if (sigma < GlobalConstants.FineTol) {
            return 0.0;
        }

        // E[(X - mu)^3] via mixture
        double third_central = 0.0;
        for (int i = 0; i < distributions.size(); i++) {
            double mi = distributions.get(i).getMean();
            double vi = distributions.get(i).getVar();
            double si = Math.sqrt(vi);
            double skewi = distributions.get(i).getSkewness();

            double delta = mi - mu;
            // Third central moment of Xi around its own mean
            double m3i = skewi * si * si * si;
            // Third central moment of Xi around global mu
            double m3_shifted = m3i + 3 * vi * delta + delta * delta * delta;

            third_central += probabilities[i] * m3_shifted;
        }

        return third_central / (sigma * sigma * sigma);
    }

    /**
     * Generates random samples from the Prior (mixture sampling).
     * Each sample is drawn from one of the alternatives selected
     * according to the prior probabilities.
     *
     * @param n number of samples
     * @param random the random number generator
     * @return array of samples
     */
    @Override
    public double[] sample(int n, Random random) {
        double[] samples = new double[n];

        // Build cumulative probability array
        double[] cumprob = new double[probabilities.length];
        cumprob[0] = probabilities[0];
        for (int i = 1; i < probabilities.length; i++) {
            cumprob[i] = cumprob[i-1] + probabilities[i];
        }

        for (int i = 0; i < n; i++) {
            // Select alternative based on probabilities
            double r = random.nextDouble();
            int idx = 0;
            for (int j = 0; j < cumprob.length; j++) {
                if (r <= cumprob[j]) {
                    idx = j;
                    break;
                }
            }
            samples[i] = distributions.get(idx).sample(1, random)[0];
        }
        return samples;
    }

    /**
     * Evaluates the mixture CDF at t.
     * F(t) = sum_i p_i * F_i(t)
     *
     * @param t the point at which to evaluate
     * @return the CDF value
     */
    @Override
    public double evalCDF(double t) {
        double cdf = 0.0;
        for (int i = 0; i < distributions.size(); i++) {
            cdf += probabilities[i] * distributions.get(i).evalCDF(t);
        }
        return cdf;
    }

    /**
     * Evaluates the mixture Laplace-Stieltjes transform.
     * L(s) = sum_i p_i * L_i(s)
     *
     * @param s the transform variable
     * @return the LST value
     */
    @Override
    public double evalLST(double s) {
        double lst = 0.0;
        for (int i = 0; i < distributions.size(); i++) {
            lst += probabilities[i] * distributions.get(i).evalLST(s);
        }
        return lst;
    }

    /**
     * Returns true indicating this is a Prior distribution.
     * Used for detection by Posterior solver.
     *
     * @return true
     */
    public boolean isPrior() {
        return true;
    }

    /**
     * Static method to check if a distribution is a Prior.
     *
     * @param dist the distribution to check
     * @return true if dist is a Prior
     */
    public static boolean isPriorDistribution(Distribution dist) {
        return dist instanceof Prior;
    }
}
