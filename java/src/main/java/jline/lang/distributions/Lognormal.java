package jline.lang.distributions;

import jline.util.Matrix;
import jline.util.Pair;
import org.apache.commons.math3.distribution.LogNormalDistribution;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class Lognormal extends ContinuousDistribution implements Serializable {
    public Lognormal(double mu, double sigma) {
        super("LogNormal", 2, new Pair<Double,Double>(0.0, Double.POSITIVE_INFINITY));
        if (sigma < 0){
            System.err.println("sigma parameter must be >= 0.0");
        }
        this.setParam(1, "mu", mu);
        this.setParam(2, "sigma", sigma);
    }

    /**
     * Gets n samples from the distribution
     *
     * @param n - the number of samples
     * @return - n samples from the distribution
     */
    @Override
    public Matrix sample(long n) {
        return this.sample(n,new Random());
    }

    @Override
    public Matrix sample(long n, Random random) {
        List<Double> samples = new ArrayList<>();
        double mu = (double) this.getParam(1).getValue();
        double sigma = (double) this.getParam(2).getValue();
        for (int i = 0; i < n; i++) {
            double z = random.nextGaussian(); // Standard Normal Distribution
            double value = Math.exp(mu + sigma * z); // Log-normal Distribution
            samples.add(value);
        }
        return new Matrix(samples);
    }

    @Override
    public double getMean() {
        double mu = (double) this.getParam(1).getValue();
        double sigma = (double) this.getParam(2).getValue();
        return Math.exp(mu + Math.pow(sigma, 2) / 2.0);
    }

    @Override
    public double getRate() {
        return 1 / getMean();
    }

    @Override
    public double getSCV() {
        double mu = (double) this.getParam(1).getValue();
        double sigma = (double) this.getParam(2).getValue();
        double ex = Math.exp(mu + Math.pow(sigma, 2) / 2);
        double var = (Math.exp(Math.pow(sigma, 2)) - 1) * Math.exp(2 * mu + Math.pow(sigma, 2));
        return var / Math.pow(ex, 2);
    }

    @Override
    public double getVar() {
        double mu = (double) this.getParam(1).getValue();
        double sigma = (double) this.getParam(2).getValue();
        return (Math.exp(Math.pow(sigma, 2)) - 1) * Math.exp(2 * mu + Math.pow(sigma, 2));
    }

    @Override
    public double getSkew() {
        double sigma = (double) this.getParam(2).getValue();
        return (Math.exp(Math.pow(sigma, 2)) + 2) * Math.sqrt(Math.exp(Math.pow(sigma, 2)) - 1);
    }

    @Override
    public double evalCDF(double t) {
        double mu = (double) this.getParam(1).getValue();
        double sigma = (double) this.getParam(2).getValue();
        LogNormalDistribution logNormalDistribution = new LogNormalDistribution(mu, sigma);
        return logNormalDistribution.cumulativeProbability(t);
//        return 0.5 + 0.5 * UTIL.erf((Math.log(t) - mu) / (Math.sqrt(2) * sigma));
    }

    @Override
    public double evalLST(double s) {
        throw new RuntimeException("Laplace-Stieltjes transform of the Lognormal distribution not available yet."); // TODO: not implemented
    }

    public static Lognormal fitMeanAndSCV(double mean, double scv) {
        double c = Math.sqrt(scv);
        double mu = Math.log(mean / Math.sqrt(Math.pow(c, 2) + 1));
        double sigma = Math.sqrt(Math.log(Math.pow(c, 2) + 1));
        return new Lognormal(mu, sigma);
    }
}

