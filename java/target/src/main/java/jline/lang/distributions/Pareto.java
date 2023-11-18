package jline.lang.distributions;

import jline.util.Pair;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.RombergIntegrator;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class Pareto extends ContinuousDistribution implements Serializable {

    public Pareto(double shape, double scale) {
        super("Pareto", 2, new Pair<Double,Double>(0.0, Double.POSITIVE_INFINITY));
        if (shape < 2){
            System.err.println("shape parameter must be >= 2.0");
        }
        this.setParam(1, "alpha", shape);
        this.setParam(2, "k", scale);
    }

    /**
     * Gets n samples from the distribution
     * @param n - the number of samples
     * @return - n samples from the distribution
     */
    @Override
    public List<Double> sample(long n) {
        return this.sample(n,new Random());
    }

    @Override
    public List<Double> sample(long n, Random random) {
        List<Double> samples = new ArrayList<>();
        double shape = (double)this.getParam(1).getValue();
        double scale = (double)this.getParam(2).getValue();

        for (long i = 0; i < n; i++) {
            double u = random.nextDouble();
            double sample = scale / Math.pow(u, 1 / shape);
            samples.add(sample);
        }

        return samples;
    }

    @Override
    public double getMean() {
        double shape = (double)this.getParam(1).getValue();
        double scale = (double)this.getParam(2).getValue();

        if (shape <= 1) {
            return Double.POSITIVE_INFINITY;
        } else {
            return shape * scale / (shape - 1);
        }
    }

    @Override
    public double getRate() {
        return 1 / getMean();
    }

    @Override
    public double getSCV() {
        double shape = (double)this.getParam(1).getValue();
        if (shape <= 2) {
            return Double.POSITIVE_INFINITY;
        } else {
            double scale = (double)this.getParam(2).getValue();
            double var = Math.pow(scale, 2) * shape / Math.pow(shape - 1, 2) / (shape - 2);
            double ex = shape * scale / (shape - 1);
            return var / Math.pow(ex, 2);
        }
    }

    @Override
    public double getVar() {
        double shape = (double)this.getParam(1).getValue();
        double scale = (double)this.getParam(2).getValue();

        if (shape <= 2) {
            return Double.POSITIVE_INFINITY;
        } else {
            return (scale * scale * shape) / ((shape - 1) * (shape - 1) * (shape - 2));
        }
    }

    @Override
    public double getSkew() {
        double shape = (double)this.getParam(1).getValue();
        if (shape <= 3) {
            return Double.POSITIVE_INFINITY;
        } else {
            return 2*(1 + shape) / (shape - 3) * Math.sqrt((shape - 2) / shape);
        }
    }

    @Override
    public double evalCDF(double t) {
        double shape = (double)this.getParam(1).getValue();
        double scale = (double)this.getParam(2).getValue();
        double k = 1/shape;
        double sigma = scale * k;
        return gpcdf(t, k, sigma, sigma/k);
    }

    public static double gpcdf(double x, double k, double sigma, double theta) {
        // Check for invalid sigma
        if (sigma <= 0) {
            return Double.NaN;
        }

        // Calculate the (x - theta) / sigma term
        double z = (x - theta) / sigma;

        // Return 0 for out-of-range values
        if (z < 0) {
            return 0.0;
        }

        // Handle the k == 0 case
        if (Math.abs(k) < Math.ulp(1.0)) {
            return -Math.expm1(-z);
        }

        // Compute the CDF value for k != 0
        double t = z * k;

        // When k < 0, the support is 0 <= x/sigma <= -1/k.
        if (t <= -1 && k < -Math.ulp(1.0)) {
            return 1.0;
        }

        return -Math.expm1((-1.0 / k) * Math.log1p(t));
    }

    @Override
    public double evalLST(double s) {
        double alpha = (double)this.getParam(1).getValue();
        double k = (double)this.getParam(2).getValue();
        RombergIntegrator integrator = new RombergIntegrator();
        UnivariateFunction function = new UnivariateFunction() {
            public double value(double t) {
                return Math.pow(t, -alpha - 1) * Math.exp(-t);
            }
        };

        double IL = integrator.integrate(10000, function, s * k, Double.POSITIVE_INFINITY);
        return alpha * Math.pow(k, alpha) * IL * Math.pow(s, 1 + alpha) / s;
    }

    public static Pareto fitMeanAndSCV(double mean, double scv){
        // Fit distribution with given mean and squared coefficient of variation (SCV=variance/mean^2)
        if (mean <= 0 || scv <= 0) {
            throw new RuntimeException("Mean and SCV should be positive.");
        }
        double shape = (scv * mean + mean * Math.sqrt(scv * (scv + 1))) / (scv * mean);
        double scale = mean + scv * mean - mean * Math.sqrt(scv * (scv + 1));
        return new Pareto(shape, scale);
    }
}

