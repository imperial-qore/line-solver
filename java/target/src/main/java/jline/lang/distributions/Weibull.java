package jline.lang.distributions;

import jline.util.Maths;
import jline.util.Pair;
import org.apache.commons.math3.distribution.WeibullDistribution;
import org.apache.commons.math3.special.Gamma;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class Weibull extends ContinuousDistribution implements Serializable {

    public Weibull(double shape, double scale) {
        super("Weibull", 2, new Pair<Double,Double>(0.0, Double.POSITIVE_INFINITY));
        if (shape < 0) {
            System.err.println("shape parameter must be >= 0.0");
        }
        this.setParam(1, "alpha", scale);
        this.setParam(2, "r", shape);
    }

    /**
     * Gets n samples from the distribution
     * @param n - the number of samples
     * @return - n samples from the distribution
     */
    @Override
    public List<Double> sample(long n) {
        return this.sample(n, new Random());
    }

    @Override
    public List<Double> sample(long n, Random random) {
        double alpha = (double) this.getParam(1).getValue();
        double r = (double) this.getParam(2).getValue();

        List<Double> samples = new ArrayList<>();
        WeibullDistribution weibullDistribution = new WeibullDistribution(alpha, r);
        for (int i = 0; i < n; i++) {
            samples.add(weibullDistribution.inverseCumulativeProbability(random.nextDouble()));
        }
        return samples;
    }

    @Override
    public double getMean() {
        double alpha = (double) this.getParam(1).getValue();
        double r = (double) this.getParam(2).getValue();

        // Mean of Weibull: scale * Gamma(1 + 1/shape)
        return alpha * Maths.gammaFunction(1 + 1 / r);
    }

    @Override
    public double getRate() {
        return 1 / this.getMean();
    }

    @Override
    public double getSCV() {
        // Get distribution squared coefficient of variation (SCV = variance / mean^2)
        return getVar() / Math.pow(getMean(), 2);
    }

    @Override
    public double getVar() {
        double alpha = (double) this.getParam(1).getValue();
        double r = (double) this.getParam(2).getValue();
        return Math.pow(alpha, 2) * (Maths.gammaFunction(1 + 2 / r) - Math.pow(Maths.gammaFunction(1 + 1 / r), 2));
    }

    @Override
    public double getSkew() {
        double r = (double) this.getParam(2).getValue();
        double gamma1 = Maths.gammaFunction(1 + 1 / r);
        double gamma2 = Maths.gammaFunction(1 + 2 / r);
        double gamma3 = Maths.gammaFunction(1 + 3 / r);

        double numerator = gamma3 - 3 * gamma2 * gamma1 + 2 * Math.pow(gamma1, 3);
        double denominator = Math.pow(gamma2 - Math.pow(gamma1, 2), 1.5);

        return numerator / denominator;
    }

    @Override
    public double evalCDF(double t) {
        // Evaluate the cumulative distribution function at t
        double alpha = (double) this.getParam(1).getValue();
        double r = (double) this.getParam(2).getValue();

        if (t <= 0.0) {
            return 0.0;
        } else {
            return 1 - Math.exp(-Math.pow(t / alpha, r));
        }
    }

    @Override
    public double evalLST(double s) {
        double alpha = (double) this.getParam(1).getValue();
        double r = (double) this.getParam(2).getValue();
        double term = Math.pow(r, alpha) * Math.pow(s, alpha);
        double incompleteGamma = Gamma.regularizedGammaQ(1 + 1 / alpha, term);
        return Math.exp(-term * incompleteGamma);
    }

    public static Weibull fitMeanAndSCV(double mean, double scv) {
        // Fit distribution with given mean and squared coefficient of variation (SCV=variance/mean^2)
        double c = Math.sqrt(scv);
        double r = Math.pow(c, -1.086); //Justus approximation (1976)
        double alpha = mean / Maths.gammaFunction(1 + 1 / r);
        return new Weibull(r, alpha);
    }
}

