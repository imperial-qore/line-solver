package jline.lang.distributions;

import jline.util.Matrix;
import jline.util.Pair;
import org.apache.commons.math3.distribution.GammaDistribution;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class Gamma extends ContinuousDistribution implements Serializable {

    public Gamma(double shape, double scale) {
        super("Gamma", 2, new Pair<Double,Double>(0.0, Double.POSITIVE_INFINITY));
        this.setParam(1, "alpha", shape);
        this.setParam(2, "beta", scale);
    }

    @Override
    public Matrix sample(long n) {
        double shape = (double) this.getParam(1).getValue();
        double scale = (double) this.getParam(2).getValue();
        GammaDistribution gammaDistribution = new GammaDistribution(shape, scale);
        double[] samples = gammaDistribution.sample((int) n);
        List<Double> sampleList = new ArrayList<>();
        for (double sample : samples) {
            sampleList.add(sample);
        }
        return new Matrix(sampleList);
    }

    public Matrix sample(long n, Random random) {
        double shape = (double) this.getParam(1).getValue();
        double scale = (double) this.getParam(2).getValue();
        GammaDistribution gammaDistribution = new GammaDistribution(shape, scale);
        List<Double> sampleList = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            sampleList.add(gammaDistribution.inverseCumulativeProbability(random.nextDouble()));
        }
        return new Matrix(sampleList);
    }

    @Override
    public double getMean() {
        double shape = (double) this.getParam(1).getValue();
        double scale = (double) this.getParam(2).getValue();
        return shape * scale;
    }

    @Override
    public double getRate() {
        return 1.0 / getMean();
    }

    @Override
    public double getSCV() {
        double shape = (double) this.getParam(1).getValue();
        return 1.0 / shape;
    }

    @Override
    public double getVar() {
        double shape = (double) this.getParam(1).getValue();
        double scale = (double) this.getParam(2).getValue();
        return shape * scale * scale;
    }

    @Override
    public double getSkew() {
        double shape = (double) this.getParam(1).getValue();
        return 2.0 / Math.sqrt(shape);
    }

    @Override
    public double evalCDF(double t) {
        double shape = (double) this.getParam(1).getValue();
        double scale = (double) this.getParam(2).getValue();
        GammaDistribution gammaDistribution = new GammaDistribution(shape, scale);
        return gammaDistribution.cumulativeProbability(t);
    }

    @Override
    public double evalLST(double s) {
        double shape = (double) this.getParam(1).getValue();
        double scale = (double) this.getParam(2).getValue();
        return Math.pow((1 / scale) / (s + 1 / scale), shape);
    }

    public static Gamma fitMeanAndSCV(double mean, double scv) {
        double shape = 1 / scv;
        double scale = mean / shape;
        return new Gamma(shape, scale);
    }
}

