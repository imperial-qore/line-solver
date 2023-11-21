package jline.lang.distributions;

import java.io.Serializable;
import java.util.Random;

import jline.util.Matrix;
import jline.util.Pair;

public class Binomial extends DiscreteDistribution implements Serializable {
    public static final int MAX_N = Integer.MAX_VALUE;
    public Binomial(double prob, int n) {
        super("Binomial", 2, new Pair<Double,Double>(0.0, (double)Math.min(n, MAX_N)));
        n = Math.min(n, MAX_N);
        this.setParam(1, "p", prob);
        this.setParam(2, "n", n);
    }

    public int getRealization(Random random) {
        int acc = 0;
        double p = (double)this.getParam(1).getValue();
        int n = (int)this.getParam(2).getValue();
        for (int i = 0; i < n; i++) {
            if (random.nextDouble() >= p) {
                acc++;
            }
        }
        return acc;
    }

    @Override
    public Matrix sample(long n) {
        return this.sample(n,null);
    }
    @Override
    public Matrix sample(long n, Random random) {
        throw new RuntimeException("Not implemented");
    }

    public double getMean() {
        return ((double)this.getParam(1).getValue())*((double)this.getParam(2).getValue());
    }
    public double getRate() {
        return 0;
    }
    public double getSCV() {
        throw new RuntimeException("Not implemented");
    }
    public double getVar() {
        double p = (double)this.getParam(1).getValue();
        double q = 1-p;
        return (p*1)*((double)this.getParam(2).getValue());
    }
    public double getSkew() {
        double p = (double)this.getParam(1).getValue();
        double q = 1-p;
        double n = (double)this.getParam(2).getValue();

        return (q-p)/(Math.sqrt(n*p*q));
    }
    public double evalCDF(double t) {
        throw new RuntimeException("Not implemented");
    }
    public double evalLST(double s) {
        throw new RuntimeException("Not implemented");
    }
}
