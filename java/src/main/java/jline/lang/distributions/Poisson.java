package jline.lang.distributions;

import java.io.Serializable;
import java.util.Random;

import jline.util.Matrix;
import jline.util.Pair;

import static org.apache.commons.math3.util.CombinatoricsUtils.factorial;

public class Poisson extends DiscreteDistribution implements Serializable {
    public Poisson(double rate) {
        super("Poisson", 1, new Pair<Double,Double>(0.0, Double.POSITIVE_INFINITY));
        this.setParam(1, "lambda", rate);
    }

    public double evalPDF(int n) {
        double lambda = this.getRate();
        double num = Math.pow(lambda, n)*Math.exp(-lambda);
        return num/((double) factorial(n));
    }

    public int getRealization(Random random) {
        int curN = 0;
        double cdfVal = evalPDF(curN);
        double cdfProb = random.nextDouble();
        while (cdfVal < cdfProb) {
            curN++;
            cdfVal += evalPDF(curN);
        }
        return curN;
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
        throw new RuntimeException("Not implemented"); // TODO: not implemented
    }


    public double getMean() {
        return this.getRate();
    }
    public double getRate() {
        return (double) this.getParam(1).getValue();
    }
    public double getSCV() {
        throw new RuntimeException("Not implemented"); // TODO: not implemented
    }
    public double getVar() {
        return this.getRate();
    }
    public double getSkew() {
        return Math.pow(this.getRate(), -0.5);
    }
    public double evalCDF(double t) {
        throw new RuntimeException("Not implemented"); // TODO: not implemented
    }
    public double evalLST(double s) {
        throw new RuntimeException("Not implemented"); // TODO: not implemented
    }
}
