package jline.lang.distributions;

import jline.lang.constant.GlobalConstants;
import jline.util.Matrix;
import jline.util.Pair;
import org.apache.commons.lang3.NotImplementedException;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * Zipf-like probability distribution
 */
public class Zipf extends DiscreteDistribution{

    public Zipf(double s){
        this(s, (int) GlobalConstants.Immediate);
    }

    /**
     * Construct a Zipf-like distribution
     * @param s - the shape
     * @param n - number of items
     */
    public Zipf(double s, int n){
        super("Zipf", 4, new Pair<Double,Double>(1.0, (double) n));
        Matrix p = new Matrix(1, n);
        double h = Zipf.genHarmonic(s, n);
        Matrix x = new Matrix(1, n);
        for(int i = 0; i < p.getNumCols(); i++){
            p.set(i, 1/Math.pow(i+1, s)/h);
            x.set(i, i+1);
        }
        setParam(1, "p", p);
        setParam(2, "x", x);
        setParam(3, "s", s);
        setParam(4, "n", n);
    }

    /**
     * Computes the distribution mean
     * @return - the mean of the distribution
     */
    @Override
    public double getMean() {
        double s = (double) this.getParam(3).getValue();
        int n = (int) this.getParam(4).getValue();
        return genHarmonic(s-1, n) / genHarmonic(s, n);
    }

    @Override
    public double getRate() {
        throw new NotImplementedException("getRate() not implemented in Zipf");
    }

    /**
     * Computes the squared coefficient of variation == variance/mean^2
     * @return - the squared coefficient of variation
     */
    @Override
    public double getSCV() {
        double s = (double) this.getParam(3).getValue();
        int n = (int) this.getParam(4).getValue();
        double ex = getMean();
        double var = genHarmonic(s-2, n) / genHarmonic(s, n) - ex * ex;
        return var/(ex * ex);
    }

    @Override
    public double getVar() {
        throw new NotImplementedException("getVar() not implemented in Zipf");
    }

    @Override
    public double getSkew() {
        throw new NotImplementedException("getSkew() not implemented in Zipf");
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

    /**
     * Evaluates the cumulative distribution function at t
     * @param t - the point where the cdf will be evaluated
     * @return - the cdf evaluated at t
     */
    @Override
    public double evalCDF(double t) {
        double s = (double) this.getParam(3).getValue();
        int n = (int) this.getParam(4).getValue();
        return genHarmonic(s, (int) t) / genHarmonic(s, n);
    }

    @Override
    public double evalLST(double s) {
        throw new NotImplementedException("evalLST() not implemented in Zipf");
    }

    public Matrix evalPMF(){
        int n = (int) this.getParam(4).getValue();
        List<Double> t = new ArrayList<>();
        for(int i = 1; i <= n; i++){
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
    public Matrix evalPMF(List<Double> t){
        double s = (double) this.getParam(3).getValue();
        int n = (int) this.getParam(4).getValue();
        List<Double> retList = new ArrayList<>();
        double Hns = genHarmonic(s, n);
        for(double d : t){
            retList.add(1/Math.pow(d, s)/Hns);
        }
        return new Matrix(retList);
    }

    /**
     * Generate harmonic numbers to normalize a Zipf-like distribution on n items with shape s
     */
    public static double genHarmonic(double s, int n){
        double Hnm = 0;
        for(int k = 1; k <= n; k++){
            Hnm += 1.0/Math.pow(k, s);
        }
        return Hnm;
    }
}
