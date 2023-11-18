package jline.lang.distributions;

import jline.lang.constant.GlobalConstants;
import jline.util.Matrix;

import java.io.Serializable;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import static java.lang.Math.exp;

public class Exp extends MarkovianDistribution  implements Serializable {
    public Exp(double lambda) {
        super("Exp", 1);
        this.setParam(1, "lambda", lambda);
    }

    @Override
    public List<Double> sample(long n) {
        return this.sample(n,null);
    }

    @Override
    public List<Double> sample(long n, Random rand)  {
        double lambda = (double)this.getParam(1).getValue();
        //return exprnd(1/lambda, n, 1);
        throw new RuntimeException("Not Implemented!");
    }

    public long getNumberOfPhases() {
        return 1;
    }

    public double evalCDF(double t) {
        double lambda = (double) this.getParam(1).getValue();
        return 1-exp(-lambda*t);
    }

    public Map<Integer, Matrix> getPH() {
        double lambda = (double) this.getParam(1).getValue();
        Matrix D0 = new Matrix(1,1,1);
        Matrix D1 = new Matrix(1,1,1);
        D0.set(0, 0, -lambda);
        D1.set(0, 0, lambda);
        Map<Integer, Matrix> res = new HashMap<Integer, Matrix>();
        res.put(0, D0);
        res.put(1, D1);
        return res;
    }

    public double evalLST(double s) {
        double lambda = (double) this.getParam(1).getValue();
        return (lambda/(lambda+s));
    }

    public void updateRate(double rate){
        this.setParam(1, "lambda", rate);
        this.mean = 1.0/rate;
        this.immediate = 1.0/rate < GlobalConstants.FineTol;
    }

    public double getSCV() {
        return 1;
    }

    public double getRate() {
        return (double) this.getParam(1).getValue();
    }

    public double getMean() {
        return 1/getRate();
    }

    public double getVar() {
        return 1/(Math.pow(getRate(),2));
    }

    public double getSkew() {
        return 2;
    }

    public String toString() {
        return String.format("jline.Exp(%f)", this.getRate());
    }

    public static Exp fitMean(double mean) {
        return new Exp(1/mean);
    }
}
