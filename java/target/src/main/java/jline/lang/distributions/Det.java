package jline.lang.distributions;

import jline.util.Matrix;
import jline.util.Pair;

import java.io.Serializable;
import java.util.*;

public class Det extends Distribution implements Serializable {

    public Det(double t) {
        super("Det", 1, new Pair<Double,Double>(t, t));
        this.setParam(1, "t", t);
    }

    public boolean isDisabled() {
        return false;
    }
    @Override
    public List<Double> sample(long n) {
        return this.sample(n,null);
    }
    @Override
    public List<Double> sample(long n, Random random) {
        List<Double> ret_list = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            ret_list.add((double)this.getParam(1).getValue());
        }
        return ret_list;
//        throw new RuntimeException("Not implemented");
    }

    public double getRate() {
        return 1 / (double) this.getParam(1).getValue();
    }

    public double getMean() {
        return (double) this.getParam(1).getValue();
    }

    public double getSCV() {
        return 0;
    }

    public double getMu() {
        return 1 / (double) this.getParam(1).getValue();
    }

    public double getPhi() {
        return 1;
    }

    public double evalCDF(double t) {
        // Evaluate the cumulative distribution function at t
        if (t < (double) this.getParam(1).getValue()) {
            return 0;
        } else {
            return 1;
        }
    }

    public double evalLST(double s) {
        // Evaluate the Laplace-Stieltjes transform of the distribution function at t
        return Math.exp(-s * (double) this.getParam(1).getValue());
    }

    public Map<Integer, Matrix> getPH() {
        double t = (double) this.getParam(1).getValue();
        Matrix D0 = new Matrix(1,1,1);
        Matrix D1 = new Matrix(1,1,1);
        D0.set(0, 0, -t);
        D1.set(0, 0, t);
        Map<Integer, Matrix> res = new HashMap<Integer, Matrix>();
        res.put(0, D0);
        res.put(1, D1);
        return res;
    }

    public boolean isImmediate() {
        return false;
    }

    public double getSkew() {
        return 0;
    }

    public double getVar() {
        return 0;
    }

    @Override
    public boolean isContinuous() {
        return true;
    }

    @Override
    public boolean isDiscrete() {
        return true;
    }


}
