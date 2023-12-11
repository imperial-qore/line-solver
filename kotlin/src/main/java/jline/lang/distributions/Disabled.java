package jline.lang.distributions;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import jline.util.Matrix;
import jline.util.Pair;

public class Disabled extends Distribution implements Serializable {

    static Disabled disabledInstance = null;

    public Disabled() {
        super("Disabled", 0, new Pair<Double,Double>(Double.NaN,Double.NaN));
    }

    public static Disabled getInstance() {
        if (disabledInstance == null) {
            disabledInstance = new Disabled();
        }
        return disabledInstance;
    }
    public boolean isDisabled() {
        return true;
    }

    @Override
    public Matrix sample(long n) {
        return this.sample(n,null);
    }
    @Override
    public Matrix sample(long n, Random random) {
        List<Double> ret_list = new ArrayList<Double>();
        for (int i = 0; i < n; i++) {
            ret_list.add(Double.NaN);
        }

        return new Matrix(ret_list);
    }

    public double getRate() {
        return Double.NaN;
    }

    public double getMean() {
        return Double.NaN;
    }

    public double getSCV() {
        return Double.NaN;
    }

    public double getMu() {
        return Double.NaN;
    }

    public double getPhi() {
        return Double.NaN;
    }

    public double evalCDF(double t) {
        return Double.NaN;
    }

    public double evalLST(double s) {
        return Double.NaN;
    }

    public boolean isImmediate() {
        return false;
    }

    public double getSkew() {
        return Double.NaN;
    }

    public double getVar() {
        return Double.NaN;
    }
}
