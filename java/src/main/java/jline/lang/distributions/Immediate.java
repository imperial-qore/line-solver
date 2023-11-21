package jline.lang.distributions;

import java.io.Serializable;
import java.util.*;

import jline.lang.constant.GlobalConstants;
import jline.util.Matrix;
import jline.util.Pair;

public class Immediate extends Distribution implements Serializable {
	
    public Immediate() {
        super("Immediate", 0, new Pair<Double,Double>(0.0,0.0));
    }

    public boolean isDisabled() {
        return false;
    }

    /**
     * Gets n samples from the distribution
     *
     * @param n - the number of samples
     * @return - n samples from the distribution
     */
    @Override
    public Matrix sample(long n) {
        return this.sample(n,null);
    }

    @Override
    public Matrix sample(long n, Random random) {
        List<Double> ret_list = new ArrayList<Double>();
        for (int i = 0; i < n; i++) {
            ret_list.add(0.0);
        }

        return new Matrix(ret_list);
    }

    public double getRate() {
        return GlobalConstants.Immediate;
    }

    public double getMean() {
        return 1 / GlobalConstants.Immediate;
    }

    public double getSCV() {
        return 0;
    }

    public double getMu() {
        return GlobalConstants.Immediate;
    }

    public double getPhi() {
        return 1;
    }

    public double evalCDF(double t) {
        return 1;
    }

    public double evalLST(double s) {
        return 1;
    }

    public Map<Integer, Matrix> getPH() {
        double t = GlobalConstants.Immediate;
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
        return true;
    }

    public double getSkew() {
        return 0;
    }

    public double getVar() {
        return 0;
    }
}
