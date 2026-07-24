/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.GlobalConstants;
import jline.util.Pair;
import jline.util.matrix.Matrix;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

/**
 * An Immediate distribution that always samples 0.
 */

public class Immediate extends Distribution implements Serializable {

    static Immediate immediateInstance = null;

    public Immediate() {
        super("Immediate", 0, new Pair<Double, Double>(GlobalConstants.Immediate, GlobalConstants.Immediate));
    }

    public static Immediate getInstance() {
        if (immediateInstance == null) {
            immediateInstance = new Immediate();
        }
        return immediateInstance;
    }

    public double evalCDF(double t) {
        return 1;
    }

    public double evalLST(double s) {
        return 1;
    }

    public double getMean() {
        // Return 0 to match MATLAB behavior (Immediate.getMean returns 0)
        return 0;
    }

    public double getMu() {
        return GlobalConstants.Immediate;
    }

    public Map<Integer, Matrix> getPH() {
        double t = GlobalConstants.Immediate;
        Matrix D0 = new Matrix(1, 1, 1);
        Matrix D1 = new Matrix(1, 1, 1);
        D0.set(0, 0, -t);
        D1.set(0, 0, t);
        Map<Integer, Matrix> res = new HashMap<Integer, Matrix>();
        res.put(0, D0);
        res.put(1, D1);
        return res;
    }

    public double getPhi() {
        return 1;
    }

    public double getRate() {
        return GlobalConstants.Immediate;
    }

    public double getSCV() {
        // Return 1 to match native Python behavior which converts Immediate to
        // Distribution(mean=0, scv=1.0) - this produces correct LQNS results
        return 1;
    }

    public double getSkewness() {
        return 0;
    }

    public double getVar() {
        return 0;
    }

    public boolean isDisabled() {
        return false;
    }

    public boolean isImmediate() {
        return true;
    }

    /**
     * Gets n samples from the distribution
     *
     * @param n - the number of samples
     * @return - n samples from the distribution
     */
    @Override
    public double[] sample(int n) {
        return this.sample(n, null);
    }

    @Override
    public double[] sample(int n, Random random) {
        double[] ret_list = new double[(int) n];
        for (int i = 0; i < n; i++) {
            ret_list[i] = 0.0;
        }

        return ret_list;
    }
}
