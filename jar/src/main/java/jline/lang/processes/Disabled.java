/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.util.Pair;

import java.io.Serializable;
import java.util.Random;

/**
 * A special distribution to denote disabled service or arrival.
 */

public class Disabled extends Distribution implements Serializable {

    static Disabled disabledInstance = null;

    public Disabled() {
        super("Disabled", 0, new Pair<Double, Double>(Double.NaN, Double.NaN));
    }

    public static Disabled instance() {
        return getInstance();
    }

    public static Disabled getInstance() {
        if (disabledInstance == null) {
            disabledInstance = new Disabled();
        }
        return disabledInstance;
    }

    public double evalCDF(double t) {
        return Double.NaN;
    }

    public double evalLST(double s) {
        return Double.NaN;
    }

    public double getMean() {
        return Double.NaN;
    }

    public double getMu() {
        return Double.NaN;
    }

    public double getPhi() {
        return Double.NaN;
    }

    public double getRate() {
        return Double.NaN;
    }

    public double getSCV() {
        return Double.NaN;
    }

    public double getSkewness() {
        return Double.NaN;
    }

    public double getVar() {
        return Double.NaN;
    }

    public boolean isDisabled() {
        return true;
    }

    public boolean isImmediate() {
        return false;
    }

    @Override
    public double[] sample(int n) {
        return this.sample(n, null);
    }

    @Override
    public double[] sample(int n, Random random) {
        double[] ret_list = new double[(int) n];
        for (int i = 0; i < n; i++) {
            ret_list[i] = Double.NaN;
        }

        return ret_list;
    }
}
