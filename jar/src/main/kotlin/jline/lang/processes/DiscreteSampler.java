/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.util.Pair;
import jline.util.RandomManager;
import jline.util.matrix.Matrix;
import org.apache.commons.math3.util.DoubleArray;
import org.apache.commons.math3.util.FastMath;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * A class for discrete distributions specified from the probability mass function
 */
public class DiscreteSampler extends DiscreteDistribution {

    public DiscreteSampler(Matrix p) {
        super("DiscreteSampler", 3, new Pair<Double, Double>(1.0, (double) p.length()));
        setParam(1, "p", p.columnMajorOrder().transpose());
        int n = p.length();
        Matrix x = new Matrix(1, n);
        for (int i = 0; i < n; i++) {
            x.set(i, i + 1);
        }
        setParam(2, "x", x.columnMajorOrder().transpose());
        Matrix f = p.columnMajorOrder().transpose().cumsumViaRow();
        double psum = p.elementSum();
        for (int i = 0; i < f.getNumCols(); i++) {
            f.set(i, f.get(i) / psum);
        }
        setParam(3, "f", f);
        setParam(4, "min", (int) x.elementMin());
        setParam(5, "max", (int) x.elementMax());
    }

    /**
     * Constructs a discrete distribution from a finite probability vector p at the points specified in vector x
     *
     * @param p - the probability of an item
     * @param x - the value of an item
     */
    public DiscreteSampler(Matrix p, Matrix x) {
        super("DiscreteSampler", 3, new Pair<Double, Double>(x.elementMin(), x.elementMax()));
        setParam(1, "p", p.columnMajorOrder().transpose());
        setParam(2, "x", x.columnMajorOrder().transpose());
        Matrix f = p.columnMajorOrder().transpose().cumsumViaRow();
        double psum = p.elementSum();
        for (int i = 0; i < f.getNumCols(); i++) {
            f.set(i, f.get(i) / psum);
        }
        setParam(3, "f", f);
    }

    public static void main(String[] args) {
        DiscreteSampler z = new DiscreteSampler(new Matrix("[0.2,0.5,0.3]"), new Matrix("[5,11,9]"));
        double[] S = z.sample(10000);
        int c1 = 0;
        for (int i = 0; i < S.length; i++) {
            if (S[i] == 9)
                c1++;
        }
        System.out.println((double) c1 / S.length);
    }

    @Override
    public double evalCDF(double t) {
        Matrix f = (Matrix) this.getParam(3).getValue();
        if (t >= 0 && t < f.length()) {
            return f.get((int) t);
        } else {
            return 0;
        }
    }

    public List<Double> evalPMF() {
        Matrix p = (Matrix) this.getParam(1).getValue();
        List<Double> retList = new ArrayList<>();
        for (int i = 0; i < p.getNumCols(); i++) {
            retList.add(p.get(i));
        }
        return retList;
    }

    @Override
    public Matrix evalPMF(List<Double> t) {
        Matrix p = (Matrix) this.getParam(1).getValue();
        Matrix x = (Matrix) this.getParam(2).getValue();
        List<Double> retList = new ArrayList<>();
        for (int i = 0; i < t.size(); i++) {
            int j = 0;
            while (j < x.length() && x.get(j) != t.get(i)) {
                j++;
            }
            retList.add(p.get(j));
        }
        return new Matrix(retList);
    }

    /**
     * Computes the distribution mean
     *
     * @return - the mean of the distribution
     */
    @Override
    public double getMean() {
        Matrix p = (Matrix) this.getParam(1).getValue();
        Matrix x = (Matrix) this.getParam(2).getValue();
        double mean = 0.0;
        for (int i = 0; i < p.length(); i++)
            mean += x.get(i) * p.get(i);
        return mean;
    }

    /**
     * Computes the distribution squared coefficient of variation (SCV = variance/mean^2)
     *
     * @return
     */
    @Override
    public double getSCV() {
        Matrix p = (Matrix) this.getParam(1).getValue();
        Matrix x = (Matrix) this.getParam(2).getValue();
        double e2 = 0.0;
        for (int k = 0; k < x.length(); k++) {
            e2 += FastMath.pow(x.get(k), 2) * p.get(k);
        }
        double ex = getMean();
        double var = e2 - ex * ex;
        return var / (ex * ex);
    }

    @Override
    public double getSkewness() {
        Matrix p = (Matrix) this.getParam(1).getValue();
        Matrix x = (Matrix) this.getParam(2).getValue();
        double e3 = 0.0;
        double e1 = getMean();
        double e2 = 0.0;
        for (int k = 0; k < x.length(); k++) {
            e2 += FastMath.pow(x.get(k), 2) * p.get(k);
        }

        for (int k = 0; k < x.length(); k++) {
            e3 += FastMath.pow(x.get(k), 3) * p.get(k);
        }

        return (e3 - 3 * e1 * e2 + 2 * e1 * e1 * e1) / FastMath.pow(getVar(), 1.5);
    }

    @Override
    public boolean isDisabled() {
        Matrix p = (Matrix) this.getParam(1).getValue();
        for (int i = 0; i < p.length(); i++) {
            if (Double.isNaN(p.get(i)))
                return true;
        }
        return false;
    }

    public double sample() {
        return this.sample(1, RandomManager.getThreadRandomAsRandom())[0];
    }

    @Override
    public double[] sample(int n) {
        return this.sample(n, RandomManager.getThreadRandomAsRandom());
    }

    public double[] sample(int n, Random random) {
        double[] x = ((Matrix) this.getParam(2).getValue()).toArray1D();
        double[] f = ((Matrix) this.getParam(3).getValue()).toArray1D();
        double[] ret = new double[n];
        for (int i = 0; i < n; i++) {
            double r = random.nextDouble();
            for (int q = 0; q < f.length; q++) {
                if (r <= f[q]) {
                    ret[i] = x[q];
                    break;
                }
            }
        }
        return ret;
    }
}
