/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.util.Pair;
import jline.util.matrix.Matrix;
import org.apache.commons.lang3.NotImplementedException;

import java.io.Serializable;
import java.util.List;

/**
 * An abstract class for discrete distributions.
 */

public abstract class DiscreteDistribution extends Distribution implements Serializable {
    public DiscreteDistribution(String name, int numParam, Pair<Double, Double> support) {
        super(name, numParam, support);
    }

    public double evalPMF(double k) {
        throw new NotImplementedException("This method must be implemented by successor classes");
    }

    public double[] evalPMF(double[] k) {
        double[] ret = new double[k.length];
        for (int i = 0; i < k.length; i++) {
            ret[i] = evalPMF(k[i]);
        }
        return ret;
    }

    public Matrix evalPMF(List<Double> t) {
        double[] array = new double[t.size()];
        for (int i = 0; i < t.size(); i++) {
            array[i] = t.get(i);
        }
        return new Matrix(evalPMF(array));
    }

    /**
     * Evaluate the Laplace-Stieltjes Transform at s
     * For discrete distributions, this is the probability generating function evaluated at e^(-s)
     * 
     * @param s the Laplace domain variable
     * @return the LST value at s
     */
    public double evalLST(double s) {
        throw new NotImplementedException("This method must be implemented by successor classes");
    }
}
