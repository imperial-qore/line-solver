/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.util.Pair;
import jline.util.matrix.MatrixCell;

import java.io.Serializable;

/**
 * An abstract class for continuous distributions
 */
abstract public class ContinuousDistribution extends Distribution implements Serializable {
    public ContinuousDistribution(String name, int numParam, Pair<Double, Double> support) {
        super(name, numParam, support);
    }

    public abstract double evalLST(double s);

    /**
     * Gets the process representation with actual distribution parameters.
     * Returns a MatrixCell containing the distribution parameters.
     *
     * @return MatrixCell with distribution-specific parameters
     */
    public abstract MatrixCell getProcess();

}
