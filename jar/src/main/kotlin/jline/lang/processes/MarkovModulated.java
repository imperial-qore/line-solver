/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.io.Serializable;

import static jline.api.mam.Map_acfcKt.map_acfc;
import static jline.api.mam.Map_gammaKt.map_gamma;

/**
 * An abstract class for a Markov-modulated point-process
 */
abstract public class MarkovModulated extends Markovian implements Serializable {
    // A Markov-modulated Markovian point process, possibly non-renewal

    public MarkovModulated(String name, int numParam) {
        super(name, numParam);
    }

    public Matrix evalACFT(int[] lags, double timescale) {
        return new Matrix(map_acfc(D(0), D(1), lags, timescale));
    }

    public double getACFDecay() {
        MatrixCell D = getProcess();
        return map_gamma(D.get(0), D.get(1));
    }

}
