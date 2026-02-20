/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.util.matrix.Matrix;

import java.io.Serializable;

/**
 * An abstract class for marked point processes
 */

abstract public class Marked extends MarkovModulated implements Serializable {
    // A Markov-modulated Markovian point process, possibly non-renewal

    public Marked(String name, int numParam) {
        super(name, numParam);
    }

    public Matrix D(int i) {
        return this.process.get(i);
    }

    public Matrix D(int i, int k) {
        if (i == 1) {
            return this.process.get(2 + k);
        } else {
            return null;
        }
    }

    public Matrix getD1k(int k) {
        return (Matrix) this.getParam(2 + k).getValue();
    }

}
