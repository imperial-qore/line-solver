/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.lib.fjcodes;

import jline.util.matrix.Matrix;

/**
 * Result of computePi.
 */
public final class PiResult {
    public final Matrix pi0;   // Steady-state probability (row vector)
    public final double En1;   // Expected number in system

    public PiResult(Matrix pi0, double En1) {
        this.pi0 = pi0;
        this.En1 = En1;
    }

    public Matrix getPi0() { return pi0; }
    public double getEn1() { return En1; }
}
