/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.mam;

import java.util.List;

import jline.util.matrix.Matrix;

/**
 * Result of LDQBD solver containing rate matrices and stationary distribution.
 */
public final class LdqbdResult {
    private final List<Matrix> R;
    private final Matrix pi;

    public LdqbdResult(List<Matrix> R, Matrix pi) {
        this.R = R;
        this.pi = pi;
    }

    public List<Matrix> getR() { return R; }
    public Matrix getPi() { return pi; }

    public List<Matrix> component1() { return R; }
    public Matrix component2() { return pi; }
}
