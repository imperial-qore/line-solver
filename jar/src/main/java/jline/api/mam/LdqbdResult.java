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
    private final List<Matrix> piCells;

    public LdqbdResult(List<Matrix> R, Matrix pi) {
        this(R, pi, null);
    }

    public LdqbdResult(List<Matrix> R, Matrix pi, List<Matrix> piCells) {
        this.R = R;
        this.pi = pi;
        this.piCells = piCells;
    }

    public List<Matrix> getR() { return R; }
    public Matrix getPi() { return pi; }

    /**
     * Per-level stationary vectors, still resolved by phase. getPi() is their
     * row sums. A caller that needs a phase-conditional quantity -- the mean of
     * something that depends on the phase and not only on the level -- cannot
     * recover it from the aggregated vector and must read these.
     *
     * @return one row vector per level, or null when the solver took the scalar
     *         fast path and there are no phases to resolve
     */
    public List<Matrix> getPiCells() { return piCells; }

    public List<Matrix> component1() { return R; }
    public Matrix component2() { return pi; }
}
