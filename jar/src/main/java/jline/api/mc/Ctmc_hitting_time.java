/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.mc;

import jline.util.matrix.Matrix;

/**
 * Mean time to reach any state in a target set from each state of a CTMC.
 *
 * <p>Continuous-time twin of {@link Dtmc_hitting_time} and the first-moment
 * special case of {@link Ctmc_passage_moments}: (-S) h = 1 on the non-target
 * block, where {@code dtmc_hitting_time} solves (I - P_NT) h = 1. Target states
 * have zero hitting time; a state that cannot reach the set has an infinite one.
 *
 * <p>Target indices are 0-based here and 1-based in the MATLAB reference.
 */
public final class Ctmc_hitting_time {

    private Ctmc_hitting_time() {
    }

    public static Matrix ctmc_hitting_time(Matrix Q, int[] targetStates) {
        int n = Q.getNumRows();
        // mall does not depend on the initial law, so a uniform one is passed
        // rather than requiring the caller to invent one.
        Matrix pi0 = new Matrix(1, n);
        for (int i = 0; i < n; i++) {
            pi0.set(0, i, 1.0 / n);
        }
        PassageMomentsResult pm =
                Ctmc_passage_moments.ctmc_passage_moments(Q, pi0, targetStates, 1);
        Matrix h = new Matrix(n, 1);
        for (int i = 0; i < n; i++) {
            h.set(i, 0, pm.mall.get(i, 0));
        }
        return h;
    }
}
