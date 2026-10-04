/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.mc;

import jline.util.matrix.Matrix;

/**
 * The phase-type form of a first passage time into a target state set.
 *
 * <p>With A the complement of the target set, S = Q(A,A) is the sub-generator
 * under which the passage has not yet completed, s0 = -S*1 is the exit vector,
 * and alpha = pi0(A). Then L(s) = alpha (sI-S)^-1 s0 + atom and
 * F(t) = 1 - alpha exp(St) 1.
 *
 * <p>ALPHA IS DELIBERATELY NOT NORMALIZED. Its mass is 1 - atom; the missing
 * mass is the ATOM AT ZERO carried by initial states already inside the target.
 * A caller that normalizes alpha and forgets the atom reports F(0) = 0 for a
 * passage that has already completed with probability atom.
 */
public class CtmcPassagePh {
    /** pi0 restricted to the non-target block, UNNORMALIZED. */
    public final Matrix alpha;
    /** Sub-generator Q(A,A). */
    public final Matrix S;
    /** Exit vector -S*1, as a column. */
    public final Matrix s0;
    /** Row of S to state index of Q. */
    public final int[] keep;
    /** Mass of pi0 already inside the target: F(0). */
    public final double atom;

    public CtmcPassagePh(Matrix alpha, Matrix S, Matrix s0, int[] keep, double atom) {
        this.alpha = alpha;
        this.S = S;
        this.s0 = s0;
        this.keep = keep;
        this.atom = atom;
    }
}
