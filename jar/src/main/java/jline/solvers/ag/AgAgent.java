/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.ag;

import jline.solvers.ag.handlers.RCATModel;
import jline.solvers.ag.handlers.Solver_ag_inap;
import jline.util.matrix.Matrix;

/**
 * The single definition of what one agent's answer is.
 *
 * <p>Every execution backend routes through here -- the serial loop, the thread
 * pool and the remote ag-worker alike -- so a remote agent's answer is the same
 * object as a local one rather than a second implementation that happens to
 * agree. Keeping one definition is what makes the identity assertions between
 * backends meaningful: if they were two code paths, the assertions would only be
 * testing that the two paths had not drifted yet.</p>
 */
public final class AgAgent {

    private AgAgent() {}

    /**
     * Agent k's generator at the reversed rates x. O(N^2); the stationary solve
     * that follows is O(N^3), which is why the cluster backend transports the
     * latter's result and recomputes this locally.
     */
    public static Matrix generator(int k, Matrix x, Matrix[] Aa, Matrix[] Pb, Matrix[] L,
                                   int[] ACT, int[] PSV, int numActions, int[] N) {
        return Solver_ag_inap.agentGenerator(k, x, Aa, Pb, L, ACT, PSV, numActions, N);
    }

    /** Agent k's stationary vector, given its generator. */
    public static Matrix stationary(Matrix Qk, RCATModel rcat, int k) {
        return Solver_ag_inap.agentStationary(Qk, rcat, k);
    }

    /**
     * Agent k's stationary vector from the QBD shape carried explicitly, for a
     * worker that holds the agent's layout without holding the whole model.
     */
    public static Matrix stationary(Matrix Qk, int mph, int nlev, int[] level) {
        return Solver_ag_inap.agentStationaryOf(Qk, mph, nlev, level);
    }
}
