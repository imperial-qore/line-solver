package jline.solvers.ctmc;

import java.util.ArrayList;
import java.util.List;

import jline.api.mc.Ctmc_solve_reducible_blkdecomp;
import jline.io.InputOutput;
import jline.lang.NetworkStruct;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;

/**
 * Single entry point for the stationary distribution of a CTMC generated from a
 * NetworkStruct.
 *
 * <p>All the stationary mass of a reducible chain lives in its bottom strongly
 * connected components, each weighted by the probability of being absorbed in it
 * from the declared initial state; every other state is transient and carries
 * zero. The block decomposition handles the irreducible case as the degenerate
 * one BSCC / no transient states, so every CTMC solve in the solver goes through
 * it and there is no dispatch that can disagree with the algorithm about whether
 * a chain is reducible.
 *
 * <p>See _kb/11-conventions-and-gotchas.md.
 */
public class CtmcStationary {

    private CtmcStationary() {
    }

    /**
     * Stationary distribution of {@code Q}, seeded from the initial state of
     * {@code sn} when that state can be located in {@code stateSpace}.
     *
     * @param Q          infinitesimal generator
     * @param stateSpace enumerated state space, rows aligned with Q
     * @param sn         network struct carrying the per-node initial state
     * @param options    solver options, used for debug reporting only
     * @return row vector of stationary probabilities, full length of Q
     */
    public static Matrix solve(Matrix Q, Matrix stateSpace, NetworkStruct sn, SolverOptions options) {
        return Ctmc_solve_reducible_blkdecomp.ctmc_solve_reducible_blkdecomp(Q, initialDistribution(Q, stateSpace, sn, options)).getLeft();
    }

    /**
     * Point mass at the initial state of {@code sn}, or null when that state is
     * absent from {@code stateSpace} (stochastic complementation may have removed
     * it, e.g. an SPN whose immediate ENABLE states were eliminated). A null seed
     * makes the block decomposition start in the SCCs with no incoming transition.
     */
    public static Matrix initialDistribution(Matrix Q, Matrix stateSpace, NetworkStruct sn, SolverOptions options) {
        if (sn == null || sn.state == null || stateSpace == null || stateSpace.getNumRows() == 0) {
            return null;
        }
        List<Double> s0parts = new ArrayList<Double>();
        for (int isf = 0; isf < sn.nstateful; isf++) {
            Matrix stateMatrix = sn.state.get(sn.stateful.get(isf));
            if (stateMatrix != null) {
                for (int ri = 0; ri < stateMatrix.getNumRows(); ri++) {
                    for (int ci = 0; ci < stateMatrix.getNumCols(); ci++) {
                        s0parts.add(stateMatrix.get(ri, ci));
                    }
                }
            }
        }
        if (s0parts.isEmpty() || s0parts.size() != stateSpace.getNumCols()) {
            return null;
        }
        Matrix s0 = new Matrix(1, s0parts.size());
        for (int i = 0; i < s0parts.size(); i++) {
            s0.set(0, i, s0parts.get(i));
        }
        int initStateIdx = Matrix.matchrow(stateSpace, s0);
        if (initStateIdx < 0) {
            if (options != null) {
                InputOutput.line_debug(options.verbose,
                        "Initial state absent from the state space, absorption seeded by the source SCCs");
            }
            return null;
        }
        Matrix pi0 = Matrix.zeros(1, Q.getNumRows());
        pi0.set(0, initStateIdx, 1.0);
        if (options != null) {
            InputOutput.line_debug(options.verbose,
                    String.format("Seeding absorption from initial state %d", initStateIdx));
        }
        return pi0;
    }
}
