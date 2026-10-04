/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.fluid.moments;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import jline.lang.NetworkStruct;
import jline.util.matrix.Matrix;

/**
 * Detects a moment-closure trajectory that has left the model.
 *
 * <p>WHY IT EXISTS. The moment-closure drift can leave the simplex: on a station
 * where {@code min(n,c)} is not the identity the Gaussian correction to the
 * per-class share can drive a coordinate negative, and since the drift is
 * conservative another grows to match. In MATLAB {@code odeset('NonNegative')}
 * projects the ACCEPTED step, so the excursion is CLAMPED rather than reported --
 * which injects mass, collapses the step size, and leaves the window never
 * returning. One MATLAB suite run sat in {@code test_CQN_Cox_CS_9} for 3h16m and
 * the 2026-08-27 run was killed after {@code test11_interlock_lqnx} had held it
 * for 100 minutes.</p>
 *
 * <p>THIS PORT CANNOT HANG THE SAME WAY, and the difference is worth stating
 * rather than papering over: the JAR never projects the clamp back into the
 * integration ({@code PassageTimeODE.computeDerivatives} evaluates the drift at
 * the integrator's own state, and {@code MethodStepHandler} clamps only what it
 * RECORDS), so the same divergence surfaces as an LSODA corrector failure or as
 * a finished window holding a state that is not a solution. The second of those
 * is silent, and this is what makes it loud. The check is therefore applied to
 * the window's recorded rows rather than as an integration-halting callback.</p>
 *
 * <p>THE TEST IS AN EXACT INVARIANT, not a heuristic bound on time or magnitude.
 * The drift conserves the population of every CLOSED CHAIN exactly, so any
 * deviation is a divergence and nothing else. The tolerance is a generous
 * fraction of that population rather than a numerical tolerance: the
 * integrator's own error is ~1e-4 relative, while the documented excursion
 * reaches 5.2e4 against a true population of 0.05. A closed model whose
 * population has moved by TOL is no longer solving the model, whatever it is
 * converging to.</p>
 *
 * <p>THE CHAIN IS THE CONSERVED UNIT, NOT THE CLASS, and the difference is the
 * whole correctness of this check. {@code sn.njobs(k)} is the population class
 * k STARTS with; class switching then moves jobs between the classes of one
 * chain, so only the chain total is invariant. Watching classes instead
 * condemns every class-switching model out of hand -- measured on
 * {@code cqn_twoclass_hyperl} (313 of 447 accepted states), on
 * {@code init_state_ps} (286 of 310) and on every one of the 162 fluid layers
 * an LQN builds under the {@code srvn.cs} encoding, where the chain sum never
 * moved at all. A cache model is the same story with the hit/miss classes.</p>
 *
 * <p>A wall-clock budget would have caught the same thing and was rejected: it
 * makes the answer depend on how busy the host is, so the same model would fall
 * back on one machine and not on another. This invariant is deterministic.</p>
 *
 * @see jline.solvers.fluid.moments.FluidNonHyperbolicException
 */
public class FluidConservationGuard {

    /** Relative population drift that counts as having left the model. */
    public static final double TOL = 0.1;

    private final int[][] idx;      // coordinate columns per watched chain
    private final double[] target;  // that chain's conserved population

    /**
     * @param sn     the network struct, for the chain membership and populations
     * @param phases (nstations,nclasses) phases per pair, as the analyzer builds
     *               it from mu; this is what fixes the coordinate blocks
     */
    public FluidConservationGuard(NetworkStruct sn, Matrix phases) {
        int M = phases.getNumRows();
        int K = phases.getNumCols();
        List<int[]> idxList = new ArrayList<int[]>();
        List<Double> targetList = new ArrayList<Double>();
        for (int c = 0; c < chainCount(sn, K); c++) {
            int[] members = chainClasses(sn, c, K);
            double total = 0;
            for (int m = 0; m < members.length; m++) {
                int k = members[m];
                total += (k < sn.njobs.length()) ? sn.njobs.get(k) : Double.POSITIVE_INFINITY;
            }
            if (!Double.isFinite(total) || total <= 0) {
                continue; // open, or absent: no conserved population to check
            }
            int width = 0;
            for (int m = 0; m < members.length; m++) {
                for (int i = 0; i < M; i++) {
                    width += (int) phases.get(i, members[m]);
                }
            }
            if (width == 0) {
                continue;
            }
            int[] cols = new int[width];
            int at = 0;
            for (int m = 0; m < members.length; m++) {
                int k = members[m];
                for (int i = 0; i < M; i++) {
                    int ph = (int) phases.get(i, k);
                    if (ph <= 0) {
                        continue;
                    }
                    int shift = (int) phases.sumSubMatrix(0, i, 0, K)
                            + (int) phases.sumSubMatrix(i, i + 1, 0, k);
                    for (int f = 0; f < ph; f++) {
                        cols[at++] = shift + f;
                    }
                }
            }
            Arrays.sort(cols);
            idxList.add(cols);
            targetList.add(Double.valueOf(total));
        }
        this.idx = idxList.toArray(new int[idxList.size()][]);
        this.target = new double[targetList.size()];
        for (int c = 0; c < this.target.length; c++) {
            this.target[c] = targetList.get(c).doubleValue();
        }
    }

    /**
     * Chains declared by {@code sn}, or one chain per class when it declares
     * none. The fallback is the safe reading and not a guess: with no chain
     * map there is no class switching to merge classes, so each class IS its
     * own conserved unit.
     */
    private static int chainCount(NetworkStruct sn, int K) {
        if (sn.chains == null || sn.chains.getNumRows() == 0) {
            return K;
        }
        return sn.chains.getNumRows();
    }

    /** The classes of chain {@code c}, as column indices into {@code phases}. */
    private static int[] chainClasses(NetworkStruct sn, int c, int K) {
        if (sn.chains == null || sn.chains.getNumRows() == 0) {
            return new int[]{c};
        }
        int n = 0;
        for (int k = 0; k < K && k < sn.chains.getNumCols(); k++) {
            if (sn.chains.get(c, k) != 0) {
                n++;
            }
        }
        int[] out = new int[n];
        int at = 0;
        for (int k = 0; k < K && k < sn.chains.getNumCols(); k++) {
            if (sn.chains.get(c, k) != 0) {
                out[at++] = k;
            }
        }
        return out;
    }

    /** True when the closure is what is being integrated, i.e. some sigma2 is nonzero. */
    public static boolean closureActive(jline.solvers.SolverOptions options) {
        if (options == null || options.config == null || options.config.moment_sigma2 == null) {
            return false;
        }
        for (int i = 0; i < options.config.moment_sigma2.length; i++) {
            if (options.config.moment_sigma2[i] != 0) {
                return true;
            }
        }
        return false;
    }

    /**
     * @param x a state vector
     * @return the chain whose conserved population has drifted past TOL, or -1
     */
    public int violation(double[] x) {
        for (int c = 0; c < idx.length; c++) {
            int[] cols = idx[c];
            if (cols.length == 0 || cols[cols.length - 1] >= x.length) {
                continue; // a caller integrating a different state vector
            }
            double mass = 0;
            for (int j = 0; j < cols.length; j++) {
                mass += x[cols[j]];
            }
            if (Math.abs(mass - target[c]) > TOL * Math.max(1.0, target[c])) {
                return c;
            }
        }
        return -1;
    }

    /**
     * Raises {@link FluidNonHyperbolicException} when the state has left the
     * model, so the fallback ladder in SolverFluid answers with 'dae' and then
     * with 'matrix'/'closing' instead of returning a state that is not a
     * solution.
     *
     * @param x the window's final state
     * @param t the time that window reached
     */
    public void assertConserved(double[] x, double t) {
        int c = violation(x);
        if (c >= 0) {
            throw new FluidNonHyperbolicException(String.format(
                    "The moment-closure drift left the model: closed chain %d lost more than %.0f%% "
                            + "of its population by t = %g, which the drift conserves exactly, so the "
                            + "excursion is a divergence rather than a solution. Falling back to a "
                            + "first-order closure.", c, 100 * TOL, t));
        }
    }
}
