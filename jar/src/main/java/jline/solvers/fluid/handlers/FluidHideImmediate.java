/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.fluid.handlers;

import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.solvers.SolverOptions;

/**
 * Whether the fluid drift of a model should be built on the stochastic complement of its
 * INSTANTANEOUS coordinates (see {@link ImmediateElimination}).
 *
 * <p>Every fluid route that builds its drift from the station/class/phase event set asks here rather
 * than reading {@code options.config.hide_immediate} directly, so that the answer is the same one
 * across matrix, closing, statedep, tbi, minnormal, refined and dae. The flag defaults to TRUE for
 * SolverFluid: a coordinate whose exit rate is GlobalConstants.Immediate is LINE's stand-in for
 * infinity, and integrating it is meaningless work no integrator does well.
 *
 * <p>THE STOCHASTIC PETRI NET ROUTE IS THE ONE EXCEPTION, and it is not a refusal. The Petri route
 * carries immediate firings as ALGEBRAIC unknowns of an index-1 DAE, which is a stronger treatment
 * than absorbing them: it keeps the firing flow itself as a solved quantity rather than folding it
 * into the timed events. That route never builds the event set this reduction acts on, so the answer
 * here is simply false and the flag is left alone rather than being turned into an error.
 */
public class FluidHideImmediate {

    private FluidHideImmediate() {
    }

    public static boolean resolve(NetworkStruct sn, SolverOptions options) {
        if (options == null || options.config == null || !options.config.hide_immediate) {
            return false;
        }
        if (sn != null && sn.nodetype != null && sn.nodetype.contains(NodeType.Transition)) {
            return false;
        }
        return true;
    }
}
