/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.fluid.moments;

/**
 * The fluid fixed point has no stationary linear noise approximation: the drift
 * Jacobian is not Hurwitz on the reachable subspace (balanced bottlenecks, a
 * saturated multiclass station, an overloaded open station).
 *
 * <p>Typed so the caller can tell it apart from any other failure.
 * {@link FluidMinNormalApplicable} cannot see this in advance -- it exists only
 * once the mean is solved -- so {@code SolverFluid} switches a RESOLVED
 * {@code minnormal} to the first-order method on this class alone. An explicit
 * {@code options.method="minnormal"} still propagates. MATLAB twin: the
 * {@code LINE:FluidNonHyperbolic} identifier raised by {@code fluid_lyapunov.m}.</p>
 */
public class FluidNonHyperbolicException extends RuntimeException {

    private static final long serialVersionUID = 1L;

    public FluidNonHyperbolicException(String message) {
        super(message);
    }
}
