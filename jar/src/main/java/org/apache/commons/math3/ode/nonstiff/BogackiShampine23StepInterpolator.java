/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package org.apache.commons.math3.ode.nonstiff;

import org.apache.commons.math3.ode.sampling.StepInterpolator;

/**
 * Dense output of {@link BogackiShampine23Integrator}.
 *
 * <p>The pair is first-same-as-last, so the drift is known at both ends of the
 * step: {@code yDotK[0]} at the start and {@code yDotK[3]} at the end. The cubic
 * Hermite polynomial through those two values and the two states is third-order
 * accurate, which is the interpolant {@code ntrp23} uses in MATLAB.</p>
 *
 * <p>The two branches below are the same polynomial written from the start and
 * from the end of the step, as the Commons Math interpolators do, so that
 * {@code theta = 0} and {@code theta = 1} return the stored states exactly. The
 * increment across the step is not stored, it is rebuilt from the stages as
 * {@code h * (2/9 k1 + 1/3 k2 + 4/9 k3)}, which is the propagation formula.</p>
 *
 * @see BogackiShampine23Integrator
 * @since 1.0
 */
class BogackiShampine23StepInterpolator extends RungeKuttaStepInterpolator {

    /** Serializable version identifier. */
    private static final long serialVersionUID = 20260731L;

    /** Propagation weight of the first stage. */
    private static final double B1 = 2.0 / 9.0;

    /** Propagation weight of the second stage. */
    private static final double B2 = 1.0 / 3.0;

    /** Propagation weight of the third stage. */
    private static final double B3 = 4.0 / 9.0;

    /** Builds an uninitialized interpolator, used as the integrator's prototype. */
    public BogackiShampine23StepInterpolator() {
        super();
    }

    /**
     * Copy constructor.
     *
     * @param interpolator interpolator to copy
     */
    BogackiShampine23StepInterpolator(final BogackiShampine23StepInterpolator interpolator) {
        super(interpolator);
    }

    /** {@inheritDoc} */
    @Override
    protected StepInterpolator doCopy() {
        return new BogackiShampine23StepInterpolator(this);
    }

    /** {@inheritDoc} */
    @Override
    protected void computeInterpolatedStateAndDerivatives(final double theta,
                                                          final double oneMinusThetaH) {

        final double theta2 = theta * theta;
        // Hermite basis: h01 weights the increment, h10 and h11 the two drifts.
        final double h01 = theta2 * (3.0 - 2.0 * theta);
        final double h10 = theta * (theta - 1.0) * (theta - 1.0);
        final double h11 = theta2 * (theta - 1.0);
        final double dh01 = 6.0 * theta * (1.0 - theta);
        final double dh10 = (3.0 * theta - 1.0) * (theta - 1.0);
        final double dh11 = theta * (3.0 * theta - 2.0);

        final double d1 = dh01 * B1 + dh10;
        final double d2 = dh01 * B2;
        final double d3 = dh01 * B3;
        final double d4 = dh11;

        if ((previousState != null) && (theta <= 0.5)) {
            final double c1 = h01 * B1 + h10;
            final double c2 = h01 * B2;
            final double c3 = h01 * B3;
            final double c4 = h11;
            for (int i = 0; i < interpolatedState.length; ++i) {
                interpolatedState[i] = previousState[i]
                        + h * (c1 * yDotK[0][i] + c2 * yDotK[1][i]
                             + c3 * yDotK[2][i] + c4 * yDotK[3][i]);
                interpolatedDerivatives[i] = d1 * yDotK[0][i] + d2 * yDotK[1][i]
                                           + d3 * yDotK[2][i] + d4 * yDotK[3][i];
            }
        } else {
            final double e = 1.0 - h01;
            final double c1 = e * B1 - h10;
            final double c2 = e * B2;
            final double c3 = e * B3;
            final double c4 = -h11;
            for (int i = 0; i < interpolatedState.length; ++i) {
                interpolatedState[i] = currentState[i]
                        - h * (c1 * yDotK[0][i] + c2 * yDotK[1][i]
                             + c3 * yDotK[2][i] + c4 * yDotK[3][i]);
                interpolatedDerivatives[i] = d1 * yDotK[0][i] + d2 * yDotK[1][i]
                                           + d3 * yDotK[2][i] + d4 * yDotK[3][i];
            }
        }
    }
}
