/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package org.apache.commons.math3.ode.nonstiff;

import org.apache.commons.math3.util.FastMath;

/**
 * Bogacki-Shampine 3(2) embedded Runge-Kutta integrator, the method MATLAB's
 * {@code ode23} implements.
 *
 * <p>LINE's fluid solver selects a fast and an accurate non-stiff integrator by
 * tolerance, mirroring {@code ode_solve.m}, whose fast handle is {@code @ode23}.
 * Apache Commons Math 3.6.1 ships no order-3 embedded Runge-Kutta method, so the
 * JAR previously stood {@code HighamHall54Integrator} in that slot and ran the
 * fast arm at order 5(4). This class closes that order mismatch.</p>
 *
 * <p>The Butcher tableau is the one published in Bogacki and Shampine, "A 3(2)
 * pair of Runge-Kutta formulas", Appl. Math. Lett. 2(4):321-325, 1989: three
 * stages advance the solution and a fourth, first-same-as-last stage evaluates
 * the drift at the new point, serving both the second-order error estimate and
 * the first stage of the next step.</p>
 *
 * <p>This class lives in the Commons Math package because
 * {@link EmbeddedRungeKuttaIntegrator}'s constructor takes a
 * {@code RungeKuttaStepInterpolator}, which is package private and therefore not
 * nameable from any other package; the code itself is original to LINE.</p>
 *
 * @see BogackiShampine23StepInterpolator
 * @since 1.0
 */
public class BogackiShampine23Integrator extends EmbeddedRungeKuttaIntegrator {

    /** Integrator name reported by {@link #getName()}. */
    private static final String METHOD_NAME = "Bogacki-Shampine 3(2)";

    /** Time steps of the internal stages, as fractions of the step size. */
    private static final double[] STATIC_C = {
        1.0 / 2.0, 3.0 / 4.0, 1.0
    };

    /** Internal weights of the stages. */
    private static final double[][] STATIC_A = {
        {1.0 / 2.0},
        {0.0, 3.0 / 4.0},
        {2.0 / 9.0, 1.0 / 3.0, 4.0 / 9.0}
    };

    /** Propagation weights of the third-order solution. */
    private static final double[] STATIC_B = {
        2.0 / 9.0, 1.0 / 3.0, 4.0 / 9.0, 0.0
    };

    /** Error weight of the first stage, third-order weight minus second-order weight. */
    private static final double E1 = 2.0 / 9.0 - 7.0 / 24.0;

    /** Error weight of the second stage. */
    private static final double E2 = 1.0 / 3.0 - 1.0 / 4.0;

    /** Error weight of the third stage. */
    private static final double E3 = 4.0 / 9.0 - 1.0 / 3.0;

    /** Error weight of the fourth, first-same-as-last stage. */
    private static final double E4 = 0.0 - 1.0 / 8.0;

    /**
     * Builds an integrator with scalar tolerances.
     *
     * @param minStep minimal step, the sign is irrelevant
     * @param maxStep maximal step, the sign is irrelevant
     * @param scalAbsoluteTolerance allowed absolute error
     * @param scalRelativeTolerance allowed relative error
     */
    public BogackiShampine23Integrator(final double minStep, final double maxStep,
                                       final double scalAbsoluteTolerance,
                                       final double scalRelativeTolerance) {
        super(METHOD_NAME, true, STATIC_C, STATIC_A, STATIC_B,
              new BogackiShampine23StepInterpolator(), minStep, maxStep,
              scalAbsoluteTolerance, scalRelativeTolerance);
    }

    /**
     * Builds an integrator with per-component tolerances.
     *
     * @param minStep minimal step, the sign is irrelevant
     * @param maxStep maximal step, the sign is irrelevant
     * @param vecAbsoluteTolerance allowed absolute error, one entry per component
     * @param vecRelativeTolerance allowed relative error, one entry per component
     */
    public BogackiShampine23Integrator(final double minStep, final double maxStep,
                                       final double[] vecAbsoluteTolerance,
                                       final double[] vecRelativeTolerance) {
        super(METHOD_NAME, true, STATIC_C, STATIC_A, STATIC_B,
              new BogackiShampine23StepInterpolator(), minStep, maxStep,
              vecAbsoluteTolerance, vecRelativeTolerance);
    }

    /**
     * Returns the order of the propagated solution, which sets the step-size
     * control exponent, matching the {@code 1/3} power {@code ode23} uses.
     *
     * @return 3
     */
    @Override
    public int getOrder() {
        return 3;
    }

    /** {@inheritDoc} */
    @Override
    protected double estimateError(final double[][] yDotK, final double[] y0, final double[] y1,
                                   final double h) {
        double error = 0;
        for (int j = 0; j < mainSetDimension; ++j) {
            final double errSum = E1 * yDotK[0][j] + E2 * yDotK[1][j]
                                + E3 * yDotK[2][j] + E4 * yDotK[3][j];
            final double yScale = FastMath.max(FastMath.abs(y0[j]), FastMath.abs(y1[j]));
            final double tol = (vecAbsoluteTolerance == null)
                    ? (scalAbsoluteTolerance + scalRelativeTolerance * yScale)
                    : (vecAbsoluteTolerance[j] + vecRelativeTolerance[j] * yScale);
            final double ratio = h * errSum / tol;
            error += ratio * ratio;
        }
        return FastMath.sqrt(error / mainSetDimension);
    }
}
