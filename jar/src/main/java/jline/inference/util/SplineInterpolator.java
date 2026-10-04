/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.inference.util;

import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

/**
 * Spline interpolation utility using Apache Commons Math3.
 */
public final class SplineInterpolator {
    private SplineInterpolator() {}

    /**
     * Interpolate data at new x-values using cubic spline interpolation.
     *
     * @param xOld original x-values (sorted, distinct)
     * @param yOld original y-values
     * @param xNew new x-values at which to evaluate
     * @return interpolated y-values at xNew
     */
    public static double[] interpolate(double[] xOld, double[] yOld, double[] xNew) {
        if (xOld.length < 2) {
            // Not enough points for interpolation, return constant
            double[] out = new double[xNew.length];
            double constVal = (yOld.length > 0) ? yOld[0] : 0.0;
            for (int i = 0; i < xNew.length; i++) {
                out[i] = constVal;
            }
            return out;
        }

        org.apache.commons.math3.analysis.interpolation.SplineInterpolator interpolator =
                new org.apache.commons.math3.analysis.interpolation.SplineInterpolator();
        PolynomialSplineFunction function = interpolator.interpolate(xOld, yOld);

        double[] out = new double[xNew.length];
        double xMin = xOld[0];
        double xMax = xOld[xOld.length - 1];
        for (int i = 0; i < xNew.length; i++) {
            double x = xNew[i];
            double xClamped = Math.max(xMin, Math.min(xMax, x));
            out[i] = function.value(xClamped);
        }
        return out;
    }
}
