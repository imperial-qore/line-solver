/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.util.matrix.Matrix;

/**
 * Fits a matrix exponential to a given mean and squared coefficient of variation.
 *
 * Mirrors MATLAB dist_fit_me.m and the native Python fit_me_mean_scv.
 */
public final class MEFit {

    private MEFit() {}

    /**
     * Fits a matrix exponential with the given mean and SCV, without a phase budget.
     *
     * @param mean target mean, positive
     * @param scv  target squared coefficient of variation, in (0, 1)
     * @return the fitted ME
     */
    public static ME fitMeanAndSCV(double mean, double scv) {
        return fitMeanAndSCV(mean, scv, 0);
    }

    /**
     * Fits a matrix exponential with the given mean and SCV.
     *
     * For {@code scv < 1} the fit is the convolution {@code X = c*Y + Z} of a scaled
     * concentrated matrix exponential Y (unit mean, minimal SCV sY for its order) with an
     * independent exponential Z. Writing {@code c + d = mean} and
     * {@code c^2*sY + d^2 = scv*mean^2},
     *
     * <pre>c = mean*(1 - sqrt(1 - (1+sY)*(1-scv)))/(1 + sY),   d = mean - c,</pre>
     *
     * so every target in {@code [sY/(1+sY), 1]} is matched EXACTLY in {@code 2n+2}
     * phases. The exponential tail is what makes the convolution reach up to SCV 1; the
     * concentrated part is what makes it reach far below the Erlang bound {@code 1/order}
     * at the same order.
     *
     * The order is the smallest tabulated one that reaches the target, capped by
     * {@code maxPhases} when given: with a phase budget an Erlang can only reach
     * {@code 1/maxPhases}, while this construction reaches {@code O(1/maxPhases^2)}, and
     * the residual SCV is then the closest achievable from below.
     *
     * {@code scv >= 1} is outside the range of a concentrated ME (its SCV never exceeds
     * 0.34), and the caller keeps its own hyperexponential fit there.
     *
     * @param mean      target mean, positive
     * @param scv       target squared coefficient of variation, in (0, 1)
     * @param maxPhases cap on the number of phases, 0 for no cap
     * @return the fitted ME
     * @throws IllegalArgumentException if the mean or the SCV is out of range
     */
    public static ME fitMeanAndSCV(double mean, double scv, int maxPhases) {
        if (Double.isNaN(mean) || Double.isInfinite(mean) || mean <= 0) {
            throw new IllegalArgumentException("MEFit mean must be a positive finite number, got " + mean);
        }
        if (Double.isNaN(scv) || scv <= 0 || scv >= 1) {
            throw new IllegalArgumentException("MEFit requires 0 < scv < 1; use a hyperexponential "
                    + "for scv >= 1 and a CME for scv = 0, got " + scv);
        }

        // Smallest tabulated order whose convolution range covers the target, subject to
        // the phase budget. reach(order) = sY/(1+sY) is the minimum SCV of the
        // convolution, which sits just below the sY of the CME alone.
        int bestOrder = -1;
        for (int order : CME.getSupportedOrders()) {
            if (maxPhases > 0 && order + 1 > maxPhases) {
                continue;
            }
            double sYcand = CME.getMinSCV(order);
            bestOrder = order; // budget-limited: keep the most concentrated one that fits
            if (sYcand / (1.0 + sYcand) <= scv) {
                break;
            }
        }
        if (bestOrder < 0) {
            throw new IllegalArgumentException("No CME order fits a budget of " + maxPhases
                    + " phases; the smallest is 3 phases plus one exponential");
        }

        Matrix[] rep = CME.representation(bestOrder);
        Matrix alphaY = rep[0];
        Matrix AY = rep[1];
        double sY = CME.getMinSCV(bestOrder);
        double reach = sY / (1.0 + sY);
        double c;
        if (scv < reach) {
            // Budget-limited: the target is below what this order can reach, so the most
            // concentrated member of the family is returned and the caller gets the
            // closest achievable SCV rather than a silent Erlang truncation.
            c = mean / (1.0 + sY);
        } else {
            c = mean * (1.0 - Math.sqrt(1.0 - (1.0 + sY) * (1.0 - scv))) / (1.0 + sY);
        }
        double d = mean - c;

        int n = alphaY.getNumElements();
        if (d <= mean * 1e-12) {
            return new CME(mean, bestOrder);
        }
        if (c <= mean * 1e-12) {
            Matrix alphaExp = new Matrix(1, 1, 1);
            alphaExp.set(0, 0, 1.0);
            Matrix Aexp = new Matrix(1, 1, 1);
            Aexp.set(0, 0, -1.0 / mean);
            return new ME(alphaExp, Aexp);
        }

        // Convolution of two matrix exponentials: the exit flow of the first block feeds
        // the entry of the second, exactly as for a phase-type.
        Matrix alpha = new Matrix(1, n + 1);
        Matrix A = new Matrix(n + 1, n + 1);
        for (int i = 0; i < n; i++) {
            alpha.set(0, i, alphaY.get(0, i));
            double rowsum = 0.0;
            for (int j = 0; j < n; j++) {
                double aij = AY.get(i, j) / c;
                A.set(i, j, aij);
                rowsum += aij;
            }
            A.set(i, n, -rowsum);
        }
        A.set(n, n, -1.0 / d);
        return new ME(alpha, A, false);
    }
}
