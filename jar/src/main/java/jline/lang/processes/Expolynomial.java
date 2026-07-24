/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.util.Pair;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.io.Serializable;
import java.util.Arrays;
import java.util.Random;

/**
 * Expolynomial distribution with density f(x) = sum ci * x^ai * exp(-li*x).
 *
 * <p>Represents an expolynomial density over a bounded domain [eft, lft],
 * matching the Sirio/ORIS GEN expolynomial format. The density is carried as a
 * Sirio expression string rather than as numeric parameters, so it is opaque to
 * this class: the class stores and transports it verbatim.
 *
 * <p>The parameter slots mirror the MATLAB reference implementation
 * (matlab/src/lang/processes/Expolynomial.m):
 * <ul>
 *   <li>param 1 - {@code density}, the Sirio density expression string</li>
 *   <li>param 2 - {@code eft}, the earliest firing time (lower bound of support)</li>
 *   <li>param 3 - {@code lft}, the latest firing time (upper bound of support,
 *       possibly {@link Double#POSITIVE_INFINITY})</li>
 * </ul>
 *
 * <p><b>Moments.</b> As in MATLAB and in the Python-native implementation, the
 * moments of an expolynomial density are not obtained by numerical integration:
 * {@link #getMean()}, {@link #getSCV()}, {@link #getSkewness()},
 * {@link #evalCDF(double)}, {@link #evalLST(double)} and {@link #sample(int, Random)}
 * all return NaN. The distribution is therefore a transport-only type in the
 * lang layer: it carries the density to engines that can interpret it (Sirio/ORIS
 * GEN), and every moment-based solver sees NaN rather than a silently wrong value.
 *
 * <p><b>Support.</b> Unlike MATLAB, which passes a NaN placeholder triple to the
 * superclass constructor, the support is reported here as the pair (eft, lft),
 * matching the {@code getSupport} of the Python-native implementation.
 */
public class Expolynomial extends ContinuousDistribution implements Serializable {

    private final String density;
    private final double eft;
    private final double lft;

    /**
     * Creates an expolynomial distribution.
     *
     * @param density density expression string in Sirio format
     * @param eft     earliest firing time (lower bound of support)
     * @param lft     latest firing time (upper bound of support, may be
     *                {@link Double#POSITIVE_INFINITY})
     */
    public Expolynomial(String density, double eft, double lft) {
        super("Expolynomial", 3, new Pair<Double, Double>(eft, lft));
        if (density == null) {
            throw new IllegalArgumentException("Expolynomial: density expression cannot be null");
        }
        this.density = density;
        this.eft = eft;
        this.lft = lft;
        setParam(1, "density", density);
        setParam(2, "eft", eft);
        setParam(3, "lft", lft);
        this.mean = Double.NaN;
        this.immediate = false;
    }

    /** Returns the density expression string in Sirio format. */
    public String getDensity() {
        return density;
    }

    /** Returns the earliest firing time, i.e. the lower bound of the support. */
    public double getEft() {
        return eft;
    }

    /**
     * Returns the latest firing time, i.e. the upper bound of the support.
     * May be {@link Double#POSITIVE_INFINITY} for an unbounded density.
     */
    public double getLft() {
        return lft;
    }

    /** Returns NaN: the mean is not obtained by numerical integration. */
    @Override
    public double getMean() {
        return Double.NaN;
    }

    /** Returns NaN: the SCV is not obtained by numerical integration. */
    @Override
    public double getSCV() {
        return Double.NaN;
    }

    /** Returns NaN: the skewness is not obtained by numerical integration. */
    @Override
    public double getSkewness() {
        return Double.NaN;
    }

    /** Returns NaN: the CDF of a Sirio density expression is not evaluated here. */
    @Override
    public double evalCDF(double t) {
        return Double.NaN;
    }

    /** Returns NaN: the LST of a Sirio density expression is not evaluated here. */
    @Override
    public double evalLST(double s) {
        return Double.NaN;
    }

    /**
     * Returns the numeric part of the process representation, {eft, lft}, as two
     * singleton matrices. The density expression is a string and hence cannot be
     * held in a MatrixCell; obtain it with {@link #getDensity()}.
     */
    @Override
    public MatrixCell getProcess() {
        MatrixCell representation = new MatrixCell();
        representation.set(0, Matrix.singleton(eft));
        representation.set(1, Matrix.singleton(lft));
        return representation;
    }

    /** Returns an array of n NaNs: sampling a Sirio density expression is not supported. */
    @Override
    public double[] sample(int n, Random random) {
        double[] samples = new double[n];
        Arrays.fill(samples, Double.NaN);
        return samples;
    }

    @Override
    public String toString() {
        return String.format("jline.Expolynomial(%s, eft=%f, lft=%f)", density, eft, lft);
    }
}
