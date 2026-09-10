/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.qsys;

/**
 * Time-varying measures of the Mt/G/infinity queue, as produced by
 * {@link Qsys_mtginf}.
 *
 * <p>Port of the struct returned by MATLAB qsys_mtginf.m.
 *
 * @since LINE 3.1.0
 */
public class QsysMtginfResult {
    /** The times at which everything below was evaluated. */
    public final double[] times;
    /** m(t), the Poisson mean of the number in system. */
    public final double[] meanNumber;
    /** Equal to meanNumber, the law being Poisson. */
    public final double[] varNumber;
    /** lambda(t). */
    public final double[] arrivalRate;
    /** delta(t) = E[lambda(t-S)], the rate of the (Poisson) departure process. */
    public final double[] departureRate;
    /** ES*lambda(t), the pointwise stationary approximation. */
    public final double[] offeredLoadPSA;
    /** E[Se] = E[S^2]/(2 ES), the time lag; NaN when ES2 was not supplied. */
    public final double meanLag;
    /** ES*lambda(t-E[Se]), the first-order lag approximation; null without ES2. */
    public final double[] lagApproximation;

    public QsysMtginfResult(double[] times, double[] meanNumber, double[] varNumber,
                            double[] arrivalRate, double[] departureRate, double[] offeredLoadPSA,
                            double meanLag, double[] lagApproximation) {
        this.times = times;
        this.meanNumber = meanNumber;
        this.varNumber = varNumber;
        this.arrivalRate = arrivalRate;
        this.departureRate = departureRate;
        this.offeredLoadPSA = offeredLoadPSA;
        this.meanLag = meanLag;
        this.lagApproximation = lagApproximation;
    }

    @Override
    public String toString() {
        return "QsysMtginfResult(points=" + (times == null ? 0 : times.length)
                + ", meanLag=" + meanLag + ")";
    }
}
