/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.qsys;

/**
 * Conditional downstream waiting time in a two-station tandem.
 *
 * <p>Returned by {@link Qsys_mm1_tandem_lindley}, one entry per supplied pair of
 * current waiting times.
 *
 * @since LINE 3.1.0
 */
public class QsysTandemLindleyResult {
    /** Conditional mean downstream waiting time of the next customer. */
    public final double[] mean;
    /** Conditional mean upstream interdeparture time, {@code 1/mu1 + q/lambda}. */
    public final double[] interdepMean;
    /** Probability the upstream server goes idle before the next arrival. */
    public final double[] idleProb;

    public QsysTandemLindleyResult(double[] mean, double[] interdepMean, double[] idleProb) {
        this.mean = mean;
        this.interdepMean = interdepMean;
        this.idleProb = idleProb;
    }
}
