/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.qsys;

/**
 * Tail bounds of a two-station tandem and the coefficients that produced them.
 *
 * <p>Returned by {@link Qsys_tandem_ub_ciucu}, one {@code S} and {@code W} entry
 * per supplied threshold.
 *
 * @since LINE 3.1.0
 */
public class QsysTandemUbResult {
    /** Upper bound on P(S &gt; x), the end-to-end sojourn time, capped at one. */
    public final double[] S;
    /** Upper bound on P(W &gt; x); NaN unless the service law is exponential. */
    public final double[] W;
    /** Tail decay rate, the positive root of E[e^{theta (Y-X)}] = 1. */
    public final double theta;
    /** E[X e^{-theta X}], the arrival functional the coefficients depend on. */
    public final double alpha;
    /** Coefficient A of the test function gamma, fixed by Lemma 4. */
    public final double A;
    /** Coefficient B of the test function gamma, fixed by Lemma 4. */
    public final double B;
    /** Coefficient C of the test function gamma, fixed by Lemma 4. */
    public final double C;
    /** Coefficient D of the test function gamma, fixed by Lemma 4. */
    public final double D;

    public QsysTandemUbResult(double[] S, double[] W, double theta, double alpha,
                              double A, double B, double C, double D) {
        this.S = S;
        this.W = W;
        this.theta = theta;
        this.alpha = alpha;
        this.A = A;
        this.B = B;
        this.C = C;
        this.D = D;
    }
}
