/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.qsys;

import jline.util.matrix.Matrix;

/**
 * Result of the matrix-analytic (M/G/1-type) analysis of a BMAP/M/1 queue.
 *
 * <p>Beyond the usual performance measures it exposes the intermediate
 * matrix-analytic quantities themselves, so that the algorithm can be inspected
 * and taught rather than only its output. Port of MATLAB qsys_bmapm1.m.</p>
 */
public class QsysBmapM1Result {
    /** Stationary vector of the BMAP phase process, generator sum_k D_k. */
    public final Matrix theta;
    /** Mean arrival rate, theta * sum_k k*D_k * e. */
    public final double lambda;
    /** Offered load lambda/mu. */
    public final double rho;
    /** Uniformization constant actually used. */
    public final double q;
    /** Randomized block A0 = (mu/q)I: a service completion, level down by one. */
    public final Matrix A0;
    /** Randomized block A1 = (1/q)(D0 - mu*I) + I: level unchanged. */
    public final Matrix A1;
    /** Boundary local block B0 = (1/q)D0 + I, used at level 0 where no service can complete. */
    public final Matrix B0;
    /** Randomized blocks Bk[k] = (1/q)D_k, raising the level by k; index 1..K. */
    public final Matrix[] Bk;
    /** A = A0 + A1 + sum_k Bk[k], the phase process of the randomized chain. */
    public final Matrix A;
    /** Stationary vector of A. */
    public final Matrix alpha;
    /** Minimal non-negative solution of G = A0 + A1*G + sum_k Bk[k]*G^(k+1). */
    public final Matrix G;
    /** alpha*(sum_k k*Bk[k])*e - alpha*A0*e; the queue is stable iff this is negative. */
    public final double drift;
    /**
     * Geometric decay rate of the level probabilities, measured as the limiting
     * ratio pi_(n+1)/pi_n. Reported rather than derived from a spectral convention
     * so that it is unambiguous.
     */
    public final double decayRate;
    /** Level probabilities pi_n as rows (level 0 first). */
    public final Matrix levelProb;
    /** Probability the system is empty (equals 1-rho exactly). */
    public final double pi0;
    /** Mean number in system. */
    public final double meanQueueLength;
    /** Server utilization, equal to rho. */
    public final double utilization;
    /** Throughput, equal to lambda. */
    public final double throughput;
    /** Level truncation used for the level distribution. */
    public final int truncLevel;
    /** Relative truncation residual at truncLevel. */
    public final double truncError;
    /** Name of the analyzer. */
    public final String analyzer;

    public QsysBmapM1Result(Matrix theta, double lambda, double rho, double q, Matrix A0, Matrix A1,
                            Matrix B0, Matrix[] Bk, Matrix A, Matrix alpha, Matrix G, double drift,
                            double decayRate, Matrix levelProb, double pi0, double meanQueueLength,
                            double utilization, double throughput, int truncLevel, double truncError,
                            String analyzer) {
        this.theta = theta;
        this.lambda = lambda;
        this.rho = rho;
        this.q = q;
        this.A0 = A0;
        this.A1 = A1;
        this.B0 = B0;
        this.Bk = Bk;
        this.A = A;
        this.alpha = alpha;
        this.G = G;
        this.drift = drift;
        this.decayRate = decayRate;
        this.levelProb = levelProb;
        this.pi0 = pi0;
        this.meanQueueLength = meanQueueLength;
        this.utilization = utilization;
        this.throughput = throughput;
        this.truncLevel = truncLevel;
        this.truncError = truncError;
        this.analyzer = analyzer;
    }
}
