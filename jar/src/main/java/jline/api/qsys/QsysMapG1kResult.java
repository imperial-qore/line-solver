/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.qsys;

import java.io.Serializable;

import jline.util.matrix.Matrix;

/**
 * Return value of {@link Qsys_mapg1k}, mirroring the MATLAB result struct of
 * {@code qsys_mapg1k.m} field for field.
 */
public final class QsysMapG1kResult implements Serializable {

    private static final long serialVersionUID = 1L;

    /** P(buffer empty). */
    public double p0;
    /** P(buffer full). */
    public double pK;
    /** Aggregate departure rate [pkts/s]. */
    public double throughput;
    /** 1 - throughput/lambda. */
    public double lossProbability;
    /** Aggregate MAP arrival rate. */
    public double lambda;
    /** E[service time]. */
    public double meanServiceTime;
    /** 1 - p0. */
    public double utilization;
    /** Offered load lambda*S. */
    public double rho;
    /** E[number in system]. */
    public double meanQueueLength;
    /** Uniformization order used. */
    public int nmax;
    /** |1 - sum_n c_n| at the truncation. */
    public double countingResidual;
    /** Stationary law of the embedded chain, 1 x K*M. */
    public Matrix sigma;
    /** 1 x M, P(level = K, phase j), summing to pK. */
    public Matrix pKvec;
    /** 1 x M, P(level = 0, phase j), summing to p0. */
    public Matrix p0vec;
    /** 1 x (K+1), P(level = l), l = 0..K. */
    public Matrix plevel;
    /** Name of the analyzer that produced this result. */
    public String analyzer = "qsys_mapg1k";
}
