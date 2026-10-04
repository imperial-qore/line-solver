/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.qsys;

import java.io.Serializable;

import jline.util.matrix.Matrix;

/**
 * Return value of {@link Qsys_mmapg1k}: the per-class throughputs and loss
 * ratios of an MMAP[K]/G/1/K buffer, mirroring the MATLAB result struct of
 * {@code qsys_mmapg1k.m}.
 */
public final class QsysMmapG1kResult implements Serializable {

    private static final long serialVersionUID = 1L;

    /** 1 x R, throughput of class k [pkts/s]. */
    public Matrix throughput;
    /** 1 x R, loss ratio of class k, in [0,1]. */
    public Matrix lossRatio;
    /** 1 x R, arrival rate of class k [pkts/s]. */
    public Matrix lambda;
    /** sum_k lambda(k). */
    public double lambdaAggregate;
    /** sum_k throughput(k). */
    public double throughputAggregate;
    /** Aggregate loss ratio. */
    public double lossAggregate;
    /** P(buffer empty). */
    public double p0;
    /** P(buffer full). */
    public double pK;
    /** 1 x M, P(buffer full, phase j). */
    public Matrix pKvec;
    /** 1 x (K+1), P(level = l). */
    public Matrix plevel;
    /** E[number in system]. */
    public double meanQueueLength;
    /** E[service time]. */
    public double meanServiceTime;
    /** 1 - p0. */
    public double utilization;
    /** Offered load. */
    public double rho;
    /** Name of the analyzer that produced this result. */
    public String analyzer = "qsys_mmapg1k";
}
