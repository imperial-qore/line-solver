/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.qsys;

import java.io.Serializable;

import jline.util.matrix.Matrix;

/**
 * Return value of {@link Qsys_mapg1k_perflow}, mirroring the MATLAB result
 * struct of {@code qsys_mapg1k_perflow.m}.
 */
public final class QsysMapG1kPerflowResult implements Serializable {

    private static final long serialVersionUID = 1L;

    /** 1 x N, throughput of flow n [pkts/s]. */
    public Matrix throughput;
    /** 1 x N, loss ratio of flow n, in [0,1]. */
    public Matrix lossRatio;
    /** 1 x N, arrival rate of flow n [pkts/s]. */
    public Matrix lambda;
    /** sum_n lambda(n). */
    public double lambdaAggregate;
    /** sum_n throughput(n). */
    public double throughputAggregate;
    /** sum_n lossRatio(n)*lambda(n)/lambdaAggregate. */
    public double lossAggregate;
    /** 1 x N, P(buffer empty) in the n-th model. */
    public Matrix p0;
    /** 1 x N, P(buffer full) in the n-th model. */
    public Matrix pK;
    /** E[service time]. */
    public double meanServiceTime;
    /** lambdaAggregate * E[S]. */
    public double rho;
    /** Name of the analyzer that produced this result. */
    public String analyzer = "qsys_mapg1k_perflow";
}
