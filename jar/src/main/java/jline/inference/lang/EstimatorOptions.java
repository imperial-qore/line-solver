/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.inference.lang;

import java.util.function.Function;

import jline.lang.Network;
import jline.solvers.NetworkSolver;

/**
 * Options for ParamEstimator.
 */
public class EstimatorOptions {
    public int verbose = 1;
    public String method = "ubr";
    public String variant = "default";
    public int iterMax = 1000;
    public double tol = 1e-3;
    public Function<Network, NetworkSolver> solverFactory = null;
    public int openPopulation = 100;
    public double[] x0 = null;

    public EstimatorOptions() {}

    public static EstimatorOptions defaultOptions() {
        return new EstimatorOptions();
    }
}
