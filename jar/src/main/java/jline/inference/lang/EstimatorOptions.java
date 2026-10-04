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

    /** Options of the variational estimator ("vi"). */
    public jline.inference.api.VariationalOptions variational =
            new jline.inference.api.VariationalOptions();
    /** Probability that a queue-length reading is faulty ("vi"). */
    public double epsilon = 0.05;
    /** Shape of the Gamma prior placed on each estimated rate ("vi"). */
    public double priorShape = 1.0;
    /** Posterior Gamma shapes left behind by the variational estimator. */
    public double[] posteriorAlpha = null;
    /** Posterior Gamma rates left behind by the variational estimator. */
    public double[] posteriorBeta = null;
    /** Evidence lower bound per iteration, left behind by the variational estimator. */
    public double[] bound = null;

    public EstimatorOptions() {}

    public static EstimatorOptions defaultOptions() {
        return new EstimatorOptions();
    }
}
