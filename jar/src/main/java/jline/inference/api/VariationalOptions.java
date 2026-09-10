/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.inference.api;

/**
 * Options of {@link Infer_variational}. A null-valued box means "derive a
 * default from the specification"; see Infer_variational for the rules.
 */
public class VariationalOptions {
    public int verbose = 0;
    public int iterMax = 20;
    public double tol = 1e-3;
    public int nsamples = 200;
    /** rate added to every feasible transition by the space expansion. */
    public double delta = 1e-3;
    public double floor = 1e-4;
    /** cap on the variational rates; derived from ymax and tmax when null. */
    public Double rateMax = null;
    public double rateCapFactor = 10.0;
    public double unifmax = 30.0;
    public double unifTol = 1e-12;
    public int unifMaxTerms = 2000;
    public Double tmax = null;
    public Double dt = null;
    public Integer ngrid = null;
    public Integer ymax = null;

    public VariationalOptions() {}

    public static VariationalOptions defaultOptions() {
        return new VariationalOptions();
    }
}
