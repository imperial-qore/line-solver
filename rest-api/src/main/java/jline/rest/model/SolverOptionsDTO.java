/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.rest.model;

import jline.solvers.SolverOptions;

/**
 * Data Transfer Object for solver options in REST API requests.
 */
public class SolverOptionsDTO {

    /**
     * Random seed for stochastic solvers.
     */
    private Integer seed;

    /**
     * Whether to output verbose information.
     */
    private Boolean verbose;

    /**
     * Number of samples for simulation-based solvers.
     */
    private Integer samples;

    /**
     * Tolerance for iterative solvers.
     */
    private Double tolerance;

    /**
     * Maximum number of iterations for iterative solvers.
     */
    private Integer maxIterations;

    /**
     * Method/variant to use for the solver.
     */
    private String method;

    /**
     * Timespan for transient analysis [start, end].
     */
    private double[] timespan;

    /**
     * Default constructor for JSON deserialization.
     */
    public SolverOptionsDTO() {
    }

    public Integer getSeed() {
        return seed;
    }

    public void setSeed(Integer seed) {
        this.seed = seed;
    }

    public Boolean getVerbose() {
        return verbose;
    }

    public void setVerbose(Boolean verbose) {
        this.verbose = verbose;
    }

    public Integer getSamples() {
        return samples;
    }

    public void setSamples(Integer samples) {
        this.samples = samples;
    }

    public Double getTolerance() {
        return tolerance;
    }

    public void setTolerance(Double tolerance) {
        this.tolerance = tolerance;
    }

    public Integer getMaxIterations() {
        return maxIterations;
    }

    public void setMaxIterations(Integer maxIterations) {
        this.maxIterations = maxIterations;
    }

    public String getMethod() {
        return method;
    }

    public void setMethod(String method) {
        this.method = method;
    }

    public double[] getTimespan() {
        return timespan;
    }

    public void setTimespan(double[] timespan) {
        this.timespan = timespan;
    }

    /**
     * Convert this DTO to a SolverOptions instance.
     * @return The corresponding SolverOptions
     */
    public SolverOptions toSolverOptions() {
        SolverOptions options = new SolverOptions();

        if (seed != null) {
            options.seed = seed;
        }
        if (verbose != null) {
            options.verbose = verbose ? jline.VerboseLevel.STD : jline.VerboseLevel.SILENT;
        }
        if (samples != null) {
            options.samples = samples;
        }
        if (tolerance != null) {
            options.iter_tol = tolerance;
        }
        if (maxIterations != null) {
            options.iter_max = maxIterations;
        }
        if (method != null) {
            options.method = method;
        }
        if (timespan != null && timespan.length == 2) {
            options.timespan = timespan;
        }

        return options;
    }
}
