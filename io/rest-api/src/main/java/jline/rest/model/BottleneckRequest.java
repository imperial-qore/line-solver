/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.rest.model;

/**
 * Request DTO for bottleneck detection analysis.
 */
public class BottleneckRequest {

    private ModelInput model;
    private String solver;
    private String analysis;
    private SolverOptionsDTO options;

    /**
     * Utilization threshold for bottleneck classification (default 0.8 = 80%).
     */
    private Double threshold;

    /**
     * Whether to include detailed metrics for all stations.
     */
    private Boolean includeAll;

    public ModelInput getModel() {
        return model;
    }

    public void setModel(ModelInput model) {
        this.model = model;
    }

    public String getSolver() {
        return solver;
    }

    public void setSolver(String solver) {
        this.solver = solver;
    }

    public String getAnalysis() {
        return analysis;
    }

    public void setAnalysis(String analysis) {
        this.analysis = analysis;
    }

    public SolverOptionsDTO getOptions() {
        return options;
    }

    public void setOptions(SolverOptionsDTO options) {
        this.options = options;
    }

    public Double getThreshold() {
        return threshold;
    }

    public void setThreshold(Double threshold) {
        this.threshold = threshold;
    }

    public Boolean getIncludeAll() {
        return includeAll;
    }

    public void setIncludeAll(Boolean includeAll) {
        this.includeAll = includeAll;
    }

    /**
     * Get effective threshold (default 0.8).
     */
    public double getEffectiveThreshold() {
        return threshold != null ? threshold : 0.8;
    }

    /**
     * Check if all stations should be included (default false).
     */
    public boolean shouldIncludeAll() {
        return includeAll != null && includeAll;
    }

    /**
     * Validate the request.
     * @return Error message if invalid, null if valid
     */
    public String validate() {
        if (model == null) {
            return "Model is required";
        }
        String modelError = model.validate();
        if (modelError != null) {
            return modelError;
        }
        if (solver == null || solver.isEmpty()) {
            return "Solver is required";
        }
        if (threshold != null && (threshold <= 0 || threshold > 1)) {
            return "Threshold must be between 0 (exclusive) and 1 (inclusive)";
        }
        return null;
    }
}
