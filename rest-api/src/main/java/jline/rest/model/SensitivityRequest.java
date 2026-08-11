/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.rest.model;

import java.util.List;

/**
 * Request DTO for sensitivity analysis.
 * Computes partial derivatives of metrics with respect to parameters.
 */
public class SensitivityRequest {

    private ModelInput model;
    private String solver;
    private String analysis;
    private List<ParameterSweep> parameters;
    private List<String> metrics;
    private SolverOptionsDTO options;

    /**
     * Perturbation size for numerical differentiation (default 0.01 = 1%).
     */
    private Double delta;

    /**
     * Whether to normalize sensitivities (elasticity).
     */
    private Boolean normalized;

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

    public List<ParameterSweep> getParameters() {
        return parameters;
    }

    public void setParameters(List<ParameterSweep> parameters) {
        this.parameters = parameters;
    }

    public List<String> getMetrics() {
        return metrics;
    }

    public void setMetrics(List<String> metrics) {
        this.metrics = metrics;
    }

    public SolverOptionsDTO getOptions() {
        return options;
    }

    public void setOptions(SolverOptionsDTO options) {
        this.options = options;
    }

    public Double getDelta() {
        return delta;
    }

    public void setDelta(Double delta) {
        this.delta = delta;
    }

    public Boolean getNormalized() {
        return normalized;
    }

    public void setNormalized(Boolean normalized) {
        this.normalized = normalized;
    }

    /**
     * Get effective delta value (default 0.01).
     */
    public double getEffectiveDelta() {
        return delta != null ? delta : 0.01;
    }

    /**
     * Check if normalization is enabled (default false).
     */
    public boolean isNormalized() {
        return normalized != null && normalized;
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
        if (parameters == null || parameters.isEmpty()) {
            return "At least one parameter is required for sensitivity analysis";
        }
        for (ParameterSweep param : parameters) {
            // For sensitivity, we only need the parameter definition, not values
            String paramError = param.validateForSensitivity();
            if (paramError != null) {
                return paramError;
            }
        }
        if (delta != null && (delta <= 0 || delta >= 1)) {
            return "Delta must be between 0 and 1 (exclusive)";
        }
        return null;
    }
}
