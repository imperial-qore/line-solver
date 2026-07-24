/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.rest.model;

import java.util.List;
import java.util.Map;

/**
 * Response DTO for sensitivity analysis results.
 */
public class SensitivityResponse {

    private String status;
    private String solver;
    private double runtime;
    private double delta;
    private boolean normalized;
    private Map<String, Object> baselineMetrics;
    private List<ParameterSensitivity> sensitivities;
    private String error;

    public String getStatus() {
        return status;
    }

    public void setStatus(String status) {
        this.status = status;
    }

    public String getSolver() {
        return solver;
    }

    public void setSolver(String solver) {
        this.solver = solver;
    }

    public double getRuntime() {
        return runtime;
    }

    public void setRuntime(double runtime) {
        this.runtime = runtime;
    }

    public double getDelta() {
        return delta;
    }

    public void setDelta(double delta) {
        this.delta = delta;
    }

    public boolean isNormalized() {
        return normalized;
    }

    public void setNormalized(boolean normalized) {
        this.normalized = normalized;
    }

    public Map<String, Object> getBaselineMetrics() {
        return baselineMetrics;
    }

    public void setBaselineMetrics(Map<String, Object> baselineMetrics) {
        this.baselineMetrics = baselineMetrics;
    }

    public List<ParameterSensitivity> getSensitivities() {
        return sensitivities;
    }

    public void setSensitivities(List<ParameterSensitivity> sensitivities) {
        this.sensitivities = sensitivities;
    }

    public String getError() {
        return error;
    }

    public void setError(String error) {
        this.error = error;
    }

    /**
     * Sensitivity of metrics with respect to a single parameter.
     */
    public static class ParameterSensitivity {
        private String parameter;
        private double baseValue;
        private Map<String, Object> derivatives;
        private String error;

        public String getParameter() {
            return parameter;
        }

        public void setParameter(String parameter) {
            this.parameter = parameter;
        }

        public double getBaseValue() {
            return baseValue;
        }

        public void setBaseValue(double baseValue) {
            this.baseValue = baseValue;
        }

        public Map<String, Object> getDerivatives() {
            return derivatives;
        }

        public void setDerivatives(Map<String, Object> derivatives) {
            this.derivatives = derivatives;
        }

        public String getError() {
            return error;
        }

        public void setError(String error) {
            this.error = error;
        }
    }
}
