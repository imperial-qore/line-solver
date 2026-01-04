/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.rest.model;

import java.util.List;
import java.util.Map;

/**
 * Response DTO for model calibration results.
 */
public class CalibrationResponse {

    private String status;
    private double runtime;
    private int iterations;
    private double finalError;
    private boolean converged;
    private List<CalibratedParameter> parameters;
    private Map<String, Object> fittedMetrics;
    private Map<String, Object> residuals;
    private String error;

    public String getStatus() {
        return status;
    }

    public void setStatus(String status) {
        this.status = status;
    }

    public double getRuntime() {
        return runtime;
    }

    public void setRuntime(double runtime) {
        this.runtime = runtime;
    }

    public int getIterations() {
        return iterations;
    }

    public void setIterations(int iterations) {
        this.iterations = iterations;
    }

    public double getFinalError() {
        return finalError;
    }

    public void setFinalError(double finalError) {
        this.finalError = finalError;
    }

    public boolean isConverged() {
        return converged;
    }

    public void setConverged(boolean converged) {
        this.converged = converged;
    }

    public List<CalibratedParameter> getParameters() {
        return parameters;
    }

    public void setParameters(List<CalibratedParameter> parameters) {
        this.parameters = parameters;
    }

    public Map<String, Object> getFittedMetrics() {
        return fittedMetrics;
    }

    public void setFittedMetrics(Map<String, Object> fittedMetrics) {
        this.fittedMetrics = fittedMetrics;
    }

    public Map<String, Object> getResiduals() {
        return residuals;
    }

    public void setResiduals(Map<String, Object> residuals) {
        this.residuals = residuals;
    }

    public String getError() {
        return error;
    }

    public void setError(String error) {
        this.error = error;
    }

    /**
     * Calibrated parameter value.
     */
    public static class CalibratedParameter {
        private String name;
        private String station;
        private String className;
        private double initialValue;
        private double fittedValue;
        private double confidence;

        public String getName() {
            return name;
        }

        public void setName(String name) {
            this.name = name;
        }

        public String getStation() {
            return station;
        }

        public void setStation(String station) {
            this.station = station;
        }

        public String getClassName() {
            return className;
        }

        public void setClassName(String className) {
            this.className = className;
        }

        public double getInitialValue() {
            return initialValue;
        }

        public void setInitialValue(double initialValue) {
            this.initialValue = initialValue;
        }

        public double getFittedValue() {
            return fittedValue;
        }

        public void setFittedValue(double fittedValue) {
            this.fittedValue = fittedValue;
        }

        public double getConfidence() {
            return confidence;
        }

        public void setConfidence(double confidence) {
            this.confidence = confidence;
        }
    }
}
