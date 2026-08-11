/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.rest.model;

import java.util.List;
import java.util.Map;

/**
 * Response DTO for what-if analysis results.
 */
public class WhatIfResponse {

    private String status;
    private String solver;
    private double runtime;
    private int totalPoints;
    private int completedPoints;
    private int failedPoints;
    private List<ParameterPoint> results;
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

    public int getTotalPoints() {
        return totalPoints;
    }

    public void setTotalPoints(int totalPoints) {
        this.totalPoints = totalPoints;
    }

    public int getCompletedPoints() {
        return completedPoints;
    }

    public void setCompletedPoints(int completedPoints) {
        this.completedPoints = completedPoints;
    }

    public int getFailedPoints() {
        return failedPoints;
    }

    public void setFailedPoints(int failedPoints) {
        this.failedPoints = failedPoints;
    }

    public List<ParameterPoint> getResults() {
        return results;
    }

    public void setResults(List<ParameterPoint> results) {
        this.results = results;
    }

    public String getError() {
        return error;
    }

    public void setError(String error) {
        this.error = error;
    }

    /**
     * Represents a single point in the parameter sweep.
     */
    public static class ParameterPoint {
        private Map<String, Double> parameters;
        private Map<String, Object> metrics;
        private String error;

        public Map<String, Double> getParameters() {
            return parameters;
        }

        public void setParameters(Map<String, Double> parameters) {
            this.parameters = parameters;
        }

        public Map<String, Object> getMetrics() {
            return metrics;
        }

        public void setMetrics(Map<String, Object> metrics) {
            this.metrics = metrics;
        }

        public String getError() {
            return error;
        }

        public void setError(String error) {
            this.error = error;
        }
    }
}
