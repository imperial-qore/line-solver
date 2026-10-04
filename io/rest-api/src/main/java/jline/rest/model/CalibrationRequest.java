/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.rest.model;

import java.util.List;
import java.util.Map;

/**
 * Request DTO for model calibration from observed metrics.
 */
public class CalibrationRequest {

    /**
     * Base model to calibrate (optional - can derive from topology).
     */
    private ModelInput model;

    /**
     * Observed metrics for each station/class.
     */
    private List<ObservedMetrics> observations;

    /**
     * Parameters to calibrate.
     */
    private List<String> parameters;

    /**
     * Calibration method: leastSquares, maxLikelihood.
     */
    private String method;

    /**
     * Solver to use for model evaluation.
     */
    private String solver;

    /**
     * Calibration options.
     */
    private CalibrationOptions options;

    public ModelInput getModel() {
        return model;
    }

    public void setModel(ModelInput model) {
        this.model = model;
    }

    public List<ObservedMetrics> getObservations() {
        return observations;
    }

    public void setObservations(List<ObservedMetrics> observations) {
        this.observations = observations;
    }

    public List<String> getParameters() {
        return parameters;
    }

    public void setParameters(List<String> parameters) {
        this.parameters = parameters;
    }

    public String getMethod() {
        return method;
    }

    public void setMethod(String method) {
        this.method = method;
    }

    public String getSolver() {
        return solver;
    }

    public void setSolver(String solver) {
        this.solver = solver;
    }

    public CalibrationOptions getOptions() {
        return options;
    }

    public void setOptions(CalibrationOptions options) {
        this.options = options;
    }

    /**
     * Validate the request.
     */
    public String validate() {
        if (observations == null || observations.isEmpty()) {
            return "Observations are required";
        }
        for (ObservedMetrics obs : observations) {
            String error = obs.validate();
            if (error != null) {
                return error;
            }
        }
        if (solver == null || solver.isEmpty()) {
            return "Solver is required";
        }
        return null;
    }

    /**
     * Observed metrics for a station/class.
     */
    public static class ObservedMetrics {
        private String station;
        private String className;
        private Double throughput;
        private Double utilization;
        private Double responseTime;
        private Double queueLength;
        private Double arrivalRate;

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

        public Double getThroughput() {
            return throughput;
        }

        public void setThroughput(Double throughput) {
            this.throughput = throughput;
        }

        public Double getUtilization() {
            return utilization;
        }

        public void setUtilization(Double utilization) {
            this.utilization = utilization;
        }

        public Double getResponseTime() {
            return responseTime;
        }

        public void setResponseTime(Double responseTime) {
            this.responseTime = responseTime;
        }

        public Double getQueueLength() {
            return queueLength;
        }

        public void setQueueLength(Double queueLength) {
            this.queueLength = queueLength;
        }

        public Double getArrivalRate() {
            return arrivalRate;
        }

        public void setArrivalRate(Double arrivalRate) {
            this.arrivalRate = arrivalRate;
        }

        public String validate() {
            if (station == null || station.isEmpty()) {
                return "Station name is required for observations";
            }
            if (throughput == null && utilization == null && responseTime == null &&
                queueLength == null && arrivalRate == null) {
                return "At least one metric is required for observation";
            }
            return null;
        }
    }

    /**
     * Calibration options.
     */
    public static class CalibrationOptions {
        private Integer maxIterations;
        private Double tolerance;
        private Double learningRate;
        private Boolean verbose;

        public Integer getMaxIterations() {
            return maxIterations;
        }

        public void setMaxIterations(Integer maxIterations) {
            this.maxIterations = maxIterations;
        }

        public Double getTolerance() {
            return tolerance;
        }

        public void setTolerance(Double tolerance) {
            this.tolerance = tolerance;
        }

        public Double getLearningRate() {
            return learningRate;
        }

        public void setLearningRate(Double learningRate) {
            this.learningRate = learningRate;
        }

        public Boolean getVerbose() {
            return verbose;
        }

        public void setVerbose(Boolean verbose) {
            this.verbose = verbose;
        }

        public int getEffectiveMaxIterations() {
            return maxIterations != null ? maxIterations : 100;
        }

        public double getEffectiveTolerance() {
            return tolerance != null ? tolerance : 1e-6;
        }
    }
}
