/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.rest.model;

import java.util.List;
import java.util.Map;

/**
 * Response DTO for trace import results.
 */
public class TraceImportResponse {

    private String status;
    private double runtime;
    private int tracesProcessed;
    private int spansProcessed;
    private int servicesDetected;
    private int operationsDetected;
    private String outputFormat;
    private String modelContent;
    private ModelSummary modelSummary;
    private List<ServiceInfo> services;
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

    public int getTracesProcessed() {
        return tracesProcessed;
    }

    public void setTracesProcessed(int tracesProcessed) {
        this.tracesProcessed = tracesProcessed;
    }

    public int getSpansProcessed() {
        return spansProcessed;
    }

    public void setSpansProcessed(int spansProcessed) {
        this.spansProcessed = spansProcessed;
    }

    public int getServicesDetected() {
        return servicesDetected;
    }

    public void setServicesDetected(int servicesDetected) {
        this.servicesDetected = servicesDetected;
    }

    public int getOperationsDetected() {
        return operationsDetected;
    }

    public void setOperationsDetected(int operationsDetected) {
        this.operationsDetected = operationsDetected;
    }

    public String getOutputFormat() {
        return outputFormat;
    }

    public void setOutputFormat(String outputFormat) {
        this.outputFormat = outputFormat;
    }

    public String getModelContent() {
        return modelContent;
    }

    public void setModelContent(String modelContent) {
        this.modelContent = modelContent;
    }

    public ModelSummary getModelSummary() {
        return modelSummary;
    }

    public void setModelSummary(ModelSummary modelSummary) {
        this.modelSummary = modelSummary;
    }

    public List<ServiceInfo> getServices() {
        return services;
    }

    public void setServices(List<ServiceInfo> services) {
        this.services = services;
    }

    public String getError() {
        return error;
    }

    public void setError(String error) {
        this.error = error;
    }

    /**
     * Summary of the generated model.
     */
    public static class ModelSummary {
        private String type;
        private int stations;
        private int classes;
        private int chains;

        public String getType() {
            return type;
        }

        public void setType(String type) {
            this.type = type;
        }

        public int getStations() {
            return stations;
        }

        public void setStations(int stations) {
            this.stations = stations;
        }

        public int getClasses() {
            return classes;
        }

        public void setClasses(int classes) {
            this.classes = classes;
        }

        public int getChains() {
            return chains;
        }

        public void setChains(int chains) {
            this.chains = chains;
        }
    }

    /**
     * Information about a detected service.
     */
    public static class ServiceInfo {
        private String name;
        private int spanCount;
        private List<String> operations;
        private double avgDurationMs;
        private double p50DurationMs;
        private double p95DurationMs;
        private double p99DurationMs;
        private double estimatedServiceRate;
        private Map<String, Integer> callerCounts;

        public String getName() {
            return name;
        }

        public void setName(String name) {
            this.name = name;
        }

        public int getSpanCount() {
            return spanCount;
        }

        public void setSpanCount(int spanCount) {
            this.spanCount = spanCount;
        }

        public List<String> getOperations() {
            return operations;
        }

        public void setOperations(List<String> operations) {
            this.operations = operations;
        }

        public double getAvgDurationMs() {
            return avgDurationMs;
        }

        public void setAvgDurationMs(double avgDurationMs) {
            this.avgDurationMs = avgDurationMs;
        }

        public double getP50DurationMs() {
            return p50DurationMs;
        }

        public void setP50DurationMs(double p50DurationMs) {
            this.p50DurationMs = p50DurationMs;
        }

        public double getP95DurationMs() {
            return p95DurationMs;
        }

        public void setP95DurationMs(double p95DurationMs) {
            this.p95DurationMs = p95DurationMs;
        }

        public double getP99DurationMs() {
            return p99DurationMs;
        }

        public void setP99DurationMs(double p99DurationMs) {
            this.p99DurationMs = p99DurationMs;
        }

        public double getEstimatedServiceRate() {
            return estimatedServiceRate;
        }

        public void setEstimatedServiceRate(double estimatedServiceRate) {
            this.estimatedServiceRate = estimatedServiceRate;
        }

        public Map<String, Integer> getCallerCounts() {
            return callerCounts;
        }

        public void setCallerCounts(Map<String, Integer> callerCounts) {
            this.callerCounts = callerCounts;
        }
    }
}
