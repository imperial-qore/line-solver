/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.rest.model;

import java.util.List;
import java.util.Map;

/**
 * Response DTO for bottleneck detection results.
 */
public class BottleneckResponse {

    private String status;
    private String solver;
    private double runtime;
    private double threshold;
    private List<Bottleneck> bottlenecks;
    private List<StationMetrics> allStations;
    private Map<String, Object> systemMetrics;
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

    public double getThreshold() {
        return threshold;
    }

    public void setThreshold(double threshold) {
        this.threshold = threshold;
    }

    public List<Bottleneck> getBottlenecks() {
        return bottlenecks;
    }

    public void setBottlenecks(List<Bottleneck> bottlenecks) {
        this.bottlenecks = bottlenecks;
    }

    public List<StationMetrics> getAllStations() {
        return allStations;
    }

    public void setAllStations(List<StationMetrics> allStations) {
        this.allStations = allStations;
    }

    public Map<String, Object> getSystemMetrics() {
        return systemMetrics;
    }

    public void setSystemMetrics(Map<String, Object> systemMetrics) {
        this.systemMetrics = systemMetrics;
    }

    public String getError() {
        return error;
    }

    public void setError(String error) {
        this.error = error;
    }

    /**
     * Represents a detected bottleneck station.
     */
    public static class Bottleneck {
        private String station;
        private double utilization;
        private double queueLength;
        private double responseTime;
        private String severity;
        private Map<String, Double> classBreakdown;

        public String getStation() {
            return station;
        }

        public void setStation(String station) {
            this.station = station;
        }

        public double getUtilization() {
            return utilization;
        }

        public void setUtilization(double utilization) {
            this.utilization = utilization;
        }

        public double getQueueLength() {
            return queueLength;
        }

        public void setQueueLength(double queueLength) {
            this.queueLength = queueLength;
        }

        public double getResponseTime() {
            return responseTime;
        }

        public void setResponseTime(double responseTime) {
            this.responseTime = responseTime;
        }

        public String getSeverity() {
            return severity;
        }

        public void setSeverity(String severity) {
            this.severity = severity;
        }

        public Map<String, Double> getClassBreakdown() {
            return classBreakdown;
        }

        public void setClassBreakdown(Map<String, Double> classBreakdown) {
            this.classBreakdown = classBreakdown;
        }
    }

    /**
     * Metrics for a single station.
     */
    public static class StationMetrics {
        private String station;
        private double utilization;
        private double queueLength;
        private double responseTime;
        private double throughput;
        private boolean isBottleneck;

        public String getStation() {
            return station;
        }

        public void setStation(String station) {
            this.station = station;
        }

        public double getUtilization() {
            return utilization;
        }

        public void setUtilization(double utilization) {
            this.utilization = utilization;
        }

        public double getQueueLength() {
            return queueLength;
        }

        public void setQueueLength(double queueLength) {
            this.queueLength = queueLength;
        }

        public double getResponseTime() {
            return responseTime;
        }

        public void setResponseTime(double responseTime) {
            this.responseTime = responseTime;
        }

        public double getThroughput() {
            return throughput;
        }

        public void setThroughput(double throughput) {
            this.throughput = throughput;
        }

        public boolean isBottleneck() {
            return isBottleneck;
        }

        public void setBottleneck(boolean bottleneck) {
            isBottleneck = bottleneck;
        }
    }
}
