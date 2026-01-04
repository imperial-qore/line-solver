/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.rest.model;

import java.util.Map;

/**
 * Response payload for the /models/solve endpoint.
 */
public class SolveResponse {

    /**
     * Status of the solve operation.
     * Values: completed, failed
     */
    private String status;

    /**
     * The solver that was used.
     */
    private String solver;

    /**
     * Runtime of the solve operation in seconds.
     */
    private double runtime;

    /**
     * The results of the solve operation.
     */
    private ResultTables results;

    /**
     * Error message if status is "failed".
     */
    private String error;

    /**
     * Default constructor.
     */
    public SolveResponse() {
    }

    /**
     * Create a successful response.
     */
    public static SolveResponse success(String solver, double runtime, ResultTables results) {
        SolveResponse response = new SolveResponse();
        response.status = "completed";
        response.solver = solver;
        response.runtime = runtime;
        response.results = results;
        return response;
    }

    /**
     * Create a failed response.
     */
    public static SolveResponse failure(String solver, String error) {
        SolveResponse response = new SolveResponse();
        response.status = "failed";
        response.solver = solver;
        response.error = error;
        return response;
    }

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

    public ResultTables getResults() {
        return results;
    }

    public void setResults(ResultTables results) {
        this.results = results;
    }

    public String getError() {
        return error;
    }

    public void setError(String error) {
        this.error = error;
    }

    /**
     * Nested class representing the result tables.
     */
    public static class ResultTables {
        private AvgTableResult avgTable;
        private AvgSysTableResult sysTable;

        public AvgTableResult getAvgTable() {
            return avgTable;
        }

        public void setAvgTable(AvgTableResult avgTable) {
            this.avgTable = avgTable;
        }

        public AvgSysTableResult getSysTable() {
            return sysTable;
        }

        public void setSysTable(AvgSysTableResult sysTable) {
            this.sysTable = sysTable;
        }
    }

    /**
     * Average table result structure.
     */
    public static class AvgTableResult {
        private String[] stations;
        private String[] classes;
        private Map<String, double[][]> metrics;

        public String[] getStations() {
            return stations;
        }

        public void setStations(String[] stations) {
            this.stations = stations;
        }

        public String[] getClasses() {
            return classes;
        }

        public void setClasses(String[] classes) {
            this.classes = classes;
        }

        public Map<String, double[][]> getMetrics() {
            return metrics;
        }

        public void setMetrics(Map<String, double[][]> metrics) {
            this.metrics = metrics;
        }
    }

    /**
     * Average system table result structure.
     */
    public static class AvgSysTableResult {
        private String[] chains;
        private Map<String, double[]> metrics;

        public String[] getChains() {
            return chains;
        }

        public void setChains(String[] chains) {
            this.chains = chains;
        }

        public Map<String, double[]> getMetrics() {
            return metrics;
        }

        public void setMetrics(Map<String, double[]> metrics) {
            this.metrics = metrics;
        }
    }
}
