/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.rest.model;

import com.google.gson.annotations.SerializedName;
import java.util.ArrayList;
import java.util.List;

/**
 * Request DTO for what-if analysis (parameter sweep).
 *
 * Supports two formats:
 * 1. Simple format: single "parameter" object with "values" array at top level
 * 2. Advanced format: "parameters" array of ParameterSweep objects
 */
public class WhatIfRequest {

    private ModelInput model;
    private String solver;
    private String analysis;
    private List<ParameterSweep> parameters;
    private List<String> metrics;
    private SolverOptionsDTO options;
    private Boolean parallel;

    // Simple format fields (for backwards compatibility with demo)
    private SimpleParameter parameter;
    private List<Double> values;

    /**
     * Simple parameter definition for backwards compatibility.
     */
    public static class SimpleParameter {
        private String type;
        private String station;
        @SerializedName("class")
        private String className;

        public String getType() { return type; }
        public void setType(String type) { this.type = type; }

        public String getStation() { return station; }
        public void setStation(String station) { this.station = station; }

        public String getClassName() { return className; }
        public void setClassName(String className) { this.className = className; }
    }

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
        // Convert simple format to parameters list if needed
        if ((parameters == null || parameters.isEmpty()) && parameter != null && values != null) {
            parameters = new ArrayList<ParameterSweep>();
            ParameterSweep sweep = new ParameterSweep();
            sweep.setStation(parameter.getStation());
            sweep.setClassName(parameter.getClassName());
            // Map type to property
            String prop = mapTypeToProperty(parameter.getType());
            sweep.setProperty(prop);
            sweep.setValues(values);
            parameters.add(sweep);
        }
        return parameters;
    }

    public void setParameters(List<ParameterSweep> parameters) {
        this.parameters = parameters;
    }

    public SimpleParameter getParameter() {
        return parameter;
    }

    public void setParameter(SimpleParameter parameter) {
        this.parameter = parameter;
    }

    public List<Double> getValues() {
        return values;
    }

    public void setValues(List<Double> values) {
        this.values = values;
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

    public Boolean getParallel() {
        return parallel;
    }

    public void setParallel(Boolean parallel) {
        this.parallel = parallel;
    }

    /**
     * Map legacy type names to property names.
     */
    private String mapTypeToProperty(String type) {
        if (type == null) return "serviceRate";
        switch (type.toLowerCase()) {
            case "service_rate":
            case "servicerate":
                return "serviceRate";
            case "arrival_rate":
            case "arrivalrate":
                return "arrivalRate";
            case "think_time":
            case "thinktime":
                return "thinkTime";
            case "servers":
            case "num_servers":
                return "servers";
            case "population":
            case "customers":
                return "population";
            default:
                return type;
        }
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

        // Check for simple format
        boolean hasSimpleFormat = parameter != null && values != null && !values.isEmpty();
        boolean hasAdvancedFormat = parameters != null && !parameters.isEmpty();

        if (!hasSimpleFormat && !hasAdvancedFormat) {
            return "Either 'parameter' with 'values' or 'parameters' array is required";
        }

        // Trigger conversion if using simple format
        List<ParameterSweep> params = getParameters();
        if (params == null || params.isEmpty()) {
            return "At least one parameter sweep is required";
        }

        for (ParameterSweep param : params) {
            String paramError = param.validate();
            if (paramError != null) {
                return paramError;
            }
        }
        return null;
    }
}
