/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.rest.model;

import com.google.gson.annotations.SerializedName;
import java.util.List;

/**
 * Defines a parameter to sweep in what-if analysis.
 *
 * Supports two field naming conventions:
 * - "type" maps to property (e.g., "service_rate" -> "serviceRate")
 * - "class" is an alias for "className"
 */
public class ParameterSweep {

    /**
     * Parameter type: service_rate, arrival_rate, think_time, servers, etc.
     * Used as alias for property.
     */
    private String type;

    /**
     * Name of the station (for station parameters).
     */
    private String station;

    /**
     * Name of the class (for class parameters).
     * The @SerializedName allows using "class" in JSON.
     */
    @SerializedName("class")
    private String className;

    /**
     * Property to vary: servers, serviceRate, arrivalRate, thinkTime, etc.
     */
    private String property;

    /**
     * List of values to sweep through.
     */
    private List<Double> values;

    /**
     * Alternative: start value for range.
     */
    private Double start;

    /**
     * Alternative: end value for range.
     */
    private Double end;

    /**
     * Alternative: step size for range.
     */
    private Double step;

    /**
     * Alternative: number of points for range.
     */
    private Integer points;

    public String getType() {
        return type;
    }

    public void setType(String type) {
        this.type = type;
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

    public String getProperty() {
        // If property is not set but type is, use type as property
        if ((property == null || property.isEmpty()) && type != null && !type.isEmpty()) {
            return mapTypeToProperty(type);
        }
        return property;
    }

    public void setProperty(String property) {
        this.property = property;
    }

    /**
     * Map legacy type names to property names.
     */
    private String mapTypeToProperty(String type) {
        if (type == null) return null;
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

    public List<Double> getValues() {
        return values;
    }

    public void setValues(List<Double> values) {
        this.values = values;
    }

    public Double getStart() {
        return start;
    }

    public void setStart(Double start) {
        this.start = start;
    }

    public Double getEnd() {
        return end;
    }

    public void setEnd(Double end) {
        this.end = end;
    }

    public Double getStep() {
        return step;
    }

    public void setStep(Double step) {
        this.step = step;
    }

    public Integer getPoints() {
        return points;
    }

    public void setPoints(Integer points) {
        this.points = points;
    }

    /**
     * Get the list of values to sweep, either from explicit values or generated from range.
     * @return List of values
     */
    public List<Double> getEffectiveValues() {
        if (values != null && !values.isEmpty()) {
            return values;
        }

        // Generate from range
        if (start != null && end != null) {
            java.util.List<Double> result = new java.util.ArrayList<Double>();

            if (step != null && step > 0) {
                for (double v = start; v <= end; v += step) {
                    result.add(v);
                }
            } else if (points != null && points > 1) {
                double stepSize = (end - start) / (points - 1);
                for (int i = 0; i < points; i++) {
                    result.add(start + i * stepSize);
                }
            } else {
                // Default to 10 points
                double stepSize = (end - start) / 9;
                for (int i = 0; i < 10; i++) {
                    result.add(start + i * stepSize);
                }
            }
            return result;
        }

        return new java.util.ArrayList<Double>();
    }

    /**
     * Validate the parameter sweep definition.
     * @return Error message if invalid, null if valid
     */
    public String validate() {
        // Use getProperty() which handles type -> property mapping
        String effectiveProperty = getProperty();
        if (effectiveProperty == null || effectiveProperty.isEmpty()) {
            return "Parameter property (or type) is required";
        }

        boolean hasValues = values != null && !values.isEmpty();
        boolean hasRange = start != null && end != null;

        if (!hasValues && !hasRange) {
            return "Parameter must have either 'values' array or 'start'/'end' range";
        }

        if (hasRange && start > end) {
            return "Parameter range 'start' must be <= 'end'";
        }

        return null;
    }

    /**
     * Validate for sensitivity analysis (values not required).
     * @return Error message if invalid, null if valid
     */
    public String validateForSensitivity() {
        String effectiveProperty = getProperty();
        if (effectiveProperty == null || effectiveProperty.isEmpty()) {
            return "Parameter property (or type) is required";
        }
        return null;
    }

    /**
     * Get a descriptive name for this parameter.
     * @return Parameter description
     */
    public String getDescription() {
        StringBuilder sb = new StringBuilder();
        if (station != null) {
            sb.append(station).append(".");
        }
        if (className != null) {
            sb.append(className).append(".");
        }
        // Use getProperty() to get the effective property
        String effectiveProperty = getProperty();
        if (effectiveProperty != null) {
            sb.append(effectiveProperty);
        }
        return sb.toString();
    }
}
