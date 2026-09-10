/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.rest.model;

import java.util.List;
import java.util.Map;

/**
 * Request DTO for importing OpenTelemetry/Jaeger traces.
 */
public class TraceImportRequest {

    /**
     * Trace format: opentelemetry, jaeger, zipkin.
     */
    private String format;

    /**
     * Trace data (JSON).
     */
    private Object traces;

    /**
     * Output model format: jsimg, lqnx.
     */
    private String outputFormat;

    /**
     * Options for trace processing.
     */
    private TraceImportOptions options;

    public String getFormat() {
        return format;
    }

    public void setFormat(String format) {
        this.format = format;
    }

    public Object getTraces() {
        return traces;
    }

    public void setTraces(Object traces) {
        this.traces = traces;
    }

    public String getOutputFormat() {
        return outputFormat;
    }

    public void setOutputFormat(String outputFormat) {
        this.outputFormat = outputFormat;
    }

    public TraceImportOptions getOptions() {
        return options;
    }

    public void setOptions(TraceImportOptions options) {
        this.options = options;
    }

    /**
     * Get effective output format (default jsimg).
     */
    public String getEffectiveOutputFormat() {
        return outputFormat != null ? outputFormat : "jsimg";
    }

    /**
     * Validate the request.
     */
    public String validate() {
        if (format == null || format.isEmpty()) {
            return "Trace format is required (opentelemetry, jaeger, zipkin)";
        }
        if (!format.equals("opentelemetry") && !format.equals("jaeger") && !format.equals("zipkin")) {
            return "Unsupported trace format: " + format + ". Supported: opentelemetry, jaeger, zipkin";
        }
        if (traces == null) {
            return "Traces data is required";
        }
        if (outputFormat != null && !outputFormat.equals("jsimg") && !outputFormat.equals("lqnx")) {
            return "Unsupported output format: " + outputFormat + ". Supported: jsimg, lqnx";
        }
        return null;
    }

    /**
     * Options for trace processing.
     */
    public static class TraceImportOptions {
        /**
         * Minimum number of spans to include a service.
         */
        private Integer minSpans;

        /**
         * Whether to aggregate multiple traces.
         */
        private Boolean aggregate;

        /**
         * Service name filter (regex).
         */
        private String serviceFilter;

        /**
         * Operation name filter (regex).
         */
        private String operationFilter;

        /**
         * Time window for aggregation (seconds).
         */
        private Integer timeWindow;

        /**
         * Percentile for service time estimation (default 50 = median).
         */
        private Integer percentile;

        public Integer getMinSpans() {
            return minSpans;
        }

        public void setMinSpans(Integer minSpans) {
            this.minSpans = minSpans;
        }

        public Boolean getAggregate() {
            return aggregate;
        }

        public void setAggregate(Boolean aggregate) {
            this.aggregate = aggregate;
        }

        public String getServiceFilter() {
            return serviceFilter;
        }

        public void setServiceFilter(String serviceFilter) {
            this.serviceFilter = serviceFilter;
        }

        public String getOperationFilter() {
            return operationFilter;
        }

        public void setOperationFilter(String operationFilter) {
            this.operationFilter = operationFilter;
        }

        public Integer getTimeWindow() {
            return timeWindow;
        }

        public void setTimeWindow(Integer timeWindow) {
            this.timeWindow = timeWindow;
        }

        public Integer getPercentile() {
            return percentile;
        }

        public void setPercentile(Integer percentile) {
            this.percentile = percentile;
        }

        public int getEffectiveMinSpans() {
            return minSpans != null ? minSpans : 1;
        }

        public int getEffectivePercentile() {
            return percentile != null ? percentile : 50;
        }

        public boolean shouldAggregate() {
            return aggregate != null && aggregate;
        }
    }
}
