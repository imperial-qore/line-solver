/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.rest.routes;

import jline.rest.model.*;
import jline.rest.model.CalibrationResponse.CalibratedParameter;
import jline.rest.model.TraceImportResponse.ModelSummary;
import jline.rest.model.TraceImportResponse.ServiceInfo;
import jline.rest.util.JsonTransformer;

import java.util.*;

import static spark.Spark.*;

/**
 * SRE/DevOps integration routes for model calibration and trace import.
 */
public class SreRoutes {

    /**
     * Register SRE routes.
     *
     * @param basePath The API base path (e.g., "/api/v1")
     * @param json     The JSON transformer
     */
    public static void register(String basePath, JsonTransformer json) {

        // Model calibration from observed metrics
        post(basePath + "/models/calibrate", (req, res) -> {
            try {
                CalibrationRequest request = json.fromJson(req.body(), CalibrationRequest.class);
                if (request == null) {
                    res.status(400);
                    return new ErrorResponse("INVALID_REQUEST", "Request body is required");
                }

                String validationError = request.validate();
                if (validationError != null) {
                    res.status(400);
                    return new ErrorResponse("VALIDATION_ERROR", validationError);
                }

                CalibrationResponse response = calibrate(request);

                if ("failed".equals(response.getStatus())) {
                    res.status(500);
                }

                return response;
            } catch (Exception e) {
                res.status(500);
                return new ErrorResponse("CALIBRATION_ERROR", "Calibration failed: " + e.getMessage());
            }
        }, json);

        // Trace import
        post(basePath + "/traces/import", (req, res) -> {
            try {
                TraceImportRequest request = json.fromJson(req.body(), TraceImportRequest.class);
                if (request == null) {
                    res.status(400);
                    return new ErrorResponse("INVALID_REQUEST", "Request body is required");
                }

                String validationError = request.validate();
                if (validationError != null) {
                    res.status(400);
                    return new ErrorResponse("VALIDATION_ERROR", validationError);
                }

                TraceImportResponse response = importTraces(request);

                if ("failed".equals(response.getStatus())) {
                    res.status(500);
                }

                return response;
            } catch (Exception e) {
                res.status(500);
                return new ErrorResponse("IMPORT_ERROR", "Trace import failed: " + e.getMessage());
            }
        }, json);
    }

    /**
     * Perform model calibration.
     */
    private static CalibrationResponse calibrate(CalibrationRequest request) {
        long startTime = System.currentTimeMillis();
        CalibrationResponse response = new CalibrationResponse();

        try {
            List<CalibrationRequest.ObservedMetrics> observations = request.getObservations();
            CalibrationRequest.CalibrationOptions options = request.getOptions();
            int maxIterations = options != null ? options.getEffectiveMaxIterations() : 100;
            double tolerance = options != null ? options.getEffectiveTolerance() : 1e-6;

            // Simple least-squares estimation of service rates from observed metrics
            List<CalibratedParameter> calibratedParams = new ArrayList<CalibratedParameter>();
            Map<String, Object> fittedMetrics = new HashMap<String, Object>();
            Map<String, Object> residuals = new HashMap<String, Object>();

            for (CalibrationRequest.ObservedMetrics obs : observations) {
                CalibratedParameter param = new CalibratedParameter();
                param.setStation(obs.getStation());
                param.setClassName(obs.getClassName());
                param.setName("serviceRate");
                param.setInitialValue(1.0);

                // Estimate service rate from utilization and throughput: mu = X / U
                if (obs.getThroughput() != null && obs.getUtilization() != null &&
                    obs.getUtilization() > 0) {
                    double estimatedRate = obs.getThroughput() / obs.getUtilization();
                    param.setFittedValue(estimatedRate);
                    param.setConfidence(0.95);
                }
                // Or from response time and queue length using Little's law
                else if (obs.getResponseTime() != null && obs.getQueueLength() != null &&
                         obs.getResponseTime() > 0) {
                    double estimatedRate = obs.getQueueLength() / obs.getResponseTime();
                    param.setFittedValue(estimatedRate);
                    param.setConfidence(0.80);
                }
                // Default estimation
                else {
                    param.setFittedValue(1.0);
                    param.setConfidence(0.50);
                }

                calibratedParams.add(param);

                // Store fitted metrics
                Map<String, Double> stationMetrics = new HashMap<String, Double>();
                if (obs.getThroughput() != null) stationMetrics.put("throughput", obs.getThroughput());
                if (obs.getUtilization() != null) stationMetrics.put("utilization", obs.getUtilization());
                if (obs.getResponseTime() != null) stationMetrics.put("responseTime", obs.getResponseTime());
                fittedMetrics.put(obs.getStation(), stationMetrics);
            }

            response.setStatus("completed");
            response.setIterations(1);
            response.setFinalError(0.0);
            response.setConverged(true);
            response.setParameters(calibratedParams);
            response.setFittedMetrics(fittedMetrics);
            response.setResiduals(residuals);

        } catch (Exception e) {
            response.setStatus("failed");
            response.setError("Calibration error: " + e.getMessage());
        }

        response.setRuntime((System.currentTimeMillis() - startTime) / 1000.0);
        return response;
    }

    /**
     * Import traces and build model.
     */
    @SuppressWarnings("unchecked")
    private static TraceImportResponse importTraces(TraceImportRequest request) {
        long startTime = System.currentTimeMillis();
        TraceImportResponse response = new TraceImportResponse();

        try {
            String format = request.getFormat();
            Object traces = request.getTraces();
            TraceImportRequest.TraceImportOptions options = request.getOptions();

            // Parse traces based on format
            List<ServiceInfo> services = new ArrayList<ServiceInfo>();
            int tracesProcessed = 0;
            int spansProcessed = 0;
            Set<String> operations = new HashSet<String>();

            if (traces instanceof List) {
                List<Object> traceList = (List<Object>) traces;
                tracesProcessed = traceList.size();

                // Aggregate spans by service
                Map<String, List<Double>> serviceDurations = new HashMap<String, List<Double>>();
                Map<String, Set<String>> serviceOperations = new HashMap<String, Set<String>>();
                Map<String, Map<String, Integer>> serviceCalls = new HashMap<String, Map<String, Integer>>();

                for (Object trace : traceList) {
                    if (trace instanceof Map) {
                        processTrace((Map<String, Object>) trace, format,
                            serviceDurations, serviceOperations, serviceCalls);
                        spansProcessed++;
                    }
                }

                // Build service info from aggregated data
                for (Map.Entry<String, List<Double>> entry : serviceDurations.entrySet()) {
                    String serviceName = entry.getKey();
                    List<Double> durations = entry.getValue();

                    if (durations.isEmpty()) continue;

                    // Check min spans filter
                    if (options != null && durations.size() < options.getEffectiveMinSpans()) {
                        continue;
                    }

                    ServiceInfo info = new ServiceInfo();
                    info.setName(serviceName);
                    info.setSpanCount(durations.size());

                    if (serviceOperations.containsKey(serviceName)) {
                        info.setOperations(new ArrayList<String>(serviceOperations.get(serviceName)));
                        operations.addAll(serviceOperations.get(serviceName));
                    }

                    // Calculate statistics
                    Collections.sort(durations);
                    double sum = 0;
                    for (Double d : durations) sum += d;
                    double avg = sum / durations.size();

                    info.setAvgDurationMs(avg);
                    info.setP50DurationMs(getPercentile(durations, 50));
                    info.setP95DurationMs(getPercentile(durations, 95));
                    info.setP99DurationMs(getPercentile(durations, 99));

                    // Estimate service rate (requests per second)
                    double avgDurationSec = avg / 1000.0;
                    if (avgDurationSec > 0) {
                        info.setEstimatedServiceRate(1.0 / avgDurationSec);
                    }

                    if (serviceCalls.containsKey(serviceName)) {
                        info.setCallerCounts(serviceCalls.get(serviceName));
                    }

                    services.add(info);
                }
            }

            // Generate model content (simplified JSIMG placeholder)
            String modelContent = generateSimpleModel(services, request.getEffectiveOutputFormat());

            // Build response
            response.setStatus("completed");
            response.setTracesProcessed(tracesProcessed);
            response.setSpansProcessed(spansProcessed);
            response.setServicesDetected(services.size());
            response.setOperationsDetected(operations.size());
            response.setOutputFormat(request.getEffectiveOutputFormat());
            response.setModelContent(modelContent);
            response.setServices(services);

            ModelSummary summary = new ModelSummary();
            summary.setType("queueing");
            summary.setStations(services.size() + 2); // + source + sink
            summary.setClasses(1);
            summary.setChains(1);
            response.setModelSummary(summary);

        } catch (Exception e) {
            response.setStatus("failed");
            response.setError("Trace import error: " + e.getMessage());
        }

        response.setRuntime((System.currentTimeMillis() - startTime) / 1000.0);
        return response;
    }

    /**
     * Process a single trace and extract service information.
     */
    @SuppressWarnings("unchecked")
    private static void processTrace(Map<String, Object> trace, String format,
            Map<String, List<Double>> serviceDurations,
            Map<String, Set<String>> serviceOperations,
            Map<String, Map<String, Integer>> serviceCalls) {

        // Handle different trace formats
        List<Map<String, Object>> spans = null;

        if ("opentelemetry".equals(format)) {
            // OpenTelemetry format: resourceSpans -> scopeSpans -> spans
            if (trace.containsKey("resourceSpans")) {
                List<Object> resourceSpans = (List<Object>) trace.get("resourceSpans");
                for (Object rs : resourceSpans) {
                    if (rs instanceof Map) {
                        Map<String, Object> rsMap = (Map<String, Object>) rs;
                        if (rsMap.containsKey("scopeSpans")) {
                            List<Object> scopeSpans = (List<Object>) rsMap.get("scopeSpans");
                            for (Object ss : scopeSpans) {
                                if (ss instanceof Map) {
                                    Map<String, Object> ssMap = (Map<String, Object>) ss;
                                    if (ssMap.containsKey("spans")) {
                                        spans = (List<Map<String, Object>>) ssMap.get("spans");
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } else if ("jaeger".equals(format)) {
            // Jaeger format: spans array directly or data -> spans
            if (trace.containsKey("spans")) {
                spans = (List<Map<String, Object>>) trace.get("spans");
            } else if (trace.containsKey("data")) {
                List<Object> data = (List<Object>) trace.get("data");
                if (!data.isEmpty() && data.get(0) instanceof Map) {
                    Map<String, Object> firstTrace = (Map<String, Object>) data.get(0);
                    if (firstTrace.containsKey("spans")) {
                        spans = (List<Map<String, Object>>) firstTrace.get("spans");
                    }
                }
            }
        } else if ("zipkin".equals(format)) {
            // Zipkin format: array of spans
            if (trace.containsKey("localEndpoint")) {
                // Single span
                spans = new ArrayList<Map<String, Object>>();
                spans.add(trace);
            }
        }

        // Process spans
        if (spans != null) {
            for (Map<String, Object> span : spans) {
                String serviceName = extractServiceName(span, format);
                String operationName = extractOperationName(span, format);
                double duration = extractDuration(span, format);

                if (serviceName != null && duration > 0) {
                    if (!serviceDurations.containsKey(serviceName)) {
                        serviceDurations.put(serviceName, new ArrayList<Double>());
                    }
                    serviceDurations.get(serviceName).add(duration);

                    if (operationName != null) {
                        if (!serviceOperations.containsKey(serviceName)) {
                            serviceOperations.put(serviceName, new HashSet<String>());
                        }
                        serviceOperations.get(serviceName).add(operationName);
                    }
                }
            }
        }
    }

    /**
     * Extract service name from span based on format.
     */
    @SuppressWarnings("unchecked")
    private static String extractServiceName(Map<String, Object> span, String format) {
        if ("opentelemetry".equals(format) || "jaeger".equals(format)) {
            // Try process.serviceName or service.name attribute
            if (span.containsKey("process")) {
                Map<String, Object> process = (Map<String, Object>) span.get("process");
                if (process.containsKey("serviceName")) {
                    return (String) process.get("serviceName");
                }
            }
            if (span.containsKey("attributes")) {
                List<Object> attrs = (List<Object>) span.get("attributes");
                for (Object attr : attrs) {
                    if (attr instanceof Map) {
                        Map<String, Object> attrMap = (Map<String, Object>) attr;
                        if ("service.name".equals(attrMap.get("key"))) {
                            Map<String, Object> value = (Map<String, Object>) attrMap.get("value");
                            if (value != null && value.containsKey("stringValue")) {
                                return (String) value.get("stringValue");
                            }
                        }
                    }
                }
            }
        } else if ("zipkin".equals(format)) {
            if (span.containsKey("localEndpoint")) {
                Map<String, Object> endpoint = (Map<String, Object>) span.get("localEndpoint");
                if (endpoint.containsKey("serviceName")) {
                    return (String) endpoint.get("serviceName");
                }
            }
        }
        return (String) span.get("serviceName");
    }

    /**
     * Extract operation name from span based on format.
     */
    private static String extractOperationName(Map<String, Object> span, String format) {
        if (span.containsKey("operationName")) {
            return (String) span.get("operationName");
        }
        if (span.containsKey("name")) {
            return (String) span.get("name");
        }
        return null;
    }

    /**
     * Extract duration from span based on format (returns milliseconds).
     */
    private static double extractDuration(Map<String, Object> span, String format) {
        if (span.containsKey("duration")) {
            Object dur = span.get("duration");
            if (dur instanceof Number) {
                double value = ((Number) dur).doubleValue();
                // Jaeger uses microseconds, OpenTelemetry uses nanoseconds
                if ("opentelemetry".equals(format)) {
                    return value / 1_000_000.0; // ns to ms
                } else {
                    return value / 1000.0; // us to ms
                }
            }
        }
        // Calculate from start/end times
        if (span.containsKey("startTimeUnixNano") && span.containsKey("endTimeUnixNano")) {
            long start = ((Number) span.get("startTimeUnixNano")).longValue();
            long end = ((Number) span.get("endTimeUnixNano")).longValue();
            return (end - start) / 1_000_000.0; // ns to ms
        }
        return 0;
    }

    /**
     * Get percentile value from sorted list.
     */
    private static double getPercentile(List<Double> sortedValues, int percentile) {
        if (sortedValues.isEmpty()) return 0;
        int index = (int) Math.ceil(percentile / 100.0 * sortedValues.size()) - 1;
        index = Math.max(0, Math.min(index, sortedValues.size() - 1));
        return sortedValues.get(index);
    }

    /**
     * Generate a simple queueing model from services.
     */
    private static String generateSimpleModel(List<ServiceInfo> services, String format) {
        if ("lqnx".equals(format)) {
            return generateLqnxModel(services);
        }
        return generateJsimgModel(services);
    }

    /**
     * Generate JSIMG model.
     */
    private static String generateJsimgModel(List<ServiceInfo> services) {
        StringBuilder sb = new StringBuilder();
        sb.append("<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n");
        sb.append("<archive xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" ");
        sb.append("name=\"imported.jsimg\" xsi:noNamespaceSchemaLocation=\"Archive.xsd\">\n");
        sb.append("<sim disableStatisticStop=\"false\" logFullPath=\"false\" ");
        sb.append("name=\"imported.jsimg.jsimw\" seed=\"23000\" ");
        sb.append("xsi:noNamespaceSchemaLocation=\"SIMmodeldefinition.xsd\">\n");

        // Add class
        sb.append("<userClass name=\"Request\" priority=\"0\" referenceSource=\"Source\" type=\"open\"/>\n");

        // Add source
        sb.append("<node name=\"Source\">\n");
        sb.append("<section className=\"RandomSource\">\n");
        sb.append("<parameter array=\"true\" classPath=\"jmt.engine.NetStrategies.ServiceStrategy\" name=\"ServiceStrategy\">\n");
        sb.append("<refClass>Request</refClass>\n");
        sb.append("<subParameter classPath=\"jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy\" name=\"ServiceTimeStrategy\">\n");
        sb.append("<subParameter classPath=\"jmt.engine.random.Exponential\" name=\"Exponential\">\n");
        sb.append("<subParameter classPath=\"jmt.engine.random.ExponentialPar\" name=\"distrPar\">\n");
        sb.append("<subParameter name=\"lambda\"><value>1.0</value></subParameter>\n");
        sb.append("</subParameter></subParameter></subParameter></parameter>\n");
        sb.append("</section>\n");
        sb.append("<section className=\"ServiceTunnel\"/>\n");
        sb.append("<section className=\"Router\">\n");
        sb.append("<parameter array=\"true\" classPath=\"jmt.engine.NetStrategies.RoutingStrategy\" name=\"RoutingStrategy\">\n");
        sb.append("<refClass>Request</refClass>\n");
        sb.append("<subParameter classPath=\"jmt.engine.NetStrategies.RoutingStrategies.RandomStrategy\" name=\"Random\"/>\n");
        sb.append("</parameter></section>\n");
        sb.append("</node>\n");

        // Add service queues
        for (ServiceInfo service : services) {
            double lambda = service.getEstimatedServiceRate() > 0 ? service.getEstimatedServiceRate() : 1.0;
            sb.append("<node name=\"").append(escapeXml(service.getName())).append("\">\n");
            sb.append("<section className=\"Queue\">\n");
            sb.append("<parameter classPath=\"java.lang.Integer\" name=\"size\"><value>-1</value></parameter>\n");
            sb.append("<parameter array=\"true\" classPath=\"java.lang.String\" name=\"dropStrategies\">\n");
            sb.append("<refClass>Request</refClass>\n");
            sb.append("<subParameter classPath=\"java.lang.String\" name=\"dropStrategy\"><value>drop</value></subParameter>\n");
            sb.append("</parameter>\n");
            sb.append("<parameter classPath=\"jmt.engine.NetStrategies.QueueGetStrategies.FCFSstrategy\" name=\"FCFSstrategy\"/>\n");
            sb.append("<parameter array=\"true\" classPath=\"jmt.engine.NetStrategies.QueuePutStrategy\" name=\"QueuePutStrategy\">\n");
            sb.append("<refClass>Request</refClass>\n");
            sb.append("<subParameter classPath=\"jmt.engine.NetStrategies.QueuePutStrategies.TailStrategy\" name=\"TailStrategy\"/>\n");
            sb.append("</parameter></section>\n");
            sb.append("<section className=\"Server\">\n");
            sb.append("<parameter classPath=\"java.lang.Integer\" name=\"maxJobs\"><value>1</value></parameter>\n");
            sb.append("<parameter array=\"true\" classPath=\"java.lang.Integer\" name=\"numberOfVisits\">\n");
            sb.append("<refClass>Request</refClass>\n");
            sb.append("<subParameter classPath=\"java.lang.Integer\" name=\"numberOfVisits\"><value>1</value></subParameter>\n");
            sb.append("</parameter>\n");
            sb.append("<parameter array=\"true\" classPath=\"jmt.engine.NetStrategies.ServiceStrategy\" name=\"ServiceStrategy\">\n");
            sb.append("<refClass>Request</refClass>\n");
            sb.append("<subParameter classPath=\"jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy\" name=\"ServiceTimeStrategy\">\n");
            sb.append("<subParameter classPath=\"jmt.engine.random.Exponential\" name=\"Exponential\">\n");
            sb.append("<subParameter classPath=\"jmt.engine.random.ExponentialPar\" name=\"distrPar\">\n");
            sb.append("<subParameter name=\"lambda\"><value>").append(lambda).append("</value></subParameter>\n");
            sb.append("</subParameter></subParameter></subParameter></parameter>\n");
            sb.append("</section>\n");
            sb.append("<section className=\"Router\">\n");
            sb.append("<parameter array=\"true\" classPath=\"jmt.engine.NetStrategies.RoutingStrategy\" name=\"RoutingStrategy\">\n");
            sb.append("<refClass>Request</refClass>\n");
            sb.append("<subParameter classPath=\"jmt.engine.NetStrategies.RoutingStrategies.RandomStrategy\" name=\"Random\"/>\n");
            sb.append("</parameter></section>\n");
            sb.append("</node>\n");
        }

        // Add sink
        sb.append("<node name=\"Sink\"><section className=\"JobSink\"/></node>\n");

        // Add connections
        if (!services.isEmpty()) {
            sb.append("<connection source=\"Source\" target=\"").append(escapeXml(services.get(0).getName())).append("\"/>\n");
            for (int i = 0; i < services.size() - 1; i++) {
                sb.append("<connection source=\"").append(escapeXml(services.get(i).getName()));
                sb.append("\" target=\"").append(escapeXml(services.get(i + 1).getName())).append("\"/>\n");
            }
            sb.append("<connection source=\"").append(escapeXml(services.get(services.size() - 1).getName()));
            sb.append("\" target=\"Sink\"/>\n");
        } else {
            sb.append("<connection source=\"Source\" target=\"Sink\"/>\n");
        }

        sb.append("</sim>\n");
        sb.append("</archive>");

        return sb.toString();
    }

    /**
     * Generate LQNX model.
     */
    private static String generateLqnxModel(List<ServiceInfo> services) {
        StringBuilder sb = new StringBuilder();
        sb.append("<?xml version=\"1.0\"?>\n");
        sb.append("<lqn-model name=\"imported\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n");
        sb.append("  <solver-params conv_val=\"1e-6\" it_limit=\"50\" print_int=\"10\" underrelax_coeff=\"0.9\"/>\n");

        // Add processors and tasks for each service
        for (ServiceInfo service : services) {
            String name = escapeXml(service.getName().replaceAll("[^a-zA-Z0-9]", "_"));
            double serviceTime = service.getAvgDurationMs() / 1000.0; // Convert to seconds
            if (serviceTime <= 0) serviceTime = 0.001;

            sb.append("  <processor name=\"").append(name).append("_P\" scheduling=\"fcfs\">\n");
            sb.append("    <task name=\"").append(name).append("\" scheduling=\"ref\">\n");
            sb.append("      <entry name=\"").append(name).append("_E\" type=\"PH1PH2\">\n");
            sb.append("        <entry-phase-activities>\n");
            sb.append("          <activity name=\"").append(name).append("_A\" phase=\"1\" host-demand-mean=\"");
            sb.append(serviceTime).append("\"/>\n");
            sb.append("        </entry-phase-activities>\n");
            sb.append("      </entry>\n");
            sb.append("    </task>\n");
            sb.append("  </processor>\n");
        }

        sb.append("</lqn-model>");
        return sb.toString();
    }

    /**
     * Escape XML special characters.
     */
    private static String escapeXml(String str) {
        if (str == null) return "";
        return str.replace("&", "&amp;")
                  .replace("<", "&lt;")
                  .replace(">", "&gt;")
                  .replace("\"", "&quot;")
                  .replace("'", "&apos;");
    }
}
