/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.streaming;

import java.io.IOException;
import java.io.OutputStream;
import java.net.HttpURLConnection;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Simple HTTP client for sending metrics to a receiver.
 * Uses Java's built-in HttpURLConnection for maximum compatibility
 * with all JVM environments including MATLAB.
 *
 * Metrics are sent as JSON via HTTP POST.
 */
public class HttpMetricsClient {

    private static final Logger logger = Logger.getLogger(HttpMetricsClient.class.getName());

    private final String endpoint;
    private final String serviceName;
    private volatile boolean isShutdown = false;
    private volatile boolean connectionFailed = false;
    private final boolean failOnConnectionError;

    /**
     * Create a new HTTP metrics client.
     * @param endpoint The HTTP receiver endpoint (e.g., "http://localhost:8080/metrics")
     * @param serviceName The service name to identify this simulation
     */
    public HttpMetricsClient(String endpoint, String serviceName) {
        this(endpoint, serviceName, false);
    }

    /**
     * Create a new HTTP metrics client.
     * @param endpoint The HTTP receiver endpoint (e.g., "http://localhost:8080/metrics")
     * @param serviceName The service name to identify this simulation
     * @param failOnConnectionError If true, throw exception on connection failure
     */
    public HttpMetricsClient(String endpoint, String serviceName, boolean failOnConnectionError) {
        // Ensure endpoint has http:// prefix
        if (!endpoint.startsWith("http://") && !endpoint.startsWith("https://")) {
            this.endpoint = "http://" + endpoint;
        } else {
            this.endpoint = endpoint;
        }
        this.serviceName = serviceName;
        this.failOnConnectionError = failOnConnectionError;

        logger.log(Level.INFO, "HttpMetricsClient initialized: endpoint={0}, serviceName={1}",
                new Object[] { this.endpoint, serviceName });
    }

    /**
     * Send a batch of metrics to the HTTP receiver.
     * @param metrics List of metric points to send
     * @return true if successful, false otherwise
     */
    public boolean sendMetrics(List<SSAMetricPoint> metrics) {
        if (isShutdown) {
            logger.warning("Cannot send metrics: client is shutdown");
            return false;
        }

        if (connectionFailed && !failOnConnectionError) {
            // Dry-run mode - silently skip
            return false;
        }

        if (metrics == null || metrics.isEmpty()) {
            return true;
        }

        try {
            String json = metricsToJson(metrics);
            return postJson(json);
        } catch (Exception e) {
            if (failOnConnectionError) {
                throw new RuntimeException("Failed to send metrics: " + e.getMessage(), e);
            }
            connectionFailed = true;
            logger.log(Level.WARNING, "Failed to send metrics - switching to dry-run mode: {0}", e.getMessage());
            return false;
        }
    }

    /**
     * Convert metrics to JSON format.
     */
    private String metricsToJson(List<SSAMetricPoint> metrics) {
        StringBuilder sb = new StringBuilder();
        sb.append("{\"serviceName\":\"").append(escapeJson(serviceName)).append("\",");
        sb.append("\"metrics\":[");

        boolean first = true;
        for (SSAMetricPoint m : metrics) {
            if (!first) {
                sb.append(",");
            }
            first = false;

            sb.append("{\"name\":\"").append(escapeJson(m.metricName)).append("\",");
            sb.append("\"value\":").append(m.value).append(",");
            sb.append("\"timestamp\":").append(m.timestampNanos).append(",");
            sb.append("\"labels\":{");

            boolean firstLabel = true;
            for (Map.Entry<String, String> entry : m.labels.entrySet()) {
                if (!firstLabel) {
                    sb.append(",");
                }
                firstLabel = false;
                sb.append("\"").append(escapeJson(entry.getKey())).append("\":\"");
                sb.append(escapeJson(entry.getValue())).append("\"");
            }
            sb.append("}}");
        }

        sb.append("]}");
        return sb.toString();
    }

    /**
     * Escape special characters for JSON.
     */
    private String escapeJson(String s) {
        if (s == null) {
            return "";
        }
        return s.replace("\\", "\\\\")
                .replace("\"", "\\\"")
                .replace("\n", "\\n")
                .replace("\r", "\\r")
                .replace("\t", "\\t");
    }

    /**
     * POST JSON data to the endpoint.
     */
    private boolean postJson(String json) throws IOException {
        URL url = new URL(endpoint);
        HttpURLConnection conn = (HttpURLConnection) url.openConnection();

        try {
            conn.setRequestMethod("POST");
            conn.setRequestProperty("Content-Type", "application/json; charset=UTF-8");
            conn.setRequestProperty("Accept", "application/json");
            conn.setDoOutput(true);
            conn.setConnectTimeout(5000);
            conn.setReadTimeout(5000);

            byte[] data = json.getBytes(StandardCharsets.UTF_8);
            conn.setFixedLengthStreamingMode(data.length);

            try (OutputStream os = conn.getOutputStream()) {
                os.write(data);
                os.flush();
            }

            int responseCode = conn.getResponseCode();
            if (responseCode >= 200 && responseCode < 300) {
                return true;
            } else {
                logger.log(Level.WARNING, "HTTP POST failed with status {0}", responseCode);
                return false;
            }
        } finally {
            conn.disconnect();
        }
    }

    /**
     * Check if connection has failed.
     * @return true if connection failed and client is in dry-run mode
     */
    public boolean isConnectionFailed() {
        return connectionFailed;
    }

    /**
     * Shutdown the client.
     */
    public void shutdown() {
        isShutdown = true;
        logger.log(Level.INFO, "HttpMetricsClient shutdown");
    }
}
