/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.streaming;

import com.sun.net.httpserver.HttpExchange;
import com.sun.net.httpserver.HttpHandler;
import com.sun.net.httpserver.HttpServer;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.net.InetSocketAddress;
import java.net.ServerSocket;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.Executors;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Mock HTTP receiver for testing streaming functionality.
 * Uses Java's built-in HttpServer for maximum compatibility.
 *
 * Receives JSON metrics via HTTP POST and captures them for verification.
 */
public class MockHttpReceiver {

    private final HttpServer server;
    private final int port;
    private final MetricsHandler handler;

    /**
     * Creates a mock HTTP receiver on a random available port.
     *
     * @throws IOException if unable to find an available port or start the server
     */
    public MockHttpReceiver() throws IOException {
        this.port = findAvailablePort();
        this.handler = new MetricsHandler();
        this.server = HttpServer.create(new InetSocketAddress(port), 0);
        this.server.createContext("/metrics", handler);
        this.server.setExecutor(Executors.newSingleThreadExecutor());
    }

    /**
     * Starts the mock receiver.
     */
    public void start() {
        server.start();
    }

    /**
     * Stops the mock receiver gracefully.
     */
    public void stop() {
        server.stop(1);
    }

    /**
     * Gets the port this receiver is listening on.
     *
     * @return the port number
     */
    public int getPort() {
        return port;
    }

    /**
     * Gets the endpoint string for this receiver.
     *
     * @return endpoint string (e.g., "http://localhost:12345/metrics")
     */
    public String getEndpoint() {
        return "http://localhost:" + port + "/metrics";
    }

    /**
     * Gets the total number of POST requests received.
     *
     * @return number of requests
     */
    public int getRequestCount() {
        return handler.getRequestCount();
    }

    /**
     * Gets the total number of metrics received across all requests.
     *
     * @return total metric count
     */
    public int getTotalMetricCount() {
        return handler.getTotalMetricCount();
    }

    /**
     * Gets all metric names received.
     *
     * @return list of metric names
     */
    public List<String> getReceivedMetricNames() {
        return handler.getReceivedMetricNames();
    }

    /**
     * Gets all received JSON payloads.
     *
     * @return list of JSON strings
     */
    public List<String> getReceivedPayloads() {
        return handler.getReceivedPayloads();
    }

    /**
     * Checks if a specific metric type was received.
     *
     * @param metricName the metric name to check
     * @return true if at least one metric with that name was received
     */
    public boolean hasReceivedMetric(String metricName) {
        return getReceivedMetricNames().contains(metricName);
    }

    /**
     * Clears all received data.
     */
    public void reset() {
        handler.reset();
    }

    /**
     * Finds an available port on localhost.
     */
    private static int findAvailablePort() throws IOException {
        try (ServerSocket socket = new ServerSocket(0)) {
            socket.setReuseAddress(true);
            return socket.getLocalPort();
        }
    }

    /**
     * HTTP handler that captures metrics from POST requests.
     */
    private static class MetricsHandler implements HttpHandler {

        private final List<String> receivedPayloads = Collections.synchronizedList(new ArrayList<String>());
        private final Set<String> metricNames = Collections.synchronizedSet(new HashSet<String>());
        private int totalMetricCount = 0;

        // Simple pattern to extract metric names from JSON
        private static final Pattern METRIC_NAME_PATTERN = Pattern.compile("\"name\"\\s*:\\s*\"([^\"]+)\"");

        @Override
        public void handle(HttpExchange exchange) throws IOException {
            if ("POST".equals(exchange.getRequestMethod())) {
                // Read request body
                StringBuilder sb = new StringBuilder();
                try (BufferedReader reader = new BufferedReader(
                        new InputStreamReader(exchange.getRequestBody(), StandardCharsets.UTF_8))) {
                    String line;
                    while ((line = reader.readLine()) != null) {
                        sb.append(line);
                    }
                }

                String json = sb.toString();
                receivedPayloads.add(json);

                // Extract metric names and count
                Matcher matcher = METRIC_NAME_PATTERN.matcher(json);
                int count = 0;
                while (matcher.find()) {
                    metricNames.add(matcher.group(1));
                    count++;
                }
                synchronized (this) {
                    totalMetricCount += count;
                }

                // Send 200 OK response
                String response = "{\"status\":\"ok\"}";
                byte[] responseBytes = response.getBytes(StandardCharsets.UTF_8);
                exchange.getResponseHeaders().set("Content-Type", "application/json");
                exchange.sendResponseHeaders(200, responseBytes.length);
                try (OutputStream os = exchange.getResponseBody()) {
                    os.write(responseBytes);
                }
            } else {
                // Method not allowed
                exchange.sendResponseHeaders(405, -1);
            }
        }

        public int getRequestCount() {
            return receivedPayloads.size();
        }

        public synchronized int getTotalMetricCount() {
            return totalMetricCount;
        }

        public List<String> getReceivedMetricNames() {
            return new ArrayList<String>(metricNames);
        }

        public List<String> getReceivedPayloads() {
            return new ArrayList<String>(receivedPayloads);
        }

        public synchronized void reset() {
            receivedPayloads.clear();
            metricNames.clear();
            totalMetricCount = 0;
        }
    }
}
