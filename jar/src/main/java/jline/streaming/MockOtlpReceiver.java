/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.streaming;

import io.grpc.Server;
import io.grpc.ServerBuilder;
import io.grpc.stub.StreamObserver;
import io.opentelemetry.proto.collector.metrics.v1.ExportMetricsServiceRequest;
import io.opentelemetry.proto.collector.metrics.v1.ExportMetricsServiceResponse;
import io.opentelemetry.proto.collector.metrics.v1.MetricsServiceGrpc;
import io.opentelemetry.proto.metrics.v1.Metric;
import io.opentelemetry.proto.metrics.v1.ResourceMetrics;
import io.opentelemetry.proto.metrics.v1.ScopeMetrics;

import java.io.IOException;
import java.net.ServerSocket;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.TimeUnit;

/**
 * Mock OTLP gRPC receiver for testing streaming functionality.
 * <p>
 * This mock server implements the OTLP MetricsService and captures all
 * received metrics for verification in tests.
 * </p>
 */
public class MockOtlpReceiver {

    private final Server server;
    private final int port;
    private final MockMetricsService metricsService;

    /**
     * Creates a mock OTLP receiver on a random available port.
     *
     * @throws IOException if unable to find an available port or start the server
     */
    public MockOtlpReceiver() throws IOException {
        this.port = findAvailablePort();
        this.metricsService = new MockMetricsService();
        this.server = ServerBuilder.forPort(port)
                .addService(metricsService)
                .build();
    }

    /**
     * Starts the mock receiver.
     *
     * @throws IOException if the server cannot be started
     */
    public void start() throws IOException {
        server.start();
        // Wait a bit for the server to be fully ready to accept connections
        waitUntilReady(1000);
    }

    /**
     * Waits until the server is ready to accept connections.
     * This is important in some JVM environments (like MATLAB) where
     * the server may not be immediately ready after start() returns.
     *
     * @param timeoutMs maximum time to wait in milliseconds
     */
    public void waitUntilReady(long timeoutMs) {
        long startTime = System.currentTimeMillis();
        while (System.currentTimeMillis() - startTime < timeoutMs) {
            try (java.net.Socket socket = new java.net.Socket()) {
                socket.connect(new java.net.InetSocketAddress("127.0.0.1", port), 100);
                // Connection successful, server is ready
                return;
            } catch (Exception e) {
                // Server not ready yet, wait a bit
                try {
                    Thread.sleep(10);
                } catch (InterruptedException ie) {
                    Thread.currentThread().interrupt();
                    return;
                }
            }
        }
        // Timeout reached, but don't throw - let the actual gRPC call handle the error
    }

    /**
     * Stops the mock receiver gracefully.
     */
    public void stop() {
        if (server != null) {
            try {
                server.shutdown();
                server.awaitTermination(5, TimeUnit.SECONDS);
            } catch (InterruptedException e) {
                server.shutdownNow();
                Thread.currentThread().interrupt();
            }
        }
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
     * Gets the endpoint string for this receiver (localhost:port).
     *
     * @return endpoint string
     */
    public String getEndpoint() {
        return "localhost:" + port;
    }

    /**
     * Gets the total number of export requests received.
     *
     * @return number of requests
     */
    public int getRequestCount() {
        return metricsService.getRequestCount();
    }

    /**
     * Gets the total number of metrics received across all requests.
     *
     * @return total metric count
     */
    public int getTotalMetricCount() {
        return metricsService.getTotalMetricCount();
    }

    /**
     * Gets all metric names received.
     *
     * @return list of metric names
     */
    public List<String> getReceivedMetricNames() {
        return metricsService.getReceivedMetricNames();
    }

    /**
     * Gets all received export requests.
     *
     * @return list of requests
     */
    public List<ExportMetricsServiceRequest> getReceivedRequests() {
        return metricsService.getReceivedRequests();
    }

    /**
     * Clears all received data.
     */
    public void reset() {
        metricsService.reset();
    }

    /**
     * Checks if a specific metric type was received.
     *
     * @param metricName the metric name to check (e.g., "queue_length", "throughput")
     * @return true if at least one metric with that name was received
     */
    public boolean hasReceivedMetric(String metricName) {
        return getReceivedMetricNames().contains(metricName);
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
     * Mock implementation of the OTLP MetricsService.
     */
    private static class MockMetricsService extends MetricsServiceGrpc.MetricsServiceImplBase {

        private final List<ExportMetricsServiceRequest> receivedRequests =
                Collections.synchronizedList(new ArrayList<ExportMetricsServiceRequest>());

        @Override
        public void export(ExportMetricsServiceRequest request,
                           StreamObserver<ExportMetricsServiceResponse> responseObserver) {
            // Store the request
            receivedRequests.add(request);

            // Send success response
            responseObserver.onNext(ExportMetricsServiceResponse.newBuilder().build());
            responseObserver.onCompleted();
        }

        public int getRequestCount() {
            return receivedRequests.size();
        }

        public int getTotalMetricCount() {
            int count = 0;
            for (ExportMetricsServiceRequest request : receivedRequests) {
                for (ResourceMetrics rm : request.getResourceMetricsList()) {
                    for (ScopeMetrics sm : rm.getScopeMetricsList()) {
                        count += sm.getMetricsCount();
                    }
                }
            }
            return count;
        }

        public List<String> getReceivedMetricNames() {
            List<String> names = new ArrayList<String>();
            for (ExportMetricsServiceRequest request : receivedRequests) {
                for (ResourceMetrics rm : request.getResourceMetricsList()) {
                    for (ScopeMetrics sm : rm.getScopeMetricsList()) {
                        for (Metric m : sm.getMetricsList()) {
                            if (!names.contains(m.getName())) {
                                names.add(m.getName());
                            }
                        }
                    }
                }
            }
            return names;
        }

        public List<ExportMetricsServiceRequest> getReceivedRequests() {
            return new ArrayList<ExportMetricsServiceRequest>(receivedRequests);
        }

        public void reset() {
            receivedRequests.clear();
        }
    }
}
