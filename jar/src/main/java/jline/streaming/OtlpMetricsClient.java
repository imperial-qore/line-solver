package jline.streaming;

import io.grpc.ManagedChannel;
import io.grpc.ManagedChannelBuilder;
import io.opentelemetry.proto.collector.metrics.v1.ExportMetricsServiceRequest;
import io.opentelemetry.proto.collector.metrics.v1.ExportMetricsServiceResponse;
import io.opentelemetry.proto.collector.metrics.v1.MetricsServiceGrpc;
import io.opentelemetry.proto.common.v1.AnyValue;
import io.opentelemetry.proto.common.v1.KeyValue;
import io.opentelemetry.proto.metrics.v1.Gauge;
import io.opentelemetry.proto.metrics.v1.Metric;
import io.opentelemetry.proto.metrics.v1.NumberDataPoint;
import io.opentelemetry.proto.metrics.v1.ResourceMetrics;
import io.opentelemetry.proto.metrics.v1.ScopeMetrics;
import io.opentelemetry.proto.resource.v1.Resource;

import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * OTLP gRPC client for sending metrics to line-est receiver.
 * Uses OpenTelemetry protobuf format for metric export.
 *
 * This client is Java 8 compatible and connects to line-est's
 * OtlpMetricsReceiver on the configured endpoint (default: localhost:4317).
 */
public class OtlpMetricsClient {

    private static final Logger logger = Logger.getLogger(OtlpMetricsClient.class.getName());

    private ManagedChannel channel;
    private MetricsServiceGrpc.MetricsServiceBlockingStub blockingStub;
    private final String serviceName;
    private volatile boolean isShutdown = false;
    private volatile boolean connectionFailed = false;
    private final boolean failOnConnectionError;

    /**
     * Create a new OTLP metrics client (dry-run mode allowed on connection failure).
     * @param endpoint The OTLP receiver endpoint in format "host:port"
     * @param serviceName The service name to identify this simulation
     */
    public OtlpMetricsClient(String endpoint, String serviceName) {
        this(endpoint, serviceName, false);
    }

    /**
     * Create a new OTLP metrics client.
     * @param endpoint The OTLP receiver endpoint in format "host:port"
     * @param serviceName The service name to identify this simulation
     * @param failOnConnectionError If true, throw exception on connection failure; if false, use dry-run mode
     */
    public OtlpMetricsClient(String endpoint, String serviceName, boolean failOnConnectionError) {
        this.serviceName = serviceName;
        this.failOnConnectionError = failOnConnectionError;

        String host;
        int port;

        String[] parts = endpoint.split(":");
        if (parts.length == 2) {
            host = parts[0];
            port = Integer.parseInt(parts[1]);
        } else {
            host = endpoint;
            port = 4317;
        }

        try {
            // Use IP address if localhost to avoid DNS issues
            String resolvedHost = "localhost".equalsIgnoreCase(host) ? "127.0.0.1" : host;
            // Use ManagedChannelBuilder for better compatibility with different JVM environments (like MATLAB)
            this.channel = ManagedChannelBuilder.forAddress(resolvedHost, port)
                    .usePlaintext()
                    .build();
            this.blockingStub = MetricsServiceGrpc.newBlockingStub(channel);

            logger.log(Level.INFO, "OtlpMetricsClient initialized: endpoint={0}:{1}, serviceName={2}",
                    new Object[] { resolvedHost, port, serviceName });
        } catch (Throwable e) {
            this.channel = null;
            this.blockingStub = null;
            this.connectionFailed = true;
            if (failOnConnectionError) {
                throw new RuntimeException("OtlpMetricsClient failed to connect to " + host + ":" + port + ": " + e.getMessage(), e);
            }
            logger.log(Level.WARNING, "OtlpMetricsClient failed to connect to {0}:{1} - running in dry-run mode: {2}",
                    new Object[] { host, port, e.getMessage() });
        }
    }

    /**
     * Send a batch of metrics to the OTLP receiver.
     * @param metrics List of metric points to send
     * @return true if successful, false otherwise
     */
    public boolean sendMetrics(List<SSAMetricPoint> metrics) {
        if (isShutdown) {
            logger.warning("Cannot send metrics: client is shutdown");
            return false;
        }

        if (connectionFailed || blockingStub == null) {
            // Dry-run mode - silently succeed without sending
            return true;
        }

        if (metrics == null || metrics.isEmpty()) {
            return true;  // Nothing to send
        }

        try {
            ExportMetricsServiceRequest request = buildRequest(metrics);
            ExportMetricsServiceResponse response = blockingStub.export(request);
            return true;
        } catch (Exception e) {
            logger.log(Level.WARNING, "Failed to send metrics to OTLP receiver: {0}", e.getMessage());
            return false;
        }
    }

    /**
     * Build an OTLP ExportMetricsServiceRequest from a list of metric points.
     * @param metrics List of metric points
     * @return OTLP request ready to send
     */
    private ExportMetricsServiceRequest buildRequest(List<SSAMetricPoint> metrics) {
        // Build resource with service name
        Resource.Builder resourceBuilder = Resource.newBuilder()
                .addAttributes(KeyValue.newBuilder()
                        .setKey("service.name")
                        .setValue(AnyValue.newBuilder()
                                .setStringValue(serviceName)
                                .build())
                        .build());

        // Build scope metrics containing all our metrics
        ScopeMetrics.Builder scopeMetricsBuilder = ScopeMetrics.newBuilder();

        for (SSAMetricPoint point : metrics) {
            Metric.Builder metricBuilder = Metric.newBuilder()
                    .setName(point.metricName);

            // Build data point
            NumberDataPoint.Builder dataPointBuilder = NumberDataPoint.newBuilder()
                    .setAsDouble(point.value)
                    .setTimeUnixNano(point.timestampNanos);

            // Add labels as attributes
            for (Map.Entry<String, String> label : point.labels.entrySet()) {
                dataPointBuilder.addAttributes(KeyValue.newBuilder()
                        .setKey(label.getKey())
                        .setValue(AnyValue.newBuilder()
                                .setStringValue(label.getValue())
                                .build())
                        .build());
            }

            // Use Gauge type for point-in-time metrics
            metricBuilder.setGauge(Gauge.newBuilder()
                    .addDataPoints(dataPointBuilder.build())
                    .build());

            scopeMetricsBuilder.addMetrics(metricBuilder.build());
        }

        // Build resource metrics
        ResourceMetrics resourceMetrics = ResourceMetrics.newBuilder()
                .setResource(resourceBuilder.build())
                .addScopeMetrics(scopeMetricsBuilder.build())
                .build();

        // Build final request
        return ExportMetricsServiceRequest.newBuilder()
                .addResourceMetrics(resourceMetrics)
                .build();
    }

    /**
     * Check if the client is connected and ready.
     * @return true if connected
     */
    public boolean isConnected() {
        if (connectionFailed || channel == null) {
            return false;
        }
        return !isShutdown && !channel.isShutdown() && !channel.isTerminated();
    }

    /**
     * Shutdown the client gracefully.
     * Waits up to 5 seconds for pending operations.
     */
    public void shutdown() {
        if (isShutdown) {
            return;
        }
        isShutdown = true;

        if (channel == null) {
            // No channel to shutdown (dry-run mode)
            return;
        }

        logger.log(Level.INFO, "Shutting down OtlpMetricsClient");
        try {
            channel.shutdown();
            if (!channel.awaitTermination(5, TimeUnit.SECONDS)) {
                logger.warning("Channel did not terminate in time, forcing shutdown");
                channel.shutdownNow();
            }
        } catch (InterruptedException e) {
            logger.warning("Interrupted during shutdown, forcing shutdown");
            channel.shutdownNow();
            Thread.currentThread().interrupt();
        }
    }

    /**
     * Get the service name configured for this client.
     * @return service name
     */
    public String getServiceName() {
        return serviceName;
    }
}
