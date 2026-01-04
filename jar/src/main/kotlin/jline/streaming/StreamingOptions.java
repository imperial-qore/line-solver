package jline.streaming;

import java.io.Serializable;

/**
 * Configuration options for SSA streaming to OTLP receiver.
 * Supports two modes: sampled (push every N events) or time-window (push averaged metrics over duration).
 */
public class StreamingOptions implements Serializable {

    /**
     * Streaming mode enumeration.
     */
    public enum StreamMode {
        /** Push metrics at fixed event frequency */
        SAMPLED,
        /** Push time-averaged metrics over windows */
        TIME_WINDOW
    }

    /** OTLP receiver endpoint (default: localhost:4317) */
    public String endpoint = "localhost:4317";

    /** Streaming mode: SAMPLED or TIME_WINDOW */
    public StreamMode mode = StreamMode.SAMPLED;

    /** For SAMPLED mode: push every N simulation events */
    public int sampleFrequency = 100;

    /** For TIME_WINDOW mode: aggregate over X seconds before pushing */
    public double timeWindowSeconds = 1.0;

    /** Service name to identify this simulation in line-est */
    public String serviceName = "line-ssa";

    /** Whether to include queue length metrics (QLen) */
    public boolean includeQueueLength = true;

    /** Whether to include utilization metrics (Util) */
    public boolean includeUtilization = true;

    /** Whether to include throughput metrics (Tput) */
    public boolean includeThroughput = true;

    /** Whether to include response time metrics (RespT) */
    public boolean includeResponseTime = true;

    /** Whether to include arrival rate metrics (ArvR) */
    public boolean includeArrivalRate = true;

    /** If true, throw an exception when connection fails instead of falling back to dry-run mode */
    public boolean failOnConnectionError = false;

    /**
     * Default constructor with default settings.
     */
    public StreamingOptions() {
    }

    /**
     * Set the OTLP receiver endpoint.
     * @param endpoint The endpoint in format "host:port"
     * @return this instance for method chaining
     */
    public StreamingOptions endpoint(String endpoint) {
        this.endpoint = endpoint;
        return this;
    }

    /**
     * Set the streaming mode.
     * @param mode SAMPLED or TIME_WINDOW
     * @return this instance for method chaining
     */
    public StreamingOptions mode(StreamMode mode) {
        this.mode = mode;
        return this;
    }

    /**
     * Set the sample frequency for SAMPLED mode.
     * @param frequency Push every N simulation events
     * @return this instance for method chaining
     */
    public StreamingOptions sampleFrequency(int frequency) {
        this.sampleFrequency = frequency;
        return this;
    }

    /**
     * Set the time window duration for TIME_WINDOW mode.
     * @param seconds Aggregate duration in seconds
     * @return this instance for method chaining
     */
    public StreamingOptions timeWindowSeconds(double seconds) {
        this.timeWindowSeconds = seconds;
        return this;
    }

    /**
     * Set the service name for identification in line-est.
     * @param name Service name
     * @return this instance for method chaining
     */
    public StreamingOptions serviceName(String name) {
        this.serviceName = name;
        return this;
    }

    /**
     * Enable or disable queue length metrics.
     * @param include true to include, false to exclude
     * @return this instance for method chaining
     */
    public StreamingOptions includeQueueLength(boolean include) {
        this.includeQueueLength = include;
        return this;
    }

    /**
     * Enable or disable utilization metrics.
     * @param include true to include, false to exclude
     * @return this instance for method chaining
     */
    public StreamingOptions includeUtilization(boolean include) {
        this.includeUtilization = include;
        return this;
    }

    /**
     * Enable or disable throughput metrics.
     * @param include true to include, false to exclude
     * @return this instance for method chaining
     */
    public StreamingOptions includeThroughput(boolean include) {
        this.includeThroughput = include;
        return this;
    }

    /**
     * Enable or disable response time metrics.
     * @param include true to include, false to exclude
     * @return this instance for method chaining
     */
    public StreamingOptions includeResponseTime(boolean include) {
        this.includeResponseTime = include;
        return this;
    }

    /**
     * Enable or disable arrival rate metrics.
     * @param include true to include, false to exclude
     * @return this instance for method chaining
     */
    public StreamingOptions includeArrivalRate(boolean include) {
        this.includeArrivalRate = include;
        return this;
    }

    /**
     * Set whether to fail on connection error or fall back to dry-run mode.
     * @param fail true to throw exception on connection error, false for dry-run mode
     * @return this instance for method chaining
     */
    public StreamingOptions failOnConnectionError(boolean fail) {
        this.failOnConnectionError = fail;
        return this;
    }

    /**
     * Parse endpoint into host and port.
     * @return array with [host, port] where port is an int as String
     */
    public String[] parseEndpoint() {
        String[] parts = endpoint.split(":");
        if (parts.length == 2) {
            return parts;
        }
        return new String[] { endpoint, "4317" };
    }

    /**
     * Get the host from the endpoint.
     * @return host part of endpoint
     */
    public String getHost() {
        return parseEndpoint()[0];
    }

    /**
     * Get the port from the endpoint.
     * @return port number
     */
    public int getPort() {
        return Integer.parseInt(parseEndpoint()[1]);
    }
}
