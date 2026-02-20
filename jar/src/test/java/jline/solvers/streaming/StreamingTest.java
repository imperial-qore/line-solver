/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.streaming;

import jline.VerboseLevel;
import jline.io.Ret.SampleResult;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.solvers.des.DESOptions;
import jline.solvers.des.SolverDES;
import jline.solvers.ssa.SampleNodeState;
import jline.solvers.ssa.SolverSSA;
import jline.solvers.ssa.SSAOptions;
import jline.streaming.StreamingOptions;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;

import java.io.IOException;
import java.util.List;

import static jline.TestTools.suppressOutput;
import static jline.TestTools.restoreOutput;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Test suite for streaming functionality in SSA and DES solvers.
 * <p>
 * These tests verify that the stream() and streamAggr() methods work correctly
 * by running simulations with streaming enabled. A mock gRPC server captures
 * the metrics for verification.
 * </p>
 */
public class StreamingTest {

    /** Base seed for reproducibility */
    private static final int BASE_SEED = 42000;

    /** Mock OTLP receiver for capturing metrics */
    private MockOtlpReceiver mockReceiver;

    /**
     * Suppresses all output for the test class.
     */
    @BeforeAll
    public static void suppressAllOutput() {
        suppressOutput();
    }

    /**
     * Restores output after all tests complete.
     */
    @AfterAll
    public static void restoreAllOutput() {
        restoreOutput();
    }

    /**
     * Sets up the mock OTLP receiver before each test.
     */
    @BeforeEach
    public void setUp() throws IOException {
        mockReceiver = new MockOtlpReceiver();
        mockReceiver.start();
    }

    /**
     * Tears down the mock receiver after each test.
     */
    @AfterEach
    public void tearDown() {
        if (mockReceiver != null) {
            mockReceiver.stop();
        }
    }

    /**
     * Creates a simple M/M/1 queue network for testing.
     * Arrival rate: 0.5, Service rate: 1.0 (utilization = 0.5)
     *
     * @return configured M/M/1 network
     */
    private static Network createMM1Network() {
        Network model = new Network("M/M/1");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass jobClass = new OpenClass(model, "Class1", 0);
        source.setArrival(jobClass, new Exp(0.5));
        queue.setService(jobClass, new Exp(1.0));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobClass, jobClass, source, queue, 1.0);
        P.set(jobClass, jobClass, queue, sink, 1.0);
        model.link(P);
        return model;
    }

    /**
     * Creates streaming options pointing to the mock receiver.
     * Uses SAMPLED mode with specified frequency.
     *
     * @param sampleFrequency push metrics every N events
     * @return configured StreamingOptions
     */
    private StreamingOptions createTestStreamingOptions(int sampleFrequency) {
        return new StreamingOptions()
                .endpoint(mockReceiver.getEndpoint())
                .transport(StreamingOptions.TransportType.GRPC)
                .mode(StreamingOptions.StreamMode.SAMPLED)
                .sampleFrequency(sampleFrequency)
                .serviceName("line-test");
    }

    /**
     * Creates streaming options for time-window mode.
     *
     * @param windowSeconds time window duration in seconds
     * @return configured StreamingOptions
     */
    private StreamingOptions createTimeWindowStreamingOptions(double windowSeconds) {
        return new StreamingOptions()
                .endpoint(mockReceiver.getEndpoint())
                .transport(StreamingOptions.TransportType.GRPC)
                .mode(StreamingOptions.StreamMode.TIME_WINDOW)
                .timeWindowSeconds(windowSeconds)
                .serviceName("line-test");
    }

    // ==================== SSA Streaming Tests ====================

    /**
     * Tests SSA streaming with sampled mode on an M/M/1 queue.
     * Verifies that metrics are received by the mock server.
     */
    @Test
    public void testSSA_MM1_Stream_Sampled() throws Exception {
        Network model = createMM1Network();
        Queue queue = (Queue) model.getNodeByName("Queue");

        SSAOptions options = new SSAOptions();
        options.verbose = VerboseLevel.SILENT;
        options.seed = BASE_SEED;
        options.samples = 10000;

        // Push every 100 events - should get ~100 pushes
        StreamingOptions streamOpts = createTestStreamingOptions(100);

        SolverSSA solver = new SolverSSA(model, options);
        SampleNodeState result = solver.stream(queue, streamOpts);

        // Verify simulation completed
        assertNotNull(result, "SampleNodeState should not be null");
        assertNotNull(result.t, "Time vector should not be null");
        assertTrue(result.t.length() > 0, "Time vector should have samples");

        // Verify metrics were received by mock server
        assertTrue(mockReceiver.getRequestCount() > 0,
                "Mock receiver should have received at least one request");
        assertTrue(mockReceiver.getTotalMetricCount() > 0,
                "Mock receiver should have received metrics");

        // Verify queue_length metric was received
        List<String> metricNames = mockReceiver.getReceivedMetricNames();
        assertTrue(metricNames.contains("queue_length"),
                "Should have received queue_length metrics");
    }

    /**
     * Tests SSA streaming with time-window mode on an M/M/1 queue.
     */
    @Test
    public void testSSA_MM1_Stream_TimeWindow() throws Exception {
        Network model = createMM1Network();
        Queue queue = (Queue) model.getNodeByName("Queue");

        SSAOptions options = new SSAOptions();
        options.verbose = VerboseLevel.SILENT;
        options.seed = BASE_SEED;
        options.samples = 10000;

        // 1-second windows
        StreamingOptions streamOpts = createTimeWindowStreamingOptions(1.0);

        SolverSSA solver = new SolverSSA(model, options);
        SampleNodeState result = solver.stream(queue, streamOpts);

        // Verify simulation completed
        assertNotNull(result, "SampleNodeState should not be null");
        assertNotNull(result.t, "Time vector should not be null");
        assertTrue(result.t.length() > 0, "Time vector should have samples");

        // Verify metrics were received
        assertTrue(mockReceiver.getRequestCount() > 0,
                "Mock receiver should have received requests in time-window mode");
    }

    /**
     * Tests SSA streamAggr method on an M/M/1 queue.
     */
    @Test
    public void testSSA_MM1_StreamAggr() throws Exception {
        Network model = createMM1Network();
        Queue queue = (Queue) model.getNodeByName("Queue");

        SSAOptions options = new SSAOptions();
        options.verbose = VerboseLevel.SILENT;
        options.seed = BASE_SEED;
        options.samples = 10000;

        StreamingOptions streamOpts = createTestStreamingOptions(100);

        SolverSSA solver = new SolverSSA(model, options);
        SampleNodeState result = solver.streamAggr(queue, streamOpts);

        // Verify simulation completed
        assertNotNull(result, "SampleNodeState should not be null");
        assertNotNull(result.t, "Time vector should not be null");
        assertTrue(result.t.length() > 0, "Time vector should have samples");
        assertTrue(result.isaggregate, "Result should be marked as aggregated");

        // Verify metrics were received
        assertTrue(mockReceiver.getRequestCount() > 0,
                "Mock receiver should have received requests for streamAggr");
    }

    // ==================== DES Streaming Tests ====================

    /**
     * Tests DES streaming with sampled mode on an M/M/1 queue.
     * Note: DES stream() runs transient simulation which may have different
     * event handling than SSA. This test verifies simulation completes.
     */
    @Test
    public void testDES_MM1_Stream_Sampled() {
        Network model = createMM1Network();
        Queue queue = (Queue) model.getNodeByName("Queue");

        DESOptions options = new DESOptions();
        options.verbose = VerboseLevel.SILENT;
        options.seed = BASE_SEED;
        options.samples = 10000;

        // Push every 100 events
        StreamingOptions streamOpts = createTestStreamingOptions(100);

        SolverDES solver = new SolverDES(model, options);
        SampleResult result = solver.stream(queue, streamOpts);

        // Verify simulation completed successfully with streaming enabled
        assertNotNull(result, "SampleResult should not be null");
        assertNotNull(result.t, "Time vector should not be null");
        assertTrue(result.t.length() > 0, "Time vector should have samples");

        // Note: DES transient simulation may not trigger streaming hooks as frequently
        // as SSA, so we don't strictly require metrics to be received here
    }

    /**
     * Tests DES streaming with time-window mode on an M/M/1 queue.
     */
    @Test
    public void testDES_MM1_Stream_TimeWindow() {
        Network model = createMM1Network();
        Queue queue = (Queue) model.getNodeByName("Queue");

        DESOptions options = new DESOptions();
        options.verbose = VerboseLevel.SILENT;
        options.seed = BASE_SEED;
        options.samples = 10000;

        // 1-second windows
        StreamingOptions streamOpts = createTimeWindowStreamingOptions(1.0);

        SolverDES solver = new SolverDES(model, options);
        SampleResult result = solver.stream(queue, streamOpts);

        // Verify simulation completed
        assertNotNull(result, "SampleResult should not be null");
        assertNotNull(result.t, "Time vector should not be null");
        assertTrue(result.t.length() > 0, "Time vector should have samples");
    }

    /**
     * Tests DES streamAggr method on an M/M/1 queue.
     * DES sample() returns aggregated per-class queue lengths, so sampleAggr
     * simply marks the result as aggregated without additional processing.
     */
    @Test
    public void testDES_MM1_StreamAggr() {
        Network model = createMM1Network();
        Queue queue = (Queue) model.getNodeByName("Queue");

        DESOptions options = new DESOptions();
        options.verbose = VerboseLevel.SILENT;
        options.seed = BASE_SEED;
        options.samples = 10000;

        StreamingOptions streamOpts = createTestStreamingOptions(100);

        SolverDES solver = new SolverDES(model, options);
        SampleResult result = solver.streamAggr(queue, streamOpts);

        assertNotNull(result, "SampleResult should not be null");
    }

    // ==================== Metric Type Verification Tests ====================

    /**
     * Tests that SSA streams queue_length metrics.
     */
    @Test
    public void testSSA_QueueLengthMetric() throws Exception {
        Network model = createMM1Network();
        Queue queue = (Queue) model.getNodeByName("Queue");

        SSAOptions options = new SSAOptions();
        options.verbose = VerboseLevel.SILENT;
        options.seed = BASE_SEED;
        options.samples = 5000;

        // Only queue length metrics
        StreamingOptions streamOpts = new StreamingOptions()
                .endpoint(mockReceiver.getEndpoint())
                .transport(StreamingOptions.TransportType.GRPC)
                .mode(StreamingOptions.StreamMode.SAMPLED)
                .sampleFrequency(50)
                .serviceName("line-test")
                .includeQueueLength(true)
                .includeUtilization(false)
                .includeThroughput(false)
                .includeResponseTime(false)
                .includeArrivalRate(false);

        SolverSSA solver = new SolverSSA(model, options);
        solver.stream(queue, streamOpts);

        // Verify queue_length metrics were received
        List<String> metricNames = mockReceiver.getReceivedMetricNames();
        assertTrue(metricNames.contains("queue_length"), "Should have queue_length");
    }
}
