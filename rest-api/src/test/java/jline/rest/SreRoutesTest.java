/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.rest;

import jline.rest.model.*;
import jline.rest.util.JsonTransformer;
import org.junit.jupiter.api.*;
import spark.Spark;

import java.io.*;
import java.net.HttpURLConnection;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.util.*;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Tests for SRE/DevOps integration endpoints (metrics, calibration, trace import).
 */
@TestMethodOrder(MethodOrderer.OrderAnnotation.class)
class SreRoutesTest {

    private static final int TEST_PORT = 4570;
    private static final String BASE_URL = "http://localhost:" + TEST_PORT + "/api/v1";
    private static LineRestServer server;
    private static JsonTransformer json;

    @BeforeAll
    static void setUp() throws Exception {
        json = new JsonTransformer();
        server = new LineRestServer(TEST_PORT);
        server.start();
        Thread.sleep(1000);
    }

    @AfterAll
    static void tearDown() {
        if (server != null) {
            server.stopServer();
        }
        Spark.awaitStop();
    }

    @Test
    @Order(1)
    void testPrometheusMetrics() throws Exception {
        HttpURLConnection conn = createConnection("/metrics", "GET");
        int status = conn.getResponseCode();
        String response = readResponse(conn);

        assertEquals(200, status);
        assertNotNull(response);
        assertTrue(response.contains("line_requests_total"));
        assertTrue(response.contains("jvm_memory_heap_bytes_used"));
        assertTrue(response.contains("line_uptime_seconds"));
    }

    @Test
    @Order(2)
    void testJsonMetrics() throws Exception {
        HttpURLConnection conn = createConnection("/metrics/json", "GET");
        int status = conn.getResponseCode();
        String response = readResponse(conn);

        assertEquals(200, status);
        assertNotNull(response);
        assertTrue(response.contains("requests"));
        assertTrue(response.contains("jvm"));
        assertTrue(response.contains("heapUsedBytes"));
    }

    @Test
    @Order(3)
    void testCalibration() throws Exception {
        CalibrationRequest request = new CalibrationRequest();
        request.setSolver("mva");

        List<CalibrationRequest.ObservedMetrics> observations = new ArrayList<CalibrationRequest.ObservedMetrics>();
        CalibrationRequest.ObservedMetrics obs = new CalibrationRequest.ObservedMetrics();
        obs.setStation("Queue1");
        obs.setThroughput(10.0);
        obs.setUtilization(0.5);
        obs.setResponseTime(0.1);
        observations.add(obs);
        request.setObservations(observations);

        HttpURLConnection conn = createConnection("/models/calibrate", "POST");
        conn.setDoOutput(true);
        conn.getOutputStream().write(json.toJson(request).getBytes(StandardCharsets.UTF_8));

        int status = conn.getResponseCode();
        String response = readResponse(conn);

        assertEquals(200, status);
        assertNotNull(response);
        assertTrue(response.contains("status"));
        assertTrue(response.contains("parameters") || response.contains("error"));
    }

    @Test
    @Order(4)
    void testCalibrationInvalidRequest() throws Exception {
        CalibrationRequest request = new CalibrationRequest();
        // Missing observations and solver

        HttpURLConnection conn = createConnection("/models/calibrate", "POST");
        conn.setDoOutput(true);
        conn.getOutputStream().write(json.toJson(request).getBytes(StandardCharsets.UTF_8));

        int status = conn.getResponseCode();
        assertEquals(400, status);
    }

    @Test
    @Order(5)
    void testTraceImportJaeger() throws Exception {
        TraceImportRequest request = new TraceImportRequest();
        request.setFormat("jaeger");
        request.setOutputFormat("jsimg");

        // Minimal Jaeger trace structure
        List<Map<String, Object>> traces = new ArrayList<Map<String, Object>>();
        Map<String, Object> trace = new HashMap<String, Object>();

        List<Map<String, Object>> spans = new ArrayList<Map<String, Object>>();
        Map<String, Object> span = new HashMap<String, Object>();
        span.put("operationName", "GET /api/users");
        span.put("serviceName", "user-service");
        span.put("duration", 50000); // 50ms in microseconds
        spans.add(span);

        trace.put("spans", spans);
        traces.add(trace);
        request.setTraces(traces);

        HttpURLConnection conn = createConnection("/traces/import", "POST");
        conn.setDoOutput(true);
        conn.getOutputStream().write(json.toJson(request).getBytes(StandardCharsets.UTF_8));

        int status = conn.getResponseCode();
        String response = readResponse(conn);

        assertEquals(200, status);
        assertNotNull(response);
        assertTrue(response.contains("status"));
    }

    @Test
    @Order(6)
    void testTraceImportOpenTelemetry() throws Exception {
        TraceImportRequest request = new TraceImportRequest();
        request.setFormat("opentelemetry");
        request.setOutputFormat("jsimg");

        // Minimal OpenTelemetry trace structure
        List<Map<String, Object>> traces = new ArrayList<Map<String, Object>>();
        Map<String, Object> trace = new HashMap<String, Object>();

        List<Map<String, Object>> resourceSpans = new ArrayList<Map<String, Object>>();
        Map<String, Object> rs = new HashMap<String, Object>();

        List<Map<String, Object>> scopeSpans = new ArrayList<Map<String, Object>>();
        Map<String, Object> ss = new HashMap<String, Object>();

        List<Map<String, Object>> spans = new ArrayList<Map<String, Object>>();
        Map<String, Object> span = new HashMap<String, Object>();
        span.put("name", "HTTP GET");
        span.put("startTimeUnixNano", 1000000000L);
        span.put("endTimeUnixNano", 1050000000L); // 50ms duration
        spans.add(span);

        ss.put("spans", spans);
        scopeSpans.add(ss);
        rs.put("scopeSpans", scopeSpans);
        resourceSpans.add(rs);
        trace.put("resourceSpans", resourceSpans);
        traces.add(trace);
        request.setTraces(traces);

        HttpURLConnection conn = createConnection("/traces/import", "POST");
        conn.setDoOutput(true);
        conn.getOutputStream().write(json.toJson(request).getBytes(StandardCharsets.UTF_8));

        int status = conn.getResponseCode();
        String response = readResponse(conn);

        assertEquals(200, status);
        assertNotNull(response);
        assertTrue(response.contains("status"));
    }

    @Test
    @Order(7)
    void testTraceImportInvalidFormat() throws Exception {
        TraceImportRequest request = new TraceImportRequest();
        request.setFormat("invalid_format");
        request.setTraces(new ArrayList<Object>());

        HttpURLConnection conn = createConnection("/traces/import", "POST");
        conn.setDoOutput(true);
        conn.getOutputStream().write(json.toJson(request).getBytes(StandardCharsets.UTF_8));

        int status = conn.getResponseCode();
        assertEquals(400, status);
    }

    @Test
    @Order(8)
    void testTraceImportMissingTraces() throws Exception {
        TraceImportRequest request = new TraceImportRequest();
        request.setFormat("jaeger");
        // Missing traces

        HttpURLConnection conn = createConnection("/traces/import", "POST");
        conn.setDoOutput(true);
        conn.getOutputStream().write(json.toJson(request).getBytes(StandardCharsets.UTF_8));

        int status = conn.getResponseCode();
        assertEquals(400, status);
    }

    @Test
    @Order(9)
    void testTraceImportLqnxOutput() throws Exception {
        TraceImportRequest request = new TraceImportRequest();
        request.setFormat("jaeger");
        request.setOutputFormat("lqnx");

        List<Map<String, Object>> traces = new ArrayList<Map<String, Object>>();
        Map<String, Object> trace = new HashMap<String, Object>();

        List<Map<String, Object>> spans = new ArrayList<Map<String, Object>>();
        Map<String, Object> span = new HashMap<String, Object>();
        span.put("operationName", "process");
        span.put("serviceName", "backend");
        span.put("duration", 100000);
        spans.add(span);

        trace.put("spans", spans);
        traces.add(trace);
        request.setTraces(traces);

        HttpURLConnection conn = createConnection("/traces/import", "POST");
        conn.setDoOutput(true);
        conn.getOutputStream().write(json.toJson(request).getBytes(StandardCharsets.UTF_8));

        int status = conn.getResponseCode();
        String response = readResponse(conn);

        assertEquals(200, status);
        assertTrue(response.contains("lqnx") || response.contains("outputFormat"));
    }

    @Test
    @Order(10)
    void testCalibrationWithOptions() throws Exception {
        CalibrationRequest request = new CalibrationRequest();
        request.setSolver("mva");

        List<CalibrationRequest.ObservedMetrics> observations = new ArrayList<CalibrationRequest.ObservedMetrics>();
        CalibrationRequest.ObservedMetrics obs = new CalibrationRequest.ObservedMetrics();
        obs.setStation("WebServer");
        obs.setClassName("RequestClass");
        obs.setThroughput(100.0);
        obs.setUtilization(0.8);
        obs.setQueueLength(4.0);
        observations.add(obs);
        request.setObservations(observations);

        CalibrationRequest.CalibrationOptions options = new CalibrationRequest.CalibrationOptions();
        options.setMaxIterations(50);
        options.setTolerance(1e-5);
        request.setOptions(options);

        HttpURLConnection conn = createConnection("/models/calibrate", "POST");
        conn.setDoOutput(true);
        conn.getOutputStream().write(json.toJson(request).getBytes(StandardCharsets.UTF_8));

        int status = conn.getResponseCode();
        String response = readResponse(conn);

        assertEquals(200, status);
        assertTrue(response.contains("fittedValue") || response.contains("parameters"));
    }

    // Helper methods

    private HttpURLConnection createConnection(String path, String method) throws Exception {
        URL url = new URL(BASE_URL + path);
        HttpURLConnection conn = (HttpURLConnection) url.openConnection();
        conn.setRequestMethod(method);
        conn.setRequestProperty("Content-Type", "application/json");
        conn.setRequestProperty("Accept", "application/json");
        return conn;
    }

    private String readResponse(HttpURLConnection conn) throws Exception {
        InputStream is;
        try {
            is = conn.getInputStream();
        } catch (IOException e) {
            is = conn.getErrorStream();
        }

        if (is == null) {
            return "";
        }

        BufferedReader reader = new BufferedReader(new InputStreamReader(is, StandardCharsets.UTF_8));
        StringBuilder sb = new StringBuilder();
        String line;
        while ((line = reader.readLine()) != null) {
            sb.append(line);
        }
        reader.close();
        return sb.toString();
    }
}
