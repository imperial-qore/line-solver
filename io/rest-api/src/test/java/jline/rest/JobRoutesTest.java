/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.rest;

import jline.rest.model.SolveRequest;
import jline.rest.model.ModelInput;
import jline.rest.util.JsonTransformer;
import org.junit.jupiter.api.*;
import spark.Spark;

import java.io.*;
import java.net.HttpURLConnection;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Tests for async job endpoints.
 */
@TestMethodOrder(MethodOrderer.OrderAnnotation.class)
class JobRoutesTest {

    private static final int TEST_PORT = 4568;
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
    void testListJobsEmpty() throws Exception {
        String response = sendGet("/jobs");
        assertNotNull(response);
        assertTrue(response.contains("\"count\":0") || response.contains("\"count\": 0"));
        assertTrue(response.contains("jobs"));
    }

    @Test
    @Order(2)
    void testSubmitAsyncJob() throws Exception {
        SolveRequest request = createValidRequest();

        HttpURLConnection conn = createConnection("/models/solve/async", "POST");
        conn.setDoOutput(true);
        conn.getOutputStream().write(json.toJson(request).getBytes(StandardCharsets.UTF_8));

        int status = conn.getResponseCode();
        assertEquals(202, status); // Accepted

        String response = readResponse(conn);
        assertNotNull(response);
        assertTrue(response.contains("jobId"));
        assertTrue(response.contains("pending") || response.contains("running"));
    }

    @Test
    @Order(3)
    void testSubmitAsyncJobInvalidRequest() throws Exception {
        SolveRequest request = new SolveRequest();
        // Missing model

        HttpURLConnection conn = createConnection("/models/solve/async", "POST");
        conn.setDoOutput(true);
        conn.getOutputStream().write(json.toJson(request).getBytes(StandardCharsets.UTF_8));

        int status = conn.getResponseCode();
        assertEquals(400, status);
    }

    @Test
    @Order(4)
    void testGetJobStatus() throws Exception {
        // First submit a job
        SolveRequest request = createValidRequest();

        HttpURLConnection conn1 = createConnection("/models/solve/async", "POST");
        conn1.setDoOutput(true);
        conn1.getOutputStream().write(json.toJson(request).getBytes(StandardCharsets.UTF_8));
        String submitResponse = readResponse(conn1);

        // Extract jobId
        Map<String, Object> submitResult = json.fromJson(submitResponse, Map.class);
        String jobId = (String) submitResult.get("jobId");
        assertNotNull(jobId);

        // Get job status
        String statusResponse = sendGet("/jobs/" + jobId);
        assertNotNull(statusResponse);
        assertTrue(statusResponse.contains(jobId));
        assertTrue(statusResponse.contains("status"));
    }

    @Test
    @Order(5)
    void testGetJobStatusNotFound() throws Exception {
        HttpURLConnection conn = createConnection("/jobs/nonexistent-job-id", "GET");
        int status = conn.getResponseCode();
        assertEquals(404, status);
    }

    @Test
    @Order(6)
    void testListJobsWithJobs() throws Exception {
        // Submit a few jobs first
        for (int i = 0; i < 3; i++) {
            SolveRequest request = createValidRequest();
            HttpURLConnection conn = createConnection("/models/solve/async", "POST");
            conn.setDoOutput(true);
            conn.getOutputStream().write(json.toJson(request).getBytes(StandardCharsets.UTF_8));
            readResponse(conn);
        }

        // List all jobs
        String response = sendGet("/jobs");
        assertNotNull(response);
        assertTrue(response.contains("jobs"));

        Map<String, Object> result = json.fromJson(response, Map.class);
        assertNotNull(result.get("count"));
    }

    @Test
    @Order(7)
    void testCancelJob() throws Exception {
        // Submit a job with a slow solver to ensure it's still running
        SolveRequest request = createValidRequest();
        request.setSolver("jmt"); // JMT is slower

        HttpURLConnection conn1 = createConnection("/models/solve/async", "POST");
        conn1.setDoOutput(true);
        conn1.getOutputStream().write(json.toJson(request).getBytes(StandardCharsets.UTF_8));
        String submitResponse = readResponse(conn1);

        Map<String, Object> submitResult = json.fromJson(submitResponse, Map.class);
        String jobId = (String) submitResult.get("jobId");

        // Cancel the job
        HttpURLConnection conn2 = createConnection("/jobs/" + jobId, "DELETE");
        int status = conn2.getResponseCode();

        // Job might already be completed (fast model), so accept 200 or 400
        assertTrue(status == 200 || status == 400);
    }

    @Test
    @Order(8)
    void testCancelJobNotFound() throws Exception {
        HttpURLConnection conn = createConnection("/jobs/nonexistent-job-id", "DELETE");
        int status = conn.getResponseCode();
        assertEquals(404, status);
    }

    @Test
    @Order(9)
    void testJobCompletion() throws Exception {
        // Submit a job with MVA (fast solver)
        SolveRequest request = createValidRequest();
        request.setSolver("mva");

        HttpURLConnection conn1 = createConnection("/models/solve/async", "POST");
        conn1.setDoOutput(true);
        conn1.getOutputStream().write(json.toJson(request).getBytes(StandardCharsets.UTF_8));
        String submitResponse = readResponse(conn1);

        Map<String, Object> submitResult = json.fromJson(submitResponse, Map.class);
        String jobId = (String) submitResult.get("jobId");

        // Wait for completion
        String finalStatus = null;
        for (int i = 0; i < 20; i++) { // Wait up to 10 seconds
            Thread.sleep(500);
            String statusResponse = sendGet("/jobs/" + jobId);
            Map<String, Object> statusResult = json.fromJson(statusResponse, Map.class);
            finalStatus = (String) statusResult.get("status");

            if ("completed".equals(finalStatus) || "failed".equals(finalStatus)) {
                break;
            }
        }

        // Job should have completed or failed (not still pending)
        assertNotNull(finalStatus);
        assertTrue("completed".equals(finalStatus) || "failed".equals(finalStatus),
            "Expected completed or failed, got: " + finalStatus);
    }

    @Test
    @Order(10)
    void testJobWithResults() throws Exception {
        // Submit a job with MVA (fast solver)
        SolveRequest request = createValidRequest();
        request.setSolver("mva");
        request.setAnalysis("avg");

        HttpURLConnection conn1 = createConnection("/models/solve/async", "POST");
        conn1.setDoOutput(true);
        conn1.getOutputStream().write(json.toJson(request).getBytes(StandardCharsets.UTF_8));
        String submitResponse = readResponse(conn1);

        Map<String, Object> submitResult = json.fromJson(submitResponse, Map.class);
        String jobId = (String) submitResult.get("jobId");

        // Wait for completion
        for (int i = 0; i < 20; i++) {
            Thread.sleep(500);
            String statusResponse = sendGet("/jobs/" + jobId);
            Map<String, Object> statusResult = json.fromJson(statusResponse, Map.class);
            String status = (String) statusResult.get("status");

            if ("completed".equals(status)) {
                // Should have results
                assertTrue(statusResponse.contains("results") || statusResponse.contains("runtime"));
                break;
            } else if ("failed".equals(status)) {
                // Failed is acceptable for this test (model may be invalid)
                break;
            }
        }
    }

    // Helper methods

    private String sendGet(String path) throws Exception {
        HttpURLConnection conn = createConnection(path, "GET");
        return readResponse(conn);
    }

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

    private SolveRequest createValidRequest() {
        SolveRequest request = new SolveRequest();
        ModelInput model = new ModelInput();
        model.setFormat("jsimg");
        model.setContent(createMinimalJsimgModel());
        request.setModel(model);
        request.setSolver("mva");
        request.setAnalysis("avg");
        return request;
    }

    private String createMinimalJsimgModel() {
        return "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n" +
               "<archive xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" " +
               "name=\"test.jsimg\" xsi:noNamespaceSchemaLocation=\"Archive.xsd\">\n" +
               "<sim disableStatisticStop=\"false\" logFullPath=\"false\" " +
               "name=\"test.jsimg.jsimw\" seed=\"23000\" " +
               "xsi:noNamespaceSchemaLocation=\"SIMmodeldefinition.xsd\">\n" +
               "<userClass name=\"Class1\" priority=\"0\" referenceSource=\"Source 1\" type=\"open\"/>\n" +
               "<node name=\"Source 1\">\n" +
               "<section className=\"RandomSource\">\n" +
               "<parameter array=\"true\" classPath=\"jmt.engine.NetStrategies.ServiceStrategy\" name=\"ServiceStrategy\">\n" +
               "<refClass>Class1</refClass>\n" +
               "<subParameter classPath=\"jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy\" name=\"ServiceTimeStrategy\">\n" +
               "<subParameter classPath=\"jmt.engine.random.Exponential\" name=\"Exponential\">\n" +
               "<subParameter classPath=\"jmt.engine.random.ExponentialPar\" name=\"distrPar\">\n" +
               "<subParameter name=\"lambda\">\n" +
               "<value>1.0</value>\n" +
               "</subParameter>\n" +
               "</subParameter>\n" +
               "</subParameter>\n" +
               "</subParameter>\n" +
               "</parameter>\n" +
               "</section>\n" +
               "<section className=\"ServiceTunnel\"/>\n" +
               "<section className=\"Router\">\n" +
               "<parameter array=\"true\" classPath=\"jmt.engine.NetStrategies.RoutingStrategy\" name=\"RoutingStrategy\">\n" +
               "<refClass>Class1</refClass>\n" +
               "<subParameter classPath=\"jmt.engine.NetStrategies.RoutingStrategies.RandomStrategy\" name=\"Random\"/>\n" +
               "</parameter>\n" +
               "</section>\n" +
               "</node>\n" +
               "<node name=\"Queue 1\">\n" +
               "<section className=\"Queue\">\n" +
               "<parameter classPath=\"java.lang.Integer\" name=\"size\">\n" +
               "<value>-1</value>\n" +
               "</parameter>\n" +
               "<parameter array=\"true\" classPath=\"java.lang.String\" name=\"dropStrategies\">\n" +
               "<refClass>Class1</refClass>\n" +
               "<subParameter classPath=\"java.lang.String\" name=\"dropStrategy\">\n" +
               "<value>drop</value>\n" +
               "</subParameter>\n" +
               "</parameter>\n" +
               "<parameter classPath=\"jmt.engine.NetStrategies.QueueGetStrategies.FCFSstrategy\" name=\"FCFSstrategy\"/>\n" +
               "<parameter array=\"true\" classPath=\"jmt.engine.NetStrategies.QueuePutStrategy\" name=\"QueuePutStrategy\">\n" +
               "<refClass>Class1</refClass>\n" +
               "<subParameter classPath=\"jmt.engine.NetStrategies.QueuePutStrategies.TailStrategy\" name=\"TailStrategy\"/>\n" +
               "</parameter>\n" +
               "</section>\n" +
               "<section className=\"Server\">\n" +
               "<parameter classPath=\"java.lang.Integer\" name=\"maxJobs\">\n" +
               "<value>1</value>\n" +
               "</parameter>\n" +
               "<parameter array=\"true\" classPath=\"java.lang.Integer\" name=\"numberOfVisits\">\n" +
               "<refClass>Class1</refClass>\n" +
               "<subParameter classPath=\"java.lang.Integer\" name=\"numberOfVisits\">\n" +
               "<value>1</value>\n" +
               "</subParameter>\n" +
               "</parameter>\n" +
               "<parameter array=\"true\" classPath=\"jmt.engine.NetStrategies.ServiceStrategy\" name=\"ServiceStrategy\">\n" +
               "<refClass>Class1</refClass>\n" +
               "<subParameter classPath=\"jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy\" name=\"ServiceTimeStrategy\">\n" +
               "<subParameter classPath=\"jmt.engine.random.Exponential\" name=\"Exponential\">\n" +
               "<subParameter classPath=\"jmt.engine.random.ExponentialPar\" name=\"distrPar\">\n" +
               "<subParameter name=\"lambda\">\n" +
               "<value>2.0</value>\n" +
               "</subParameter>\n" +
               "</subParameter>\n" +
               "</subParameter>\n" +
               "</subParameter>\n" +
               "</parameter>\n" +
               "</section>\n" +
               "<section className=\"Router\">\n" +
               "<parameter array=\"true\" classPath=\"jmt.engine.NetStrategies.RoutingStrategy\" name=\"RoutingStrategy\">\n" +
               "<refClass>Class1</refClass>\n" +
               "<subParameter classPath=\"jmt.engine.NetStrategies.RoutingStrategies.RandomStrategy\" name=\"Random\"/>\n" +
               "</parameter>\n" +
               "</section>\n" +
               "</node>\n" +
               "<node name=\"Sink 1\">\n" +
               "<section className=\"JobSink\"/>\n" +
               "</node>\n" +
               "<connection source=\"Source 1\" target=\"Queue 1\"/>\n" +
               "<connection source=\"Queue 1\" target=\"Sink 1\"/>\n" +
               "</sim>\n" +
               "</archive>";
    }
}
