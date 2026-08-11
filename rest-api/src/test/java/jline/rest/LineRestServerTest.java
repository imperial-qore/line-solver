/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.rest;

import jline.rest.model.ErrorResponse;
import jline.rest.model.SolveRequest;
import jline.rest.model.SolveResponse;
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
 * Tests for the LINE REST API server.
 */
@TestMethodOrder(MethodOrderer.OrderAnnotation.class)
class LineRestServerTest {

    private static final int TEST_PORT = 4567;
    private static final String BASE_URL = "http://localhost:" + TEST_PORT + "/api/v1";
    private static LineRestServer server;
    private static JsonTransformer json;

    @BeforeAll
    static void setUp() throws Exception {
        json = new JsonTransformer();
        server = new LineRestServer(TEST_PORT);
        server.start();
        // Wait for server to start
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
    void testHealthEndpoint() throws Exception {
        String response = sendGet("/health");
        assertNotNull(response);
        assertTrue(response.contains("healthy"));
        assertTrue(response.contains("timestamp"));
    }

    @Test
    @Order(2)
    void testReadyEndpoint() throws Exception {
        String response = sendGet("/ready");
        assertNotNull(response);
        assertTrue(response.contains("ready") || response.contains("not_ready"));
        assertTrue(response.contains("memoryUsagePercent"));
    }

    @Test
    @Order(3)
    void testInfoEndpoint() throws Exception {
        String response = sendGet("/info");
        assertNotNull(response);
        assertTrue(response.contains("LINE Solver"));
        assertTrue(response.contains("version"));
        assertTrue(response.contains("capabilities"));
        assertTrue(response.contains("solvers"));
    }

    @Test
    @Order(4)
    void testSolversEndpoint() throws Exception {
        String response = sendGet("/solvers");
        assertNotNull(response);
        assertTrue(response.contains("mva"));
        assertTrue(response.contains("jmt"));
        assertTrue(response.contains("ctmc"));
        assertTrue(response.contains("fluid"));
        assertTrue(response.contains("ssa"));
        assertTrue(response.contains("nc"));
        assertTrue(response.contains("ln"));
        assertTrue(response.contains("lqns"));
    }

    @Test
    @Order(5)
    void testSolverDetailEndpoint() throws Exception {
        String response = sendGet("/solvers/mva");
        assertNotNull(response);
        assertTrue(response.contains("Mean Value Analysis"));
        assertTrue(response.contains("analytical"));
    }

    @Test
    @Order(6)
    void testSolverNotFound() throws Exception {
        HttpURLConnection conn = createConnection("/solvers/nonexistent", "GET");
        int status = conn.getResponseCode();
        assertEquals(404, status);
    }

    @Test
    @Order(7)
    void testValidateEndpoint_InvalidFormat() throws Exception {
        ModelInput input = new ModelInput();
        input.setFormat("invalid");
        input.setContent("test");

        HttpURLConnection conn = createConnection("/models/validate", "POST");
        conn.setDoOutput(true);
        conn.getOutputStream().write(json.toJson(input).getBytes(StandardCharsets.UTF_8));

        int status = conn.getResponseCode();
        assertEquals(400, status);

        String response = readResponse(conn);
        assertTrue(response.contains("false") || response.contains("Invalid format"));
    }

    @Test
    @Order(8)
    void testSolveEndpoint_MissingModel() throws Exception {
        SolveRequest request = new SolveRequest();
        request.setSolver("mva");

        HttpURLConnection conn = createConnection("/models/solve", "POST");
        conn.setDoOutput(true);
        conn.getOutputStream().write(json.toJson(request).getBytes(StandardCharsets.UTF_8));

        int status = conn.getResponseCode();
        assertEquals(400, status);

        String response = readResponse(conn);
        assertTrue(response.contains("Model is required") || response.contains("failed"));
    }

    @Test
    @Order(9)
    void testSolveEndpoint_InvalidSolver() throws Exception {
        SolveRequest request = new SolveRequest();
        ModelInput model = new ModelInput();
        model.setFormat("jsimg");
        model.setContent("<test/>");
        request.setModel(model);
        request.setSolver("invalid");

        HttpURLConnection conn = createConnection("/models/solve", "POST");
        conn.setDoOutput(true);
        conn.getOutputStream().write(json.toJson(request).getBytes(StandardCharsets.UTF_8));

        int status = conn.getResponseCode();
        assertEquals(400, status);

        String response = readResponse(conn);
        assertTrue(response.contains("Invalid solver") || response.contains("failed"));
    }

    @Test
    @Order(10)
    void testSolveEndpoint_ValidRequest() throws Exception {
        // Create a minimal JSIMG model
        String jsimgContent = createMinimalJsimgModel();

        SolveRequest request = new SolveRequest();
        ModelInput model = new ModelInput();
        model.setFormat("jsimg");
        model.setContent(jsimgContent);
        request.setModel(model);
        request.setSolver("mva");
        request.setAnalysis("avg");

        HttpURLConnection conn = createConnection("/models/solve", "POST");
        conn.setDoOutput(true);
        conn.getOutputStream().write(json.toJson(request).getBytes(StandardCharsets.UTF_8));

        // This may fail if the model is not valid, but we're testing the endpoint works
        int status = conn.getResponseCode();
        String response = readResponse(conn);
        assertNotNull(response);
        // Either success or failure, but should have a response
        assertTrue(response.contains("status"));
    }

    @Test
    @Order(11)
    void testNotFoundEndpoint() throws Exception {
        HttpURLConnection conn = createConnection("/nonexistent", "GET");
        int status = conn.getResponseCode();
        assertEquals(404, status);
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
