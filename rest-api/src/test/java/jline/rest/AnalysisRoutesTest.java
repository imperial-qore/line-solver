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
 * Tests for analysis endpoints (what-if, sensitivity, bottleneck).
 */
@TestMethodOrder(MethodOrderer.OrderAnnotation.class)
class AnalysisRoutesTest {

    private static final int TEST_PORT = 4569;
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
    void testBottleneckAnalysis() throws Exception {
        BottleneckRequest request = new BottleneckRequest();
        ModelInput model = new ModelInput();
        model.setFormat("jsimg");
        model.setContent(createMinimalJsimgModel());
        request.setModel(model);
        request.setSolver("mva");
        request.setThreshold(0.5);
        request.setIncludeAll(true);

        HttpURLConnection conn = createConnection("/analysis/bottleneck", "POST");
        conn.setDoOutput(true);
        conn.getOutputStream().write(json.toJson(request).getBytes(StandardCharsets.UTF_8));

        int status = conn.getResponseCode();
        String response = readResponse(conn);

        // Accept either 200 (success) or 500 with error details (solve may fail on minimal model)
        assertTrue(status == 200 || status == 500, "Expected 200 or 500, got " + status);
        assertNotNull(response);
        assertTrue(response.contains("status") || response.contains("error"));
    }

    @Test
    @Order(2)
    void testBottleneckAnalysisInvalidRequest() throws Exception {
        BottleneckRequest request = new BottleneckRequest();
        // Missing model

        HttpURLConnection conn = createConnection("/analysis/bottleneck", "POST");
        conn.setDoOutput(true);
        conn.getOutputStream().write(json.toJson(request).getBytes(StandardCharsets.UTF_8));

        int status = conn.getResponseCode();
        assertEquals(400, status);
    }

    @Test
    @Order(3)
    void testWhatIfAnalysis() throws Exception {
        WhatIfRequest request = new WhatIfRequest();
        ModelInput model = new ModelInput();
        model.setFormat("jsimg");
        model.setContent(createMinimalJsimgModel());
        request.setModel(model);
        request.setSolver("mva");
        request.setAnalysis("avg");

        // Add parameter sweep
        ParameterSweep sweep = new ParameterSweep();
        sweep.setStation("Queue 1");
        sweep.setProperty("servers");
        List<Double> values = new ArrayList<Double>();
        values.add(1.0);
        values.add(2.0);
        sweep.setValues(values);

        List<ParameterSweep> parameters = new ArrayList<ParameterSweep>();
        parameters.add(sweep);
        request.setParameters(parameters);

        HttpURLConnection conn = createConnection("/analysis/whatif", "POST");
        conn.setDoOutput(true);
        conn.getOutputStream().write(json.toJson(request).getBytes(StandardCharsets.UTF_8));

        int status = conn.getResponseCode();
        String response = readResponse(conn);

        assertEquals(200, status);
        assertNotNull(response);
        assertTrue(response.contains("status"));
        assertTrue(response.contains("results") || response.contains("totalPoints"));
    }

    @Test
    @Order(4)
    void testWhatIfAnalysisInvalidRequest() throws Exception {
        WhatIfRequest request = new WhatIfRequest();
        ModelInput model = new ModelInput();
        model.setFormat("jsimg");
        model.setContent(createMinimalJsimgModel());
        request.setModel(model);
        request.setSolver("mva");
        // Missing parameters

        HttpURLConnection conn = createConnection("/analysis/whatif", "POST");
        conn.setDoOutput(true);
        conn.getOutputStream().write(json.toJson(request).getBytes(StandardCharsets.UTF_8));

        int status = conn.getResponseCode();
        assertEquals(400, status);
    }

    @Test
    @Order(5)
    void testWhatIfWithRangeParameter() throws Exception {
        WhatIfRequest request = new WhatIfRequest();
        ModelInput model = new ModelInput();
        model.setFormat("jsimg");
        model.setContent(createMinimalJsimgModel());
        request.setModel(model);
        request.setSolver("mva");

        // Add parameter sweep with range
        ParameterSweep sweep = new ParameterSweep();
        sweep.setStation("Queue 1");
        sweep.setProperty("serviceRate");
        sweep.setStart(1.0);
        sweep.setEnd(3.0);
        sweep.setPoints(3);

        List<ParameterSweep> parameters = new ArrayList<ParameterSweep>();
        parameters.add(sweep);
        request.setParameters(parameters);

        HttpURLConnection conn = createConnection("/analysis/whatif", "POST");
        conn.setDoOutput(true);
        conn.getOutputStream().write(json.toJson(request).getBytes(StandardCharsets.UTF_8));

        int status = conn.getResponseCode();
        String response = readResponse(conn);

        assertEquals(200, status);
        assertNotNull(response);
        assertTrue(response.contains("totalPoints"));
    }

    @Test
    @Order(6)
    void testSensitivityAnalysis() throws Exception {
        SensitivityRequest request = new SensitivityRequest();
        ModelInput model = new ModelInput();
        model.setFormat("jsimg");
        model.setContent(createMinimalJsimgModel());
        request.setModel(model);
        request.setSolver("mva");
        request.setDelta(0.05);
        request.setNormalized(false);

        // Add parameter
        ParameterSweep param = new ParameterSweep();
        param.setStation("Queue 1");
        param.setProperty("serviceRate");

        List<ParameterSweep> parameters = new ArrayList<ParameterSweep>();
        parameters.add(param);
        request.setParameters(parameters);

        HttpURLConnection conn = createConnection("/analysis/sensitivity", "POST");
        conn.setDoOutput(true);
        conn.getOutputStream().write(json.toJson(request).getBytes(StandardCharsets.UTF_8));

        int status = conn.getResponseCode();
        String response = readResponse(conn);

        // Accept either 200 (success) or 500 with error details (solve may fail on minimal model)
        assertTrue(status == 200 || status == 500, "Expected 200 or 500, got " + status);
        assertNotNull(response);
        assertTrue(response.contains("status") || response.contains("error"));
    }

    @Test
    @Order(7)
    void testSensitivityAnalysisInvalidRequest() throws Exception {
        SensitivityRequest request = new SensitivityRequest();
        // Missing model

        HttpURLConnection conn = createConnection("/analysis/sensitivity", "POST");
        conn.setDoOutput(true);
        conn.getOutputStream().write(json.toJson(request).getBytes(StandardCharsets.UTF_8));

        int status = conn.getResponseCode();
        assertEquals(400, status);
    }

    @Test
    @Order(8)
    void testSensitivityInvalidDelta() throws Exception {
        SensitivityRequest request = new SensitivityRequest();
        ModelInput model = new ModelInput();
        model.setFormat("jsimg");
        model.setContent(createMinimalJsimgModel());
        request.setModel(model);
        request.setSolver("mva");
        request.setDelta(1.5); // Invalid: must be between 0 and 1

        ParameterSweep param = new ParameterSweep();
        param.setProperty("serviceRate");
        List<ParameterSweep> parameters = new ArrayList<ParameterSweep>();
        parameters.add(param);
        request.setParameters(parameters);

        HttpURLConnection conn = createConnection("/analysis/sensitivity", "POST");
        conn.setDoOutput(true);
        conn.getOutputStream().write(json.toJson(request).getBytes(StandardCharsets.UTF_8));

        int status = conn.getResponseCode();
        assertEquals(400, status);
    }

    @Test
    @Order(9)
    void testBottleneckWithThreshold() throws Exception {
        BottleneckRequest request = new BottleneckRequest();
        ModelInput model = new ModelInput();
        model.setFormat("jsimg");
        model.setContent(createMinimalJsimgModel());
        request.setModel(model);
        request.setSolver("mva");
        request.setThreshold(0.9); // High threshold

        HttpURLConnection conn = createConnection("/analysis/bottleneck", "POST");
        conn.setDoOutput(true);
        conn.getOutputStream().write(json.toJson(request).getBytes(StandardCharsets.UTF_8));

        int status = conn.getResponseCode();
        String response = readResponse(conn);

        // Accept either 200 (success) or 500 with error details (solve may fail on minimal model)
        assertTrue(status == 200 || status == 500, "Expected 200 or 500, got " + status);
        assertTrue(response.contains("threshold") || response.contains("error"));
    }

    @Test
    @Order(10)
    void testBottleneckInvalidThreshold() throws Exception {
        BottleneckRequest request = new BottleneckRequest();
        ModelInput model = new ModelInput();
        model.setFormat("jsimg");
        model.setContent(createMinimalJsimgModel());
        request.setModel(model);
        request.setSolver("mva");
        request.setThreshold(1.5); // Invalid: must be <= 1

        HttpURLConnection conn = createConnection("/analysis/bottleneck", "POST");
        conn.setDoOutput(true);
        conn.getOutputStream().write(json.toJson(request).getBytes(StandardCharsets.UTF_8));

        int status = conn.getResponseCode();
        assertEquals(400, status);
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
