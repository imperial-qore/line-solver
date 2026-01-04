/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.io;

import jline.lang.processes.APH;
import jline.lang.workflow.Workflow;
import jline.lang.workflow.WorkflowActivity;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Test suite for WfCommonsLoader.
 */
public class WfCommonsLoaderTest {

    @TempDir
    Path tempDir;

    @Test
    public void testSimpleLinearWorkflow() throws IOException {
        String json = "{\n" +
                "  \"name\": \"LinearWorkflow\",\n" +
                "  \"schemaVersion\": \"1.5\",\n" +
                "  \"workflow\": {\n" +
                "    \"specification\": {\n" +
                "      \"tasks\": [\n" +
                "        {\"id\": \"task1\", \"parents\": [], \"children\": [\"task2\"]},\n" +
                "        {\"id\": \"task2\", \"parents\": [\"task1\"], \"children\": [\"task3\"]},\n" +
                "        {\"id\": \"task3\", \"parents\": [\"task2\"], \"children\": []}\n" +
                "      ]\n" +
                "    },\n" +
                "    \"execution\": {\n" +
                "      \"tasks\": [\n" +
                "        {\"id\": \"task1\", \"runtimeInSeconds\": 1.0},\n" +
                "        {\"id\": \"task2\", \"runtimeInSeconds\": 2.0},\n" +
                "        {\"id\": \"task3\", \"runtimeInSeconds\": 3.0}\n" +
                "      ]\n" +
                "    }\n" +
                "  }\n" +
                "}";

        Workflow wf = WfCommonsLoader.loadFromString(json);
        APH ph = wf.toPH();

        double expectedMean = 1.0 + 2.0 + 3.0;
        assertEquals(expectedMean, ph.getMean(), 0.1);
        assertEquals(3, wf.getActivities().size());
    }

    @Test
    public void testForkJoinWorkflow() throws IOException {
        String json = "{\n" +
                "  \"name\": \"ForkJoinWorkflow\",\n" +
                "  \"schemaVersion\": \"1.5\",\n" +
                "  \"workflow\": {\n" +
                "    \"specification\": {\n" +
                "      \"tasks\": [\n" +
                "        {\"id\": \"start\", \"parents\": [], \"children\": [\"branch1\", \"branch2\"]},\n" +
                "        {\"id\": \"branch1\", \"parents\": [\"start\"], \"children\": [\"end\"]},\n" +
                "        {\"id\": \"branch2\", \"parents\": [\"start\"], \"children\": [\"end\"]},\n" +
                "        {\"id\": \"end\", \"parents\": [\"branch1\", \"branch2\"], \"children\": []}\n" +
                "      ]\n" +
                "    },\n" +
                "    \"execution\": {\n" +
                "      \"tasks\": [\n" +
                "        {\"id\": \"start\", \"runtimeInSeconds\": 1.0},\n" +
                "        {\"id\": \"branch1\", \"runtimeInSeconds\": 2.0},\n" +
                "        {\"id\": \"branch2\", \"runtimeInSeconds\": 3.0},\n" +
                "        {\"id\": \"end\", \"runtimeInSeconds\": 0.5}\n" +
                "      ]\n" +
                "    }\n" +
                "  }\n" +
                "}";

        Workflow wf = WfCommonsLoader.loadFromString(json);
        APH ph = wf.toPH();

        // Fork-join should have mean > sum of serial parts
        assertTrue(ph.getMean() > 4.0);
        assertEquals(4, wf.getActivities().size());
    }

    @Test
    public void testMissingExecutionData() throws IOException {
        String json = "{\n" +
                "  \"name\": \"NoExecutionData\",\n" +
                "  \"schemaVersion\": \"1.5\",\n" +
                "  \"workflow\": {\n" +
                "    \"specification\": {\n" +
                "      \"tasks\": [\n" +
                "        {\"id\": \"task1\", \"parents\": [], \"children\": []}\n" +
                "      ]\n" +
                "    }\n" +
                "  }\n" +
                "}";

        WfCommonsOptions options = new WfCommonsOptions().setDefaultRuntime(5.0);
        Workflow wf = WfCommonsLoader.loadFromString(json, options);
        APH ph = wf.toPH();

        assertEquals(5.0, ph.getMean(), 0.1);
    }

    @Test
    public void testMetadataStorage() throws IOException {
        String json = "{\n" +
                "  \"name\": \"MetadataTest\",\n" +
                "  \"schemaVersion\": \"1.5\",\n" +
                "  \"workflow\": {\n" +
                "    \"specification\": {\n" +
                "      \"tasks\": [\n" +
                "        {\"id\": \"task1\", \"name\": \"My Task\", \"inputFiles\": [\"file1.txt\"], \"parents\": [], \"children\": []}\n" +
                "      ]\n" +
                "    },\n" +
                "    \"execution\": {\n" +
                "      \"tasks\": [\n" +
                "        {\"id\": \"task1\", \"runtimeInSeconds\": 2.0, \"avgCPU\": 75.5}\n" +
                "      ]\n" +
                "    }\n" +
                "  }\n" +
                "}";

        WfCommonsOptions options = new WfCommonsOptions().setStoreMetadata(true);
        Workflow wf = WfCommonsLoader.loadFromString(json, options);

        WorkflowActivity act = wf.getActivity("task1");
        assertTrue(act.hasMetadata());
        assertEquals("task1", act.getMetadataValue("taskId"));
        assertEquals(75.5, act.getMetadataValue("avgCPU"));
    }

    @Test
    public void testSchemaValidationMissingVersion() {
        String json = "{\n" +
                "  \"name\": \"InvalidWorkflow\",\n" +
                "  \"workflow\": {\n" +
                "    \"specification\": {\n" +
                "      \"tasks\": [{\"id\": \"task1\"}]\n" +
                "    }\n" +
                "  }\n" +
                "}";

        assertThrows(IllegalArgumentException.class, () -> {
            WfCommonsLoader.loadFromString(json);
        });
    }

    @Test
    public void testSchemaValidationMissingWorkflow() {
        String json = "{\n" +
                "  \"name\": \"InvalidWorkflow\",\n" +
                "  \"schemaVersion\": \"1.5\"\n" +
                "}";

        assertThrows(IllegalArgumentException.class, () -> {
            WfCommonsLoader.loadFromString(json);
        });
    }

    @Test
    public void testDeterministicDistribution() throws IOException {
        String json = "{\n" +
                "  \"name\": \"DeterministicTest\",\n" +
                "  \"schemaVersion\": \"1.5\",\n" +
                "  \"workflow\": {\n" +
                "    \"specification\": {\n" +
                "      \"tasks\": [\n" +
                "        {\"id\": \"task1\", \"parents\": [], \"children\": []}\n" +
                "      ]\n" +
                "    },\n" +
                "    \"execution\": {\n" +
                "      \"tasks\": [\n" +
                "        {\"id\": \"task1\", \"runtimeInSeconds\": 3.0}\n" +
                "      ]\n" +
                "    }\n" +
                "  }\n" +
                "}";

        WfCommonsOptions options = WfCommonsOptions.deterministic();
        Workflow wf = WfCommonsLoader.loadFromString(json, options);
        APH ph = wf.toPH();

        assertEquals(3.0, ph.getMean(), 0.1);
    }

    @Test
    public void testFileLoading() throws IOException {
        String json = "{\n" +
                "  \"name\": \"FileTest\",\n" +
                "  \"schemaVersion\": \"1.5\",\n" +
                "  \"workflow\": {\n" +
                "    \"specification\": {\n" +
                "      \"tasks\": [\n" +
                "        {\"id\": \"A\", \"parents\": [], \"children\": [\"B\"]},\n" +
                "        {\"id\": \"B\", \"parents\": [\"A\"], \"children\": []}\n" +
                "      ]\n" +
                "    },\n" +
                "    \"execution\": {\n" +
                "      \"tasks\": [\n" +
                "        {\"id\": \"A\", \"runtimeInSeconds\": 1.0},\n" +
                "        {\"id\": \"B\", \"runtimeInSeconds\": 2.0}\n" +
                "      ]\n" +
                "    }\n" +
                "  }\n" +
                "}";

        File tempFile = tempDir.resolve("test_workflow.json").toFile();
        FileWriter writer = new FileWriter(tempFile);
        writer.write(json);
        writer.close();

        Workflow wf = WfCommonsLoader.load(tempFile.getAbsolutePath());
        APH ph = wf.toPH();

        assertEquals(3.0, ph.getMean(), 0.1);
    }

    @Test
    public void testValidateFile() throws IOException {
        // Valid file
        String validJson = "{\n" +
                "  \"name\": \"ValidWorkflow\",\n" +
                "  \"schemaVersion\": \"1.5\",\n" +
                "  \"workflow\": {\n" +
                "    \"specification\": {\n" +
                "      \"tasks\": [{\"id\": \"task1\", \"parents\": [], \"children\": []}]\n" +
                "    }\n" +
                "  }\n" +
                "}";

        File validFile = tempDir.resolve("valid.json").toFile();
        FileWriter validWriter = new FileWriter(validFile);
        validWriter.write(validJson);
        validWriter.close();

        assertTrue(WfCommonsLoader.validateFile(validFile.getAbsolutePath()));

        // Invalid file
        String invalidJson = "{\"name\": \"invalid\"}";
        File invalidFile = tempDir.resolve("invalid.json").toFile();
        FileWriter invalidWriter = new FileWriter(invalidFile);
        invalidWriter.write(invalidJson);
        invalidWriter.close();

        assertFalse(WfCommonsLoader.validateFile(invalidFile.getAbsolutePath()));
    }

    @Test
    public void testWorkflowFromWfCommons() throws IOException {
        String json = "{\n" +
                "  \"name\": \"StaticMethodTest\",\n" +
                "  \"schemaVersion\": \"1.5\",\n" +
                "  \"workflow\": {\n" +
                "    \"specification\": {\n" +
                "      \"tasks\": [\n" +
                "        {\"id\": \"A\", \"parents\": [], \"children\": [\"B\"]},\n" +
                "        {\"id\": \"B\", \"parents\": [\"A\"], \"children\": []}\n" +
                "      ]\n" +
                "    },\n" +
                "    \"execution\": {\n" +
                "      \"tasks\": [\n" +
                "        {\"id\": \"A\", \"runtimeInSeconds\": 1.0},\n" +
                "        {\"id\": \"B\", \"runtimeInSeconds\": 2.0}\n" +
                "      ]\n" +
                "    }\n" +
                "  }\n" +
                "}";

        File tempFile = tempDir.resolve("static_test.json").toFile();
        FileWriter writer = new FileWriter(tempFile);
        writer.write(json);
        writer.close();

        Workflow wf = Workflow.fromWfCommons(tempFile.getAbsolutePath());
        APH ph = wf.toPH();

        assertEquals(3.0, ph.getMean(), 0.1);
    }

    @Test
    public void testComplexDAG() throws IOException {
        // DAG: task1 -> [task2, task3] -> task4 -> task5
        String json = "{\n" +
                "  \"name\": \"ComplexDAG\",\n" +
                "  \"schemaVersion\": \"1.5\",\n" +
                "  \"workflow\": {\n" +
                "    \"specification\": {\n" +
                "      \"tasks\": [\n" +
                "        {\"id\": \"task1\", \"parents\": [], \"children\": [\"task2\", \"task3\"]},\n" +
                "        {\"id\": \"task2\", \"parents\": [\"task1\"], \"children\": [\"task4\"]},\n" +
                "        {\"id\": \"task3\", \"parents\": [\"task1\"], \"children\": [\"task4\"]},\n" +
                "        {\"id\": \"task4\", \"parents\": [\"task2\", \"task3\"], \"children\": [\"task5\"]},\n" +
                "        {\"id\": \"task5\", \"parents\": [\"task4\"], \"children\": []}\n" +
                "      ]\n" +
                "    },\n" +
                "    \"execution\": {\n" +
                "      \"tasks\": [\n" +
                "        {\"id\": \"task1\", \"runtimeInSeconds\": 1.0},\n" +
                "        {\"id\": \"task2\", \"runtimeInSeconds\": 2.0},\n" +
                "        {\"id\": \"task3\", \"runtimeInSeconds\": 1.5},\n" +
                "        {\"id\": \"task4\", \"runtimeInSeconds\": 0.5},\n" +
                "        {\"id\": \"task5\", \"runtimeInSeconds\": 1.0}\n" +
                "      ]\n" +
                "    }\n" +
                "  }\n" +
                "}";

        Workflow wf = WfCommonsLoader.loadFromString(json);
        APH ph = wf.toPH();

        assertEquals(5, wf.getActivities().size());
        assertTrue(ph.getMean() > 0);
    }
}
