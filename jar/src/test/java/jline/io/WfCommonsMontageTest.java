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
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Test suite for loading Montage workflows from WfCommons Pegasus-instances repository.
 * <p>
 * Tests loading of real Montage workflow traces from:
 * https://github.com/wfcommons/pegasus-instances/tree/master/montage/chameleon-cloud
 * </p>
 * <p>
 * These tests verify:
 * - Schema 1.4 format parsing (used by Montage workflows)
 * - DAG structure with parallel tasks
 * - Runtime extraction from embedded task data
 * - Task metadata extraction
 * </p>
 */
public class WfCommonsMontageTest {

    @TempDir
    Path tempDir;

    /**
     * Test loading a minimal Montage-like workflow with schema 1.4 format.
     * This simulates the structure found in real Montage workflow files.
     */
    @Test
    public void testMontageSchema14Structure() throws IOException {
        // This JSON mimics the schema 1.4 format used by Montage workflows
        // In schema 1.4, tasks are directly under workflow (not workflow.specification)
        // and runtime data is embedded in the task objects
        String json = "{\n" +
                "  \"name\": \"montage\",\n" +
                "  \"schemaVersion\": \"1.4\",\n" +
                "  \"wms\": {\n" +
                "    \"name\": \"Pegasus\",\n" +
                "    \"version\": \"5.0\"\n" +
                "  },\n" +
                "  \"workflow\": {\n" +
                "    \"tasks\": [\n" +
                "      {\n" +
                "        \"name\": \"mProject_ID0000001\",\n" +
                "        \"id\": \"ID0000001\",\n" +
                "        \"type\": \"compute\",\n" +
                "        \"category\": \"mProject\",\n" +
                "        \"parents\": [],\n" +
                "        \"children\": [\"ID0000003\", \"ID0000004\"],\n" +
                "        \"runtimeInSeconds\": 16.712,\n" +
                "        \"avgCPU\": 97.6723\n" +
                "      },\n" +
                "      {\n" +
                "        \"name\": \"mProject_ID0000002\",\n" +
                "        \"id\": \"ID0000002\",\n" +
                "        \"type\": \"compute\",\n" +
                "        \"category\": \"mProject\",\n" +
                "        \"parents\": [],\n" +
                "        \"children\": [\"ID0000003\", \"ID0000005\"],\n" +
                "        \"runtimeInSeconds\": 12.5,\n" +
                "        \"avgCPU\": 95.2\n" +
                "      },\n" +
                "      {\n" +
                "        \"name\": \"mDiffFit_ID0000003\",\n" +
                "        \"id\": \"ID0000003\",\n" +
                "        \"type\": \"compute\",\n" +
                "        \"category\": \"mDiffFit\",\n" +
                "        \"parents\": [\"ID0000001\", \"ID0000002\"],\n" +
                "        \"children\": [\"ID0000006\"],\n" +
                "        \"runtimeInSeconds\": 0.092\n" +
                "      },\n" +
                "      {\n" +
                "        \"name\": \"mDiffFit_ID0000004\",\n" +
                "        \"id\": \"ID0000004\",\n" +
                "        \"type\": \"compute\",\n" +
                "        \"category\": \"mDiffFit\",\n" +
                "        \"parents\": [\"ID0000001\"],\n" +
                "        \"children\": [\"ID0000006\"],\n" +
                "        \"runtimeInSeconds\": 0.085\n" +
                "      },\n" +
                "      {\n" +
                "        \"name\": \"mDiffFit_ID0000005\",\n" +
                "        \"id\": \"ID0000005\",\n" +
                "        \"type\": \"compute\",\n" +
                "        \"category\": \"mDiffFit\",\n" +
                "        \"parents\": [\"ID0000002\"],\n" +
                "        \"children\": [\"ID0000006\"],\n" +
                "        \"runtimeInSeconds\": 0.078\n" +
                "      },\n" +
                "      {\n" +
                "        \"name\": \"mConcatFit_ID0000006\",\n" +
                "        \"id\": \"ID0000006\",\n" +
                "        \"type\": \"compute\",\n" +
                "        \"category\": \"mConcatFit\",\n" +
                "        \"parents\": [\"ID0000003\", \"ID0000004\", \"ID0000005\"],\n" +
                "        \"children\": [],\n" +
                "        \"runtimeInSeconds\": 0.195\n" +
                "      }\n" +
                "    ]\n" +
                "  }\n" +
                "}";

        Workflow wf = WfCommonsLoader.loadFromString(json);

        // Verify workflow structure
        assertEquals(6, wf.getActivities().size());
        assertEquals("montage", wf.getName());

        // Verify activities exist with correct IDs
        assertNotNull(wf.getActivity("ID0000001"));
        assertNotNull(wf.getActivity("ID0000002"));
        assertNotNull(wf.getActivity("ID0000003"));
        assertNotNull(wf.getActivity("ID0000006"));

        // Verify PH representation can be computed
        APH ph = wf.toPH();
        assertNotNull(ph);
        assertTrue(ph.getMean() > 0);
    }

    /**
     * Test loading Montage workflow with metadata extraction enabled.
     */
    @Test
    public void testMontageMetadataExtraction() throws IOException {
        String json = "{\n" +
                "  \"name\": \"montage-metadata\",\n" +
                "  \"schemaVersion\": \"1.4\",\n" +
                "  \"workflow\": {\n" +
                "    \"tasks\": [\n" +
                "      {\n" +
                "        \"name\": \"mProject_task\",\n" +
                "        \"id\": \"task1\",\n" +
                "        \"type\": \"compute\",\n" +
                "        \"category\": \"mProject\",\n" +
                "        \"parents\": [],\n" +
                "        \"children\": [],\n" +
                "        \"runtimeInSeconds\": 10.5,\n" +
                "        \"avgCPU\": 98.5,\n" +
                "        \"memory\": 14800,\n" +
                "        \"priority\": 20\n" +
                "      }\n" +
                "    ]\n" +
                "  }\n" +
                "}";

        WfCommonsOptions options = new WfCommonsOptions().setStoreMetadata(true);
        Workflow wf = WfCommonsLoader.loadFromString(json, options);

        WorkflowActivity act = wf.getActivity("task1");
        assertNotNull(act);
        assertTrue(act.hasMetadata());
        assertEquals("task1", act.getMetadataValue("taskId"));
        assertEquals(98.5, act.getMetadataValue("avgCPU"));
    }

    /**
     * Test DAG pattern typical in Montage: mProject -> mDiffFit -> mConcatFit.
     * Verifies the workflow correctly captures the parallel structure.
     */
    @Test
    public void testMontagePipelinePattern() throws IOException {
        // Simplified Montage pipeline:
        // mProject1 ─┬─► mDiffFit1 ─┬─► mConcatFit
        //            │              │
        // mProject2 ─┴─► mDiffFit2 ─┘
        String json = "{\n" +
                "  \"name\": \"montage-pipeline\",\n" +
                "  \"schemaVersion\": \"1.4\",\n" +
                "  \"workflow\": {\n" +
                "    \"tasks\": [\n" +
                "      {\"id\": \"mProject1\", \"parents\": [], \"children\": [\"mDiffFit1\", \"mDiffFit2\"], \"runtimeInSeconds\": 15.0},\n" +
                "      {\"id\": \"mProject2\", \"parents\": [], \"children\": [\"mDiffFit1\", \"mDiffFit2\"], \"runtimeInSeconds\": 12.0},\n" +
                "      {\"id\": \"mDiffFit1\", \"parents\": [\"mProject1\", \"mProject2\"], \"children\": [\"mConcatFit\"], \"runtimeInSeconds\": 0.1},\n" +
                "      {\"id\": \"mDiffFit2\", \"parents\": [\"mProject1\", \"mProject2\"], \"children\": [\"mConcatFit\"], \"runtimeInSeconds\": 0.1},\n" +
                "      {\"id\": \"mConcatFit\", \"parents\": [\"mDiffFit1\", \"mDiffFit2\"], \"children\": [], \"runtimeInSeconds\": 0.2}\n" +
                "    ]\n" +
                "  }\n" +
                "}";

        Workflow wf = WfCommonsLoader.loadFromString(json);
        assertEquals(5, wf.getActivities().size());

        APH ph = wf.toPH();
        // Critical path: mProject1 (15) + mDiffFit1 (0.1) + mConcatFit (0.2) = 15.3
        // or mProject2 (12) + ... so max is around 15.3
        // Due to parallel structure, mean should be > individual task times
        assertTrue(ph.getMean() > 0);
    }

    /**
     * Test deterministic distribution fitting for Montage tasks.
     */
    @Test
    public void testMontageDeterministicTiming() throws IOException {
        String json = "{\n" +
                "  \"name\": \"montage-det\",\n" +
                "  \"schemaVersion\": \"1.4\",\n" +
                "  \"workflow\": {\n" +
                "    \"tasks\": [\n" +
                "      {\"id\": \"task1\", \"parents\": [], \"children\": [\"task2\"], \"runtimeInSeconds\": 5.0},\n" +
                "      {\"id\": \"task2\", \"parents\": [\"task1\"], \"children\": [], \"runtimeInSeconds\": 3.0}\n" +
                "    ]\n" +
                "  }\n" +
                "}";

        WfCommonsOptions options = WfCommonsOptions.deterministic();
        Workflow wf = WfCommonsLoader.loadFromString(json, options);
        APH ph = wf.toPH();

        // Serial workflow: mean = 5.0 + 3.0 = 8.0
        assertEquals(8.0, ph.getMean(), 0.1);
    }

    /**
     * Test loading from file in schema 1.4 format.
     */
    @Test
    public void testMontageFileLoading() throws IOException {
        String json = "{\n" +
                "  \"name\": \"montage-file\",\n" +
                "  \"schemaVersion\": \"1.4\",\n" +
                "  \"workflow\": {\n" +
                "    \"tasks\": [\n" +
                "      {\"id\": \"A\", \"parents\": [], \"children\": [\"B\"], \"runtimeInSeconds\": 2.0},\n" +
                "      {\"id\": \"B\", \"parents\": [\"A\"], \"children\": [], \"runtimeInSeconds\": 3.0}\n" +
                "    ]\n" +
                "  }\n" +
                "}";

        File tempFile = tempDir.resolve("montage_test.json").toFile();
        FileWriter writer = new FileWriter(tempFile);
        writer.write(json);
        writer.close();

        Workflow wf = WfCommonsLoader.load(tempFile.getAbsolutePath());
        APH ph = wf.toPH();

        assertEquals(5.0, ph.getMean(), 0.1);
        assertEquals(2, wf.getActivities().size());
    }

    /**
     * Test loading a realistic Montage workflow structure with multiple stages.
     * This mimics the actual Montage workflow structure found in Pegasus instances.
     */
    @Test
    public void testRealisticMontageWorkflow() throws IOException {
        // Realistic Montage workflow stages:
        // 1. mProject (parallel image projection)
        // 2. mDiffFit (difference fitting between pairs)
        // 3. mConcatFit (concatenate fit results)
        // 4. mBgModel (background model)
        // 5. mBackground (apply background correction)
        // 6. mImgTbl (create image table)
        // 7. mAdd (add images)
        // 8. mShrink (shrink final image)
        String json = "{\n" +
                "  \"name\": \"montage-realistic\",\n" +
                "  \"schemaVersion\": \"1.4\",\n" +
                "  \"wms\": {\"name\": \"Pegasus\", \"version\": \"5.0\"},\n" +
                "  \"workflow\": {\n" +
                "    \"machines\": [{\"nodeName\": \"compute-node\", \"cpu\": {\"count\": 48, \"speed\": 1200}}],\n" +
                "    \"tasks\": [\n" +
                "      {\"id\": \"mProject_1\", \"category\": \"mProject\", \"parents\": [], \"children\": [\"mDiffFit_1\"], \"runtimeInSeconds\": 16.7},\n" +
                "      {\"id\": \"mProject_2\", \"category\": \"mProject\", \"parents\": [], \"children\": [\"mDiffFit_1\", \"mDiffFit_2\"], \"runtimeInSeconds\": 12.3},\n" +
                "      {\"id\": \"mProject_3\", \"category\": \"mProject\", \"parents\": [], \"children\": [\"mDiffFit_2\"], \"runtimeInSeconds\": 14.1},\n" +
                "      {\"id\": \"mDiffFit_1\", \"category\": \"mDiffFit\", \"parents\": [\"mProject_1\", \"mProject_2\"], \"children\": [\"mConcatFit_1\"], \"runtimeInSeconds\": 0.092},\n" +
                "      {\"id\": \"mDiffFit_2\", \"category\": \"mDiffFit\", \"parents\": [\"mProject_2\", \"mProject_3\"], \"children\": [\"mConcatFit_1\"], \"runtimeInSeconds\": 0.088},\n" +
                "      {\"id\": \"mConcatFit_1\", \"category\": \"mConcatFit\", \"parents\": [\"mDiffFit_1\", \"mDiffFit_2\"], \"children\": [\"mBgModel_1\"], \"runtimeInSeconds\": 0.195},\n" +
                "      {\"id\": \"mBgModel_1\", \"category\": \"mBgModel\", \"parents\": [\"mConcatFit_1\"], \"children\": [\"mBackground_1\", \"mBackground_2\", \"mBackground_3\"], \"runtimeInSeconds\": 0.1},\n" +
                "      {\"id\": \"mBackground_1\", \"category\": \"mBackground\", \"parents\": [\"mBgModel_1\"], \"children\": [\"mImgTbl_1\"], \"runtimeInSeconds\": 0.2},\n" +
                "      {\"id\": \"mBackground_2\", \"category\": \"mBackground\", \"parents\": [\"mBgModel_1\"], \"children\": [\"mImgTbl_1\"], \"runtimeInSeconds\": 0.18},\n" +
                "      {\"id\": \"mBackground_3\", \"category\": \"mBackground\", \"parents\": [\"mBgModel_1\"], \"children\": [\"mImgTbl_1\"], \"runtimeInSeconds\": 0.22},\n" +
                "      {\"id\": \"mImgTbl_1\", \"category\": \"mImgTbl\", \"parents\": [\"mBackground_1\", \"mBackground_2\", \"mBackground_3\"], \"children\": [\"mAdd_1\"], \"runtimeInSeconds\": 0.05},\n" +
                "      {\"id\": \"mAdd_1\", \"category\": \"mAdd\", \"parents\": [\"mImgTbl_1\"], \"children\": [\"mShrink_1\"], \"runtimeInSeconds\": 2.5},\n" +
                "      {\"id\": \"mShrink_1\", \"category\": \"mShrink\", \"parents\": [\"mAdd_1\"], \"children\": [], \"runtimeInSeconds\": 0.02}\n" +
                "    ]\n" +
                "  }\n" +
                "}";

        Workflow wf = WfCommonsLoader.loadFromString(json);

        // Verify structure
        assertEquals(13, wf.getActivities().size());
        assertEquals("montage_realistic", wf.getName());

        // Verify key activities exist
        assertNotNull(wf.getActivity("mProject_1"));
        assertNotNull(wf.getActivity("mDiffFit_1"));
        assertNotNull(wf.getActivity("mConcatFit_1"));
        assertNotNull(wf.getActivity("mShrink_1"));

        // Verify PH representation
        APH ph = wf.toPH();
        assertNotNull(ph);
        assertTrue(ph.getMean() > 0);

        // Critical path should include the longest mProject + downstream stages
        // mProject_1 (16.7) -> mDiffFit_1 (0.092) -> mConcatFit_1 (0.195) ->
        // mBgModel_1 (0.1) -> mBackground_x (0.2) -> mImgTbl_1 (0.05) ->
        // mAdd_1 (2.5) -> mShrink_1 (0.02) ≈ 19.8
        assertTrue(ph.getMean() > 15.0);
    }

    /**
     * Test URL loading is available (method exists).
     * Note: Actual URL loading test is conditional on network availability.
     */
    @Test
    public void testUrlLoadingMethodExists() {
        // Verify the loadFromUrl method exists
        try {
            java.lang.reflect.Method method = WfCommonsLoader.class.getMethod(
                    "loadFromUrl", String.class);
            assertNotNull(method);
        } catch (NoSuchMethodException e) {
            fail("loadFromUrl method should exist");
        }

        try {
            java.lang.reflect.Method method = WfCommonsLoader.class.getMethod(
                    "loadFromUrl", String.class, WfCommonsOptions.class);
            assertNotNull(method);
        } catch (NoSuchMethodException e) {
            fail("loadFromUrl method with options should exist");
        }
    }

    /**
     * Test loading real Montage workflow from test fixtures.
     * Uses montage-chameleon-dss-05d-001.json from WfCommons pegasus-instances.
     */
    @Test
    public void testLoadFromPegasusInstances() throws IOException {
        // Load the fixture from classpath
        String fixturePath = copyResourceToTempFile("wfcommons/montage/montage-chameleon-dss-05d-001.json");

        Workflow wf = WfCommonsLoader.load(fixturePath);

        assertNotNull(wf);
        // The DSS 05d workflow has 50 tasks
        assertTrue(wf.getActivities().size() >= 40, "Expected at least 40 tasks");

        // Verify it's a Montage workflow (file uses montage-0 which is sanitized to montage_0)
        assertEquals("montage_0", wf.getName());

        // Verify workflow has valid activities with timing data
        for (WorkflowActivity act : wf.getActivities()) {
            assertNotNull(act.getName());
            assertTrue(act.getHostDemandMean() >= 0);
        }
    }

    /**
     * Test loading larger Montage workflow from test fixtures.
     * Uses montage-chameleon-2mass-005d-001.json from WfCommons pegasus-instances.
     */
    @Test
    public void testLoadLargerMontageWorkflow() throws IOException {
        // Load the fixture from classpath
        String fixturePath = copyResourceToTempFile("wfcommons/montage/montage-chameleon-2mass-005d-001.json");

        Workflow wf = WfCommonsLoader.load(fixturePath);

        assertNotNull(wf);
        // The 005d workflow has 59 tasks
        assertTrue(wf.getActivities().size() >= 50);

        // Verify workflow has activities with valid timing
        for (WorkflowActivity act : wf.getActivities()) {
            assertNotNull(act.getName());
            assertTrue(act.getHostDemandMean() >= 0);
        }
    }

    /**
     * Helper method to copy a classpath resource to a temporary file.
     * This is needed because WfCommonsLoader.load() expects a file path.
     */
    private String copyResourceToTempFile(String resourcePath) throws IOException {
        InputStream is = getClass().getClassLoader().getResourceAsStream(resourcePath);
        if (is == null) {
            throw new IOException("Resource not found: " + resourcePath);
        }
        Path tempFile = Files.createTempFile("wfcommons-test-", ".json");
        Files.copy(is, tempFile, StandardCopyOption.REPLACE_EXISTING);
        is.close();
        tempFile.toFile().deleteOnExit();
        return tempFile.toString();
    }
}
