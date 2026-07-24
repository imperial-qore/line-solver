/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.examples.java.GettingStarted;
import jline.lang.Environment;
import jline.lang.Network;
import jline.lang.layered.LayeredNetwork;
import jline.solvers.NetworkAvgTable;
import jline.solvers.mva.SolverMVA;
import jline.solvers.nc.SolverNC;
import jline.util.Maths;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Test;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;

import static org.junit.jupiter.api.Assertions.*;
import static jline.TestTools.withSuppressedOutput;

/**
 * Test class for Getting Started tutorial examples.
 *
 * Tests that all tutorial methods execute without errors and produce valid models.
 */
public class GettingStartedTest {

    private static VerboseLevel originalVerboseLevel;
    private static PrintStream originalOut;
    private static PrintStream originalErr;
    private static ByteArrayOutputStream suppressedOut;
    private static ByteArrayOutputStream suppressedErr;

    @BeforeAll
    public static void setUpClass() {
        // Save original verbose level and set to SILENT for tests
        originalVerboseLevel = GlobalConstants.getVerbose();
        GlobalConstants.setVerbose(VerboseLevel.SILENT);

        // Redirect System.out and System.err to suppress output
        originalOut = System.out;
        originalErr = System.err;
        suppressedOut = new ByteArrayOutputStream();
        suppressedErr = new ByteArrayOutputStream();
        System.setOut(new PrintStream(suppressedOut));
        System.setErr(new PrintStream(suppressedErr));

        // Ensure MATLAB-compatible random number generation
        Maths.setRandomNumbersMatlab(true);
    }

    @AfterAll
    public static void tearDownClass() {
        // Restore original verbose level
        GlobalConstants.setVerbose(originalVerboseLevel);

        // Restore System.out and System.err
        System.setOut(originalOut);
        System.setErr(originalErr);
    }

    @Test
    public void testTut01Mm1Basics() {
        withSuppressedOutput(() -> {
            Network model = GettingStarted.tut01_mm1_basics();
            assertNotNull(model, "tut01 should return a valid Network");
            assertEquals("M/M/1", model.getName(), "Model name should be 'M/M/1'");
            assertEquals(3, model.getNumberOfNodes(), "Model should have 3 nodes (Source, Queue, Sink)");
            assertEquals(1, model.getNumberOfClasses(), "Model should have 1 job class");
        });
    }

    @Test
    public void testTut02Mg1MulticlassSolvers() {
        withSuppressedOutput(() -> {
            Network model = GettingStarted.tut02_mg1_multiclass_solvers();
            assertNotNull(model, "tut02 should return a valid Network");
            assertEquals("M/G/1", model.getName(), "Model name should be 'M/G/1'");
            assertEquals(3, model.getNumberOfNodes(), "Model should have 3 nodes (Source, Queue, Sink)");
            assertEquals(2, model.getNumberOfClasses(), "Model should have 2 job classes");
        });
    }

    @Test
    public void testTut03Repairmen() {
        withSuppressedOutput(() -> {
            Network model = GettingStarted.tut03_repairmen();
            assertNotNull(model, "tut03 should return a valid Network");
            assertEquals("MRP", model.getName(), "Model name should be 'MRP'");
            assertEquals(2, model.getNumberOfNodes(), "Model should have 2 nodes (Delay, Queue)");
            assertEquals(1, model.getNumberOfClasses(), "Model should have 1 job class");
        });
    }

    @Test
    public void testTut04LbRouting() {
        withSuppressedOutput(() -> {
            Network model = GettingStarted.tut04_lb_routing();
            assertNotNull(model, "tut04 should return a valid Network");
            assertEquals("RRLB", model.getName(), "Model name should be 'RRLB'");
            assertEquals(5, model.getNumberOfNodes(), "Model should have 5 nodes (Source, LB, Queue1, Queue2, Sink)");
            assertEquals(1, model.getNumberOfClasses(), "Model should have 1 job class");
        });
    }

    @Test
    public void testTut05CompletesFlag() {
        withSuppressedOutput(() -> {
            Network model = GettingStarted.tut05_completes_flag();
            assertNotNull(model, "tut05 should return a valid Network");
            assertEquals("RL", model.getName(), "Model name should be 'RL'");
            assertEquals(2, model.getNumberOfNodes(), "Model should have 2 nodes (Queue + auto-generated ClassSwitch)");
            assertEquals(3, model.getNumberOfClasses(), "Model should have 3 job classes");
        });
    }

    @Test
    public void testTut06CacheLruZipf() {
        withSuppressedOutput(() -> {
            Network model = GettingStarted.tut06_cache_lru_zipf();
            assertNotNull(model, "tut06 should return a valid Network");
            assertEquals("Model", model.getName(), "Model name should be 'Model'");
            assertTrue(model.getNumberOfNodes() >= 3, "Model should have at least 3 nodes");
            assertEquals(3, model.getNumberOfClasses(), "Model should have 3 job classes");
        });
    }

    @Test
    public void testTut07ResptCdf() {
        withSuppressedOutput(() -> {
            Network model = GettingStarted.tut07_respt_cdf();
            assertNotNull(model, "tut07 should return a valid Network");
            assertEquals("Model", model.getName(), "Model name should be 'Model'");
            assertEquals(2, model.getNumberOfNodes(), "Model should have 2 nodes (Delay, Queue)");
            assertEquals(1, model.getNumberOfClasses(), "Model should have 1 job class");
        });
    }

    @Test
    public void testTut08OptLoadBalancing() {
        withSuppressedOutput(() -> {
            Network model = GettingStarted.tut08_opt_load_balancing();
            assertNotNull(model, "tut08 should return a valid Network");
            assertEquals("LoadBalCQN", model.getName(), "Model name should be 'LoadBalCQN'");
            assertEquals(3, model.getNumberOfNodes(), "Model should have 3 nodes (Delay, Queue1, Queue2)");
            assertEquals(1, model.getNumberOfClasses(), "Model should have 1 job class");
        });
    }

    @Test
    public void testTut09DepProcessAnalysis() {
        withSuppressedOutput(() -> {
            // tut09 is a void method that runs departure process analysis
            // Just verify it executes without throwing exceptions
            assertDoesNotThrow(() -> GettingStarted.tut09_dep_process_analysis(),
                "tut09 departure process analysis should execute without errors");
        });
    }

    @Test
    public void testTut10LqnBasics() {
        withSuppressedOutput(() -> {
            LayeredNetwork model = GettingStarted.tut10_lqn_basics();
            assertNotNull(model, "tut10 should return a valid LayeredNetwork");
            assertEquals("ClientDBSystem", model.getName(), "Model name should be 'ClientDBSystem'");
            // Verify LQN structure (hosts are processors in LQN terminology)
            assertNotNull(model.getHosts(), "LQN should have hosts (processors)");
            assertFalse(model.getHosts().isEmpty(), "LQN should have at least one host (processor)");
            assertNotNull(model.getTasks(), "LQN should have tasks");
            assertFalse(model.getTasks().isEmpty(), "LQN should have at least one task");
        });
    }

    @Test
    public void testTut11RandomEnv() throws Exception {
        withSuppressedOutput(() -> {
            try {
                Environment env = GettingStarted.tut11_random_env();
                assertNotNull(env, "tut11 should return a valid Environment");
                assertEquals("ServerModes", env.getName(), "Environment name should be 'ServerModes'");
                assertEquals(2, env.getNumberOfStages(), "Environment should have 2 stages (Fast, Slow)");

                // Verify each stage has a valid model
                for (int i = 0; i < env.getNumberOfStages(); i++) {
                    Network stageModel = env.getModel(i);
                    assertNotNull(stageModel, "Stage " + i + " should have a valid model");
                }
            } catch (Exception e) {
                throw new RuntimeException(e);
            }
        });
    }

    @Test
    public void testGalleryMerl1() {
        withSuppressedOutput(() -> {
            Object[] components = GettingStarted.gallery_merl1();
            assertNotNull(components, "gallery_merl1 should return components array");
            assertEquals(5, components.length, "Should return 5 components [model, source, queue, sink, oclass]");

            Network model = (Network) components[0];
            assertNotNull(model, "Model should not be null");
            assertEquals("M/E/1", model.getName(), "Model name should be 'M/E/1'");
        });
    }
}
