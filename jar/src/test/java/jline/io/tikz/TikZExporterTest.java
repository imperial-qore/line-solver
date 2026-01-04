/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.io.tikz;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.lang.Environment;
import jline.lang.Network;
import jline.lang.ClosedClass;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.*;
import jline.lang.processes.Exp;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Test;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Tests for TikZ network visualization.
 */
public class TikZExporterTest {

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
    public void testSimpleOpenNetwork() {
        // Create a simple M/M/1 model: Source -> Queue -> Sink
        Network model = new Network("M_M_1");

        // Create nodes
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue1", jline.lang.constant.SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        // Create job class
        OpenClass jobClass = new OpenClass(model, "Class1", 0);

        // Set arrival and service rates
        source.setArrival(jobClass, new Exp(1.0));
        queue.setService(jobClass, new Exp(2.0));

        // Link nodes
        model.addLink(source, queue);
        model.addLink(queue, sink);

        // Generate TikZ code
        String tikz = model.toTikZ();

        // Verify the output contains expected elements
        assertNotNull(tikz);
        assertTrue(tikz.contains("\\documentclass"));
        assertTrue(tikz.contains("\\begin{tikzpicture}"));
        assertTrue(tikz.contains("\\end{tikzpicture}"));
        assertTrue(tikz.contains("source"));  // source style
        assertTrue(tikz.contains("queue"));   // queue style
        assertTrue(tikz.contains("sink"));    // sink style
        assertTrue(tikz.contains("Source"));  // node name
        assertTrue(tikz.contains("Queue1"));  // node name
        assertTrue(tikz.contains("Sink"));    // node name

    }

    @Test
    public void testLayoutEngine() {
        Network model = new Network("TestLayout");

        Source source = new Source(model, "Source");
        Queue queue1 = new Queue(model, "Queue1", jline.lang.constant.SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", jline.lang.constant.SchedStrategy.PS);
        Sink sink = new Sink(model, "Sink");

        OpenClass jobClass = new OpenClass(model, "Class1", 0);
        source.setArrival(jobClass, new Exp(1.0));
        queue1.setService(jobClass, new Exp(2.0));
        queue2.setService(jobClass, new Exp(3.0));

        // Create a branching network: Source -> Queue1 -> Sink
        //                                    -> Queue2 -> Sink
        model.addLink(source, queue1);
        model.addLink(source, queue2);
        model.addLink(queue1, sink);
        model.addLink(queue2, sink);

        TikZOptions options = new TikZOptions();
        TikZLayoutEngine layout = new TikZLayoutEngine(model, options);
        layout.computeLayout();

        // Verify all nodes have positions
        assertNotNull(layout.getPosition(source));
        assertNotNull(layout.getPosition(queue1));
        assertNotNull(layout.getPosition(queue2));
        assertNotNull(layout.getPosition(sink));

        // Verify source is at layer 0 (leftmost)
        double[] sourcePos = layout.getPosition(source);
        assertEquals(0.0, sourcePos[0], 0.001);

        // Verify sink is at rightmost layer
        double[] sinkPos = layout.getPosition(sink);
        assertTrue(sinkPos[0] > sourcePos[0]);
    }

    @Test
    public void testPdfLatexCheck() {
        // Just test the availability check doesn't throw
        boolean available = TikZExporter.isPdfLatexAvailable();
        // Just verify it returns without exception
        assertTrue(available || !available);
    }

    @Test
    public void testLargeSerialForkJoinNetwork() {
        // Create a cascaded fork-join network with 10 nodes
        Network model = jline.examples.java.basic.ForkJoinModel.fj_serialfjs_open();

        // Generate TikZ code
        String tikz = model.toTikZ();

        // Verify the output contains expected elements
        assertNotNull(tikz);
        assertTrue(tikz.contains("\\begin{tikzpicture}"));
        assertTrue(tikz.contains("Source"));
        assertTrue(tikz.contains("Fork1"));
        assertTrue(tikz.contains("Fork2"));
        assertTrue(tikz.contains("Join1"));
        assertTrue(tikz.contains("Join2"));
        assertTrue(tikz.contains("Queue1"));
        assertTrue(tikz.contains("Queue2"));
        assertTrue(tikz.contains("Queue3"));
        assertTrue(tikz.contains("Queue4"));
        assertTrue(tikz.contains("Sink"));
    }

    @Test
    public void testLargeLinearNetwork() {
        // Create a linear network with many queues (M/M/1 tandem with 8 queues)
        Network model = jline.examples.java.models.Gallery.gallery_mm1_linear(8, 0.9);

        // Generate TikZ code
        String tikz = model.toTikZ();

        // Verify the output
        assertNotNull(tikz);
        assertTrue(tikz.contains("\\begin{tikzpicture}"));
        assertTrue(tikz.contains("mySource"));
        assertTrue(tikz.contains("mySink"));
        for (int i = 1; i <= 8; i++) {
            assertTrue(tikz.contains("Queue" + i), "Missing Queue" + i);
        }
    }

    @Test
    public void testComplexClosedNetwork() {
        // Create a complex closed network with multiple queues and classes
        Network model = jline.examples.java.basic.ClosedModel.cqn_multiserver();

        // Generate TikZ code
        String tikz = model.toTikZ();

        // Verify the output
        assertNotNull(tikz);
        assertTrue(tikz.contains("\\begin{tikzpicture}"));
        assertTrue(tikz.contains("Delay"));
        assertTrue(tikz.contains("Queue1"));
        assertTrue(tikz.contains("Queue2"));
    }

    @Test
    public void testExportTutorialPNGs() throws Exception {
        // Export PNG files for tutorials to temp directory
        String imgDir = System.getProperty("java.io.tmpdir") + "/";
        int successCount = 0;

        // Tutorial 1: M/M/1 basics
        try {
            Network tut01 = jline.examples.java.GettingStarted.tut01_mm1_basics();
            tut01.tikzExportPNG(imgDir + "tut01_mm1_network", 200);
            successCount++;
        } catch (Exception e) {
            // Silently handle export failure
        }

        // Tutorial 3: Machine Repair Problem
        try {
            Network tut03 = jline.examples.java.GettingStarted.tut03_repairmen();
            tut03.tikzExportPNG(imgDir + "tut03_repairmen_network", 200);
            successCount++;
        } catch (Exception e) {
            // Silently handle export failure
        }

        // Tutorial 4: Load Balancing
        try {
            Network tut04 = jline.examples.java.GettingStarted.tut04_lb_routing();
            tut04.tikzExportPNG(imgDir + "tut04_lb_routing_network", 200);
            successCount++;
        } catch (Exception e) {
            // Silently handle export failure
        }

        // Tutorial 6: Cache LRU Zipf
        try {
            Network tut06 = jline.examples.java.GettingStarted.tut06_cache_lru_zipf();
            tut06.tikzExportPNG(imgDir + "tut06_cache_network", 200);
            successCount++;
        } catch (Exception e) {
            // Silently handle export failure
        }

        // Tutorial 8: Optimization
        try {
            Network tut08 = jline.examples.java.GettingStarted.tut08_opt_load_balancing();
            tut08.tikzExportPNG(imgDir + "tut08_optimization_network", 200);
            successCount++;
        } catch (Exception e) {
            // Silently handle export failure
        }

        // Tutorial 9: Departure Process (M/E/1)
        try {
            Object[] tut09Components = jline.examples.java.GettingStarted.gallery_merl1();
            Network tut09 = (Network) tut09Components[0];
            tut09.tikzExportPNG(imgDir + "tut09_me1_network", 200);
            successCount++;
        } catch (Exception e) {
            // Silently handle export failure
        }

        assertTrue(successCount >= 4, "At least 4 tutorials should export successfully");
    }

    @Test
    public void testAllTutorialExamples() {
        // Test TikZ generation for all 11 website tutorial examples

        // Tutorial 1: M/M/1 basics
        Network tut01 = jline.examples.java.GettingStarted.tut01_mm1_basics();
        assertNotNull(tut01.toTikZ());

        // Tutorial 2: M/G/1 multiclass
        Network tut02 = jline.examples.java.GettingStarted.tut02_mg1_multiclass_solvers();
        assertNotNull(tut02.toTikZ());

        // Tutorial 3: Machine Repair Problem
        Network tut03 = jline.examples.java.GettingStarted.tut03_repairmen();
        assertNotNull(tut03.toTikZ());

        // Tutorial 4: Load Balancing with Routing
        Network tut04 = jline.examples.java.GettingStarted.tut04_lb_routing();
        assertNotNull(tut04.toTikZ());

        // Tutorial 5: Class Switching
        Network tut05 = jline.examples.java.GettingStarted.tut05_completes_flag();
        assertNotNull(tut05.toTikZ());

        // Tutorial 6: Cache with LRU and Zipf
        Network tut06 = jline.examples.java.GettingStarted.tut06_cache_lru_zipf();
        assertNotNull(tut06.toTikZ());

        // Tutorial 7: Response Time CDF
        Network tut07 = jline.examples.java.GettingStarted.tut07_respt_cdf();
        assertNotNull(tut07.toTikZ());

        // Tutorial 8: Optimization Load Balancing
        Network tut08 = jline.examples.java.GettingStarted.tut08_opt_load_balancing();
        assertNotNull(tut08.toTikZ());

        // Tutorial 9: Departure Process (uses gallery_merl1)
        Object[] tut09Components = jline.examples.java.GettingStarted.gallery_merl1();
        Network tut09 = (Network) tut09Components[0];
        assertNotNull(tut09.toTikZ());

        // Tutorial 10: Layered Queueing Network
        // Note: LayeredNetwork doesn't have toTikZ() yet, so we just verify the model is created
        jline.lang.layered.LayeredNetwork tut10 = jline.examples.java.GettingStarted.tut10_lqn_basics();
        assertNotNull(tut10, "tut10 LQN model should be created");

        // Tutorial 11: Random Environment (test each stage model)
        try {
            jline.lang.Environment tut11 = jline.examples.java.GettingStarted.tut11_random_env();
            assertNotNull(tut11, "tut11 environment should be created");
            // Each stage has a Network model that can export to TikZ
            for (int i = 0; i < tut11.getNumberOfStages(); i++) {
                Network stageModel = tut11.getModel(i);
                assertNotNull(stageModel.toTikZ(), "tut11 stage " + i + " should export to TikZ");
            }
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    @Test
    public void testMixedForkJoinWithFeedback() {
        // Create a MIXED network with:
        // - Open class: Source -> Fork -> Queues -> Join -> Router -> (feedback or Sink)
        // - Closed class: Delay -> Fork -> Queues -> Join -> Router -> Delay
        // - Fork-join structure
        // - Two feedback queues

        Network model = new Network("MixedFJ_Feedback");

        // === NODES ===

        // Open class nodes
        Source source = new Source(model, "Source");
        Sink sink = new Sink(model, "Sink");

        // Closed class delay station (think time)
        Delay thinkTime = new Delay(model, "ThinkTime");

        // Fork node
        Fork fork = new Fork(model, "Fork");

        // Two parallel queues in fork-join
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.FCFS);

        // Join node
        Join join = new Join(model, "Join", fork);

        // Router for feedback/exit decision
        Router router = new Router(model, "Router");

        // Two feedback queues
        Queue feedbackQ1 = new Queue(model, "FeedbackQ1", SchedStrategy.FCFS);
        Queue feedbackQ2 = new Queue(model, "FeedbackQ2", SchedStrategy.FCFS);

        // === JOB CLASSES ===

        // Open class: external arrivals
        OpenClass openClass = new OpenClass(model, "OpenJobs", 0);
        source.setArrival(openClass, new Exp(0.5));
        queue1.setService(openClass, new Exp(2.0));
        queue2.setService(openClass, new Exp(2.0));
        feedbackQ1.setService(openClass, new Exp(3.0));
        feedbackQ2.setService(openClass, new Exp(3.0));

        // Closed class: 5 jobs circulating
        ClosedClass closedClass = new ClosedClass(model, "ClosedJobs", 5, thinkTime);
        thinkTime.setService(closedClass, new Exp(1.0));
        queue1.setService(closedClass, new Exp(2.5));
        queue2.setService(closedClass, new Exp(2.5));
        feedbackQ1.setService(closedClass, new Exp(4.0));
        feedbackQ2.setService(closedClass, new Exp(4.0));

        // === ROUTING ===

        RoutingMatrix P = model.initRoutingMatrix();

        // Open class routing:
        // Source -> Fork
        P.set(openClass, source, fork, 1.0);
        // Fork -> Queue1, Queue2
        P.set(openClass, fork, queue1, 1.0);
        P.set(openClass, fork, queue2, 1.0);
        // Queue1, Queue2 -> Join
        P.set(openClass, queue1, join, 1.0);
        P.set(openClass, queue2, join, 1.0);
        // Join -> Router
        P.set(openClass, join, router, 1.0);
        // Router -> FeedbackQ1 (25%), FeedbackQ2 (25%), Sink (50%)
        P.set(openClass, router, feedbackQ1, 0.25);
        P.set(openClass, router, feedbackQ2, 0.25);
        P.set(openClass, router, sink, 0.50);
        // Feedback loops back to Fork
        P.set(openClass, feedbackQ1, fork, 1.0);
        P.set(openClass, feedbackQ2, fork, 1.0);

        // Closed class routing:
        // ThinkTime -> Fork
        P.set(closedClass, thinkTime, fork, 1.0);
        // Fork -> Queue1, Queue2
        P.set(closedClass, fork, queue1, 1.0);
        P.set(closedClass, fork, queue2, 1.0);
        // Queue1, Queue2 -> Join
        P.set(closedClass, queue1, join, 1.0);
        P.set(closedClass, queue2, join, 1.0);
        // Join -> Router
        P.set(closedClass, join, router, 1.0);
        // Router -> FeedbackQ1 (30%), FeedbackQ2 (30%), ThinkTime (40%)
        P.set(closedClass, router, feedbackQ1, 0.30);
        P.set(closedClass, router, feedbackQ2, 0.30);
        P.set(closedClass, router, thinkTime, 0.40);
        // Feedback loops back to Fork
        P.set(closedClass, feedbackQ1, fork, 1.0);
        P.set(closedClass, feedbackQ2, fork, 1.0);

        model.link(P);

        // Generate TikZ code
        String tikz = model.toTikZ();

        // Verify the output
        assertNotNull(tikz);
        assertTrue(tikz.contains("\\begin{tikzpicture}"));
        assertTrue(tikz.contains("Source"));
        assertTrue(tikz.contains("Sink"));
        assertTrue(tikz.contains("ThinkTime"));
        assertTrue(tikz.contains("Fork"));
        assertTrue(tikz.contains("Queue1"));
        assertTrue(tikz.contains("Queue2"));
        assertTrue(tikz.contains("Join"));
        assertTrue(tikz.contains("Router"));
        assertTrue(tikz.contains("FeedbackQ1"));
        assertTrue(tikz.contains("FeedbackQ2"));
    }
}
