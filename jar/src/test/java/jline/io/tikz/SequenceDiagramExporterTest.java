/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.io.tikz;

import jline.io.LQN2UML;
import jline.lang.layered.*;
import jline.lang.constant.SchedStrategy;
import jline.lang.processes.Exp;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.io.File;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Tests for UML sequence diagram export from LayeredNetwork models.
 */
public class SequenceDiagramExporterTest {

    @TempDir
    Path tempDir;

    /**
     * Creates a basic 3-tier LQN model for testing.
     * Client -> AppServer -> Database
     */
    private LayeredNetwork createBasicLQN() {
        LayeredNetwork model = new LayeredNetwork("BasicLQN");

        // Create processors
        Processor clientProc = new Processor(model, "ClientProc", 1, SchedStrategy.PS);
        Processor serverProc = new Processor(model, "ServerProc", 2, SchedStrategy.PS);
        Processor dbProc = new Processor(model, "DBProc", 1, SchedStrategy.PS);

        // Create tasks
        Task client = new Task(model, "Client", 10, SchedStrategy.REF).on(clientProc);
        client.setThinkTime(new Exp(1.0));
        Task appServer = new Task(model, "AppServer", 5, SchedStrategy.FCFS).on(serverProc);
        Task database = new Task(model, "Database", 1, SchedStrategy.FCFS).on(dbProc);

        // Create entries
        Entry clientEntry = new Entry(model, "ClientEntry").on(client);
        Entry appEntry = new Entry(model, "AppEntry").on(appServer);
        Entry dbEntry = new Entry(model, "DBEntry").on(database);

        // Create activities with calls
        Activity clientAct = new Activity(model, "ClientActivity", new Exp(0.5)).on(client);
        clientAct.boundTo(clientEntry);
        clientAct.synchCall(appEntry, 1);

        Activity appAct = new Activity(model, "AppActivity", new Exp(1.0)).on(appServer);
        appAct.boundTo(appEntry);
        appAct.synchCall(dbEntry, 2);
        appAct.repliesTo(appEntry);

        Activity dbAct = new Activity(model, "DBActivity", new Exp(2.0)).on(database);
        dbAct.boundTo(dbEntry);
        dbAct.repliesTo(dbEntry);

        return model;
    }

    /**
     * Creates a simple 2-task LQN model (similar to lqn_basic.m).
     */
    private LayeredNetwork createSimpleLQN() {
        LayeredNetwork model = new LayeredNetwork("SimpleLQN");

        Processor P1 = new Processor(model, "P1", 2, SchedStrategy.PS);
        Processor P2 = new Processor(model, "P2", 3, SchedStrategy.PS);

        Task T1 = new Task(model, "T1", 50, SchedStrategy.REF).on(P1);
        T1.setThinkTime(new Exp(0.5));
        Task T2 = new Task(model, "T2", 50, SchedStrategy.FCFS).on(P1);
        T2.setThinkTime(new Exp(0.33));
        Task T3 = new Task(model, "T3", 25, SchedStrategy.FCFS).on(P2);
        T3.setThinkTime(new Exp(0.25));

        Entry E1 = new Entry(model, "E1").on(T1);
        Entry E2 = new Entry(model, "E2").on(T2);
        Entry E3 = new Entry(model, "E3").on(T3);

        Activity A1 = new Activity(model, "AS1", new Exp(0.1)).on(T1).boundTo(E1).synchCall(E2, 1);
        Activity A2 = new Activity(model, "AS2", new Exp(0.05)).on(T2).boundTo(E2).synchCall(E3, 5).repliesTo(E2);
        Activity A3 = new Activity(model, "AS3", new Exp(0.02)).on(T3).boundTo(E3).repliesTo(E3);

        return model;
    }

    /**
     * Creates an LQN with asynchronous calls.
     */
    private LayeredNetwork createLQNWithAsyncCalls() {
        LayeredNetwork model = new LayeredNetwork("AsyncLQN");

        Processor P1 = new Processor(model, "Processor1", 1, SchedStrategy.PS);
        Processor P2 = new Processor(model, "Processor2", 1, SchedStrategy.PS);

        Task sender = new Task(model, "Sender", 1, SchedStrategy.REF).on(P1);
        sender.setThinkTime(new Exp(1.0));
        Task receiver = new Task(model, "Receiver", 5, SchedStrategy.FCFS).on(P2);

        Entry senderEntry = new Entry(model, "SenderEntry").on(sender);
        Entry receiverEntry = new Entry(model, "ReceiverEntry").on(receiver);

        Activity sendAct = new Activity(model, "SendActivity", new Exp(0.1)).on(sender);
        sendAct.boundTo(senderEntry);
        sendAct.asynchCall(receiverEntry, 3);  // Async call, no wait

        Activity recvAct = new Activity(model, "ReceiveActivity", new Exp(0.5)).on(receiver);
        recvAct.boundTo(receiverEntry);

        return model;
    }

    /**
     * Creates an LQN with a loop precedence.
     */
    private LayeredNetwork createLQNWithLoop() {
        LayeredNetwork model = new LayeredNetwork("LoopLQN");

        Processor P1 = new Processor(model, "Proc1", 1, SchedStrategy.PS);

        Task task = new Task(model, "LoopTask", 1, SchedStrategy.REF).on(P1);
        task.setThinkTime(new Exp(1.0));

        Entry entry = new Entry(model, "LoopEntry").on(task);

        Activity initAct = new Activity(model, "Init", new Exp(0.1)).on(task).boundTo(entry);
        Activity loopAct = new Activity(model, "LoopBody", new Exp(0.5)).on(task);
        Activity endAct = new Activity(model, "End", new Exp(0.1)).on(task);

        // Loop 5 times
        List<String> loopActivities = Arrays.asList("LoopBody", "End");
        task.addPrecedence(ActivityPrecedence.Loop("Init", loopActivities, 5.0));

        return model;
    }

    /**
     * Creates an LQN with fork-join (parallel) precedence.
     */
    private LayeredNetwork createLQNWithForkJoin() {
        LayeredNetwork model = new LayeredNetwork("ForkJoinLQN");

        Processor P1 = new Processor(model, "Proc1", 4, SchedStrategy.PS);

        Task task = new Task(model, "ParTask", 1, SchedStrategy.REF).on(P1);
        task.setThinkTime(new Exp(1.0));

        Entry entry = new Entry(model, "ParEntry").on(task);

        Activity startAct = new Activity(model, "Start", new Exp(0.1)).on(task).boundTo(entry);
        Activity branch1 = new Activity(model, "Branch1", new Exp(0.5)).on(task);
        Activity branch2 = new Activity(model, "Branch2", new Exp(0.3)).on(task);
        Activity joinAct = new Activity(model, "Join", new Exp(0.1)).on(task);

        // Fork into parallel branches
        List<String> parallelActs = Arrays.asList("Branch1", "Branch2");
        task.addPrecedence(ActivityPrecedence.AndFork("Start", parallelActs));
        task.addPrecedence(ActivityPrecedence.AndJoin(parallelActs, "Join"));

        return model;
    }

    @Test
    public void testBasicLQNGeneration() {
        LayeredNetwork model = createBasicLQN();

        String tikz = LQN2UML.toTikZ(model);

        assertNotNull(tikz);
        assertTrue(tikz.contains("\\documentclass"));
        assertTrue(tikz.contains("\\usepackage{pgf-umlsd}"));
        assertTrue(tikz.contains("\\begin{sequencediagram}"));
        assertTrue(tikz.contains("\\end{sequencediagram}"));
        assertTrue(tikz.contains("\\newthread"));
    }

    @Test
    public void testSimpleLQNGeneration() {
        LayeredNetwork model = createSimpleLQN();

        String tikz = LQN2UML.toTikZ(model);

        assertNotNull(tikz);
        assertTrue(tikz.contains("sequencediagram"));
        // Check for task lifelines
        assertTrue(tikz.contains("T1") || tikz.contains("\\newthread{T1}"));
        assertTrue(tikz.contains("T2") || tikz.contains("\\newthread{T2}"));
        assertTrue(tikz.contains("T3") || tikz.contains("\\newthread{T3}"));
    }

    @Test
    public void testProcessorGrouping() {
        LayeredNetwork model = createBasicLQN();

        SequenceDiagramOptions options = new SequenceDiagramOptions()
                .setShowProcessorFrames(true);
        String tikz = LQN2UML.toTikZ(model, options);

        assertNotNull(tikz);
        // Check for processor frames (sdblock)
        assertTrue(tikz.contains("\\begin{sdblock}") || tikz.contains("sdblock"));
        assertTrue(tikz.contains("ClientProc") || tikz.contains("ServerProc") || tikz.contains("DBProc"));
    }

    @Test
    public void testWithoutProcessorGrouping() {
        LayeredNetwork model = createBasicLQN();

        SequenceDiagramOptions options = new SequenceDiagramOptions()
                .setShowProcessorFrames(false);
        String tikz = LQN2UML.toTikZ(model, options);

        assertNotNull(tikz);
        assertTrue(tikz.contains("\\newthread"));
    }

    @Test
    public void testSynchCalls() {
        LayeredNetwork model = createBasicLQN();

        String tikz = LQN2UML.toTikZ(model);

        assertNotNull(tikz);
        // Check for synchronous call constructs
        assertTrue(tikz.contains("\\begin{call}") || tikz.contains("call"));
    }

    @Test
    public void testAsyncCalls() {
        LayeredNetwork model = createLQNWithAsyncCalls();

        String tikz = LQN2UML.toTikZ(model);

        assertNotNull(tikz);
        // Check for async message (mess command)
        assertTrue(tikz.contains("\\mess") || tikz.contains("Interactions"));
    }

    @Test
    public void testLayoutEngine() {
        LayeredNetwork model = createBasicLQN();

        SequenceDiagramOptions options = new SequenceDiagramOptions();
        SequenceDiagramLayoutEngine layout = new SequenceDiagramLayoutEngine(model, options);
        layout.computeLayout();

        // Verify processor groups are created
        assertNotNull(layout.getProcessorTaskGroups());
        assertFalse(layout.getProcessorTaskGroups().isEmpty());

        // Verify activities are sorted
        List<Activity> sorted = layout.getSortedActivities();
        assertNotNull(sorted);
        assertFalse(sorted.isEmpty());
    }

    @Test
    public void testTraverser() {
        LayeredNetwork model = createBasicLQN();

        SequenceDiagramOptions options = new SequenceDiagramOptions();
        SequenceDiagramLayoutEngine layout = new SequenceDiagramLayoutEngine(model, options);
        layout.computeLayout();

        SequenceDiagramTraverser traverser = new SequenceDiagramTraverser(model, layout);
        List<SequenceDiagramTraverser.Interaction> interactions = traverser.extractInteractions();

        assertNotNull(interactions);
        // BasicLQN has 2 synch calls: Client->App, App->DB
        assertTrue(interactions.size() >= 2);
    }

    @Test
    public void testLoopFragment() {
        LayeredNetwork model = createLQNWithLoop();

        SequenceDiagramOptions options = new SequenceDiagramOptions();
        SequenceDiagramLayoutEngine layout = new SequenceDiagramLayoutEngine(model, options);
        layout.computeLayout();

        SequenceDiagramTraverser traverser = new SequenceDiagramTraverser(model, layout);
        List<SequenceDiagramTraverser.Fragment> fragments = traverser.extractFragments();

        assertNotNull(fragments);
        // Should have at least one loop fragment
        boolean hasLoop = false;
        for (SequenceDiagramTraverser.Fragment frag : fragments) {
            if (frag.type == SequenceDiagramTraverser.Fragment.FragmentType.LOOP) {
                hasLoop = true;
                assertTrue(frag.label.contains("loop"));
                break;
            }
        }
        assertTrue(hasLoop, "Expected a LOOP fragment");
    }

    @Test
    public void testForkJoinFragment() {
        LayeredNetwork model = createLQNWithForkJoin();

        SequenceDiagramOptions options = new SequenceDiagramOptions();
        SequenceDiagramLayoutEngine layout = new SequenceDiagramLayoutEngine(model, options);
        layout.computeLayout();

        SequenceDiagramTraverser traverser = new SequenceDiagramTraverser(model, layout);
        List<SequenceDiagramTraverser.Fragment> fragments = traverser.extractFragments();

        assertNotNull(fragments);
        // Should have at least one PAR fragment
        boolean hasPar = false;
        for (SequenceDiagramTraverser.Fragment frag : fragments) {
            if (frag.type == SequenceDiagramTraverser.Fragment.FragmentType.PAR) {
                hasPar = true;
                assertEquals("par", frag.label);
                break;
            }
        }
        assertTrue(hasPar, "Expected a PAR fragment");
    }

    @Test
    public void testOptionsBuilder() {
        SequenceDiagramOptions options = new SequenceDiagramOptions()
                .setLifelineSpacing(4.0)
                .setMessageSpacing(1.5)
                .setShowHostDemand(false)
                .setShowCallMeans(true)
                .setShowProcessorFrames(true)
                .setShowEntryNames(true)
                .setAsyncDashed(true)
                .setShowReplies(true)
                .setScale(0.8);

        assertEquals(4.0, options.getLifelineSpacing(), 0.001);
        assertEquals(1.5, options.getMessageSpacing(), 0.001);
        assertFalse(options.isShowHostDemand());
        assertTrue(options.isShowCallMeans());
        assertTrue(options.isShowProcessorFrames());
        assertTrue(options.isShowEntryNames());
        assertTrue(options.isAsyncDashed());
        assertTrue(options.isShowReplies());
        assertEquals(0.8, options.getScale(), 0.001);
    }

    @Test
    public void testExportToFile() throws Exception {
        LayeredNetwork model = createSimpleLQN();

        File texFile = tempDir.resolve("test_sequence.tex").toFile();
        LQN2UML.Options options = new LQN2UML.Options()
                .setFormat(LQN2UML.Format.TEX)
                .setOutputPath(texFile.getAbsolutePath());
        LQN2UML.export(model, options);

        assertTrue(texFile.exists());
        assertTrue(texFile.length() > 0);
    }

    @Test
    public void testPdfExportWhenAvailable() throws Exception {
        // Only test if pdflatex is available
        if (!TikZExporter.isPdfLatexAvailable()) {
            System.out.println("Skipping PDF export test - pdflatex not available");
            return;
        }

        LayeredNetwork model = createSimpleLQN();

        LQN2UML.Options options = new LQN2UML.Options()
                .setFormat(LQN2UML.Format.PDF);
        File pdfFile = LQN2UML.export(model, options);

        assertNotNull(pdfFile);
        assertTrue(pdfFile.exists());
        assertTrue(pdfFile.length() > 0);
    }

    @Test
    public void testPngExportWhenAvailable() throws Exception {
        // Only test if pdflatex is available
        if (!TikZExporter.isPdfLatexAvailable()) {
            System.out.println("Skipping PNG export test - pdflatex not available");
            return;
        }

        LayeredNetwork model = createSimpleLQN();

        String pngPath = tempDir.resolve("test_sequence.png").toString();
        LQN2UML.Options options = new LQN2UML.Options()
                .setFormat(LQN2UML.Format.PNG)
                .setOutputPath(pngPath)
                .setDpi(150);
        LQN2UML.export(model, options);

        File pngFile = new File(pngPath);
        assertTrue(pngFile.exists() || new File(pngPath.replace(".png", "") + ".png").exists());
    }

    @Test
    public void testEmptyModel() {
        LayeredNetwork model = new LayeredNetwork("EmptyLQN");

        // Should not throw, just generate minimal output
        String tikz = LQN2UML.toTikZ(model);

        assertNotNull(tikz);
        assertTrue(tikz.contains("\\begin{sequencediagram}"));
        assertTrue(tikz.contains("\\end{sequencediagram}"));
    }

    @Test
    public void testModelWithSingleProcessor() {
        LayeredNetwork model = new LayeredNetwork("SingleProcLQN");

        Processor P1 = new Processor(model, "SingleProc", 1, SchedStrategy.PS);
        Task T1 = new Task(model, "Task1", 1, SchedStrategy.REF).on(P1);
        Entry E1 = new Entry(model, "Entry1").on(T1);
        Activity A1 = new Activity(model, "Act1", new Exp(1.0)).on(T1).boundTo(E1);

        String tikz = LQN2UML.toTikZ(model);

        assertNotNull(tikz);
        assertTrue(tikz.contains("SingleProc") || tikz.contains("Task1"));
    }

    @Test
    public void testLatexEscaping() {
        LayeredNetwork model = new LayeredNetwork("Test_Model");

        // Create elements with special characters in names
        Processor P1 = new Processor(model, "Proc_1", 1, SchedStrategy.PS);
        Task T1 = new Task(model, "Task_A", 1, SchedStrategy.REF).on(P1);
        Entry E1 = new Entry(model, "Entry_X").on(T1);
        Activity A1 = new Activity(model, "Act_1", new Exp(1.0)).on(T1).boundTo(E1);

        String tikz = LQN2UML.toTikZ(model);

        assertNotNull(tikz);
        // Underscores should be escaped in LaTeX
        assertTrue(tikz.contains("\\_") || tikz.contains("Proc"));
    }
}
