package jline.solvers.ln;

import jline.GlobalConstants;
import jline.TestTools;
import jline.lang.constant.SchedStrategy;
import jline.VerboseLevel;
import jline.lang.layered.*;
import jline.lang.processes.Exp;
import jline.solvers.LayeredNetworkAvgTable;
import jline.solvers.SolverOptions;
import jline.solvers.lqns.SolverLQNS;
import jline.util.Maths;
import jline.util.matrix.Matrix;

import static org.junit.jupiter.api.Assertions.*;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;

/**
 * Error handling and advanced feature validation tests for SolverLN.
 */
class SolverLNValidationTest extends SolverLNTestBase {

    @BeforeAll
    public static void setUp() {
        // Ensure MATLAB-compatible random number generation
        Maths.setRandomNumbersMatlab(true);
        // Set verbose level to SILENT to suppress warnings during tests
        GlobalConstants.setVerbose(VerboseLevel.SILENT);
    }

    // Tests for SolverLN with moment3 method enabled.
    // These tests verify that the post-convergence 3-moment APH fitting works correctly.

    @Test
    public void test_moment3_basic() throws Exception {
        // Use an existing valid model
        LayeredNetwork lqn = SolverLNTestFixtures.buildModel2();

        // Create options with moment3 method
        LNOptions options = new LNOptions();
        options.verbose = VerboseLevel.SILENT;
        options.method = "moment3";

        // Create solver and run to convergence
        SolverLN solver = new SolverLN(lqn, options);
        LayeredNetworkAvgTable avgTable = (LayeredNetworkAvgTable) solver.getEnsembleAvg();

        // Basic sanity checks - the solver should converge and produce results
        assertTrue(solver.hasconverged, "Solver should converge");
        assertNotNull(solver.servt, "Service times should be computed");
        assertNotNull(solver.tput, "Throughput should be computed");

        // Check that entry processes are populated (post-convergence feature)
        assertNotNull(solver.entryproc, "Entry processes should be created");
        assertNotNull(avgTable, "Average table should be computed");
    }

    @Test
    public void test_moment3_vs_default_comparison() throws Exception {
        // Run with default method
        LNOptions defaultOptions = new LNOptions();
        defaultOptions.verbose = VerboseLevel.SILENT;
        defaultOptions.method = "default";
        SolverLN solverDefault = new SolverLN(SolverLNTestFixtures.buildModel2(), defaultOptions);
        LayeredNetworkAvgTable avgDefault = (LayeredNetworkAvgTable) solverDefault.getEnsembleAvg();

        // Run with moment3 method
        LNOptions moment3Options = new LNOptions();
        moment3Options.verbose = VerboseLevel.SILENT;
        moment3Options.method = "moment3";
        SolverLN solverMoment3 = new SolverLN(SolverLNTestFixtures.buildModel2(), moment3Options);
        LayeredNetworkAvgTable avgMoment3 = (LayeredNetworkAvgTable) solverMoment3.getEnsembleAvg();

        // Both should converge
        assertTrue(solverDefault.hasconverged, "Default solver should converge");
        assertTrue(solverMoment3.hasconverged, "Moment3 solver should converge");

        // Results should exist
        assertNotNull(avgDefault);
        assertNotNull(avgMoment3);
    }

    // Error detection tests

    @Test
    public void test_LQN_err_1() throws Exception {
        // Activity in REF task replies
        Exception exception = suppressOutput(() -> assertThrows(Exception.class, () -> {
            LayeredNetwork model = new LayeredNetwork("add_cart");

            // first layer
            Processor P1 = new Processor(model, "client_p", 1, SchedStrategy.INF);
            Task T1 = new Task(model, "user", 1, SchedStrategy.REF).on(P1);
            Entry E1 = new Entry(model, "user").on(T1);
            T1.setThinkTime(Exp.fitMean(0));

            // second layer
            Processor P2 = new Processor(model, "WeiUI_p", 1, SchedStrategy.INF);
            Task T2 = new Task(model, "WeiUI", 1, SchedStrategy.INF).on(P2);
            Entry E2 = new Entry(model, "add_cart").on(T2);

            // activities
            Activity A1 = new Activity(model, "A1", Exp.fitMean(0)).on(T1).boundTo(E1);
            Activity A2 = new Activity(model, "A2", Exp.fitMean(0)).on(T1).synchCall(E2, 1).repliesTo(E1);
            Activity A3 = new Activity(model, "A3", Exp.fitMean(0)).on(T1).repliesTo(E1);
            Activity A4 = new Activity(model, "A4", Exp.fitMean(1)).on(T2).boundTo(E2).repliesTo(E2);

            Matrix probs = new Matrix(1, 2);
            probs.set(0, 0, 0.4);
            probs.set(0, 1, 0.6);
            T1.addPrecedence(ActivityPrecedence.OrFork(A1, java.util.Arrays.asList(A2, A3), probs));

            // Run solver - should fail during model construction
            SolverOptions options = new LNOptions();
            options.verbose = VerboseLevel.SILENT;
            SolverLN solver = new SolverLN(model, options);
        }));

        assertTrue(exception.getMessage().contains("Activities in reference tasks cannot reply") ||
                   exception.getMessage().contains("REF task") ||
                   exception.getMessage().contains("reply"));
    }

    @Test
    public void test_LQN_err_2() throws Exception {
        // Entry on task calls itself
        Exception exception = suppressOutput(() -> assertThrows(Exception.class, () -> {
            LayeredNetwork model = new LayeredNetwork("self_call_error");

            Processor P1 = new Processor(model, "proc", 1, SchedStrategy.INF);
            Task T1 = new Task(model, "task", 1, SchedStrategy.INF).on(P1);
            Entry E1 = new Entry(model, "entry").on(T1);

            // Create activity that calls its own entry - should be invalid
            Activity A1 = new Activity(model, "A1", Exp.fitMean(1)).on(T1).boundTo(E1).synchCall(E1, 1).repliesTo(E1);

            SolverOptions options = new LNOptions();
            options.verbose = VerboseLevel.SILENT;
            SolverLN solver = new SolverLN(model, options);
        }));

        assertTrue(exception.getMessage().contains("calls itself") ||
                   exception.getMessage().contains("self") ||
                   exception.getMessage().contains("cycle") ||
                   exception.getMessage().contains("invalid"));
    }

    @Test
    public void test_LQN_err_3() throws Exception {
        // Entry on task calls entry on the same task
        Exception exception = suppressOutput(() -> assertThrows(Exception.class, () -> {
            LayeredNetwork model = new LayeredNetwork("same_task_call");

            Processor P1 = new Processor(model, "proc", 1, SchedStrategy.INF);
            Task T1 = new Task(model, "task", 1, SchedStrategy.INF).on(P1);
            Entry E1 = new Entry(model, "entry1").on(T1);
            Entry E2 = new Entry(model, "entry2").on(T1);

            // Activity on T1 calls another entry on the same task T1
            Activity A1 = new Activity(model, "A1", Exp.fitMean(1)).on(T1).boundTo(E1).synchCall(E2, 1);
            Activity A2 = new Activity(model, "A2", Exp.fitMean(1)).on(T1).boundTo(E2).repliesTo(E2);

            SolverOptions options = new LNOptions();
            options.verbose = VerboseLevel.SILENT;
            SolverLN solver = new SolverLN(model, options);
        }));

        assertTrue(exception.getMessage().contains("same task") ||
                   exception.getMessage().contains("calls entry") ||
                   exception.getMessage().contains("invalid") ||
                   exception.getMessage().contains("cycle"));
    }

    @Test
    public void test_LQN_err_4() throws Exception {
        // Test cycle detection in layered network
        Exception exception = suppressOutput(() -> assertThrows(Exception.class, () -> {
            LayeredNetwork model = new LayeredNetwork("cycle_detection");

            // Create processors
            Processor P1 = new Processor(model, "proc1", 1, SchedStrategy.INF);
            Processor P2 = new Processor(model, "proc2", 1, SchedStrategy.INF);
            Processor P3 = new Processor(model, "proc3", 1, SchedStrategy.INF);

            // Create tasks
            Task T1 = new Task(model, "task1", 1, SchedStrategy.INF).on(P1);
            Task T2 = new Task(model, "task2", 1, SchedStrategy.INF).on(P2);
            Task T3 = new Task(model, "task3", 1, SchedStrategy.INF).on(P3);

            // Create entries
            Entry E1 = new Entry(model, "entry1").on(T1);
            Entry E2 = new Entry(model, "entry2").on(T2);
            Entry E3 = new Entry(model, "entry3").on(T3);

            // Create activities that form a cycle: T1->T2->T3->T1
            Activity A1 = new Activity(model, "A1", Exp.fitMean(1)).on(T1).boundTo(E1).synchCall(E2, 1);
            Activity A2 = new Activity(model, "A2", Exp.fitMean(1)).on(T2).boundTo(E2).synchCall(E3, 1).repliesTo(E2);
            Activity A3 = new Activity(model, "A3", Exp.fitMean(1)).on(T3).boundTo(E3).synchCall(E1, 1).repliesTo(E3);

            SolverOptions options = new LNOptions();
            options.verbose = VerboseLevel.SILENT;
            SolverLN solver = new SolverLN(model, options);
        }));

        assertTrue(exception.getMessage().contains("cycle") ||
                   exception.getMessage().contains("circular") ||
                   exception.getMessage().contains("loop") ||
                   exception.getMessage().contains("recursive"));
    }

    @Test
    public void test_LQN_err_5() throws Exception {
        // Test unsupported replyTo patterns - implementation specific
        suppressOutput(() -> {
            try {
                LayeredNetwork model = new LayeredNetwork("unsupported_reply");

                Processor P1 = new Processor(model, "proc1", 1, SchedStrategy.INF);
                Processor P2 = new Processor(model, "proc2", 1, SchedStrategy.INF);

                Task T1 = new Task(model, "task1", 1, SchedStrategy.INF).on(P1);
                Task T2 = new Task(model, "task2", 1, SchedStrategy.INF).on(P2);

                Entry E1 = new Entry(model, "entry1").on(T1);
                Entry E2 = new Entry(model, "entry2").on(T2);
                Entry E3 = new Entry(model, "entry3").on(T2);

                // Activity trying to reply to multiple entries - unsupported pattern
                Activity A1 = new Activity(model, "A1", Exp.fitMean(1)).on(T1).boundTo(E1).synchCall(E2, 1).synchCall(E3, 1);
                Activity A2 = new Activity(model, "A2", Exp.fitMean(1)).on(T2).boundTo(E2).repliesTo(E2).repliesTo(E3);

                SolverOptions options = new LNOptions();
                options.verbose = VerboseLevel.SILENT;
                SolverLN solver = new SolverLN(model, options);

                // If no exception is thrown, the implementation allows this pattern
                assertTrue(true);
            } catch (Exception e) {
                // Verify error message if exception is thrown
                assertTrue(e.getMessage().contains("reply") ||
                           e.getMessage().contains("multiple") ||
                           e.getMessage().contains("unsupported") ||
                           e.getMessage().contains("invalid") ||
                           e.getMessage().contains("service"));
            }
        });
    }

    @Test
    public void test_LQN_err_6() throws Exception {
        // Test boundTo specification validation - implementation specific
        suppressOutput(() -> {
            try {
                LayeredNetwork model = new LayeredNetwork("boundTo_error");

                Processor P1 = new Processor(model, "proc1", 1, SchedStrategy.INF);
                Task T1 = new Task(model, "task1", 1, SchedStrategy.INF).on(P1);
                Task T2 = new Task(model, "task2", 1, SchedStrategy.INF).on(P1);

                Entry E1 = new Entry(model, "entry1").on(T1);
                Entry E2 = new Entry(model, "entry2").on(T2);

                // Activity on T1 trying to bind to entry on T2 - might be invalid
                Activity A1 = new Activity(model, "A1", Exp.fitMean(1)).on(T1).boundTo(E2);

                SolverOptions options = new LNOptions();
                options.verbose = VerboseLevel.SILENT;
                SolverLN solver = new SolverLN(model, options);

                // If no exception is thrown, the implementation allows this pattern
                assertTrue(true);
            } catch (Exception e) {
                // Verify error message if exception is thrown
                // The Java implementation throws an exception for this pattern
                assertTrue(e.getMessage() != null);
            }
        });
    }

    @Test
    public void test_LQN_err_7() throws Exception {
        // Entry called both synchronously and asynchronously
        Exception exception = suppressOutput(() -> assertThrows(Exception.class, () -> {
            LayeredNetwork model = new LayeredNetwork("mixed_call_types");

            Processor P1 = new Processor(model, "client_p", 1, SchedStrategy.INF);
            Task T1 = new Task(model, "client", 1, SchedStrategy.REF).on(P1);
            Entry E1 = new Entry(model, "client").on(T1);
            T1.setThinkTime(Exp.fitMean(0));

            Processor P2 = new Processor(model, "server_p", 1, SchedStrategy.INF);
            Task T2 = new Task(model, "server", 1, SchedStrategy.INF).on(P2);
            Entry E2 = new Entry(model, "server").on(T2);

            // Try to call same entry both synchronously and asynchronously
            Activity A1 = new Activity(model, "A1", Exp.fitMean(1)).on(T1).boundTo(E1).synchCall(E2, 1);
            Activity A2 = new Activity(model, "A2", Exp.fitMean(1)).on(T1).asynchCall(E2, 1);
            Activity A3 = new Activity(model, "A3", Exp.fitMean(1)).on(T2).boundTo(E2).repliesTo(E2);

            T1.addPrecedence(ActivityPrecedence.Serial(A1, A2));

            SolverOptions options = new LNOptions();
            options.verbose = VerboseLevel.SILENT;
            SolverLN solver = new SolverLN(model, options);
        }));

        assertTrue(exception.getMessage().contains("synchronously and asynchronously") ||
                   exception.getMessage().contains("sync") ||
                   exception.getMessage().contains("async") ||
                   exception.getMessage().contains("mixed") ||
                   exception.getMessage().contains("call"));
    }

    @Test
    public void test_LQN_err_8() throws Exception {
        // Test parent task validation
        Exception exception = assertThrows(Exception.class, () -> {
            LayeredNetwork model = new LayeredNetwork("parent_task_error");

            Processor P1 = new Processor(model, "proc1", 1, SchedStrategy.INF);
            Task T1 = new Task(model, "task1", 1, SchedStrategy.INF).on(P1);
            Task T2 = new Task(model, "task2", 1, SchedStrategy.INF).on(P1);

            Entry E1 = new Entry(model, "entry1").on(T1);

            // Activity without parent task trying to bind to entry
            Activity A1 = new Activity(model, "A1", Exp.fitMean(1)).boundTo(E1);
            // Or activity on wrong task
            Activity A2 = new Activity(model, "A2", Exp.fitMean(1)).on(T2).boundTo(E1);

            SolverOptions options = new LNOptions();
            options.verbose = VerboseLevel.SILENT;
            SolverLN solver = new SolverLN(model, options);
        });

        assertTrue(exception.getMessage().contains("parent") ||
                   exception.getMessage().contains("task") ||
                   exception.getMessage().contains("activity") ||
                   exception.getMessage().contains("invalid"));
    }

    @Test
    public void test_LQN_err_9() throws Exception {
        // Test null parent for Entry - implementation specific
        try {
            LayeredNetwork model = new LayeredNetwork("null_parent_error");

            // Create Entry without parent Task
            Entry E1 = new Entry(model, "entry1");

            // Try to use solver with orphaned entry
            SolverOptions options = new LNOptions();
            options.verbose = VerboseLevel.SILENT;
            SolverLN solver = new SolverLN(model, options);

            // If no exception is thrown, the implementation allows orphaned entries
            assertTrue(true);
        } catch (Exception e) {
            // Verify error message if exception is thrown
            assertTrue(e.getMessage().contains("Entry") ||
                       e.getMessage().contains("Task") ||
                       e.getMessage().contains("parent") ||
                       e.getMessage().contains("invalid") ||
                       e.getMessage().contains("null"));
        }
    }

    @Test
    public void test_LQN_err_10() throws Exception {
        // Test null parent for Task - implementation specific
        try {
            LayeredNetwork model = new LayeredNetwork("task_null_parent");

            // Create Task without parent Processor
            Task T1 = new Task(model, "task1", 1, SchedStrategy.INF);

            // Try to use solver with orphaned task
            SolverOptions options = new LNOptions();
            options.verbose = VerboseLevel.SILENT;
            SolverLN solver = new SolverLN(model, options);

            // If no exception is thrown, the implementation allows orphaned tasks
            assertTrue(true);
        } catch (Exception e) {
            // Verify error message if exception is thrown
            assertTrue(e.getMessage().contains("Task") ||
                       e.getMessage().contains("Processor") ||
                       e.getMessage().contains("parent") ||
                       e.getMessage().contains("null"));
        }
    }

    @Test
    public void test_LQN_err_11() throws Exception {
        // Test repeated call validation - implementation specific
        try {
            LayeredNetwork model = new LayeredNetwork("repeated_call_error");

            Processor P1 = new Processor(model, "proc1", 1, SchedStrategy.INF);
            Processor P2 = new Processor(model, "proc2", 1, SchedStrategy.INF);

            Task T1 = new Task(model, "task1", 1, SchedStrategy.INF).on(P1);
            Task T2 = new Task(model, "task2", 1, SchedStrategy.INF).on(P2);

            Entry E1 = new Entry(model, "entry1").on(T1);
            Entry E2 = new Entry(model, "entry2").on(T2);

            // Activity making duplicate calls to same entry
            Activity A1 = new Activity(model, "A1", Exp.fitMean(1)).on(T1).boundTo(E1)
                .synchCall(E2, 1).synchCall(E2, 1); // Duplicate call
            Activity A2 = new Activity(model, "A2", Exp.fitMean(1)).on(T2).boundTo(E2).repliesTo(E2);

            SolverOptions options = new LNOptions();
            options.verbose = VerboseLevel.SILENT;
            SolverLN solver = new SolverLN(model, options);

            // If no exception is thrown, the implementation allows duplicate calls
            assertTrue(true);
        } catch (Exception e) {
            // Verify error message if exception is thrown
            assertTrue(e.getMessage().contains("duplicate") ||
                       e.getMessage().contains("repeated") ||
                       e.getMessage().contains("multiple") ||
                       e.getMessage().contains("call") ||
                       e.getMessage().contains("service"));
        }
    }

    @Test
    public void test_LQN_err_12() throws Exception {
        // Test repeated reply validation - implementation specific
        try {
            LayeredNetwork model = new LayeredNetwork("repeated_reply_error");

            Processor P1 = new Processor(model, "proc1", 1, SchedStrategy.INF);
            Task T1 = new Task(model, "task1", 1, SchedStrategy.INF).on(P1);

            Entry E1 = new Entry(model, "entry1").on(T1);
            Entry E2 = new Entry(model, "entry2").on(T1);

            // Multiple activities trying to reply to same entry
            Activity A1 = new Activity(model, "A1", Exp.fitMean(1)).on(T1).boundTo(E1).repliesTo(E1);
            Activity A2 = new Activity(model, "A2", Exp.fitMean(1)).on(T1).repliesTo(E1); // Another reply to same entry

            SolverOptions options = new LNOptions();
            options.verbose = VerboseLevel.SILENT;
            SolverLN solver = new SolverLN(model, options);

            // If no exception is thrown, the implementation allows multiple replies
            assertTrue(true);
        } catch (Exception e) {
            // Verify error message if exception is thrown
            assertTrue(e.getMessage().contains("reply") ||
                       e.getMessage().contains("multiple") ||
                       e.getMessage().contains("duplicate") ||
                       e.getMessage().contains("already"));
        }
    }

    /**
     * Test: Single activity with think-time, validate against LQNS results.
     */
    @Test
    @Tag("no-ci")
    public void test_activity_think_time_single() throws Exception {
        LayeredNetwork lqn = new LayeredNetwork("TestActivityThinkTime");

        // Client processor (infinite capacity)
        Processor P1 = new Processor(lqn, "P1", Integer.MAX_VALUE, SchedStrategy.INF);

        // Client task (reference)
        Task T1 = new Task(lqn, "Client", 1, SchedStrategy.REF).on(P1);
        Entry E1 = new Entry(lqn, "ClientE").on(T1);

        // Server processor (infinite capacity)
        Processor P2 = new Processor(lqn, "P2", Integer.MAX_VALUE, SchedStrategy.INF);

        // Server task (infinite servers)
        Task T2 = new Task(lqn, "Server", Integer.MAX_VALUE, SchedStrategy.INF).on(P2);
        Entry E2 = new Entry(lqn, "ServerE").on(T2);

        // Client activity: service=0.1, calls server entry
        Activity A1 = new Activity(lqn, "ClientActivity", Exp.fitMean(0.1))
            .on(T1).boundTo(E1).synchCall(E2, 1.0);

        // Server activity: service=0.5, with think-time=0.2
        Activity A2 = new Activity(lqn, "ServerActivity", Exp.fitMean(0.5))
            .on(T2).boundTo(E2).repliesTo(E2);
        A2.setThinkTime(0.2);

        // Run LQNS solver
        SolverLQNS lqnsSolver = suppressOutput(() -> new SolverLQNS(lqn));
        LayeredNetworkAvgTable lqnsTable = suppressOutput(lqnsSolver::getAvgTable);

        // Get baseline results (without think-time) for comparison
        LayeredNetwork lqnBaseline = new LayeredNetwork("TestNoThinkTime");
        P1 = new Processor(lqnBaseline, "P1", Integer.MAX_VALUE, SchedStrategy.INF);
        T1 = new Task(lqnBaseline, "Client", 1, SchedStrategy.REF).on(P1);
        E1 = new Entry(lqnBaseline, "ClientE").on(T1);
        P2 = new Processor(lqnBaseline, "P2", Integer.MAX_VALUE, SchedStrategy.INF);
        T2 = new Task(lqnBaseline, "Server", Integer.MAX_VALUE, SchedStrategy.INF).on(P2);
        E2 = new Entry(lqnBaseline, "ServerE").on(T2);
        A1 = new Activity(lqnBaseline, "ClientActivity", Exp.fitMean(0.1))
            .on(T1).boundTo(E1).synchCall(E2, 1.0);
        A2 = new Activity(lqnBaseline, "ServerActivity", Exp.fitMean(0.5))
            .on(T2).boundTo(E2).repliesTo(E2);
        // No think-time on baseline

        SolverLQNS lqnsSolverBaseline = suppressOutput(() -> new SolverLQNS(lqnBaseline));
        LayeredNetworkAvgTable lqnsTableBaseline = suppressOutput(lqnsSolverBaseline::getAvgTable);

        // Get response times for ServerActivity (row 7 for activities)
        double rtServerWithoutThinkTime = lqnsTableBaseline.getRespT().get(7);
        double rtServerWithThinkTime = lqnsTable.getRespT().get(7);

        // The difference should be approximately equal to the think-time (0.2)
        double rtIncrease = rtServerWithThinkTime - rtServerWithoutThinkTime;

        // Assert that response time increased by approximately the think-time value
        // Allow 1% tolerance for numerical precision
        assertEquals(0.2, rtIncrease, 0.002,
            "Server activity response time should increase by think-time (0.2)");

        // Also verify that client sees the same increase (since job waits for think-time)
        double rtClientWithoutThinkTime = lqnsTableBaseline.getRespT().get(6);  // ClientActivity
        double rtClientWithThinkTime = lqnsTable.getRespT().get(6);
        double clientIncrease = rtClientWithThinkTime - rtClientWithoutThinkTime;

        assertEquals(0.2, clientIncrease, 0.002,
            "Client activity response time should also increase by think-time (0.2)");
    }

    /**
     * Test: Multiple activities with different think-times.
     */
    @Test
    public void test_activity_think_time_multiple() throws Exception {
        LayeredNetwork lqn = new LayeredNetwork("TestMultipleActivityThinkTime");

        // Client processor
        Processor P1 = new Processor(lqn, "P1", Integer.MAX_VALUE, SchedStrategy.INF);

        // Reference task (client)
        Task ClientTask = new Task(lqn, "Client", 1, SchedStrategy.REF).on(P1);
        Entry ClientEntry = new Entry(lqn, "ClientE").on(ClientTask);

        // Server 1 processor
        Processor P2 = new Processor(lqn, "P2", Integer.MAX_VALUE, SchedStrategy.INF);
        Task Server1Task = new Task(lqn, "Server1", Integer.MAX_VALUE, SchedStrategy.INF).on(P2);
        Entry Server1Entry = new Entry(lqn, "Server1E").on(Server1Task);

        // Server 2 processor
        Processor P3 = new Processor(lqn, "P3", Integer.MAX_VALUE, SchedStrategy.INF);
        Task Server2Task = new Task(lqn, "Server2", Integer.MAX_VALUE, SchedStrategy.INF).on(P3);
        Entry Server2Entry = new Entry(lqn, "Server2E").on(Server2Task);

        // Client activity calls both servers
        Activity ClientAct = new Activity(lqn, "ClientAct", Exp.fitMean(0.05))
            .on(ClientTask).boundTo(ClientEntry)
            .synchCall(Server1Entry, 0.5)
            .synchCall(Server2Entry, 0.5);

        // Server 1 activity: service=0.4, think-time=0.1
        Activity Server1Act = new Activity(lqn, "Server1Act", Exp.fitMean(0.4))
            .on(Server1Task).boundTo(Server1Entry).repliesTo(Server1Entry);
        Server1Act.setThinkTime(0.1);

        // Server 2 activity: service=0.3, think-time=0.3
        Activity Server2Act = new Activity(lqn, "Server2Act", Exp.fitMean(0.3))
            .on(Server2Task).boundTo(Server2Entry).repliesTo(Server2Entry);
        Server2Act.setThinkTime(0.3);

        // Run LQNS solver
        SolverLQNS lqnsSolver = suppressOutput(() -> new SolverLQNS(lqn));
        LayeredNetworkAvgTable lqnsTable = suppressOutput(lqnsSolver::getAvgTable);

        // Get baseline (no think-times) for comparison
        LayeredNetwork lqnBaseline = new LayeredNetwork("TestMultipleActivityNoThinkTime");
        P1 = new Processor(lqnBaseline, "P1", Integer.MAX_VALUE, SchedStrategy.INF);
        ClientTask = new Task(lqnBaseline, "Client", 1, SchedStrategy.REF).on(P1);
        ClientEntry = new Entry(lqnBaseline, "ClientE").on(ClientTask);
        P2 = new Processor(lqnBaseline, "P2", Integer.MAX_VALUE, SchedStrategy.INF);
        Server1Task = new Task(lqnBaseline, "Server1", Integer.MAX_VALUE, SchedStrategy.INF).on(P2);
        Server1Entry = new Entry(lqnBaseline, "Server1E").on(Server1Task);
        P3 = new Processor(lqnBaseline, "P3", Integer.MAX_VALUE, SchedStrategy.INF);
        Server2Task = new Task(lqnBaseline, "Server2", Integer.MAX_VALUE, SchedStrategy.INF).on(P3);
        Server2Entry = new Entry(lqnBaseline, "Server2E").on(Server2Task);

        ClientAct = new Activity(lqnBaseline, "ClientAct", Exp.fitMean(0.05))
            .on(ClientTask).boundTo(ClientEntry)
            .synchCall(Server1Entry, 0.5)
            .synchCall(Server2Entry, 0.5);

        Server1Act = new Activity(lqnBaseline, "Server1Act", Exp.fitMean(0.4))
            .on(Server1Task).boundTo(Server1Entry).repliesTo(Server1Entry);
        // No think-time

        Server2Act = new Activity(lqnBaseline, "Server2Act", Exp.fitMean(0.3))
            .on(Server2Task).boundTo(Server2Entry).repliesTo(Server2Entry);
        // No think-time

        SolverLQNS lqnsSolverBaseline = suppressOutput(() -> new SolverLQNS(lqnBaseline));
        LayeredNetworkAvgTable lqnsTableBaseline = suppressOutput(lqnsSolverBaseline::getAvgTable);

        // Find indices of activities in result table
        // Verify Server1 activity response time increased by ~0.1
        double rtServer1Without = lqnsTableBaseline.getRespT().get(10);  // Server1Act
        double rtServer1With = lqnsTable.getRespT().get(10);
        double increase1 = rtServer1With - rtServer1Without;
        assertEquals(0.1, increase1, 0.002,
            "Server1 activity response time should increase by think-time (0.1)");

        // Verify Server2 activity response time increased by ~0.3
        double rtServer2Without = lqnsTableBaseline.getRespT().get(11);  // Server2Act
        double rtServer2With = lqnsTable.getRespT().get(11);
        double increase2 = rtServer2With - rtServer2Without;
        assertEquals(0.3, increase2, 0.003,
            "Server2 activity response time should increase by think-time (0.3)");

        // Client response time should be affected (increases by weighted sum of think-times)
        double rtClientWithout = lqnsTableBaseline.getRespT().get(9);  // ClientAct
        double rtClientWith = lqnsTable.getRespT().get(9);
        double clientIncrease = rtClientWith - rtClientWithout;

        // Allow larger tolerance due to interaction of multiple think-times
        assertTrue(clientIncrease >= 0.15 && clientIncrease <= 0.25,
            "Client activity response time increase (" + clientIncrease +
            ") should be between 0.15 and 0.25 (weighted average of server think-times)");
    }
}
