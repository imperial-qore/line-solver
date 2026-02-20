package jline.solvers.qns;

import static org.junit.jupiter.api.Assertions.*;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.lang.processes.HyperExp;
import jline.solvers.NetworkAvgTable;
import jline.solvers.SolverOptions;
import jline.solvers.jmt.SolverJMT;
import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.util.Maths;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.condition.EnabledIf;

/**
 * Test suite for SolverQNS - Tests Open and Mixed Queueing Networks
 *
 * DISABLED: SolverQNS implementation is incomplete - qnsolver output does not include
 * Source/Sink stations but the implementation expects them. Requires refactoring to
 * properly handle open queueing networks with Source/Sink nodes.
 */
public class SolverQNSTest {

    @BeforeAll
    public static void setUp() {
        // Ensure MATLAB-compatible random number generation
        Maths.setRandomNumbersMatlab(true);
        // Set verbose level to SILENT to suppress warnings during tests
        GlobalConstants.setVerbose(VerboseLevel.SILENT);
    }

    double tolerance = 1e-2;  // QNS is an approximation method, use relaxed tolerance
    double relaxedTolerance = 5e-2;  // For multiserver approximations

    /**
     * Check if QNS solver is available on the system
     */
    static boolean isQNSAvailable() {
        return SolverQNS.isAvailable();
    }

    // ===== Open Queueing Network (OQN) Tests =====

    /**
     * Test simple M/M/1 queue: Source → Queue → Sink
     * Validates basic open network functionality
     */
    @Test
    @EnabledIf("isQNSAvailable")
    public void testOQN_MM1_Basic() {
        Network model = new Network("M/M/1");

        double lambda = 0.8;  // Arrival rate
        double mu = 1.0;      // Service rate

        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        OpenClass jobClass = new OpenClass(model, "Class1", 0);
        source.setArrival(jobClass, new Exp(lambda));
        queue.setService(jobClass, new Exp(mu));

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobClass, jobClass, source, queue, 1.0);
        P.set(jobClass, jobClass, queue, sink, 1.0);
        model.link(P);

        SolverQNS solver = new SolverQNS(model);
        NetworkAvgTable avgTable = solver.getAvgTable();

        assertNotNull(avgTable, "avgTable should not be null");
        
        if (avgTable.getTput().isEmpty()) {
            // If qnsolver returned 0s or failed to compute valid metrics
            fail("SolverQNS returned empty results (likely qnsolver returned 0s)");
        }

        // Get results for Queue station by name (qnsolver doesn't output Source/Sink)
        NetworkAvgTable queueResults = avgTable.tget("Queue");
        assertNotNull(queueResults, "Queue results should not be null");

        if (queueResults.getTput().isEmpty()) {
            // If qnsolver returned 0s or failed to compute valid metrics
            fail("SolverQNS returned empty results (likely qnsolver returned 0s)");
        }

        // Throughput should equal arrival rate in stable open system
        double tput = queueResults.getTput().get(0);
        assertEquals(lambda, tput, tolerance, "Throughput should equal arrival rate");

        // Utilization check: rho = lambda/mu
        double expectedUtil = lambda / mu;
        double actualUtil = queueResults.getUtil().get(0);
        assertEquals(expectedUtil, actualUtil, tolerance, "Utilization should be lambda/mu");
    }

    /**
     * Test M/M/c multiserver queue with 3 servers
     * Validates multiserver approximation methods
     */
    @Test
    @EnabledIf("isQNSAvailable")
    public void testOQN_MMc_Multiserver() {
        Network model = new Network("M/M/3");

        double lambda = 2.0;
        double mu = 1.0;
        int c = 3;

        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        queue.setNumberOfServers(c);
        Sink sink = new Sink(model, "Sink");

        OpenClass jobClass = new OpenClass(model, "Class1", 0);
        source.setArrival(jobClass, new Exp(lambda));
        queue.setService(jobClass, new Exp(mu));

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobClass, jobClass, source, queue, 1.0);
        P.set(jobClass, jobClass, queue, sink, 1.0);
        model.link(P);

        SolverQNS solver = new SolverQNS(model);
        NetworkAvgTable avgTable = solver.getAvgTable();

        assertNotNull(avgTable, "avgTable should not be null");

        // Get results by station name (qnsolver doesn't output Source/Sink)
        NetworkAvgTable queueResults = avgTable.tget("Queue");
        assertNotNull(queueResults, "Queue results should not be null");

        // Throughput should equal arrival rate
        double tput = queueResults.getTput().get(0);
        assertEquals(lambda, tput, tolerance, "Throughput should equal arrival rate");

        // Utilization: rho = lambda/(c*mu)
        double expectedUtil = lambda / (c * mu);
        double actualUtil = queueResults.getUtil().get(0);
        assertEquals(expectedUtil, actualUtil, relaxedTolerance, "Utilization should be lambda/(c*mu)");
    }

    /**
     * Test tandem network: Source → Queue1 → Queue2 → Sink
     * Validates flow conservation across multiple stations
     */
    @Test
    @EnabledIf("isQNSAvailable")
    public void testOQN_Tandem() {
        Network model = new Network("Tandem");

        double lambda = 0.5;
        double mu1 = 1.0;
        double mu2 = 0.8;

        Source source = new Source(model, "Source");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        OpenClass jobClass = new OpenClass(model, "Class1", 0);
        source.setArrival(jobClass, new Exp(lambda));
        queue1.setService(jobClass, new Exp(mu1));
        queue2.setService(jobClass, new Exp(mu2));

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobClass, jobClass, source, queue1, 1.0);
        P.set(jobClass, jobClass, queue1, queue2, 1.0);
        P.set(jobClass, jobClass, queue2, sink, 1.0);
        model.link(P);

        SolverQNS solver = new SolverQNS(model);
        NetworkAvgTable avgTable = solver.getAvgTable();

        assertNotNull(avgTable, "avgTable should not be null");

        // Get results by station name (qnsolver doesn't output Source/Sink)
        NetworkAvgTable queue1Results = avgTable.tget("Queue1");
        NetworkAvgTable queue2Results = avgTable.tget("Queue2");
        assertNotNull(queue1Results, "Queue1 results should not be null");
        assertNotNull(queue2Results, "Queue2 results should not be null");

        // Both queues should have same throughput (flow conservation)
        double tput1 = queue1Results.getTput().get(0);
        double tput2 = queue2Results.getTput().get(0);
        assertEquals(lambda, tput1, tolerance, "Queue1 throughput should equal arrival rate");
        assertEquals(lambda, tput2, tolerance, "Queue2 throughput should equal arrival rate");

        // Utilization checks
        assertEquals(lambda / mu1, queue1Results.getUtil().get(0), tolerance, "Queue1 utilization");
        assertEquals(lambda / mu2, queue2Results.getUtil().get(0), tolerance, "Queue2 utilization");
    }

    /**
     * Test open network with delay station
     * Validates infinite server (delay) node handling
     */
    @Test
    @EnabledIf("isQNSAvailable")
    public void testOQN_WithDelay() {
        Network model = new Network("OQN_Delay");

        double lambda = 0.6;
        double delayRate = 2.0;
        double queueRate = 1.5;

        Source source = new Source(model, "Source");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue", SchedStrategy.PS);
        Sink sink = new Sink(model, "Sink");

        OpenClass jobClass = new OpenClass(model, "Class1", 0);
        source.setArrival(jobClass, new Exp(lambda));
        delay.setService(jobClass, new Exp(delayRate));
        queue.setService(jobClass, new Exp(queueRate));

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobClass, jobClass, source, delay, 1.0);
        P.set(jobClass, jobClass, delay, queue, 1.0);
        P.set(jobClass, jobClass, queue, sink, 1.0);
        model.link(P);

        SolverQNS solver = new SolverQNS(model);
        NetworkAvgTable avgTable = solver.getAvgTable();

        assertNotNull(avgTable, "avgTable should not be null");

        // Get results by station name (qnsolver doesn't output Source/Sink)
        NetworkAvgTable delayResults = avgTable.tget("Delay");
        NetworkAvgTable queueResults = avgTable.tget("Queue");
        assertNotNull(delayResults, "Delay results should not be null");
        assertNotNull(queueResults, "Queue results should not be null");

        // Throughput checks
        double tputDelay = delayResults.getTput().get(0);
        double tputQueue = queueResults.getTput().get(0);
        assertEquals(lambda, tputDelay, tolerance, "Delay throughput should equal arrival rate");
        assertEquals(lambda, tputQueue, tolerance, "Queue throughput should equal arrival rate");

        // Utilization check for queue
        assertEquals(lambda / queueRate, queueResults.getUtil().get(0), tolerance, "Queue utilization");
    }

    /**
     * Test open network with multiclass
     * Validates handling of multiple job classes
     */
    @Test
    @EnabledIf("isQNSAvailable")
    public void testOQN_Multiclass() {
        Network model = new Network("OQN_Multiclass");

        double lambda1 = 0.3;
        double lambda2 = 0.4;
        double mu1_class1 = 1.0;
        double mu1_class2 = 1.2;

        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.PS);
        Sink sink = new Sink(model, "Sink");

        OpenClass class1 = new OpenClass(model, "Class1", 0);
        OpenClass class2 = new OpenClass(model, "Class2", 0);

        source.setArrival(class1, new Exp(lambda1));
        source.setArrival(class2, new Exp(lambda2));
        queue.setService(class1, new Exp(mu1_class1));
        queue.setService(class2, new Exp(mu1_class2));

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(class1, class1, source, queue, 1.0);
        P.set(class1, class1, queue, sink, 1.0);
        P.set(class2, class2, source, queue, 1.0);
        P.set(class2, class2, queue, sink, 1.0);
        model.link(P);

        SolverQNS solver = new SolverQNS(model);
        NetworkAvgTable avgTable = solver.getAvgTable();

        assertNotNull(avgTable, "avgTable should not be null");

        // Throughput checks for both classes (use named access since qnsolver omits Source/Sink)
        NetworkAvgTable queueResults = avgTable.tget("Queue");
        double tput_class1 = queueResults.getTput().get(0);  // Queue, Class1
        double tput_class2 = queueResults.getTput().get(1);  // Queue, Class2
        assertEquals(lambda1, tput_class1, tolerance, "Class1 throughput");
        assertEquals(lambda2, tput_class2, tolerance, "Class2 throughput");
    }

    // ===== Closed Queueing Network (CQN) Tests =====

    /**
     * Test simple closed network with single class
     * Validates basic closed network functionality
     */
    @Test
    @EnabledIf("isQNSAvailable")
    public void testCQN_Basic() {
        Network model = new Network("CQN_Basic");

        int N = 5;  // Number of jobs

        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);

        ClosedClass jobClass = new ClosedClass(model, "Class1", N, delay, 0);

        delay.setService(jobClass, new Exp(2.0));
        queue.setService(jobClass, new Exp(1.0));

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobClass, jobClass, delay, queue, 1.0);
        P.set(jobClass, jobClass, queue, delay, 1.0);
        model.link(P);

        SolverQNS solver = new SolverQNS(model);
        NetworkAvgTable avgTable = solver.getAvgTable();

        assertNotNull(avgTable, "avgTable should not be null");

        // Throughput should be the same at both stations
        double tputDelay = avgTable.getTput().get(0);
        double tputQueue = avgTable.getTput().get(1);
        assertEquals(tputDelay, tputQueue, tolerance, "Throughput should be equal at both stations");

        // Total queue length should equal N
        double totalQLen = avgTable.getQLen().get(0) + avgTable.getQLen().get(1);
        assertEquals(N, totalQLen, tolerance, "Total queue length should equal N");

        // Throughput should be positive
        assertTrue(tputDelay > 0, "Throughput should be positive");
    }

    /**
     * Test closed network with multiserver queue
     * Validates multiserver approximation in closed networks
     */
    @Test
    @EnabledIf("isQNSAvailable")
    public void testCQN_Multiserver() {
        Network model = new Network("CQN_Multiserver");

        int N = 10;
        int c = 3;

        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        queue.setNumberOfServers(c);

        ClosedClass jobClass = new ClosedClass(model, "Class1", N, delay, 0);

        delay.setService(jobClass, new Exp(1.5));
        queue.setService(jobClass, new Exp(1.0));

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobClass, jobClass, delay, queue, 1.0);
        P.set(jobClass, jobClass, queue, delay, 1.0);
        model.link(P);

        SolverQNS solver = new SolverQNS(model);
        NetworkAvgTable avgTable = solver.getAvgTable();

        assertNotNull(avgTable, "avgTable should not be null");

        // Flow conservation: throughput should be equal at both stations
        double tputDelay = avgTable.getTput().get(0);
        double tputQueue = avgTable.getTput().get(1);
        assertEquals(tputDelay, tputQueue, tolerance, "Throughput conservation");

        // Total jobs should equal N
        double totalQLen = avgTable.getQLen().get(0) + avgTable.getQLen().get(1);
        assertEquals(N, totalQLen, tolerance, "Total jobs should equal N");
    }

    /**
     * Test closed network with multiple queues (tandem)
     * Validates flow conservation in closed tandem networks
     */
    @Test
    @EnabledIf("isQNSAvailable")
    public void testCQN_Tandem() {
        Network model = new Network("CQN_Tandem");

        int N = 8;

        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);
        Queue queue3 = new Queue(model, "Queue3", SchedStrategy.PS);

        ClosedClass jobClass = new ClosedClass(model, "Class1", N, queue1, 0);

        queue1.setService(jobClass, new Exp(2.0));
        queue2.setService(jobClass, new Exp(3.0));
        queue3.setService(jobClass, new Exp(1.5));

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobClass, jobClass, queue1, queue2, 1.0);
        P.set(jobClass, jobClass, queue2, queue3, 1.0);
        P.set(jobClass, jobClass, queue3, queue1, 1.0);
        model.link(P);

        SolverQNS solver = new SolverQNS(model);
        NetworkAvgTable avgTable = solver.getAvgTable();

        assertNotNull(avgTable, "avgTable should not be null");

        // All stations should have the same throughput
        double tput1 = avgTable.getTput().get(0);
        double tput2 = avgTable.getTput().get(1);
        double tput3 = avgTable.getTput().get(2);
        assertEquals(tput1, tput2, tolerance, "Queue1 and Queue2 throughput");
        assertEquals(tput1, tput3, tolerance, "Queue1 and Queue3 throughput");

        // Total queue length should equal N
        double totalQLen = avgTable.getQLen().get(0) + avgTable.getQLen().get(1) + avgTable.getQLen().get(2);
        assertEquals(N, totalQLen, tolerance, "Total jobs should equal N");
    }

    /**
     * Test closed network with multiclass
     * Validates multiple closed classes
     */
    @Test
    @EnabledIf("isQNSAvailable")
    public void testCQN_Multiclass() {
        Network model = new Network("CQN_Multiclass");

        int N1 = 3;
        int N2 = 5;

        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue", SchedStrategy.PS);

        ClosedClass class1 = new ClosedClass(model, "Class1", N1, delay, 0);
        ClosedClass class2 = new ClosedClass(model, "Class2", N2, delay, 0);

        delay.setService(class1, new Exp(1.0));
        delay.setService(class2, new Exp(1.5));
        queue.setService(class1, new Exp(2.0));
        queue.setService(class2, new Exp(2.5));

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(class1, class1, delay, queue, 1.0);
        P.set(class1, class1, queue, delay, 1.0);
        P.set(class2, class2, delay, queue, 1.0);
        P.set(class2, class2, queue, delay, 1.0);
        model.link(P);

        SolverQNS solver = new SolverQNS(model);
        NetworkAvgTable avgTable = solver.getAvgTable();

        assertNotNull(avgTable, "avgTable should not be null");

        // Use named access instead of indices for clarity
        NetworkAvgTable delayClass1 = avgTable.tget("Delay", "Class1");
        NetworkAvgTable delayClass2 = avgTable.tget("Delay", "Class2");
        NetworkAvgTable queueClass1 = avgTable.tget("Queue", "Class1");
        NetworkAvgTable queueClass2 = avgTable.tget("Queue", "Class2");

        double tput1_delay = delayClass1.getTput().get(0);
        double tput1_queue = queueClass1.getTput().get(0);
        double tput2_delay = delayClass2.getTput().get(0);
        double tput2_queue = queueClass2.getTput().get(0);

        assertEquals(tput1_delay, tput1_queue, tolerance, "Class1 throughput conservation");
        assertEquals(tput2_delay, tput2_queue, tolerance, "Class2 throughput conservation");

        // Total jobs per class using named access
        double qlen1_delay = delayClass1.getQLen().get(0);
        double qlen1_queue = queueClass1.getQLen().get(0);
        double qlen2_delay = delayClass2.getQLen().get(0);
        double qlen2_queue = queueClass2.getQLen().get(0);

        double totalQ1 = qlen1_delay + qlen1_queue;
        double totalQ2 = qlen2_delay + qlen2_queue;
        assertEquals(N1, totalQ1, tolerance, "Total jobs for Class1");
        assertEquals(N2, totalQ2, tolerance, "Total jobs for Class2");
    }

    /**
     * Test closed network with complex routing (central server model)
     * Validates probabilistic routing in closed networks
     */
    @Test
    @EnabledIf("isQNSAvailable")
    public void testCQN_CentralServer() {
        Network model = new Network("CQN_CentralServer");

        int N = 6;

        Delay delay = new Delay(model, "Delay");  // Think time
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.FCFS);

        ClosedClass jobClass = new ClosedClass(model, "Class1", N, delay, 0);

        delay.setService(jobClass, new Exp(1.0));
        queue1.setService(jobClass, new Exp(2.0));
        queue2.setService(jobClass, new Exp(3.0));

        // Central server routing: Delay → Queue1 (60%) or Queue2 (40%)
        // Both queues return to Delay
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobClass, jobClass, delay, queue1, 0.6);
        P.set(jobClass, jobClass, delay, queue2, 0.4);
        P.set(jobClass, jobClass, queue1, delay, 1.0);
        P.set(jobClass, jobClass, queue2, delay, 1.0);
        model.link(P);

        SolverQNS solver = new SolverQNS(model);
        NetworkAvgTable avgTable = solver.getAvgTable();

        assertNotNull(avgTable, "avgTable should not be null");

        // Visit ratios
        double tputDelay = avgTable.getTput().get(0);
        double tputQueue1 = avgTable.getTput().get(1);
        double tputQueue2 = avgTable.getTput().get(2);

        // Queue1 should have 60% of delay throughput, Queue2 should have 40%
        assertEquals(0.6 * tputDelay, tputQueue1, relaxedTolerance, "Queue1 throughput ratio");
        assertEquals(0.4 * tputDelay, tputQueue2, relaxedTolerance, "Queue2 throughput ratio");

        // Total jobs should equal N
        double totalQLen = avgTable.getQLen().get(0) + avgTable.getQLen().get(1) + avgTable.getQLen().get(2);
        assertEquals(N, totalQLen, tolerance, "Total jobs should equal N");
    }

    // ===== Mixed Queueing Network (MQN) Tests =====

    /**
     * Test simple mixed network with one closed and one open class
     * Validates basic mixed network functionality
     */
    @Test
    @EnabledIf("isQNSAvailable")
    public void testMQN_Basic() {
        Network model = new Network("MQN_Basic");

        // Nodes
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue", SchedStrategy.PS);
        Source source = new Source(model, "Source");
        Sink sink = new Sink(model, "Sink");

        // Classes
        ClosedClass closedClass = new ClosedClass(model, "ClosedClass", 2, delay, 0);
        OpenClass openClass = new OpenClass(model, "OpenClass", 0);

        // Service definitions
        delay.setService(closedClass, new Exp(2.0));
        delay.setService(openClass, new Exp(3.0));
        queue.setService(closedClass, new Exp(1.5));
        queue.setService(openClass, new Exp(1.0));
        source.setArrival(openClass, new Exp(0.1));

        // Routing
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(closedClass, closedClass, delay, queue, 1.0);
        P.set(closedClass, closedClass, queue, delay, 1.0);
        P.set(openClass, openClass, source, delay, 1.0);
        P.set(openClass, openClass, delay, queue, 1.0);
        P.set(openClass, openClass, queue, sink, 1.0);
        model.link(P);

        SolverQNS solver = new SolverQNS(model);
        NetworkAvgTable avgTable = solver.getAvgTable();

        assertNotNull(avgTable, "avgTable should not be null");

        // Use named access for clarity
        NetworkAvgTable delayOpen = avgTable.tget("Delay", "OpenClass");
        NetworkAvgTable delayClosed = avgTable.tget("Delay", "ClosedClass");

        // Open class throughput should equal arrival rate
        double openTput = delayOpen.getTput().get(0);
        assertEquals(0.1, openTput, tolerance, "Open class throughput should equal arrival rate");

        // Closed class should have positive throughput
        double closedTput = delayClosed.getTput().get(0);
        assertTrue(closedTput > 0, "Closed class should have positive throughput");
    }

    /**
     * Test mixed network with multiserver queues
     * Validates multiserver handling in mixed models
     */
    @Test
    @EnabledIf("isQNSAvailable")
    public void testMQN_Multiserver() {
        Network model = new Network("MQN_Multiserver");

        // Nodes
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.PS);
        queue1.setNumberOfServers(2);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);
        queue2.setNumberOfServers(3);
        Source source = new Source(model, "Source");
        Sink sink = new Sink(model, "Sink");

        // Classes
        ClosedClass closedClass = new ClosedClass(model, "ClosedClass", 3, queue1, 0);
        OpenClass openClass = new OpenClass(model, "OpenClass", 0);

        // Service definitions
        queue1.setService(closedClass, new Exp(1.0));
        queue1.setService(openClass, new Exp(1.5));
        queue2.setService(closedClass, new Exp(2.0));
        queue2.setService(openClass, new Exp(2.5));
        source.setArrival(openClass, new Exp(0.5));

        // Routing
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(closedClass, closedClass, queue1, queue2, 1.0);
        P.set(closedClass, closedClass, queue2, queue1, 1.0);
        P.set(openClass, openClass, source, queue1, 1.0);
        P.set(openClass, openClass, queue1, queue2, 1.0);
        P.set(openClass, openClass, queue2, sink, 1.0);
        model.link(P);

        SolverQNS solver = new SolverQNS(model);
        NetworkAvgTable avgTable = solver.getAvgTable();
        assertNotNull(avgTable, "avgTable should not be null");

        // Open class throughput checks
        double openTput1 = avgTable.getTput().get(1);  // Open class at queue1
        double openTput2 = avgTable.getTput().get(3);  // Open class at queue2
        assertEquals(0.5, openTput1, tolerance, "Open class throughput at Queue1");
        assertEquals(0.5, openTput2, tolerance, "Open class throughput at Queue2");

        // Both classes should have positive throughputs
        assertTrue(avgTable.getTput().get(0) > 0, "Closed class at Queue1");
        assertTrue(avgTable.getTput().get(2) > 0, "Closed class at Queue2");
    }

    /**
     * Test mixed network with complex routing
     * Uses example from MixedModel.mqn_basic()
     */
    @Test
    @EnabledIf("isQNSAvailable")
    public void testMQN_ComplexRouting() {
        Network model = new Network("MQN_Complex");

        // Nodes
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue", SchedStrategy.PS);
        Source source = new Source(model, "Source");
        Sink sink = new Sink(model, "Sink");

        // Classes
        ClosedClass closedClass = new ClosedClass(model, "ClosedClass", 2, delay, 0);
        OpenClass openClass = new OpenClass(model, "OpenClass", 0);

        // Service definitions with different distributions
        delay.setService(closedClass, new Erlang(3, 2));
        delay.setService(openClass, new HyperExp(0.5, 3.0, 10.0));
        queue.setService(closedClass, new HyperExp(0.1, 1.0, 10.0));
        queue.setService(openClass, new Exp(1));
        source.setArrival(openClass, new Exp(0.1));

        // Routing
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(closedClass, closedClass, delay, queue, 1.0);
        P.set(closedClass, closedClass, queue, delay, 1.0);
        P.set(openClass, openClass, source, delay, 1.0);
        P.set(openClass, openClass, delay, queue, 1.0);
        P.set(openClass, openClass, queue, sink, 1.0);
        model.link(P);

        SolverQNS solver = new SolverQNS(model);
        NetworkAvgTable avgTable = solver.getAvgTable();

        assertNotNull(avgTable, "avgTable should not be null");

        // Basic sanity checks
        assertTrue(avgTable.getTput().get(0) > 0, "Closed class throughput at delay");
        assertTrue(avgTable.getTput().get(1) > 0, "Closed class throughput at queue");
        assertTrue(avgTable.getTput().get(2) > 0, "Open class throughput at delay");
        assertTrue(avgTable.getTput().get(3) > 0, "Open class throughput at queue");

        // Queue lengths should be non-negative
        assertTrue(avgTable.getQLen().get(0) >= 0, "Queue length at delay (closed)");
        assertTrue(avgTable.getQLen().get(1) >= 0, "Queue length at queue (closed)");
        assertTrue(avgTable.getQLen().get(2) >= 0, "Queue length at delay (open)");
        assertTrue(avgTable.getQLen().get(3) >= 0, "Queue length at queue (open)");
    }

    /**
     * Test that solver options are properly handled
     */
    @Test
    @EnabledIf("isQNSAvailable")
    public void testSolverOptions() {
        Network model = new Network("Options_Test");

        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        queue.setNumberOfServers(2);
        Sink sink = new Sink(model, "Sink");

        OpenClass jobClass = new OpenClass(model, "Class1", 0);
        source.setArrival(jobClass, new Exp(0.8));
        queue.setService(jobClass, new Exp(1.0));

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobClass, jobClass, source, queue, 1.0);
        P.set(jobClass, jobClass, queue, sink, 1.0);
        model.link(P);

        // Test with different multiserver methods
        SolverOptions options = new SolverOptions();
        options.verbose = VerboseLevel.SILENT;
        options.config.multiserver = "conway";

        SolverQNS solver = new SolverQNS(model, options);
        NetworkAvgTable avgTable = solver.getAvgTable();

        assertNotNull(avgTable, "avgTable should not be null with Conway method");

        // Get results by station name (qnsolver doesn't output Source/Sink)
        NetworkAvgTable queueResults = avgTable.tget("Queue");
        assertNotNull(queueResults, "Queue results should not be null with Conway method");
        assertEquals(0.8, queueResults.getTput().get(0), tolerance, "Throughput with Conway method");
    }
}
