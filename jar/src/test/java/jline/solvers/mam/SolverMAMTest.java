package jline.solvers.mam;

import jline.api.qsys.QsysMapPhResult;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SolverType;
import jline.solvers.NetworkAvgTable;
import jline.solvers.NetworkSolver;
import jline.solvers.SolverOptions;
import jline.lang.OpenClass;
import jline.lang.ClosedClass;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Source;
import jline.lang.nodes.Sink;
import jline.lang.constant.SchedStrategy;
import jline.lang.processes.Distribution;
import jline.lang.processes.Exp;
import jline.lang.processes.MAP;
import jline.lang.processes.ME;
import jline.lang.processes.PH;
import jline.lang.processes.RAP;
import jline.util.matrix.Matrix;
import jline.GlobalConstants;
import jline.VerboseLevel;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

import static jline.api.qsys.Qsys_mapph1Kt.qsys_mapph1;
import static jline.api.qsys.Qsys_mapmap1Kt.qsys_mapmap1;
import static jline.api.qsys.Qsys_mapg1Kt.qsys_mapg1;
import static jline.api.qsys.Qsys_phph1Kt.qsys_phph1;
import static jline.api.qsys.Qsys_mapmcKt.qsys_mapmc;
import static jline.api.qsys.Qsys_mapmcKt.qsys_mapm1;
import static jline.examples.java.basic.ClosedModel.*;
import static jline.examples.java.basic.MixedModel.*;
import static jline.examples.java.basic.OpenModel.*;

/**
 * Manual test demonstrating the usage of the Matrix Analytic Methods (MAM) solver.
 * 
 * <p>This example shows how to:
 * <ul>
 *   <li>Use SolverMAM with different network types (open, closed, mixed)</li>
 *   <li>Apply the MNA (Matrix Normalizing Algorithm) method</li>
 *   <li>Solve various queueing network models</li>
 *   <li>Extract and display performance metrics</li>
 * </ul>
 * 
 * @see SolverMAM
 * @see MAMOptions
 * @see NetworkAvgTable
 */
public class SolverMAMTest {

    @BeforeAll
    public static void setUpVerbosity() {
        GlobalConstants.setVerbose(VerboseLevel.SILENT);
    }

    /**
     * Main method to run all MAM solver demonstrations.
     *
     * @param args command line arguments (not used)
     */
    public static void main(String[] args) {
        runAll();
    }

    /**
     * Runs all test cases demonstrating different MAM solver applications.
     * 
     * <p>Executes tests for:
     * <ul>
     *   <li>MNA method on closed and open models</li>
     *   <li>Mixed open/closed networks</li>
     *   <li>Various closed network configurations</li>
     *   <li>Open network examples</li>
     * </ul>
     */
    public static void runAll() {
        test_mna_1();
        test_mna_2();
        test_mam_mixed_1();
        test_mam_closed_1();
        test_mam_open_1();
    }

    /**
     * Demonstrates MNA (Matrix Normalizing Algorithm) on a closed model.
     * 
     * <p>Uses MAMOptions to specify the MNA method for solving
     * a closed queueing network example.
     */
    public static void test_mna_1() {
        Network model1 = cqn_repairmen();

        NetworkSolver solver1 = new SolverMAM(model1, new MAMOptions().method("mna"));
        //solver1.getAvgTable().print();
    }

    /**
     * Demonstrates MNA method on an open queueing network.
     * 
     * <p>Shows how to:
     * <ul>
     *   <li>Configure SolverOptions for MAM solver type</li>
     *   <li>Specify MNA as the solution method</li>
     *   <li>Retrieve and print average performance metrics</li>
     * </ul>
     */
    public static void test_mna_2() {

        Network model1 = oqn_cs_routing();

        SolverOptions options = new SolverOptions(SolverType.MAM);
        options.method = "mna";

        NetworkSolver solver1 = new SolverMAM(model1, options);
        NetworkAvgTable t1 = solver1.getAvgTable();
        //t1.print(options);
    }

    /**
     * Demonstrates MAM solver on mixed open/closed networks.
     * 
     * <p>Solves multiple mixed network examples and displays:
     * <ul>
     *   <li>Average performance metrics for each model</li>
     *   <li>Runtime statistics for performance comparison</li>
     * </ul>
     */
    public static void test_mam_mixed_1() {


        Network model1 = mqn_basic();
        Network model2 = mqn_multiserver_ps();
        Network model3 = mqn_multiserver_fcfs();
        Network model4 = mqn_singleserver_fcfs();
        Network model5 = mqn_singleserver_ps();


        SolverMAM solver1 = new SolverMAM(model1);
        //solver1.getAvgTable().print(new SolverOptions());

        SolverMAM solver2 = new SolverMAM(model2);
        //solver2.getAvgTable().print(new SolverOptions());

        SolverMAM solver3 = new SolverMAM(model3);
        //solver3.getAvgTable().print(new SolverOptions());

        SolverMAM solver4 = new SolverMAM(model4);
        //solver4.getAvgTable().print(new SolverOptions());

        SolverMAM solver5 = new SolverMAM(model5);
        //solver5.getAvgTable().print(new SolverOptions());
    }

    /**
     * Demonstrates MAM solver on various closed network configurations.
     * 
     * <p>Tests include:
     * <ul>
     *   <li>Simple closed networks</li>
     *   <li>Multi-class closed networks</li>
     *   <li>Networks with complex routing</li>
     * </ul>
     * 
     * <p>Prints performance metrics and runtime for each model.
     */
    public static void test_mam_closed_1() {

        Network model1 = cqn_repairmen();
        SolverMAM solver1 = new SolverMAM(model1);
        //solver1.getAvgTable().print(new SolverOptions());

        Network model2 = cqn_twoclass_hyperl();
        SolverMAM solver2 = new SolverMAM(model2);
        //solver2.getAvgTable().print(new SolverOptions());

        Network model3 = cqn_threeclass_hyperl();
        SolverMAM solver3 = new SolverMAM(model3);
        //solver3.getAvgTable().print(new SolverOptions());

        Network model4 = cqn_oneline();
        SolverMAM solver4 = new SolverMAM(model4);
        //solver4.getAvgTable().print(new SolverOptions());

        Network model5 = cqn_repairmen_multi();
        SolverMAM solver5 = new SolverMAM(model5);
        //solver5.getAvgTable().print(new SolverOptions());

        Network model6 = cqn_twoqueues_multi();
        SolverMAM solver6 = new SolverMAM(model6);
        //solver6.getAvgTable().print(new SolverOptions());

    }

    /**
     * Demonstrates MAM solver on various open network models.
     * 
     * <p>Includes examples of:
     * <ul>
     *   <li>Single-class open networks</li>
     *   <li>Multi-class open networks</li>
     *   <li>Networks with different scheduling strategies</li>
     * </ul>
     * 
     * <p>Displays average table results and solver runtime.
     */
    public static void test_mam_open_1() {

        Network model1 = oqn_basic();
        Network model2 = oqn_oneline();
        Network model3 = oqn_cs_routing();
        Network model4 = oqn_vsinks();
        Network model5 = oqn_fourqueues();

        SolverMAM solver1 = new SolverMAM(model1);
        //solver1.getAvgTable().print(new SolverOptions());

        SolverMAM solver2 = new SolverMAM(model2);
        //solver2.getAvgTable().print(new SolverOptions());

        SolverMAM solver3 = new SolverMAM(model3);
        //solver3.getAvgTable().print(new SolverOptions());

        SolverMAM solver4 = new SolverMAM(model4);
        //solver4.getAvgTable().print(new SolverOptions());

        SolverMAM solver5 = new SolverMAM(model5);
        //solver5.getAvgTable().print(new SolverOptions());
    }

    @Test
    public void testMAMAcceptsMixedModel() {
        // Create a mixed model with both open and closed classes
        Network model = new Network("MixedModel");

        // Add nodes
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        Delay delay = new Delay(model, "Delay");

        // Add an open class
        OpenClass openClass = new OpenClass(model, "OpenClass");
        source.setArrival(openClass, new Exp(1.0));
        queue.setService(openClass, new Exp(2.0));

        // Add a closed class
        ClosedClass closedClass = new ClosedClass(model, "ClosedClass", 5, delay);
        queue.setService(closedClass, new Exp(3.0));
        delay.setService(closedClass, new Exp(1.0));

        // Link nodes using a single combined routing matrix
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(openClass, openClass, source, queue, 1.0);
        P.set(openClass, openClass, queue, sink, 1.0);
        P.set(closedClass, closedClass, delay, queue, 1.0);
        P.set(closedClass, closedClass, queue, delay, 1.0);
        model.link(P);

        // MAM now supports mixed models via the dec.source method
        SolverMAM solver = new SolverMAM(model);
        NetworkAvgTable avgTable = solver.getAvgTable();
        assertNotNull(avgTable, "MAM should produce results for mixed models");
    }
    
    @Test
    public void testMAMAcceptsOpenModel() {
        // Create an open-only model
        Network model = new Network("OpenModel");

        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        OpenClass openClass = new OpenClass(model, "OpenClass");
        source.setArrival(openClass, new Exp(1.0));
        queue.setService(openClass, new Exp(2.0));

        model.link(model.serialRouting(source, queue, sink));

        // Create solver - should not throw for open-only model
        SolverMAM solver = new SolverMAM(model);

        // This should not throw an exception for open-only models
        assertDoesNotThrow(() -> {
            // Note: The actual analysis might fail for other reasons,
            // but it should not fail due to mixed model check
            try {
                solver.runAnalyzer();
            } catch (Exception e) {
                // Make sure it's not the mixed model error
                assertFalse(e.getMessage().contains("mixed models"),
                           "Should not fail due to mixed model check: " + e.getMessage());
                // Re-throw to show the actual error (if any)
                throw e;
            }
        });
    }

    // ===== Qsys MAP/PH Queue Analyzer Tests =====

    private static final double QSYS_TOLERANCE = 1e-4;

    /**
     * Creates exponential (1-phase PH) representation.
     */
    private Matrix[] createExponentialPH(double rate) {
        Matrix alpha = new Matrix(1, 1);
        alpha.set(0, 0, 1.0);
        Matrix T = new Matrix(1, 1);
        T.set(0, 0, -rate);
        return new Matrix[]{alpha, T};
    }

    /**
     * Creates Poisson (exponential inter-arrival) MAP representation.
     */
    private Matrix[] createPoissonMAP(double lambda) {
        Matrix D0 = new Matrix(1, 1);
        D0.set(0, 0, -lambda);
        Matrix D1 = new Matrix(1, 1);
        D1.set(0, 0, lambda);
        return new Matrix[]{D0, D1};
    }

    @Test
    public void testQsysMapPh1_MM1() {
        // Test MAP/PH/1 with exponential arrival and service (= M/M/1)
        double lambda = 0.8;
        double mu = 1.0;
        double rho = lambda / mu;

        Matrix[] arrival = createPoissonMAP(lambda);
        Matrix[] service = createExponentialPH(mu);

        QsysMapPhResult result = qsys_mapph1(arrival[0], arrival[1], service[0], service[1]);

        // M/M/1 analytical formula for mean queue length
        double expectedQL = rho / (1 - rho);

        assertEquals(expectedQL, result.getMeanQueueLength(), QSYS_TOLERANCE * expectedQL,
                "Mean queue length should match M/M/1 formula");
        assertEquals(rho, result.getUtilization(), QSYS_TOLERANCE,
                "Utilization should match rho");
        assertEquals("BUTools:MMAPPH1FCFS", result.getAnalyzer(),
                "Should use BUTools analyzer");
    }

    @Test
    public void testQsysMapMap1_MM1() {
        // Test MAP/MAP/1 with exponential arrival and service
        double lambda = 0.6;
        double mu = 1.0;
        double rho = lambda / mu;

        Matrix[] arrival = createPoissonMAP(lambda);
        Matrix[] service = createPoissonMAP(mu);

        QsysMapPhResult result = qsys_mapmap1(arrival[0], arrival[1], service[0], service[1]);

        double expectedQL = rho / (1 - rho);

        assertEquals(expectedQL, result.getMeanQueueLength(), QSYS_TOLERANCE * expectedQL,
                "Mean queue length should match M/M/1 formula");
    }

    @Test
    public void testQsysPhPh1_MM1() {
        // Test PH/PH/1 with exponential arrival and service
        double lambda = 0.75;
        double mu = 1.0;
        double rho = lambda / mu;

        Matrix[] arrival = createExponentialPH(lambda);
        Matrix[] service = createExponentialPH(mu);

        QsysMapPhResult result = qsys_phph1(arrival[0], arrival[1], service[0], service[1]);

        double expectedQL = rho / (1 - rho);

        assertEquals(expectedQL, result.getMeanQueueLength(), QSYS_TOLERANCE * expectedQL,
                "Mean queue length should match M/M/1 formula");
    }

    @Test
    public void testQsysMapMc_MM2() {
        // Test MAP/M/2 with Poisson arrivals
        double lambda = 1.5;
        double mu = 1.0;
        int c = 2;
        double rho = lambda / (mu * c);

        Matrix[] arrival = createPoissonMAP(lambda);

        QsysMapPhResult result = qsys_mapmc(arrival[0], arrival[1], mu, c);

        assertTrue(result.getUtilization() > 0 && result.getUtilization() < 1,
                "Utilization should be between 0 and 1: " + result.getUtilization());
        assertEquals(rho, result.getUtilization(), QSYS_TOLERANCE,
                "Utilization should match rho");
        assertEquals("Q-MAM:MAP/M/2", result.getAnalyzer(),
                "Should use Q-MAM analyzer for multi-server");
    }

    @Test
    public void testQsysMapG1_MM1() {
        // Test MAP/G/1 with exponential service (cv = 1)
        double lambda = 0.7;
        double mu = 1.0;
        double rho = lambda / mu;

        Matrix[] arrival = createPoissonMAP(lambda);
        double meanService = 1.0 / mu;
        double cv = 1.0;

        QsysMapPhResult result = qsys_mapg1(arrival[0], arrival[1], meanService, cv);

        double expectedQL = rho / (1 - rho);

        assertEquals(expectedQL, result.getMeanQueueLength(), 2 * QSYS_TOLERANCE * expectedQL,
                "Mean queue length should approximately match M/M/1 formula");
    }

    @Test
    public void testQsysMapPh1_TwoPhaseMAP() {
        // Test with a 2-phase MMPP
        Matrix D0 = new Matrix(2, 2);
        D0.set(0, 0, -2.0);
        D0.set(0, 1, 1.0);
        D0.set(1, 0, 1.0);
        D0.set(1, 1, -3.0);

        Matrix D1 = new Matrix(2, 2);
        D1.set(0, 0, 1.0);
        D1.set(0, 1, 0.0);
        D1.set(1, 0, 0.0);
        D1.set(1, 1, 2.0);

        Matrix[] service = createExponentialPH(2.0);

        QsysMapPhResult result = qsys_mapph1(D0, D1, service[0], service[1]);

        assertNotNull(result);
        assertTrue(result.getMeanQueueLength() > 0,
                "Mean queue length should be positive");
        assertTrue(result.getUtilization() > 0 && result.getUtilization() < 1,
                "Utilization should be between 0 and 1");
    }

    // ==================== ME/RAP Tests (from SolverMAMMETest) ====================

    private static final double TOLERANCE = 1e-6;
    private static final double LOOSE_TOLERANCE = 1e-3;

    // ==================== MAP/ME/1 Tests ====================

    @Test
    public void test_map_me_exp_1() {
        // Test MAP/ME/1 where ME is equivalent to Exp
        // This should work seamlessly with MAM solver

        Network model = new Network("MAP/ME(Exp)/1");

        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        OpenClass jobClass = new OpenClass(model, "Class1", 0);

        // MAP arrivals (simple 2-state MAP)
        Matrix D0 = new Matrix(new double[][]{
                {-2.0, 1.0},
                {0.5, -1.5}
        });
        Matrix D1 = new Matrix(new double[][]{
                {0.5, 0.5},
                {0.5, 0.5}
        });
        MAP mapArrival = new MAP(D0, D1);
        source.setArrival(jobClass, mapArrival);

        // ME service from Exp(2.0)
        ME meService = ME.fromExp(2.0);
        queue.setService(jobClass, meService);

        // Routing
        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        routingMatrix.addConnection(jobClass, jobClass, source, queue, 1.0);
        routingMatrix.addConnection(jobClass, jobClass, queue, sink, 1.0);
        model.link(routingMatrix);

        // Solve with MAM
        SolverOptions options = new SolverOptions();
        options.verbose = VerboseLevel.SILENT;

        SolverMAM solverMAM = new SolverMAM(model, options);
        NetworkAvgTable avgTableMAM = solverMAM.getAvgTable();

        assertNotNull(avgTableMAM);

        // Verify ME service mean
        assertEquals(0.5, meService.getMean(), TOLERANCE);

        // Verify solver completed successfully
        double qLen = getQueueLengthForClass(avgTableMAM, queue, jobClass);
        assertTrue(qLen > 0, "Queue length should be positive");
    }

    @Test
    public void test_map_me_erlang_1() {
        // Test MAP/ME/1 where ME is equivalent to Erlang-2

        Network model = new Network("MAP/ME(Erlang)/1");

        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        OpenClass jobClass = new OpenClass(model, "Class1", 0);

        // Poisson arrivals (special case of MAP)
        Matrix D0_poisson = new Matrix(new double[][]{{-0.8}});
        Matrix D1_poisson = new Matrix(new double[][]{{0.8}});
        MAP mapArrival = new MAP(D0_poisson, D1_poisson);
        source.setArrival(jobClass, mapArrival);

        // ME service from Erlang-2 with rate 2.0
        // Mean = 2/2 = 1.0, SCV = 1/2 = 0.5
        ME meService = ME.fromErlang(2, 2.0);
        queue.setService(jobClass, meService);

        // Routing
        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        routingMatrix.addConnection(jobClass, jobClass, source, queue, 1.0);
        routingMatrix.addConnection(jobClass, jobClass, queue, sink, 1.0);
        model.link(routingMatrix);

        // Solve with MAM
        SolverOptions options = new SolverOptions();
        options.verbose = VerboseLevel.SILENT;

        SolverMAM solverMAM = new SolverMAM(model, options);
        NetworkAvgTable avgTableMAM = solverMAM.getAvgTable();

        assertNotNull(avgTableMAM);

        // Verify ME moments
        assertEquals(1.0, meService.getMean(), TOLERANCE);
        assertEquals(0.5, meService.getSCV(), TOLERANCE);

        // Verify utilization: ρ = λ * E[S] = 0.8 * 1.0 = 0.8
        double util = getUtilizationForClass(avgTableMAM, queue, jobClass);
        assertEquals(0.8, util, LOOSE_TOLERANCE);
    }

    // ==================== RAP/PH/1 Tests ====================

    @Test
    public void test_rap_poisson_ph_1() {
        // Test RAP/PH/1 where RAP is equivalent to Poisson

        Network model = new Network("RAP(Poisson)/PH/1");

        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        OpenClass jobClass = new OpenClass(model, "Class1", 0);

        // RAP arrivals from Poisson with rate 0.8
        RAP rapArrival = RAP.fromPoisson(0.8);
        source.setArrival(jobClass, rapArrival);

        // PH service (Erlang-2 with rate 2.0)
        Matrix alpha = new Matrix(new double[]{1.0, 0.0});
        Matrix A = new Matrix(new double[][]{
                {-2.0, 2.0},
                {0.0, -2.0}
        });
        PH phService = new PH(alpha, A);
        queue.setService(jobClass, phService);

        // Routing
        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        routingMatrix.addConnection(jobClass, jobClass, source, queue, 1.0);
        routingMatrix.addConnection(jobClass, jobClass, queue, sink, 1.0);
        model.link(routingMatrix);

        // Solve with MAM
        SolverOptions options = new SolverOptions();
        options.verbose = VerboseLevel.SILENT;

        SolverMAM solverMAM = new SolverMAM(model, options);
        NetworkAvgTable avgTableMAM = solverMAM.getAvgTable();

        assertNotNull(avgTableMAM);

        // Verify RAP mean
        assertEquals(1.25, rapArrival.getMean(), TOLERANCE);

        // Verify solver completed successfully
        double qLen = getQueueLengthForClass(avgTableMAM, queue, jobClass);
        assertTrue(qLen > 0, "Queue length should be positive");
    }

    @Test
    public void test_rap_erlang_ph_1() {
        // Test RAP/PH/1 where RAP is equivalent to Erlang renewal

        Network model = new Network("RAP(Erlang)/PH/1");

        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        OpenClass jobClass = new OpenClass(model, "Class1", 0);

        // RAP arrivals from Erlang-2 with rate 1.0
        RAP rapArrival = RAP.fromErlang(2, 1.0);
        source.setArrival(jobClass, rapArrival);

        // Exponential service (mean = 1.0)
        queue.setService(jobClass, Exp.fitRate(1.0));

        // Routing
        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        routingMatrix.addConnection(jobClass, jobClass, source, queue, 1.0);
        routingMatrix.addConnection(jobClass, jobClass, queue, sink, 1.0);
        model.link(routingMatrix);

        // Solve with MAM
        SolverOptions options = new SolverOptions();
        options.verbose = VerboseLevel.SILENT;

        SolverMAM solverMAM = new SolverMAM(model, options);
        NetworkAvgTable avgTableMAM = solverMAM.getAvgTable();

        assertNotNull(avgTableMAM);

        // Verify RAP moments
        assertEquals(2.0, rapArrival.getMean(), TOLERANCE);
        assertEquals(0.5, rapArrival.getSCV(), TOLERANCE);

        // Arrival rate = 0.5, service rate = 1.0
        // Utilization = 0.5
        double util = getUtilizationForClass(avgTableMAM, queue, jobClass);
        assertEquals(0.5, util, LOOSE_TOLERANCE);
    }

    // ==================== ME/ME/1 Test ====================

    @Test
    public void test_me_me_1() {
        // Test ME/ME/1 queue with MAM solver

        Network model = new Network("ME/ME/1");

        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        OpenClass jobClass = new OpenClass(model, "Class1", 0);

        // ME arrivals (Exp, rate 0.8)
        ME meArrival = ME.fromExp(0.8);
        source.setArrival(jobClass, meArrival);

        // ME service (Exp, rate 2.0)
        ME meService = ME.fromExp(2.0);
        queue.setService(jobClass, meService);

        // Routing
        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        routingMatrix.addConnection(jobClass, jobClass, source, queue, 1.0);
        routingMatrix.addConnection(jobClass, jobClass, queue, sink, 1.0);
        model.link(routingMatrix);

        // Solve with MAM
        SolverOptions options = new SolverOptions();
        options.verbose = VerboseLevel.SILENT;

        SolverMAM solverMAM = new SolverMAM(model, options);
        NetworkAvgTable avgTableMAM = solverMAM.getAvgTable();

        assertNotNull(avgTableMAM);

        // Verify moments
        assertEquals(1.25, meArrival.getMean(), TOLERANCE);
        assertEquals(0.5, meService.getMean(), TOLERANCE);

        // Utilization = 0.8 / 2.0 = 0.4
        double util = getUtilizationForClass(avgTableMAM, queue, jobClass);
        assertEquals(0.4, util, LOOSE_TOLERANCE);
    }

    // ==================== Validation Against qsys_mapph1 ====================

    @Test
    public void test_me_vs_qsys_mapph1() {
        // Validate ME/PH/1 queue against qsys_mapph1 API function

        // Create MAP arrivals (Poisson)
        Matrix D0 = new Matrix(new double[][]{{-0.8}});
        Matrix D1 = new Matrix(new double[][]{{0.8}});

        // Create ME service (Erlang-2)
        Matrix alpha = new Matrix(new double[]{1.0, 0.0});
        Matrix A = new Matrix(new double[][]{
                {-2.0, 2.0},
                {0.0, -2.0}
        });

        // Compute using qsys_mapph1
        QsysMapPhResult qsysResult = qsys_mapph1(D0, D1, alpha, A);

        // Create network model
        Network model = new Network("MAP/ME/1");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass jobClass = new OpenClass(model, "Class1", 0);

        source.setArrival(jobClass, new MAP(D0, D1));
        queue.setService(jobClass, new ME(alpha, A));

        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        routingMatrix.addConnection(jobClass, jobClass, source, queue, 1.0);
        routingMatrix.addConnection(jobClass, jobClass, queue, sink, 1.0);
        model.link(routingMatrix);

        // Solve with MAM
        SolverMAM solverMAM = new SolverMAM(model, new SolverOptions());
        NetworkAvgTable avgTableMAM = solverMAM.getAvgTable();

        // Compare results
        double qLenQsys = qsysResult.getMeanQueueLength();
        double qLenMAM = getQueueLengthForClass(avgTableMAM, queue, jobClass);

        assertEquals(qLenQsys, qLenMAM, LOOSE_TOLERANCE);
    }

    @Test
    public void test_rap_from_map_consistency() {
        // Verify RAP.fromMAP produces identical results to original MAP

        // Create a MAP
        Matrix D0 = new Matrix(new double[][]{
                {-2.0, 1.0},
                {0.5, -1.5}
        });
        Matrix D1 = new Matrix(new double[][]{
                {0.5, 0.5},
                {0.5, 0.5}
        });

        MAP mapArrival = new MAP(D0, D1);
        RAP rapArrival = RAP.fromMAP(mapArrival);

        // Create two identical networks, one with MAP, one with RAP
        Network model1 = createMapPhModelForME(mapArrival);
        Network model2 = createMapPhModelForME(rapArrival);

        // Solve both
        SolverMAM solver1 = new SolverMAM(model1, new SolverOptions());
        SolverMAM solver2 = new SolverMAM(model2, new SolverOptions());

        NetworkAvgTable table1 = solver1.getAvgTable();
        NetworkAvgTable table2 = solver2.getAvgTable();

        // Results should be identical
        Queue queue1 = (Queue) model1.getNodes().get(1);
        Queue queue2 = (Queue) model2.getNodes().get(1);
        OpenClass class1 = (OpenClass) model1.getClasses().get(0);
        OpenClass class2 = (OpenClass) model2.getClasses().get(0);

        double qLen1 = getQueueLengthForClass(table1, queue1, class1);
        double qLen2 = getQueueLengthForClass(table2, queue2, class2);

        assertEquals(qLen1, qLen2, TOLERANCE);
    }

    // ==================== Helper Methods ====================

    private Network createMapPhModelForME(Distribution arrivalDist) {
        Network model = new Network("MAP/PH/1");

        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        OpenClass jobClass = new OpenClass(model, "Class1", 0);

        source.setArrival(jobClass, arrivalDist);
        queue.setService(jobClass, Exp.fitRate(2.0));

        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        routingMatrix.addConnection(jobClass, jobClass, source, queue, 1.0);
        routingMatrix.addConnection(jobClass, jobClass, queue, sink, 1.0);
        model.link(routingMatrix);

        return model;
    }

    private double getQueueLengthForClass(NetworkAvgTable table, Queue queue, OpenClass jobClass) {
        for (int i = 0; i < table.getStationNames().size(); i++) {
            if (table.getStationNames().get(i).equals(queue.getName()) &&
                table.getClassNames().get(i).equals(jobClass.getName())) {
                return table.getQLen().get(i);
            }
        }
        return 0.0;
    }

    private double getUtilizationForClass(NetworkAvgTable table, Queue queue, OpenClass jobClass) {
        for (int i = 0; i < table.getStationNames().size(); i++) {
            if (table.getStationNames().get(i).equals(queue.getName()) &&
                table.getClassNames().get(i).equals(jobClass.getName())) {
                return table.getUtil().get(i);
            }
        }
        return 0.0;
    }
}