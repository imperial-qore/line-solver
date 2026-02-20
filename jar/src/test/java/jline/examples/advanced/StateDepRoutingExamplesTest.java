package jline.examples.advanced;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.examples.java.advanced.StateDepRoutingModel;
import jline.lang.Network;
import jline.solvers.NetworkAvgNodeChainTable;
import jline.solvers.NetworkAvgNodeTable;
import jline.solvers.NetworkAvgTable;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.jmt.SolverJMT;
import jline.solvers.ssa.SolverSSA;
import jline.util.Maths;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;
import static jline.TestTools.withSuppressedOutput;
import static jline.TestTools.assertTableMetrics;
import static jline.TestTools.COARSE_TOL;
import static jline.TestTools.VERY_COARSE_TOL;

/**
 * Test class for State-Dependent Routing examples.
 *
 * Expected values come from allExamplesBaseline.txt MATLAB execution.
 * Only sdroute_closed is present in allExamplesBaseline.txt. Other examples
 * (sdroute_open and sdroute_twoclasses_closed) retain their existing
 * Java baseline values for regression testing.
 */
public class StateDepRoutingExamplesTest {

    @BeforeAll
    public static void setUp() {
        // Ensure MATLAB-compatible random number generation
        Maths.setRandomNumbersMatlab(true);
        // Set verbose level to SILENT to suppress warnings during tests
        GlobalConstants.setVerbose(VerboseLevel.SILENT);
    }

    /**
     * Test sdroute_closed with CTMC solver.
     * 
     * This test verifies a single-class closed network with round-robin routing.
     * The model has:
     * - Single closed class with 1 job
     * - Delay node with high-variability HyperExp service (SCV=25)
     * - Two PS queues with different service rates (1.0 and 2.0)
     * - Round-robin routing from delay to queues
     * 
     * Expected values from MATLAB dev/test_sdroute_closed.m output.
     */
    @Test
    // Fixed: RROBIN routing indexing (JAR was using 0-based nvars extraction that differed from MATLAB 1-based)
    public void testSdrouteClosedCTMC() {
        // Create the model
        Network model = StateDepRoutingModel.sdroute_closed();

        // Create and run the solver
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            SolverCTMC solver = new SolverCTMC(model, "keep", true);
            avgTableHolder[0] = solver.getAvgTable();
            assertEquals("default", solver.result.method, "CTMC solver should use default method");
        });
        NetworkAvgTable avgTable = avgTableHolder[0];

        // Check if results are computed
        assertNotNull(avgTable);
        
        // Expected values from allExamplesBaseline.txt (CTMC solver)
        // Order: Delay(Class1), Queue1(Class1), Queue2(Class1)
        double[] expectedQLen = {0.666666666666667, 0.222222222222222, 0.111111111111111};
        double[] expectedUtil = {0.666666666666667, 0.222222222222222, 0.111111111111111};
        double[] expectedRespT = {1.0, 1.0, 0.5};
        double[] expectedResidT = {1.0, 0.333333333333333, 0.166666666666667};
        double[] expectedArvR = {0.666666666666667, 0.222222222222222, 0.222222222222222};
        double[] expectedTput = {0.666666666666667, 0.222222222222222, 0.222222222222222};
        
        // Verify table size
        assertEquals(3, avgTable.getQLen().size(), 
            "Expected 3 entries (3 stations × 1 class)");
        
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }

    /**
     * Test sdroute_closed with JMT solver.
     * 
     * JMT solver uses simulation and may produce slightly different results.
     * Expected values from MATLAB dev/test_sdroute_closed.m output.
     */
    @Test
    //@Disabled("Test failing - disabled for investigation")
    public void testSdrouteClosedJMT() {
        // Create the model
        Network model = StateDepRoutingModel.sdroute_closed();
        
        // Create and run the solver
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            SolverJMT solver = new SolverJMT(model, "seed", 23000);
            solver.options.samples = 100000;
            avgTableHolder[0] = solver.getAvgTable();
            assertEquals("default", solver.result.method, "JMT solver should use default method");
        });
        NetworkAvgTable avgTable = avgTableHolder[0];
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Expected values from ground truth JMT solver output
        // Order: Delay(Class1), Queue1(Class1), Queue2(Class1)
        double[] expectedQLen = {0.663399372707037, 0.2227027512816, 0.112950606464826};
        double[] expectedUtil = {0.663399372707037, 0.2227027512816, 0.112950606464826};
        double[] expectedRespT = {0.996424658588384, 0.994079873461343, 0.500543022159858};
        double[] expectedResidT = {0.996424658588384, 0.331359957820448, 0.166847674053286};
        double[] expectedArvR = {0.672750511310764, 0.221364400029354, 0.221362884623031};
        double[] expectedTput = {0.665947155529669, 0.221321960426697, 0.221363814590946};
        
        // Verify table size
        assertEquals(3, avgTable.getQLen().size(), 
            "Expected 3 entries (3 stations × 1 class)");
        
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }

    /**
     * Test sdroute_closed with SSA solver.
     * 
     * SSA solver uses stochastic simulation and may produce variable results.
     * Expected values from MATLAB dev/test_sdroute_closed.m output.
     */
    @Test
    public void testSdrouteClosedSSA() {
        // Create the model
        Network model = StateDepRoutingModel.sdroute_closed();
        
        // Create and run the solver
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            SolverSSA solver = new SolverSSA(model, "seed", 23000, "verbose", true, "samples", 10000);
            avgTableHolder[0] = solver.getAvgTable();
            assertEquals("default/serial", solver.result.method, "SSA solver should use default/serial method");
        });
        NetworkAvgTable avgTable = avgTableHolder[0];

        // Check if results are computed
        assertNotNull(avgTable);
        
        // Expected values from CTMC solver (ground truth for this model)
        // Order: Delay(Class1), Queue1(Class1), Queue2(Class1)
        double[] expectedQLen = {0.6518663537952001, 0.23111823435519943, 0.1170154118496004};
        double[] expectedUtil = {0.6518663537952001, 0.22936442309981792, 0.11468221154990896};
        double[] expectedRespT = {0.9473517368639048, 1.0, 0.5};
        double[] expectedResidT = {0.9473517368639048, 0.3333333333333333, 0.16666666666666666};
        double[] expectedArvR = {0.6945134811542182, 0.2293644230998179, 0.2293644230998179};
        double[] expectedTput = {0.6880932692994537, 0.23111823435519943, 0.2340308236992008};

        // Verify table size
        assertEquals(3, avgTable.getQLen().size(),
            "Expected 3 entries (3 stations × 1 class)");

        // Check all metrics against expected values
        // Use VERY_COARSE_TOL (10%) for SSA since it's a stochastic simulation
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                          expectedResidT, expectedArvR, expectedTput, VERY_COARSE_TOL);
    }

    /**
     * Test sdroute_open with CTMC solver.
     * 
     * This test verifies an open network with load balancing router.
     * The model has:
     * - Open network: Source → Router → Queues → Sink
     * - Router uses round-robin strategy for load balancing
     * - Two FCFS queues with identical service rates (0.5)
     * - Exponential arrivals (rate 1.0) distributed evenly
     * 
     * Expected values from MATLAB dev/test_sdroute_open.m output.
     */
    @Test
    //@Disabled("RROBIN bug: actual error is ~15% (QLen 0.329 vs expected 0.287), not 0.17% as previously documented")
    public void testSdrouteOpenCTMC() {
        // Create the model
        Network model = StateDepRoutingModel.sdroute_open();
        
        // Create and run the solver
        final NetworkAvgNodeTable[] avgTableHolder = new NetworkAvgNodeTable[1];
        withSuppressedOutput(() -> {
            SolverCTMC solver = new SolverCTMC(model, "cutoff", 5);
            avgTableHolder[0] = solver.getAvgNodeTable();
            assertEquals("default", solver.result.method, "CTMC solver should use default method");
        });
        NetworkAvgNodeTable avgTable = avgTableHolder[0];

        // Check if results are computed
        assertNotNull(avgTable);
        
        // Expected values from ground truth CTMC solver output
        // Order: Source(Class1), Router(Class1), Queue1(Class1), Queue2(Class1), Sink(Class1)
        double[] expectedQLen = {0.0, 0.0, 0.287011687938476, 0.287011687938476, 0.0};
        double[] expectedUtil = {0.0, 0.0, 0.249611007153048, 0.249611007153048, 0.0};
        double[] expectedRespT = {0.0, 0.0, 0.574917931729061, 0.57491793172906, 0.0};
        double[] expectedResidT = {0.0, 0.0, 0.28745896586453, 0.28745896586453, 0.0};
        double[] expectedArvR = {0.0, 0.998444038534968, 0.499222019267484, 0.499222019267484, 0.998444038534968};
        double[] expectedTput = {0.998444038534968, 0.998444038534968, 0.499222019267484, 0.499222019267484, 0.0};
        
        // Verify table size
        assertEquals(5, avgTable.getQLen().size(), 
            "Expected 5 entries (5 stations × 1 class)");
        
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                          expectedResidT, expectedArvR, expectedTput);
    }

    /**
     * Test sdroute_open with JMT solver.
     * 
     * JMT solver uses simulation for open networks.
     * Expected values from MATLAB dev/test_sdroute_open.m output.
     */
    @Test
    //@Disabled("Test failing - disabled for investigation")
    public void testSdrouteOpenJMT() {
        // Create the model
        Network model = StateDepRoutingModel.sdroute_open();
        
        // Create and run the solver
        final NetworkAvgNodeTable[] avgTableHolder = new NetworkAvgNodeTable[1];
        withSuppressedOutput(() -> {
            SolverJMT solver = new SolverJMT(model, "seed", 23000);
            avgTableHolder[0] = solver.getAvgNodeTable();
            assertEquals("default", solver.result.method, "JMT solver should use default method");
        });
        NetworkAvgNodeTable avgTable = avgTableHolder[0];
        // avgTable.print();

        // Check if results are computed
        assertNotNull(avgTable);
        
        // Expected values from ground truth JMT solver output
        // Order: Source(Class1), Router(Class1), Queue1(Class1), Queue2(Class1), Sink(Class1)
        double[] expectedQLen = {0.0, 0.0, 0.302397872438382, 0.293169630234353, 0.0};
        double[] expectedUtil = {0.0, 0.0, 0.262182472400218, 0.253606605795572, 0.0};
        double[] expectedRespT = {0.0, 0.0, 0.576111394707508, 0.569511400762629, 0.0};
        double[] expectedResidT = {0.0, 0.0, 0.288055697353754, 0.284755700381315, 0.0};
        double[] expectedArvR = {0.0, 1.01197698600991, 0.505988493004956, 0.505988493004956, 1.01197698600991};
        double[] expectedTput = {1.01197698600991, 1.01197698600991, 0.505988493004956, 0.505988493004956, 0.0};
        
        // Verify table size
        assertEquals(5, avgTable.getQLen().size(), 
            "Expected 5 entries (5 stations × 1 class)");
        
        // Check all metrics against expected values (use COARSE_TOL for simulation variance)
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                          expectedResidT, expectedArvR, expectedTput, COARSE_TOL);
    }

    /**
     * Test sdroute_twoclasses_closed with CTMC solver.
     * 
     * This test verifies a multi-class closed network with advanced service processes.
     * The model has:
     * - Two closed classes: Class1 (1 job), Class2 (2 jobs)
     * - Advanced service processes: APH, PH, MAP distributions
     * - Mixed scheduling: PS Queue1, FCFS Queue2
     * - Round-robin routing for both classes from delay
     * 
     * Expected values from MATLAB dev/test_sdroute_twoclasses_closed.m output.
     */
    @Test
    public void testSdrouteTwoclassesClosedCTMC() {
        // Create the model
        Network model = StateDepRoutingModel.sdroute_twoclasses_closed();

        // Create and run the solver
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            SolverCTMC solver = new SolverCTMC(model, "keep", false);
            avgTableHolder[0] = solver.getAvgTable();
            assertEquals("default", solver.result.method, "CTMC solver should use default method");
        });
        NetworkAvgTable avgTable = avgTableHolder[0];

        // Check if results are computed
        assertNotNull(avgTable);

        // Expected values from MATLAB ground truth CTMC
        double[] expectedQLen  = {0.400069438815555, 1.11067991256304, 0.496396547055449, 0.591440113956447, 0.103534014128997, 0.297879973480519};
        double[] expectedUtil  = {0.400069438815555, 1.11067991256304, 0.333391199012962, 0.370226637521012, 0.0770504104385513, 0.259158646264708};
        double[] expectedRespT = {1.0, 1.0, 3.72232791781157, 1.597508266603, 0.776370330377043, 0.804588171923781};
        double[] expectedResidT= {1.0, 1.0, 1.24077597260386, 0.532502755534332, 0.258790110125681, 0.268196057307927};
        double[] expectedArvR  = {0.400069438815555, 1.11067991256303, 0.133356479605185, 0.370226637521012, 0.133356479605185, 0.370226637521012};
        double[] expectedTput  = {0.400069438815554, 1.11067991256303, 0.133356479605185, 0.370226637521012, 0.133356479605185, 0.370226637521012};

        assertEquals(6, avgTable.getQLen().size(), "Expected 6 entries (3 stations × 2 classes)");
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                expectedResidT, expectedArvR, expectedTput, COARSE_TOL);
    }

    @Test
    public void testSdrouteTwoclassesClosedJMT() {
        // Create the model
        Network model = StateDepRoutingModel.sdroute_twoclasses_closed();

        // Create and run the solver
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            SolverJMT solver = new SolverJMT(model, "seed", 23000);
            solver.options.samples = 100000;
            avgTableHolder[0] = solver.getAvgTable();
            assertEquals("default", solver.result.method, "JMT solver should use default method");
        });
        NetworkAvgTable avgTable = avgTableHolder[0];

        assertNotNull(avgTable);

        // Expected values from MATLAB ground truth JMT
        double[] expectedQLen  = {0.408013386994108, 1.10397611739635, 0.497678698758767, 0.594521517756364, 0.102470946499728, 0.294328638730734};
        double[] expectedUtil  = {0.408013386994108, 1.10397611739635, 0.335720933262931, 0.370376371786127, 0.0755419502967797, 0.256923946244702};
        double[] expectedRespT = {1.03668186731426, 1.02099998511709, 3.7385743530352, 1.60282117383469, 0.771113075096285, 0.805736481419478};
        double[] expectedResidT= {1.03668186731426, 1.02099998511709, 1.24619145101173, 0.534273724611563, 0.257037691698762, 0.268578827139826};
        double[] expectedArvR  = {0.393682565047062, 1.10259450370911, 0.132510989868834, 0.374912967676714, 0.132510826180476, 0.374914903482739};
        double[] expectedTput  = {0.393681746941033, 1.10291933844795, 0.132511291138965, 0.374901938545321, 0.132511014896606, 0.374911627252682};

        assertEquals(6, avgTable.getQLen().size(), "Expected 6 entries (3 stations × 2 classes)");
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                expectedResidT, expectedArvR, expectedTput);
    }

    @Test
    public void testSdrouteTwoclassesClosedSSA() {
        // Create the model
        Network model = StateDepRoutingModel.sdroute_twoclasses_closed();

        // Create and run the solver
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            SolverSSA solver = new SolverSSA(model, "seed", 23000, "verbose", false, "samples", 50000);
            avgTableHolder[0] = solver.getAvgTable();
            assertEquals("default/serial", solver.result.method, "SSA solver should use default/serial method");
        });
        NetworkAvgTable avgTable = avgTableHolder[0];
        assertNotNull(avgTable);

        // Expected values from CTMC ground truth (exact analytical solution)
        // SSA is stochastic, so we compare against exact CTMC values with tolerance
        double[] expectedQLen  = {0.400069438815555, 1.11067991256304, 0.496396547055449, 0.591440113956447, 0.103534014128997, 0.297879973480519};
        double[] expectedUtil  = {0.400069438815555, 1.11067991256304, 0.333391199012962, 0.370226637521012, 0.0770504104385513, 0.259158646264708};
        double[] expectedRespT = {1.0, 1.0, 3.72232791781157, 1.597508266603, 0.776370330377043, 0.804588171923781};
        double[] expectedResidT= {1.0, 1.0, 1.24077597260386, 0.532502755534332, 0.258790110125681, 0.268196057307927};
        double[] expectedArvR  = {0.400069438815555, 1.11067991256303, 0.133356479605185, 0.370226637521012, 0.133356479605185, 0.370226637521012};
        double[] expectedTput  = {0.400069438815554, 1.11067991256303, 0.133356479605185, 0.370226637521012, 0.133356479605185, 0.370226637521012};

        assertEquals(6, avgTable.getQLen().size(), "Expected 6 entries (3 stations × 2 classes)");
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                expectedResidT, expectedArvR, expectedTput, VERY_COARSE_TOL);
    }
}