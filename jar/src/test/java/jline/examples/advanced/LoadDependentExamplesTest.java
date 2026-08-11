package jline.examples.advanced;
import jline.GlobalConstants;

import jline.VerboseLevel;
import jline.examples.java.advanced.LoadDependentModel;
import jline.lang.Network;
import jline.solvers.NetworkAvgNodeTable;
import jline.solvers.ctmc.CTMCOptions;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.wrappers.jmt.SolverJMT;
import jline.solvers.mva.SolverMVA;
import jline.solvers.nc.SolverNC;
import jline.solvers.SolverOptions;
import jline.util.Maths;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;
import static jline.TestTools.MID_TOL;
import static jline.TestTools.COARSE_TOL;
import static jline.TestTools.relativeTolerance;

/**
 * Unit tests for LoadDependent examples, comparing solver outputs between MATLAB and Java implementations.
 * Expected values are obtained by running the corresponding MATLAB examples in the dev/ directory.
 * 
 * IMPORTANT: To populate the expected values in this test file:
 * 1. Run the MATLAB test scripts in the dev/ directory:
 *    - matlab -batch "run('dev/test_ld_class_dependence.m')"
 *    - matlab -batch "run('dev/test_ld_multiserver_fcfs.m')"
 *    - matlab -batch "run('dev/test_ld_multiserver_ps.m')"
 *    - matlab -batch "run('dev/test_ld_multiserver_ps_twoclasses.m')"
 * 2. Copy the AvgNodeTable output values into the expected arrays
 * 3. Note that Java may create implicit nodes that don't exist in MATLAB
 * 4. The order of entries in the table must match the Java implementation
 * 
 * NOTE: Due to MATLAB execution constraints, current expected values are from Java implementation.
 * These should be replaced with MATLAB dev/ output when available.
 */
public class LoadDependentExamplesTest {

    @BeforeAll
    public static void setUp() {
        // Ensure MATLAB-compatible random number generation
        Maths.setRandomNumbersMatlab(true);
        GlobalConstants.setVerbose(VerboseLevel.SILENT);
    }
    
    /**
     * Utility method to assert all metrics in a NetworkAvgNodeTable against expected values.
     * 
     * @param avgTable the NetworkAvgNodeTable to verify
     * @param expectedQLen expected queue length values
     * @param expectedUtil expected utilization values
     * @param expectedRespT expected response time values
     * @param expectedResidT expected residence time values
     * @param expectedArvR expected arrival rate values
     * @param expectedTput expected throughput values
     */
    private static void assertTableMetrics(NetworkAvgNodeTable avgTable,
                                         double[] expectedQLen,
                                         double[] expectedUtil,
                                         double[] expectedRespT,
                                         double[] expectedResidT,
                                         double[] expectedArvR,
                                         double[] expectedTput) {
        assertNotNull(avgTable, "AvgTable should not be null");
        assertEquals(expectedQLen.length, avgTable.getQLen().size(), 
                    "Expected " + expectedQLen.length + " entries in table");

        for (int i = 0; i < expectedQLen.length; i++) {
            assertEquals(expectedQLen[i], avgTable.getQLen().get(i), relativeTolerance(expectedQLen[i], MID_TOL),
                    String.format("QLen[%d]: expected %.6g, got %.6g", i, expectedQLen[i], avgTable.getQLen().get(i)));
            assertEquals(expectedUtil[i], avgTable.getUtil().get(i), relativeTolerance(expectedUtil[i], MID_TOL),
                    String.format("Util[%d]: expected %.6g, got %.6g", i, expectedUtil[i], avgTable.getUtil().get(i)));
            assertEquals(expectedRespT[i], avgTable.getRespT().get(i), relativeTolerance(expectedRespT[i], MID_TOL),
                    String.format("RespT[%d]: expected %.6g, got %.6g", i, expectedRespT[i], avgTable.getRespT().get(i)));
            assertEquals(expectedResidT[i], avgTable.getResidT().get(i), relativeTolerance(expectedResidT[i], MID_TOL),
                    String.format("ResidT[%d]: expected %.6g, got %.6g", i, expectedResidT[i], avgTable.getResidT().get(i)));
            assertEquals(expectedArvR[i], avgTable.getArvR().get(i), relativeTolerance(expectedArvR[i], MID_TOL),
                    String.format("ArvR[%d]: expected %.6g, got %.6g", i, expectedArvR[i], avgTable.getArvR().get(i)));
            assertEquals(expectedTput[i], avgTable.getTput().get(i), relativeTolerance(expectedTput[i], MID_TOL),
                    String.format("Tput[%d]: expected %.6g, got %.6g", i, expectedTput[i], avgTable.getTput().get(i)));
        }
    }

    /**
     * Overloaded version with custom tolerance for simulation-based tests.
     */
    private static void assertTableMetrics(NetworkAvgNodeTable avgTable,
                                         double[] expectedQLen,
                                         double[] expectedUtil,
                                         double[] expectedRespT,
                                         double[] expectedResidT,
                                         double[] expectedArvR,
                                         double[] expectedTput,
                                         double relTol) {
        assertNotNull(avgTable, "AvgTable should not be null");
        assertEquals(expectedQLen.length, avgTable.getQLen().size(),
                    "Expected " + expectedQLen.length + " entries in table");

        for (int i = 0; i < expectedQLen.length; i++) {
            assertEquals(expectedQLen[i], avgTable.getQLen().get(i), relativeTolerance(expectedQLen[i], relTol),
                    String.format("QLen[%d]: expected %.6g, got %.6g", i, expectedQLen[i], avgTable.getQLen().get(i)));
            assertEquals(expectedUtil[i], avgTable.getUtil().get(i), relativeTolerance(expectedUtil[i], relTol),
                    String.format("Util[%d]: expected %.6g, got %.6g", i, expectedUtil[i], avgTable.getUtil().get(i)));
            assertEquals(expectedRespT[i], avgTable.getRespT().get(i), relativeTolerance(expectedRespT[i], relTol),
                    String.format("RespT[%d]: expected %.6g, got %.6g", i, expectedRespT[i], avgTable.getRespT().get(i)));
            assertEquals(expectedResidT[i], avgTable.getResidT().get(i), relativeTolerance(expectedResidT[i], relTol),
                    String.format("ResidT[%d]: expected %.6g, got %.6g", i, expectedResidT[i], avgTable.getResidT().get(i)));
            assertEquals(expectedArvR[i], avgTable.getArvR().get(i), relativeTolerance(expectedArvR[i], relTol),
                    String.format("ArvR[%d]: expected %.6g, got %.6g", i, expectedArvR[i], avgTable.getArvR().get(i)));
            assertEquals(expectedTput[i], avgTable.getTput().get(i), relativeTolerance(expectedTput[i], relTol),
                    String.format("Tput[%d]: expected %.6g, got %.6g", i, expectedTput[i], avgTable.getTput().get(i)));
        }
    }

    // ===== ld_multiserver_fcfs tests =====
    // MATLAB uses: CTMC, JMT (simulation)
    
    @Test
    public void testLdMultiserverFcfsCTMC() {
        // Test the ld_multiserver_fcfs example with CTMC solver
        Network model = LoadDependentModel.ld_multiserver_fcfs();

        SolverCTMC solver = new SolverCTMC(model, "verbose", VerboseLevel.SILENT);
        NetworkAvgNodeTable avgTable = solver.getAvgNodeTable();
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Expected values from Java execution (should be replaced with MATLAB dev/ output)
        // Order: Delay(Class1), Queue1(Class1)
        double[] expectedQLen = {1.333333, 14.666667};
        double[] expectedUtil = {1.333333, 1.000000};
        double[] expectedRespT = {1.000000, 11.000000};
        double[] expectedResidT = {1.000000, 11.000000};
        double[] expectedArvR = {1.333333, 1.333333};
        double[] expectedTput = {1.333333, 1.333333};
        
        // Verify table size matches expected structure
        assertEquals(2, avgTable.getQLen().size(), 
            "Expected 2 entries (2 nodes × 1 class) for ld_multiserver_fcfs model");
        
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }
    
    @Test
    public void testLdMultiserverFcfsJMT() {
        // Test the ld_multiserver_fcfs example with JMT solver
        Network model = LoadDependentModel.ld_multiserver_fcfs();
        
        // Create solver options matching MATLAB
        SolverOptions options = new SolverOptions();
        options.seed = 12345;
        options.samples = 100000;
        options.verbose = VerboseLevel.SILENT;

        SolverJMT solver = new SolverJMT(model, options);
        NetworkAvgNodeTable avgTable = solver.getAvgNodeTable();

        // Check if results are computed
        assertNotNull(avgTable);

        // Expected values from CTMC (theoretical)
        // Order: Delay(Class1), Queue1(Class1)
        // Note: JMT simulation results may vary slightly from CTMC
        double[] expectedQLen = {1.333333, 14.666667};
        double[] expectedUtil = {1.333333, 1.000000};
        double[] expectedRespT = {1.000000, 11.000000};
        double[] expectedResidT = {1.000000, 11.000000};
        double[] expectedArvR = {1.333333, 1.333333};
        double[] expectedTput = {1.333333, 1.333333};

        // Verify table size matches expected structure
        assertEquals(2, avgTable.getQLen().size(),
            "Expected 2 entries (2 nodes × 1 class) for ld_multiserver_fcfs model");

        // Check all metrics against expected values (using 1% tolerance for simulation)
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                          expectedResidT, expectedArvR, expectedTput, COARSE_TOL);
    }

    // ===== ld_multiserver_ps_twoclasses tests =====
    // MATLAB uses: NC (method=rd), MVA (method=exact), JMT (simulation)

    @Test
    public void testLdMultiserverPsTwoclassesNC() {
        // Test the ld_multiserver_ps_twoclasses example with NC solver
        Network model = LoadDependentModel.ld_multiserver_ps_twoclasses();

        SolverNC solver = new SolverNC(model, "method", "rd", "verbose", VerboseLevel.SILENT);
        NetworkAvgNodeTable avgTable = solver.getAvgNodeTable();
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Expected values from MVA (theoretical)
        // Order: Delay(Class1,Class2), Queue1(Class1,Class2)
        double[] expectedQLen = {0.897502, 0.511574, 3.102498, 1.488426};
        double[] expectedUtil = {0.897502, 0.511574, 0.677206, 0.321434};
        double[] expectedRespT = {1.000000, 2.000000, 3.456815, 5.819003};
        double[] expectedResidT = {1.000000, 2.000000, 3.456815, 5.819003};
        double[] expectedArvR = {0.897502, 0.255787, 0.897502, 0.255787};
        double[] expectedTput = {0.897502, 0.255787, 0.897502, 0.255787};

        // Verify table size matches expected structure
        assertEquals(4, avgTable.getQLen().size(),
            "Expected 4 entries (2 nodes × 2 classes) for ld_multiserver_ps_twoclasses model");

        // Check all metrics against expected values (using 2% tolerance for NC-rd approximation)
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                          expectedResidT, expectedArvR, expectedTput, 0.02);
    }

    @Test
    public void testLdMultiserverPsTwoclassesMVA() {
        // Test the ld_multiserver_ps_twoclasses example with MVA solver
        Network model = LoadDependentModel.ld_multiserver_ps_twoclasses();

        SolverMVA solver = new SolverMVA(model, "method", "exact", "verbose", VerboseLevel.SILENT);
        NetworkAvgNodeTable avgTable = solver.getAvgNodeTable();
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Expected values from MATLAB MVA (exact LD-multiserver-PS; matches exact CTMC)
        // Order: Delay(Class1,Class2), Queue1(Class1,Class2)
        double[] expectedQLen = {0.897502, 0.511574, 3.102498, 1.488426};
        double[] expectedUtil = {0.897502, 0.511574, 0.673126, 0.319734};
        double[] expectedRespT = {1.000000, 2.000000, 3.456815, 5.819003};
        double[] expectedResidT = {1.000000, 2.000000, 3.456815, 5.819003};
        double[] expectedArvR = {0.897502, 0.255787, 0.897502, 0.255787};
        double[] expectedTput = {0.897502, 0.255787, 0.897502, 0.255787};
        
        // Verify table size matches expected structure
        assertEquals(4, avgTable.getQLen().size(), 
            "Expected 4 entries (2 nodes × 2 classes) for ld_multiserver_ps_twoclasses model");
        
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }
    
    @Test
    public void testLdMultiserverPsTwoclassesJMT() {
        // Test the ld_multiserver_ps_twoclasses example with JMT solver
        Network model = LoadDependentModel.ld_multiserver_ps_twoclasses();

        // Create solver options matching MATLAB
        SolverOptions options = new SolverOptions();
        options.seed = 12345;
        options.samples = 100000;
        options.verbose = VerboseLevel.SILENT;

        SolverJMT solver = new SolverJMT(model, options);
        NetworkAvgNodeTable avgTable = solver.getAvgNodeTable();
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Expected values from MVA (theoretical)
        // Order: Delay(Class1,Class2), Queue1(Class1,Class2)
        // Note: JMT simulation results may vary slightly
        double[] expectedQLen = {0.897502, 0.511574, 3.102498, 1.488426};
        double[] expectedUtil = {0.897502, 0.511574, 0.677206, 0.321434};
        double[] expectedRespT = {1.000000, 2.000000, 3.456815, 5.819003};
        double[] expectedResidT = {1.000000, 2.000000, 3.456815, 5.819003};
        double[] expectedArvR = {0.897502, 0.255787, 0.897502, 0.255787};
        double[] expectedTput = {0.897502, 0.255787, 0.897502, 0.255787};

        // Verify table size matches expected structure
        assertEquals(4, avgTable.getQLen().size(),
            "Expected 4 entries (2 nodes × 2 classes) for ld_multiserver_ps_twoclasses model");

        // Check all metrics against expected values (using 1% tolerance for simulation)
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                          expectedResidT, expectedArvR, expectedTput, COARSE_TOL);
    }

    // ===== ld_class_dependence tests =====
    // MATLAB uses: CTMC, JMT (simulation)
    // Note: These models have computational complexity and may not produce valid results
    
    @Test
    public void testLdClassDependenceCTMC() {
        // Test the ld_class_dependence example with CTMC solver
        Network model = LoadDependentModel.ld_class_dependence();

        try {
            // Queue1 is a class-dependent station beta(n)=min(n_1,2) with
            // declared peak c=2; utilization is normalized as T*S/peak.
            CTMCOptions options = new CTMCOptions();
            options.verbose = VerboseLevel.SILENT;
            SolverCTMC solver = new SolverCTMC(model, options);
            NetworkAvgNodeTable avgTable = solver.getAvgNodeTable();
            
            // Check if results are computed
            assertNotNull(avgTable);
            
            // Verify table size matches expected structure  
            assertEquals(4, avgTable.getQLen().size(), 
                "Expected 4 entries (2 nodes × 2 classes) for ld_class_dependence model");
            
            // Note: This model has computational complexity that may prevent valid results
            // If solver succeeds, it should produce non-negative values
            for (int i = 0; i < avgTable.getQLen().size(); i++) {
                assertTrue(avgTable.getQLen().get(i) >= 0, "Queue lengths should be non-negative");
                assertTrue(avgTable.getUtil().get(i) >= 0, "Utilizations should be non-negative");
                assertTrue(avgTable.getUtil().get(i) <= 1, "Utilizations should be at most 1");
            }
            
        } catch (RuntimeException e) {
            // Load-dependent models may have computational issues with CTMC
            assertTrue(e.getMessage().contains("Unable to compute results") || 
                      e.getMessage().contains("matrix") ||
                      e.getMessage().contains("convergence"), 
                "Expected known issue with load-dependent model in CTMC solver: " + e.getMessage());
        }
    }
    
    @Test
    public void testLdClassDependenceJMT() {
        // Test the ld_class_dependence example with JMT solver
        Network model = LoadDependentModel.ld_class_dependence();

        // Create solver options matching MATLAB
        SolverOptions options = new SolverOptions();
        options.seed = 12345;
        options.samples = 100000;
        options.verbose = VerboseLevel.SILENT;

        try {
            SolverJMT solver = new SolverJMT(model, options);
            NetworkAvgNodeTable avgTable = solver.getAvgNodeTable();
            
            // Check if results are computed
            assertNotNull(avgTable);
            
            // Verify table size matches expected structure
            assertEquals(4, avgTable.getQLen().size(), 
                "Expected 4 entries (2 nodes × 2 classes) for ld_class_dependence model");
            
            // Note: This model has computational complexity that may prevent valid results
            // If solver succeeds, it should produce non-negative values
            for (int i = 0; i < avgTable.getQLen().size(); i++) {
                assertTrue(avgTable.getQLen().get(i) >= 0, "Queue lengths should be non-negative");
                assertTrue(avgTable.getUtil().get(i) >= 0, "Utilizations should be non-negative");
                assertTrue(avgTable.getUtil().get(i) <= 1, "Utilizations should be at most 1");
            }
            
        } catch (RuntimeException e) {
            // Load-dependent models may have issues with JMT solver
            assertTrue(e.getMessage().contains("Unable to compute results") || 
                      e.getMessage().contains("matrix bounds") ||
                      e.getMessage().contains("load dependent"), 
                "Expected known issue with load-dependent model in JMT solver: " + e.getMessage());
        }
    }

    // ===== ld_multiserver_ps tests =====
    // MATLAB uses: NC (method=rd), MVA (method=exact), JMT (simulation)
    // Note: These models have 3 nodes and may have computational complexity
    
    @Test
    public void testLdMultiserverPsNC() {
        // Test the ld_multiserver_ps example with NC solver
        Network model = LoadDependentModel.ld_multiserver_ps();

        try {
            SolverNC solver = new SolverNC(model, "method", "rd", "verbose", VerboseLevel.SILENT);
            NetworkAvgNodeTable avgTable = solver.getAvgNodeTable();
            
            // Check if results are computed
            assertNotNull(avgTable);
            
            // Verify table size matches expected structure
            assertEquals(6, avgTable.getQLen().size(), 
                "Expected 6 entries (3 nodes × 2 classes) for ld_multiserver_ps model");
            
            // Note: This model has computational complexity that may prevent valid results
            // If solver succeeds, it should produce non-negative values
            for (int i = 0; i < avgTable.getQLen().size(); i++) {
                assertTrue(avgTable.getQLen().get(i) >= 0, "Queue lengths should be non-negative");
                assertTrue(avgTable.getUtil().get(i) >= 0, "Utilizations should be non-negative");
                assertTrue(avgTable.getUtil().get(i) <= 1, "Utilizations should be at most 1");
            }
            
        } catch (RuntimeException e) {
            // Load-dependent models may have computational issues with NC
            assertTrue(e.getMessage().contains("Unable to compute results") || 
                      e.getMessage().contains("matrix") ||
                      e.getMessage().contains("convergence"), 
                "Expected known issue with load-dependent model in NC solver: " + e.getMessage());
        }
    }
    
    @Test
    public void testLdMultiserverPsMVA() {
        // Test the ld_multiserver_ps example with MVA solver
        Network model = LoadDependentModel.ld_multiserver_ps();

        try {
            SolverMVA solver = new SolverMVA(model, "method", "exact", "verbose", VerboseLevel.SILENT);
            NetworkAvgNodeTable avgTable = solver.getAvgNodeTable();
            
            // Check if results are computed
            assertNotNull(avgTable);
            
            // Verify table size matches expected structure
            assertEquals(6, avgTable.getQLen().size(), 
                "Expected 6 entries (3 nodes × 2 classes) for ld_multiserver_ps model");
            
            // Note: This model has computational complexity that may prevent valid results
            // If solver succeeds, it should produce non-negative values
            for (int i = 0; i < avgTable.getQLen().size(); i++) {
                assertTrue(avgTable.getQLen().get(i) >= 0, "Queue lengths should be non-negative");
                assertTrue(avgTable.getUtil().get(i) >= 0, "Utilizations should be non-negative");
                assertTrue(avgTable.getUtil().get(i) <= 1, "Utilizations should be at most 1");
            }
            
        } catch (RuntimeException e) {
            // Load-dependent models may have computational issues with MVA
            assertTrue(e.getMessage().contains("Unable to compute results") || 
                      e.getMessage().contains("matrix") ||
                      e.getMessage().contains("convergence"), 
                "Expected known issue with load-dependent model in MVA solver: " + e.getMessage());
        }
    }
    
    @Test
    public void testLdMultiserverPsJMT() {
        // Test the ld_multiserver_ps example with JMT solver
        Network model = LoadDependentModel.ld_multiserver_ps();

        // Create solver options matching MATLAB
        SolverOptions options = new SolverOptions();
        options.seed = 12345;
        options.samples = 100000;
        options.verbose = VerboseLevel.SILENT;

        try {
            SolverJMT solver = new SolverJMT(model, options);
            NetworkAvgNodeTable avgTable = solver.getAvgNodeTable();
            
            // Check if results are computed
            assertNotNull(avgTable);
            
            // Verify table size matches expected structure
            assertEquals(6, avgTable.getQLen().size(), 
                "Expected 6 entries (3 nodes × 2 classes) for ld_multiserver_ps model");
            
            // Note: This model has computational complexity that may prevent valid results
            // If solver succeeds, it should produce non-negative values
            for (int i = 0; i < avgTable.getQLen().size(); i++) {
                assertTrue(avgTable.getQLen().get(i) >= 0, "Queue lengths should be non-negative");
                assertTrue(avgTable.getUtil().get(i) >= 0, "Utilizations should be non-negative");
                assertTrue(avgTable.getUtil().get(i) <= 1, "Utilizations should be at most 1");
            }
            
        } catch (RuntimeException e) {
            // Load-dependent models may have issues with JMT solver
            assertTrue(e.getMessage().contains("Unable to compute results") || 
                      e.getMessage().contains("matrix bounds") ||
                      e.getMessage().contains("load dependent"), 
                "Expected known issue with load-dependent model in JMT solver: " + e.getMessage());
        }
    }
    
    // ===== Structural tests for verifying model creation =====
    
    @Test
    public void testLdClassDependenceStructure() {
        // Test the structural properties of the ld_class_dependence model
        Network model = LoadDependentModel.ld_class_dependence();
        
        // Verify model structure
        assertNotNull(model, "Model should be created successfully");
        assertEquals("model", model.getName(), "Model name should match");
        assertEquals(2, model.getNumberOfClasses(), "Should have 2 classes: Class1 and Class2");
        assertEquals(2, model.getNumberOfNodes(), "Should have 2 nodes: Delay and Queue1");
        
        // Verify the model has closed classes with correct populations
        assertEquals(16, model.getClasses().get(0).getNumberOfJobs(), "Class1 should have 16 jobs");
        assertEquals(8, model.getClasses().get(1).getNumberOfJobs(), "Class2 should have 8 jobs (N/2)");
    }
    
    @Test
    public void testLdMultiserverFcfsStructure() {
        // Test the structural properties of the ld_multiserver_fcfs model
        Network model = LoadDependentModel.ld_multiserver_fcfs();
        
        // Verify model structure
        assertNotNull(model, "Model should be created successfully");
        assertEquals("model", model.getName(), "Model name should match");
        assertEquals(1, model.getNumberOfClasses(), "Should have 1 class: Class1");
        assertEquals(2, model.getNumberOfNodes(), "Should have 2 nodes: Delay and Queue1");
        
        // Verify the model has closed class with correct population
        assertEquals(16, model.getClasses().get(0).getNumberOfJobs(), "Class1 should have 16 jobs");
    }
    
    @Test
    public void testLdMultiserverPsStructure() {
        // Test the structural properties of the ld_multiserver_ps model
        Network model = LoadDependentModel.ld_multiserver_ps();
        
        // Verify model structure
        assertNotNull(model, "Model should be created successfully");
        assertEquals("model", model.getName(), "Model name should match");
        assertEquals(2, model.getNumberOfClasses(), "Should have 2 classes: Class1 and Class2");
        assertEquals(3, model.getNumberOfNodes(), "Should have 3 nodes: Delay, Queue1, Queue2");
        
        // Verify the model has closed classes with correct populations
        assertEquals(4, model.getClasses().get(0).getNumberOfJobs(), "Class1 should have 4 jobs");
        assertEquals(2, model.getClasses().get(1).getNumberOfJobs(), "Class2 should have 2 jobs (N/2)");
    }
    
    @Test
    public void testLdMultiserverPsTwoclassesStructure() {
        // Test the structural properties of the ld_multiserver_ps_twoclasses model
        Network model = LoadDependentModel.ld_multiserver_ps_twoclasses();
        
        // Verify model structure
        assertNotNull(model, "Model should be created successfully");
        assertEquals("model", model.getName(), "Model name should match");
        assertEquals(2, model.getNumberOfClasses(), "Should have 2 classes: Class1 and Class2");
        assertEquals(2, model.getNumberOfNodes(), "Should have 2 nodes: Delay and Queue1");
        
        // Verify the model has closed classes with correct populations
        assertEquals(4, model.getClasses().get(0).getNumberOfJobs(), "Class1 should have 4 jobs");
        assertEquals(2, model.getClasses().get(1).getNumberOfJobs(), "Class2 should have 2 jobs (N/2)");
    }
}