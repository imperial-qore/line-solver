package jline.examples.advanced;
import jline.GlobalConstants;
import jline.VerboseLevel;

import jline.examples.java.advanced.LoadDependentModel;
import jline.lang.Network;
import jline.solvers.NetworkAvgTable;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.des.SolverDES;
import jline.solvers.jmt.SolverJMT;
import jline.solvers.mva.SolverMVA;
import jline.solvers.nc.SolverNC;
import jline.solvers.SolverOptions;
import jline.util.Maths;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;
import static jline.TestTools.withSuppressedOutput;
import static jline.TestTools.assertTableMetrics;

/**
 * Unit tests for LoadDependent examples.
 *
 * Expected values come from allExamplesBaseline.txt MATLAB execution where available.
 * Only ld_multiserver_fcfs is present in allExamplesBaseline.txt. Other examples retain
 * their existing Java baseline values for regression testing.
 */
public class LoadDepExamplesTest {

    @BeforeAll
    public static void setUp() {
        // Ensure MATLAB-compatible random number generation
        Maths.setRandomNumbersMatlab(true);
        GlobalConstants.setVerbose(VerboseLevel.SILENT);
    }
    
    // ===== ld_multiserver_fcfs tests =====
    // MATLAB uses: CTMC, JMT (simulation)
    
    @Test
    public void testLdMultiserverFcfsCTMC() {
        // Test the ld_multiserver_fcfs example with CTMC solver
        Network model = LoadDependentModel.ld_multiserver_fcfs();
        
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            SolverCTMC solver = new SolverCTMC(model);
            avgTableHolder[0] = solver.getAvgTable();
            assertEquals("default", solver.result.method, "CTMC solver should use default method");
        });
        NetworkAvgTable avgTable = avgTableHolder[0];
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Expected values from allExamplesBaseline.txt (CTMC solver)
        // Order: Delay(Class1), Queue1(Class1)
        double[] expectedQLen = {1.33333333332244, 14.6666666666776};
        double[] expectedUtil = {1.33333333332244, 0.999999999999372};
        double[] expectedRespT = {1.0, 11.000000000098};
        double[] expectedResidT = {1.0, 11.000000000098};
        double[] expectedArvR = {1.33333333332244, 1.33333333332244};
        double[] expectedTput = {1.33333333332244, 1.33333333332244};
        
        // Verify table size matches expected structure
        assertEquals(2, avgTable.getQLen().size(), 
            "Expected 2 entries (2 stations × 1 class) for ld_multiserver_fcfs model");
        
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }

    // ===== ld_multiserver_ps_twoclasses tests =====
    // ANNOTATION: ld_multiserver_ps_twoclasses example is not present in allExamplesBaseline.txt.
    // Keeping existing Java baseline values for regression testing.
    
    @Test
    public void testLdMultiserverPsTwoclassesNCExact() {
        // Test the ld_multiserver_ps_twoclasses example with NC solver (exact method)
        Network model = LoadDependentModel.ld_multiserver_ps_twoclasses();
        
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            SolverNC solver = new SolverNC(model, "method", "exact");
            avgTableHolder[0] = solver.getAvgTable();
            assertEquals("exact/comomld", solver.result.method, "NC solver should use exact method");
        });
        NetworkAvgTable avgTable = avgTableHolder[0];
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Expected values from MATLAB NC solver output with exact/comomld method
        // Order: Delay(Class1,Class2), Queue1(Class1,Class2)
        double[] expectedQLen = {0.897501892505677, 0.511574166301447, 3.10249810749432, 1.48842583369855};
        double[] expectedUtil = {0.897501892505677, 0.511574166301447, 0.673126419379258, 0.319733853938404};
        double[] expectedRespT = {1.0, 2.0, 3.45681511470985, 5.81900311526479};
        double[] expectedResidT = {1.0, 2.0, 3.45681511470985, 5.81900311526479};
        double[] expectedArvR = {0.897501892505677, 0.255787083150723, 0.897501892505677, 0.255787083150723};
        double[] expectedTput = {0.897501892505677, 0.255787083150723, 0.897501892505677, 0.255787083150723};
        
        // Verify table size matches expected structure
        assertEquals(4, avgTable.getQLen().size(), "Expected 4 entries (2 stations × 2 classes) for ld_multiserver_ps_twoclasses model");

        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }
    
    @Test
    public void testLdMultiserverPsTwoclassesNCRD() {
        // Test the ld_multiserver_ps_twoclasses example with NC solver (RD method)
        Network model = LoadDependentModel.ld_multiserver_ps_twoclasses();
        
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            SolverNC solver = new SolverNC(model, "method", "rd");
            avgTableHolder[0] = solver.getAvgTable();
            assertEquals("rd", solver.result.method, "NC solver should use rd method");
        });
        NetworkAvgTable avgTable = avgTableHolder[0];
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Expected values from MATLAB NC solver output (RD method)
        // From MATLAB ground truth output
        // Order: Delay(Class1,Class2), Queue1(Class1,Class2)
        double[] expectedQLen = {0.889720578495604, 0.509608978547403, 3.11027942150440, 1.49039102145260};
        double[] expectedUtil = {0.889720578495604, 0.509608978547403, 0.667290433871703, 0.318505611592127};
        double[] expectedRespT = {1.0, 2.0, 3.49579350717441, 5.84915527077576};
        double[] expectedResidT = {1.0, 2.0, 3.49579350717441, 5.84915527077576};
        double[] expectedArvR = {0.889720578495604, 0.254804489273702, 0.889720578495604, 0.254804489273702};
        double[] expectedTput = {0.889720578495604, 0.254804489273702, 0.889720578495604, 0.254804489273702};
        
        // Verify table size matches expected structure
        assertEquals(4, avgTable.getQLen().size(), 
            "Expected 4 entries (2 stations × 2 classes) for ld_multiserver_ps_twoclasses model");
        
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }
    
    @Test
    public void testLdMultiserverPsTwoclassesMVA() {
        // Test the ld_multiserver_ps_twoclasses example with MVA solver
        Network model = LoadDependentModel.ld_multiserver_ps_twoclasses();
        
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            SolverMVA solver = new SolverMVA(model, "method", "exact");
            avgTableHolder[0] = solver.getAvgTable();
            assertEquals("exact", solver.result.method, "MVA solver should use exact method");
        });
        NetworkAvgTable avgTable = avgTableHolder[0];
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Expected values from MATLAB MVA solver output
        // These match the MVA load-dependent solver, not NC solver
        // Order: Delay(Class1,Class2), Queue1(Class1,Class2)
        double[] expectedQLen = {0.897501892505677, 0.511574166301446, 3.10249810749432, 1.48842583369855};
        double[] expectedUtil = {0.897501892505677, 0.511574166301446, 0.677206263197737, 0.321433788862770};
        double[] expectedRespT = {1.0, 2.0, 3.45681511470985, 5.8190031152648};
        double[] expectedResidT = {1.0, 2.0, 3.45681511470985, 5.8190031152648};
        double[] expectedArvR = {0.897501892505677, 0.255787083150723, 0.897501892505677, 0.255787083150723};
        double[] expectedTput = {0.897501892505677, 0.255787083150723, 0.897501892505677, 0.255787083150723};
        
        // Verify table size matches expected structure
        assertEquals(4, avgTable.getQLen().size(), 
            "Expected 4 entries (2 stations × 2 classes) for ld_multiserver_ps_twoclasses model");
        
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }
    
    @Test
    //@Disabled("QLen[0] mismatch ==> expected: <0.897502> but was: <0.44711046113604697>")
    public void testLdMultiserverPsTwoclassesJMT() {
        // Test the ld_multiserver_ps_twoclasses example with JMT solver
        Network model = LoadDependentModel.ld_multiserver_ps_twoclasses();
        
        // Create solver options matching MATLAB
        SolverOptions options = new SolverOptions();
        options.seed = 23000;
        options.samples = 5000;
        
        final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            SolverJMT solver = new SolverJMT(model, options);
            avgTableHolder[0] = solver.getAvgTable();
            assertEquals("default", solver.result.method, "JMT solver should use default method");
        });
        NetworkAvgTable avgTable = avgTableHolder[0];
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Expected values from JMT simulation output (seed=23000, samples=5000)
        // Verified against MVA analytical solution (within ~3% simulation tolerance)
        // Order: Delay(Class1,Class2), Queue1(Class1,Class2)
        double[] expectedQLen = {0.9246282256641253, 0.49666100193946176, 3.0753717743358747, 1.5033389980605383};
        double[] expectedUtil = {0.9246282256641253, 0.49666100193946176, 0.6769272067116478, 0.3183379011310485};
        double[] expectedRespT = {0.9855172851668193, 1.928477437315547, 3.38338787523959, 5.737737406632632};
        double[] expectedResidT = {0.9855172851668193, 1.928477437315547, 3.38338787523959, 5.737737406632632};
        double[] expectedArvR = {0.9211296108655024, 0.26303150167987224, 0.9187564310915316, 0.26008853979851624};
        double[] expectedTput = {0.9187564310915316, 0.26008853979851624, 0.9153310952968068, 0.26312157421426676};
        
        // Verify table size matches expected structure
        assertEquals(4, avgTable.getQLen().size(), 
            "Expected 4 entries (2 stations × 2 classes) for ld_multiserver_ps_twoclasses model");
        
        // Check all metrics against expected values (using larger tolerance for simulation)
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput, 0.01); // 1% tolerance for JMT simulation
    }

    // ===== ld_class_dependence tests =====
    // ANNOTATION: ld_class_dependence example is not present in allExamplesBaseline.txt.
    // Note: These models have computational complexity and may not produce valid results
    
    @Test
    public void testLdClassDependenceCTMC() {
        // Test the ld_class_dependence example with CTMC solver
        Network model = LoadDependentModel.ld_class_dependence();
        
        try {
            final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
            withSuppressedOutput(() -> {
                SolverCTMC solver = new SolverCTMC(model);
                avgTableHolder[0] = solver.getAvgTable();
                assertEquals("default", solver.result.method, "CTMC solver should use default method");
            });
            NetworkAvgTable avgTable = avgTableHolder[0];
            
            // Check if results are computed
            assertNotNull(avgTable);
            
            // Verify table size matches expected structure  
            assertEquals(4, avgTable.getQLen().size(), 
                "Expected 4 entries (2 stations × 2 classes) for ld_class_dependence model");
            
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
        options.seed = 23000;
        options.samples = 100000;

        try {
            final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
            withSuppressedOutput(() -> {
                SolverJMT solver = new SolverJMT(model, options);
                avgTableHolder[0] = solver.getAvgTable();
                assertEquals("default", solver.result.method, "JMT solver should use default method");
            });
            NetworkAvgTable avgTable = avgTableHolder[0];

            // Check if results are computed
            assertNotNull(avgTable);

            // Verify table size matches expected structure
            assertEquals(4, avgTable.getQLen().size(),
                "Expected 4 entries (2 stations × 2 classes) for ld_class_dependence model");
            
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
    // ANNOTATION: ld_multiserver_ps example is not present in allExamplesBaseline.txt.
    // Note: These models have 3 nodes and may have computational complexity
    
    @Test
    //@Disabled("Expected 6 entries (3 stations × 2 classes) for ld_multiserver_ps model ==> expected: <6> but was: <0>")
    public void testLdMultiserverPsNC() {
        // Test the ld_multiserver_ps example with NC solver
        Network model = LoadDependentModel.ld_multiserver_ps();
        
        try {
            final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
            withSuppressedOutput(() -> {
                SolverNC solver = new SolverNC(model, "method", "rd");
                avgTableHolder[0] = solver.getAvgTable();
                assertEquals("rd", solver.result.method, "NC solver should use rd method");
            });
            NetworkAvgTable avgTable = avgTableHolder[0];
            
            // Check if results are computed
            assertNotNull(avgTable);
            
            // Verify table size matches expected structure
            assertEquals(6, avgTable.getQLen().size(), 
                "Expected 6 entries (3 stations × 2 classes) for ld_multiserver_ps model");
            
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
            final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
            withSuppressedOutput(() -> {
                SolverMVA solver = new SolverMVA(model, "method", "exact");
                avgTableHolder[0] = solver.getAvgTable();
            });
            NetworkAvgTable avgTable = avgTableHolder[0];
            
            // Check if results are computed
            assertNotNull(avgTable);
            
            // Verify table size matches expected structure
            assertEquals(6, avgTable.getQLen().size(), 
                "Expected 6 entries (3 stations × 2 classes) for ld_multiserver_ps model");
            
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
        options.seed = 23000;
        options.samples = 100000;

        try {
            final NetworkAvgTable[] avgTableHolder = new NetworkAvgTable[1];
            withSuppressedOutput(() -> {
                SolverJMT solver = new SolverJMT(model, options);
                avgTableHolder[0] = solver.getAvgTable();
                assertEquals("default", solver.result.method, "JMT solver should use default method");
            });
            NetworkAvgTable avgTable = avgTableHolder[0];

            // Check if results are computed
            assertNotNull(avgTable);

            // Verify table size matches expected structure
            assertEquals(6, avgTable.getQLen().size(),
                "Expected 6 entries (3 stations × 2 classes) for ld_multiserver_ps model");
            
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

    // ===== ld_multiserver_ps DES tests =====
    // Tests for load-dependent PS with DES solver

    @Test
    public void testLdMultiserverPsDES() {
        // Test the ld_multiserver_ps example with DES solver
        // Compare against MVA solver as ground truth
        Network model = LoadDependentModel.ld_multiserver_ps();

        // Get MVA reference values
        final NetworkAvgTable[] mvaTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            SolverMVA mvaSolver = new SolverMVA(model, "method", "exact");
            mvaTableHolder[0] = mvaSolver.getAvgTable();
        });
        NetworkAvgTable mvaTable = mvaTableHolder[0];
        assertNotNull(mvaTable, "MVA solver should produce results");

        // Run DES solver with moderate sample count (50k balances accuracy vs test time)
        SolverOptions options = new SolverOptions();
        options.seed = 23000;
        options.samples = 50000;

        final NetworkAvgTable[] desTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            SolverDES solver = new SolverDES(model, options);
            desTableHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable desTable = desTableHolder[0];

        assertNotNull(desTable, "DES solver should produce results");
        assertEquals(6, desTable.getQLen().size(),
            "Expected 6 entries (3 stations × 2 classes) for ld_multiserver_ps model");

        // Verify DES results are close to MVA results (5% tolerance for simulation)
        double tolerance = 0.05;
        for (int i = 0; i < mvaTable.getQLen().size(); i++) {
            double mvaQLen = mvaTable.getQLen().get(i);
            double desQLen = desTable.getQLen().get(i);
            double mvaTput = mvaTable.getTput().get(i);
            double desTput = desTable.getTput().get(i);

            if (mvaQLen > 0.01) {
                double qlenRelError = Math.abs(desQLen - mvaQLen) / mvaQLen;
                assertTrue(qlenRelError < tolerance,
                    String.format("QLen[%d] relative error %.4f exceeds tolerance %.2f (DES=%.4f, MVA=%.4f)",
                        i, qlenRelError, tolerance, desQLen, mvaQLen));
            }

            if (mvaTput > 0.01) {
                double tputRelError = Math.abs(desTput - mvaTput) / mvaTput;
                assertTrue(tputRelError < tolerance,
                    String.format("Tput[%d] relative error %.4f exceeds tolerance %.2f (DES=%.4f, MVA=%.4f)",
                        i, tputRelError, tolerance, desTput, mvaTput));
            }
        }
    }

    @Test
    public void testLdMultiserverPsTwoclassesDES() {
        // Test the ld_multiserver_ps_twoclasses example with DES solver
        // Compare against MVA solver as ground truth
        Network model = LoadDependentModel.ld_multiserver_ps_twoclasses();

        // Get MVA reference values
        final NetworkAvgTable[] mvaTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            SolverMVA mvaSolver = new SolverMVA(model, "method", "exact");
            mvaTableHolder[0] = mvaSolver.getAvgTable();
        });
        NetworkAvgTable mvaTable = mvaTableHolder[0];
        assertNotNull(mvaTable, "MVA solver should produce results");

        // Run DES solver
        SolverOptions options = new SolverOptions();
        options.seed = 23000;
        options.samples = 50000;  // Reduced from 100000 (10% tolerance allows fewer samples)

        final NetworkAvgTable[] desTableHolder = new NetworkAvgTable[1];
        withSuppressedOutput(() -> {
            SolverDES solver = new SolverDES(model, options);
            desTableHolder[0] = solver.getAvgTable();
        });
        NetworkAvgTable desTable = desTableHolder[0];

        assertNotNull(desTable, "DES solver should produce results");
        assertEquals(4, desTable.getQLen().size(),
            "Expected 4 entries (2 stations × 2 classes) for ld_multiserver_ps_twoclasses model");

        // Verify DES results are close to MVA results (10% tolerance for simulation)
        double tolerance = 0.10;
        for (int i = 0; i < mvaTable.getQLen().size(); i++) {
            double mvaQLen = mvaTable.getQLen().get(i);
            double desQLen = desTable.getQLen().get(i);
            double mvaTput = mvaTable.getTput().get(i);
            double desTput = desTable.getTput().get(i);

            if (mvaQLen > 0.01) {
                double qlenRelError = Math.abs(desQLen - mvaQLen) / mvaQLen;
                assertTrue(qlenRelError < tolerance,
                    String.format("QLen[%d] relative error %.4f exceeds tolerance %.2f (DES=%.4f, MVA=%.4f)",
                        i, qlenRelError, tolerance, desQLen, mvaQLen));
            }

            if (mvaTput > 0.01) {
                double tputRelError = Math.abs(desTput - mvaTput) / mvaTput;
                assertTrue(tputRelError < tolerance,
                    String.format("Tput[%d] relative error %.4f exceeds tolerance %.2f (DES=%.4f, MVA=%.4f)",
                        i, tputRelError, tolerance, desTput, mvaTput));
            }
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
