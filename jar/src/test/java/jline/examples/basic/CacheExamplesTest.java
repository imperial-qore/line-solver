package jline.examples.basic;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.examples.java.basic.CacheModel;
import jline.lang.Network;
import jline.lang.layered.LayeredNetwork;
import jline.solvers.NetworkAvgNodeTable;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.jmt.SolverJMT;
import jline.solvers.mva.SolverMVA;
import jline.solvers.nc.SolverNC;
import jline.solvers.ssa.SolverSSA;
import jline.util.Maths;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.*;
import static jline.TestTools.assertTableMetrics;
import static jline.TestTools.withSuppressedOutput;

public class CacheExamplesTest {

    @BeforeAll
    public static void setUp() {
        // Ensure MATLAB-compatible random number generation
        Maths.setRandomNumbersMatlab(true);
        // Set verbose level to SILENT to suppress warnings during tests
        GlobalConstants.setVerbose(VerboseLevel.SILENT);
    }
    
    
    // ======================== cache_replc_rr tests ========================
    // MATLAB uses: CTMC (cutoff=1), SSA, MVA, NC
    
    @Test
    //@Disabled("Test failing - disabled for investigation")
    public void testCacheReplcRRCTMC() {
        // Test the cache_replc_rr example with CTMC solver
        Network model = CacheModel.cache_replc_rr();
        
        SolverCTMC solver = new SolverCTMC(model, "cutoff", 1, "verbose", VerboseLevel.SILENT);
        NetworkAvgNodeTable avgTable = solver.getAvgNodeTable();

        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify the executed method
        assertNotNull(solver.result, "Solver result should not be null");
        assertEquals("default", solver.result.method, 
            "CTMC solver should use default method");
        
        // Expected values based on baseline output (6 non-zero rows)
        // The baseline filters out rows where all metrics are zero
        // Order: Source(InitClass), Cache(InitClass), Cache(HitClass), Cache(MissClass), Sink(HitClass), Sink(MissClass)
        double[] expectedQLen = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedUtil = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedRespT = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedResidT = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedArvR = {0, 0, 0, 1.99999992, 0, 0, 0, 1.1807919439343, 0.819207976065707};
        double[] expectedTput = {1.99999992, 0, 0, 0, 1.1807919439343, 0.819207976065707, 0, 0, 0};
        
        // Verify table size matches MATLAB output
        assertEquals(9, avgTable.getQLen().size(), "Expected 9 entries (3 nodes × 3 classes)");
        
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, expectedResidT, expectedArvR, expectedTput);
    }
    
    @Test
    //@Disabled("Test failing - disabled for investigation")
    public void testCacheReplcRRSSA() {
        // Test the cache_replc_rr example with SSA solver
        Network model = CacheModel.cache_replc_rr();
        
        SolverSSA solver = new SolverSSA(model, "samples", 10000, "verbose", VerboseLevel.SILENT, "method", "serial", "seed", 23000);
        NetworkAvgNodeTable avgTable = solver.getAvgNodeTable();
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify the executed method
        assertNotNull(solver.result, "Solver result should not be null");
        assertEquals("serial", solver.result.method, 
            "SSA solver should use serial method");
        
        // Expected values based on baseline output (6 non-zero rows)
        // Note: SSA results may vary slightly due to sampling
        double[] expectedQLen = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedUtil = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedRespT = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedResidT = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedArvR = {0, 0, 0, 2, 0, 0, 0, 1.16946678088066, 0.83053321911934};
        double[] expectedTput = {2, 0, 0, 0, 1.16946678088066, 0.83053321911934, 0, 0, 0};
        
        // Verify table size matches MATLAB output
        assertEquals(9, avgTable.getQLen().size(), "Expected 9 entries (3 nodes × 3 classes)");
        
        // Check all metrics against expected values (using larger tolerance for SSA)
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }

    //@Disabled
    @Test
    public void testCacheReplcRRMVA() {
        // Test the cache_replc_rr example with MVA solver
        Network model = CacheModel.cache_replc_rr();
        
        SolverMVA solver = new SolverMVA(model, "verbose", VerboseLevel.SILENT);
        NetworkAvgNodeTable avgTable = solver.getAvgNodeTable();
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify the executed method
        assertNotNull(solver.result, "Solver result should not be null");
        assertEquals("default/fpi", solver.result.method, 
            "MVA solver should use default/fpi method for cache models");
        
        // Expected values based on baseline output (6 non-zero rows)
        double[] expectedQLen = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedUtil = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedRespT = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedResidT = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedArvR = {0, 0, 0, 2, 0, 0, 0, 1.11270096048447, 0.887299039515532};
        double[] expectedTput = {2, 0, 0, 0, 1.11270096048447, 0.887299039515532, 0, 0, 0};
        
        // Verify table size matches MATLAB output
        assertEquals(9, avgTable.getQLen().size(), "Expected 9 entries (3 nodes × 3 classes)");
        
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }
    
    @Test
    public void testCacheReplcRRNC() {
        // Test the cache_replc_rr example with NC solver
        Network model = CacheModel.cache_replc_rr();
        
        SolverNC solver = new SolverNC(model, "verbose", VerboseLevel.SILENT);
        NetworkAvgNodeTable avgTable = solver.getAvgNodeTable();
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify the executed method
        assertNotNull(solver.result, "Solver result should not be null");
        assertEquals("default/rayint", solver.result.method, 
            "NC solver should use default/rayint method for cache models");
        
        // Expected values based on baseline output (6 non-zero rows)
        double[] expectedQLen = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedUtil = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedRespT = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedResidT = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedArvR = {0, 0, 0, 2, 0, 0, 0, 1.17218326985538, 0.827816730144616};
        double[] expectedTput = {2, 0, 0, 0, 1.17218326985538, 0.827816730144616, 0, 0, 0};
        
        // Verify table size matches MATLAB output
        assertEquals(9, avgTable.getQLen().size(), "Expected 9 entries (3 nodes × 3 classes)");
        
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }


    @Test
    public void testCacheReplcRRJMT() {
        // Test the cache_replc_rr example with JMT solver
        Network model = CacheModel.cache_replc_rr();

        SolverJMT solver = new SolverJMT(model, "seed", 23000, "verbose", VerboseLevel.SILENT, "samples", 100000);
        NetworkAvgNodeTable avgTable = solver.getAvgNodeTable();

        // Check if results are computed
        assertNotNull(avgTable);

        // Verify the executed method
        assertNotNull(solver.result, "Solver result should not be null");
        assertEquals("default", solver.result.method,
            "JMT solver should use default method");

        // Expected values based on MATLAB JMT output for cache models (with tolerance for simulation)
        double[] expectedQLen = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedUtil = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedRespT = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedResidT = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedArvR = {0, 0, 0, 2, 0, 0, 0, 1.11270096048447, 0.887299039515532};
        double[] expectedTput = {2, 0, 0, 0, 1.11270096048447, 0.887299039515532, 0, 0, 0};

        // Verify table size
        assertEquals(9, avgTable.getQLen().size(), "Expected 9 entries (3 nodes × 3 classes)");

        // Check all metrics against expected values (using 10% tolerance for simulation results)
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                          expectedResidT, expectedArvR, expectedTput, 0.10);
    }
    
    
    // ======================== cache_replc_fifo tests ========================
    // MATLAB uses: CTMC (no cutoff), SSA, MVA
    
    @Test
    public void testCacheReplcFIFOCTMC() {
        // Test the cache_replc_fifo example with CTMC solver
        Network model = CacheModel.cache_replc_fifo();
        
        SolverCTMC solver = new SolverCTMC(model, "keep", false, "verbose", VerboseLevel.SILENT);
        NetworkAvgNodeTable avgTable = solver.getAvgNodeTable();
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify the executed method
        assertNotNull(solver.result, "Solver result should not be null");
        assertEquals("default", solver.result.method, 
            "CTMC solver should use default method");
        
        // Expected values from MATLAB CTMC solver (dev/ directory output)
        // Java produces 9 entries (3 nodes × 3 classes) while MATLAB aggregates to 7
        double[] expectedQLen = {0.99999998, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedUtil = {0.99999998, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedRespT = {1, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedResidT = {1, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedArvR = {0.99999998, 0, 0, 0.99999998, 0, 0, 0, 0.399999992, 0.599999988};
        double[] expectedTput = {0.99999998, 0, 0, 0, 0.399999992, 0.599999988, 0.99999998, 0, 0};
        
        // Verify table size matches Java implementation
        assertEquals(9, avgTable.getQLen().size(), "Expected 9 entries (3 nodes × 3 classes)");
        
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }
    
    @Test
    //@Disabled("Test failing - disabled for investigation")
    public void testCacheReplcFIFOSSA() {
        // Test the cache_replc_fifo example with SSA solver
        Network model = CacheModel.cache_replc_fifo();
        
        SolverSSA solver = new SolverSSA(model, "samples", 10000, "verbose", VerboseLevel.SILENT, "method", "serial", "seed", 23000);
        NetworkAvgNodeTable avgTable = solver.getAvgNodeTable();
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify the executed method
        assertNotNull(solver.result, "Solver result should not be null");
        assertEquals("serial", solver.result.method, 
            "SSA solver should use serial method");
        
        // Expected values from MATLAB SSA solver (dev/ directory output)
        // Note: SSA results may vary slightly due to sampling
        double[] expectedQLen = {0.999999980265041, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedUtil = {0.999999980265041, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedRespT = {1, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedResidT = {1, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedArvR = {0.999999980265041, 0, 0, 0.999999980265041, 0, 0, 0, 0.408351410489633, 0.591648569775407};
        double[] expectedTput = {0.999999980265041, 0, 0, 0, 0.408351410489634, 0.591648569775407, 0.999999980265041, 0, 0};
        
        // Verify table size matches Java implementation
        assertEquals(9, avgTable.getQLen().size(), "Expected 9 entries (3 nodes × 3 classes)");
        
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }
    
    @Test
    public void testCacheReplcFIFOMVA() {
        // Test the cache_replc_fifo example with MVA solver
        Network model = CacheModel.cache_replc_fifo();
        
        SolverMVA solver = new SolverMVA(model, "verbose", VerboseLevel.SILENT);
        NetworkAvgNodeTable avgTable = solver.getAvgNodeTable();
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify the executed method
        assertNotNull(solver.result, "Solver result should not be null");
        assertEquals("default/fpi", solver.result.method, 
            "MVA solver should use default/fpi method for cache models");
        
        // Expected values from MATLAB MVA solver (dev/ directory output)
        double[] expectedQLen = {1, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedUtil = {1, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedRespT = {1, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedResidT = {1, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedArvR = {1, 0, 0, 1, 0, 0, 0, 0.4, 0.6};
        double[] expectedTput = {1, 0, 0, 0, 0.4, 0.6, 1, 0, 0};
        
        // Verify table size matches Java implementation
        assertEquals(9, avgTable.getQLen().size(), "Expected 9 entries (3 nodes × 3 classes)");
        
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }


    @Test
    public void testCacheReplcFIFOJMT() {
        // Test the cache_replc_fifo example with JMT solver
        Network model = CacheModel.cache_replc_fifo();
        
        SolverJMT solver = new SolverJMT(model, "seed", 23000, "verbose", VerboseLevel.SILENT, "samples", 100000);
        NetworkAvgNodeTable avgTable = solver.getAvgNodeTable();
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify the executed method
        assertNotNull(solver.result, "Solver result should not be null");
        assertEquals("default", solver.result.method, 
            "JMT solver should use default method");
        
        // Expected values based on MATLAB JMT output (with tolerance for stochastic variation)
        double[] expectedQLen = {1, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedUtil = {1, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedRespT = {1, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedResidT = {1, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedArvR = {1, 0, 0, 1, 0, 0, 0, 0.4, 0.6};
        double[] expectedTput = {1, 0, 0, 0, 0.4, 0.6, 1, 0, 0};

        // Verify table size
        assertEquals(9, avgTable.getQLen().size(), "Expected 9 entries (3 nodes × 3 classes)");

        // Check all metrics against expected values (using 10% tolerance for simulation results)
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                          expectedResidT, expectedArvR, expectedTput, 0.10);
    }


    // ======================== cache_replc_lru tests ========================
    // MATLAB uses: CTMC (no cutoff), SSA, MVA
    
    @Test
    public void testCacheReplcLRUCTMC() {
        // Test the cache_replc_lru example with CTMC solver
        Network model = CacheModel.cache_replc_lru();
        
        SolverCTMC solver = new SolverCTMC(model, "keep", false, "verbose", VerboseLevel.SILENT);
        NetworkAvgNodeTable avgTable = solver.getAvgNodeTable();
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify the executed method
        assertNotNull(solver.result, "Solver result should not be null");
        assertEquals("default", solver.result.method, 
            "CTMC solver should use default method");
        
        // Expected values from MATLAB CTMC solver (dev/ directory output)
        // Java produces 9 entries (3 nodes × 3 classes) while MATLAB aggregates
        double[] expectedQLen = {0.99999998, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedUtil = {0.99999998, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedRespT = {1, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedResidT = {1, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedArvR = {0.99999998, 0, 0, 0.99999998, 0, 0, 0, 0.399999992, 0.599999988};
        double[] expectedTput = {0.99999998, 0, 0, 0, 0.399999992, 0.599999988, 0.99999998, 0, 0};
        
        // Verify table size matches Java implementation
        assertEquals(9, avgTable.getQLen().size(), "Expected 9 entries (3 nodes × 3 classes)");
        
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }
    
    @Test
    // @Disabled - MAPE was 0.0580%, Max APE was 0.9087%
    public void testCacheReplcLRUSSA() {
        // Test the cache_replc_lru example with SSA solver
        Network model = CacheModel.cache_replc_lru();
        
        SolverSSA solver = new SolverSSA(model, "samples", 10000, "verbose", VerboseLevel.SILENT, "method", "serial", "seed", 23000);
        NetworkAvgNodeTable avgTable = solver.getAvgNodeTable();
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify the executed method
        assertNotNull(solver.result, "Solver result should not be null");
        assertEquals("serial", solver.result.method, 
            "SSA solver should use serial method");
        
        // Previous MAPE: 0.0580%, Max APE: 0.9087%
        // Expected values from MATLAB SSA solver (dev/ directory output)
        // Note: SSA results may vary slightly due to sampling, using MATLAB baseline
        double[] expectedQLen = {0.9999999802650407, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedUtil = {0.9999999802650407, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedRespT = {1, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedResidT = {1, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedArvR = {0.9999999802650407, 0, 0, 0.9999999802650416, 0, 0, 5.551115013574754E-17, 0.4046406673133472, 0.5953593129516938};
        double[] expectedTput = {0.9999999802650407, 0, 0, 0, 0.40464066731334725, 0.5953593129516941, 0.999999980265041, 0, 0};
        
        // Verify table size matches Java implementation
        assertEquals(9, avgTable.getQLen().size(), "Expected 9 entries (3 nodes × 3 classes)");
        
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }
    
    @Test
    public void testCacheReplcLRUMVA() {
        // Test the cache_replc_lru example with MVA solver
        Network model = CacheModel.cache_replc_lru();
        
        SolverMVA solver = new SolverMVA(model, "verbose", VerboseLevel.SILENT);
        NetworkAvgNodeTable avgTable = solver.getAvgNodeTable();
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify the executed method
        assertNotNull(solver.result, "Solver result should not be null");
        assertEquals("default/fpi", solver.result.method, 
            "MVA solver should use default/fpi method for cache models");
        
        // Expected values from MATLAB MVA solver (dev/ directory output)
        double[] expectedQLen = {1, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedUtil = {1, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedRespT = {1, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedResidT = {1, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedArvR = {1, 0, 0, 1, 0, 0, 0, 0.4, 0.6};
        double[] expectedTput = {1, 0, 0, 0, 0.4, 0.6, 1, 0, 0};
        
        // Verify table size matches Java implementation
        assertEquals(9, avgTable.getQLen().size(), "Expected 9 entries (3 nodes × 3 classes)");
        
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }
    
    @Test
    public void testCacheReplcLRUJMT() {
        // Test the cache_replc_lru example with JMT solver
        Network model = CacheModel.cache_replc_lru();

        SolverJMT solver = new SolverJMT(model, "seed", 23000, "verbose", VerboseLevel.SILENT, "samples", 100000);
        NetworkAvgNodeTable avgTable = solver.getAvgNodeTable();

        // Check if results are computed
        assertNotNull(avgTable);

        // Verify the executed method
        assertNotNull(solver.result, "Solver result should not be null");
        assertEquals("default", solver.result.method,
            "JMT solver should use default method");

        // Expected values based on MATLAB JMT output (with tolerance for stochastic variation)
        double[] expectedQLen = {1, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedUtil = {1, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedRespT = {1, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedResidT = {1, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedArvR = {1, 0, 0, 1, 0, 0, 0, 0.4, 0.6};
        double[] expectedTput = {1, 0, 0, 0, 0.4, 0.6, 1, 0, 0};

        // Verify table size
        assertEquals(9, avgTable.getQLen().size(), "Expected 9 entries (3 nodes × 3 classes)");

        // Check all metrics against expected values (using 10% tolerance for simulation results)
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                          expectedResidT, expectedArvR, expectedTput, 0.10);
    }
    
    
    // ======================== cache_compare_replc tests ========================
    // MATLAB uses: CTMC (cutoff=1), MVA, NC (if supported)
    
    @Test
    //@Disabled("Test failing - disabled for investigation")
    public void testCacheCompareReplcCTMC() {
        // Test the cache_compare_replc example with CTMC solver
        Network model = CacheModel.cache_compare_replc();
        
        SolverCTMC solver = new SolverCTMC(model, "keep", false, "cutoff", 1, "verbose", VerboseLevel.SILENT);
        NetworkAvgNodeTable avgTable = solver.getAvgNodeTable();

        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify the executed method
        assertNotNull(solver.result, "Solver result should not be null");
        assertEquals("default", solver.result.method, 
            "CTMC solver should use default method");
        
        // Expected values based on Java CTMC solver output (3 nodes × 3 classes = 9 entries)
        double[] expectedQLen = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedUtil = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedRespT = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedResidT = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedArvR = {0, 0, 0, 1.99999992, 0, 0, 0, 1.426764538, 0.5732353819};
        double[] expectedTput = {1.99999992, 0, 0, 0, 1.426764538, 0.5732353819, 0, 0, 0};
        
        // Verify table size matches MATLAB output
        assertEquals(9, avgTable.getQLen().size(), "Expected 9 entries (3 nodes × 3 classes)");
        
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }

    //@Disabled
    @Test
    public void testCacheCompareReplcMVA() {
        // Test the cache_compare_replc example with MVA solver
        Network model = CacheModel.cache_compare_replc();
        
        SolverMVA solver = new SolverMVA(model, "verbose", VerboseLevel.SILENT);
        NetworkAvgNodeTable avgTable = solver.getAvgNodeTable();
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify the executed method
        assertNotNull(solver.result, "Solver result should not be null");
        assertEquals("default/fpi", solver.result.method, 
            "MVA solver should use default/fpi method for cache models");
        
        // Expected values based on Java MVA solver output (3 nodes × 3 classes = 9 entries)
        double[] expectedQLen = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedUtil = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedRespT = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedResidT = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedArvR = {0, 0, 0, 2, 0, 0, 0, 1.39346155, 0.6065384502};
        double[] expectedTput = {2, 0, 0, 0, 1.39346155, 0.6065384502, 0, 0, 0};
        
        // Verify table size matches MATLAB output
        assertEquals(9, avgTable.getQLen().size(), "Expected 9 entries (3 nodes × 3 classes)");
        
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }
    
    @Test
    public void testCacheCompareReplcNC() {
        // Test the cache_compare_replc example with NC solver
        Network model = CacheModel.cache_compare_replc();
        
        SolverNC solver = new SolverNC(model, "verbose", VerboseLevel.SILENT);
        NetworkAvgNodeTable avgTable = solver.getAvgNodeTable();
        
        // Check if results are computed
        assertNotNull(avgTable);
        
        // Verify the executed method
        assertNotNull(solver.result, "Solver result should not be null");
        assertEquals("default/rayint", solver.result.method, 
            "NC solver should use default/rayint method for cache models");
        
        // Expected values based on Java NC solver output (3 nodes × 3 classes = 9 entries)
        double[] expectedQLen = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedUtil = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedRespT = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedResidT = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedArvR = {0, 0, 0, 2, 0, 0, 0, 1.406113146, 0.5938868541};
        double[] expectedTput = {2, 0, 0, 0, 1.406113146, 0.5938868541, 0, 0, 0};
        
        // Verify table size matches Java implementation
        assertEquals(9, avgTable.getQLen().size(), "Expected 9 entries (3 nodes × 3 classes)");
        
        // Check all metrics against expected values
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT, 
                          expectedResidT, expectedArvR, expectedTput);
    }

    @Test
    public void testCacheCompareReplcJMT() {
        // Test the cache_compare_replc example with JMT solver
        Network model = CacheModel.cache_compare_replc();

        SolverJMT solver = new SolverJMT(model, "seed", 23000, "verbose", VerboseLevel.SILENT, "samples", 100000);
        NetworkAvgNodeTable avgTable = solver.getAvgNodeTable();

        // Check if results are computed
        assertNotNull(avgTable);

        // Verify the executed method
        assertNotNull(solver.result, "Solver result should not be null");
        assertEquals("default", solver.result.method,
            "JMT solver should use default method");

        // Expected values based on MATLAB JMT output (with tolerance for stochastic variation)
        double[] expectedQLen = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedUtil = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedRespT = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedResidT = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        double[] expectedArvR = {0, 0, 0, 2, 0, 0, 0, 1.406113146, 0.5938868541};
        double[] expectedTput = {2, 0, 0, 0, 1.406113146, 0.5938868541, 0, 0, 0};

        // Verify table size
        assertEquals(9, avgTable.getQLen().size(), "Expected 9 entries (3 nodes × 3 classes)");

        // Check all metrics against expected values (using 10% tolerance for simulation results)
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                          expectedResidT, expectedArvR, expectedTput, 0.10);
    }
    
    
    // ======================== cache_replc_routing tests ========================
    // MATLAB uses: CTMC, SSA, MVA, NC


    @Test
    public void testCacheReplcRoutingJMT() {
        // Test the cache_replc_routing example with JMT solver
        Network model = CacheModel.cache_replc_routing();

        SolverJMT solver = new SolverJMT(model, "seed", 23000, "verbose", VerboseLevel.SILENT, "samples", 100000);
        NetworkAvgNodeTable avgTable = solver.getAvgNodeTable();

        // Check if results are computed
        assertNotNull(avgTable);

        // Verify the executed method
        assertNotNull(solver.result, "Solver result should not be null");
        assertEquals("default", solver.result.method,
            "JMT solver should use default method");

        // Expected values based on Java JMT solver output (seed=23000, samples=100000)
        // Note: This model has 6 nodes × 3 classes = 18 entries
        // Indices: Node0(0-2), Node1(3-5), Node2(6-8), Node3(9-11), Node4(12-14), Node5(15-17)
        // Classes: InitClass(0), HitClass(1), MissClass(2)
        double[] expectedQLen = new double[18];
        double[] expectedUtil = new double[18];
        double[] expectedRespT = new double[18];
        double[] expectedResidT = new double[18];
        double[] expectedArvR = new double[18];
        double[] expectedTput = new double[18];

        // Source (node 0): InitClass throughput = 2
        expectedTput[0] = 2;

        // Cache (node 1): InitClass arrival = 2, HitClass/MissClass throughput
        expectedArvR[3] = 2;   // Cache/InitClass arrival
        expectedTput[4] = 0.8; // Cache/HitClass throughput (hit prob ~0.4)
        expectedTput[5] = 1.2; // Cache/MissClass throughput (miss prob ~0.6)

        // Router (node 2): HitClass/MissClass pass through
        expectedArvR[7] = 0.8;  // Router/HitClass arrival
        expectedTput[7] = 0.8;  // Router/HitClass throughput
        expectedArvR[8] = 1.2;  // Router/MissClass arrival
        expectedTput[8] = 1.2;  // Router/MissClass throughput

        // Delay1 (node 3): HitClass and MissClass with random routing (50%)
        // Expected values based on Java JMT solver output (seed=23000, samples=100000)
        expectedQLen[10] = 0.04;    // Delay1/HitClass queue length
        expectedUtil[10] = 0.04;    // Delay1/HitClass utilization
        expectedRespT[10] = 0.1;    // Delay1/HitClass response time
        expectedResidT[10] = 0.0201; // Delay1/HitClass residence time
        expectedArvR[10] = 0.4;     // Delay1/HitClass arrival (0.8 * 0.5)
        expectedTput[10] = 0.4;     // Delay1/HitClass throughput
        expectedQLen[11] = 0.6;     // Delay1/MissClass queue length
        expectedUtil[11] = 0.6;     // Delay1/MissClass utilization
        expectedRespT[11] = 1.0;    // Delay1/MissClass response time
        expectedResidT[11] = 0.298; // Delay1/MissClass residence time
        expectedArvR[11] = 0.6;     // Delay1/MissClass arrival (1.2 * 0.5)
        expectedTput[11] = 0.6;     // Delay1/MissClass throughput

        // Delay2 (node 4): HitClass and MissClass with random routing (50%)
        // Expected values based on Java JMT solver output (seed=23000, samples=100000)
        expectedQLen[13] = 0.02;    // Delay2/HitClass queue length
        expectedUtil[13] = 0.02;    // Delay2/HitClass utilization
        expectedRespT[13] = 0.05;   // Delay2/HitClass response time
        expectedResidT[13] = 0.01;  // Delay2/HitClass residence time
        expectedArvR[13] = 0.4;     // Delay2/HitClass arrival
        expectedTput[13] = 0.4;     // Delay2/HitClass throughput
        expectedQLen[14] = 0.3;     // Delay2/MissClass queue length
        expectedUtil[14] = 0.3;     // Delay2/MissClass utilization
        expectedRespT[14] = 0.5;    // Delay2/MissClass response time
        expectedResidT[14] = 0.149; // Delay2/MissClass residence time
        expectedArvR[14] = 0.6;     // Delay2/MissClass arrival
        expectedTput[14] = 0.6;     // Delay2/MissClass throughput

        // Sink (node 5): Arrivals only
        expectedArvR[16] = 0.8;  // Sink/HitClass arrival
        expectedArvR[17] = 1.2;  // Sink/MissClass arrival

        // Verify table size
        assertEquals(18, avgTable.getQLen().size(), "Expected 18 entries (6 nodes × 3 classes)");

        // Check all metrics against expected values (using 10% tolerance for simulation results)
        assertTableMetrics(avgTable, expectedQLen, expectedUtil, expectedRespT,
                          expectedResidT, expectedArvR, expectedTput, 0.10);
    }
    
    
    // ======================== Structural tests for complex models ========================
    
    @Test
    public void testCacheReplcRouting() {
        // Test the cache_replc_routing example - model creation and structural validation
        // Note: This model has complex routing that requires specific solver handling
        Network model = CacheModel.cache_replc_routing();

        // Verify model was created successfully with correct structure
        assertNotNull(model, "Model should be created successfully");
        assertEquals("model", model.getName(), "Model name should match");
        assertEquals(6, model.getNumberOfNodes(), "Should have 6 nodes: Source, Cache, Router, Delay1, Delay2, Sink");
        assertEquals(3, model.getNumberOfClasses(), "Should have 3 classes: InitClass, HitClass, MissClass");

        // Validate model structure is consistent
        assertTrue(model.getNumberOfNodes() > 0, "Model should have nodes");
        assertTrue(model.getNumberOfClasses() > 0, "Model should have classes");
    }

    @Test
    public void testCacheReplcRoutingCTMC() {
        // Test the cache_replc_routing example with CTMC solver
        // Validates against MATLAB CTMC ground truth values
        Network model = CacheModel.cache_replc_routing();

        final NetworkAvgNodeTable[] avgTableHolder = new NetworkAvgNodeTable[1];
        withSuppressedOutput(() -> {
            SolverCTMC solver = new SolverCTMC(model, "cutoff", 1);
            avgTableHolder[0] = solver.getAvgNodeTable();
        });
        NetworkAvgNodeTable avgTable = avgTableHolder[0];

        assertNotNull(avgTable, "CTMC solver should produce results");

        // JAR AvgNodeTable: all (node x class) combos, 6 nodes × 3 classes = 18 entries
        // Index layout: Source(Init=0,Hit=1,Miss=2), Cache(Init=3,Hit=4,Miss=5),
        //   Router(Init=6,Hit=7,Miss=8), Delay1(Init=9,Hit=10,Miss=11),
        //   Delay2(Init=12,Hit=13,Miss=14), Sink(Init=15,Hit=16,Miss=17)
        List<Double> tput = avgTable.getTput();

        double tol = 1e-4;
        // MATLAB ground truth values (from CTMC solver with cutoff=1)
        assertEquals(1.980789069644373, tput.get(0), tol, "Source InitClass Tput");
        assertEquals(0.792315627857749, tput.get(4), tol, "Cache HitClass Tput");
        assertEquals(1.188473441786624, tput.get(5), tol, "Cache MissClass Tput");
        assertEquals(0.792315627857749, tput.get(7), tol, "Router HitClass Tput");
        assertEquals(1.188473441786624, tput.get(8), tol, "Router MissClass Tput");
        assertEquals(0.389554071016588, tput.get(10), tol, "Delay1 HitClass Tput");
        assertEquals(0.424640622506156, tput.get(11), tol, "Delay1 MissClass Tput");
        assertEquals(0.402227919095865, tput.get(13), tol, "Delay2 HitClass Tput");
        assertEquals(0.583730769969780, tput.get(14), tol, "Delay2 MissClass Tput");

        // Verify queue lengths at Delay nodes
        List<Double> qLen = avgTable.getQLen();
        assertEquals(0.038955407101659, qLen.get(10), tol, "Delay1 HitClass QLen");
        assertEquals(0.424640622506156, qLen.get(11), tol, "Delay1 MissClass QLen");
        assertEquals(0.020111395954793, qLen.get(13), tol, "Delay2 HitClass QLen");
        assertEquals(0.291865384984890, qLen.get(14), tol, "Delay2 MissClass QLen");
    }
    
    @Test
    public void testLcqSinglehost() {
        // Test the lcq_singlehost example - LayeredNetwork model structural validation
        LayeredNetwork model = CacheModel.lcq_singlehost();
        
        // Verify model was created successfully with correct layered structure
        assertNotNull(model, "LayeredNetwork model should be created successfully");
        assertEquals("cacheInLayeredNetwork", model.getName(), "Model name should match");
        assertNotNull(model.getStruct(), "LayeredNetwork struct should be initialized");
        
        // Verify layered network structure
        assertEquals(2, model.getStruct().nhosts, "Should have 2 hosts: P1, PC");
        assertEquals(2, model.getStruct().ntasks, "Should have 2 tasks: T1, C2 (CacheTask)");
        assertEquals(2, model.getStruct().nentries, "Should have 2 entries: E1, I2 (ItemEntry)");
        
        // Validate structural consistency
        assertTrue(model.getStruct().nhosts > 0, "Should have at least one host");
        assertTrue(model.getStruct().ntasks > 0, "Should have at least one task");
        assertTrue(model.getStruct().nentries > 0, "Should have at least one entry");
    }
    
    @Test
    public void testLcqThreehosts() {
        // Test the lcq_threehosts example - LayeredNetwork model structural validation
        LayeredNetwork model = CacheModel.lcq_threehosts();
        
        // Verify model was created successfully with correct layered structure
        assertNotNull(model, "LayeredNetwork model should be created successfully");
        assertEquals("LQNwithCaching", model.getName(), "Model name should match");
        assertNotNull(model.getStruct(), "LayeredNetwork struct should be initialized");
        
        // Verify layered network structure
        assertEquals(3, model.getStruct().nhosts, "Should have 3 hosts: P1, Pc, P2");
        assertEquals(3, model.getStruct().ntasks, "Should have 3 tasks: T1, CT (CacheTask), T2");
        assertEquals(3, model.getStruct().nentries, "Should have 3 entries: E1, IE (ItemEntry), E2");
        
        // Validate structural consistency for multi-tier layered network
        assertTrue(model.getStruct().nhosts > 0, "Should have at least one host");
        assertTrue(model.getStruct().ntasks > 0, "Should have at least one task");
        assertTrue(model.getStruct().nentries > 0, "Should have at least one entry");
        assertTrue(model.getStruct().nhosts >= model.getStruct().ntasks, "Hosts should accommodate tasks");
    }
    
}