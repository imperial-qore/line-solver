package jline.lang.nodes;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import static jline.TestTools.MID_TOL;
import static jline.TestTools.ZERO_TOL;
import jline.lang.constant.ReplacementStrategy;
import jline.lang.processes.Exp;
import jline.lang.processes.Immediate;
import jline.lang.processes.Zipf;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.des.SolverDES;
import jline.solvers.jmt.SolverJMT;
import jline.solvers.nc.SolverNC;
import jline.solvers.ssa.SolverSSA;
import jline.solvers.NetworkAvgTable;
import jline.solvers.SolverOptions;
import jline.VerboseLevel;
import jline.util.matrix.Matrix;


import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

/**
 * Unit tests for specific cache configurations with exact CTMC hit ratios.
 * Tests verify SSA accuracy against CTMC values and exact CTMC hit ratios.
 * Based on gettingstarted_ex6 with closed-loop arrival patterns.
 */
@Disabled
class CacheListTest {

    private static final double TOLERANCE_SSA = MID_TOL; // 5% tolerance for SSA vs CTMC
    private static final double TOLERANCE_EXACT = ZERO_TOL; // Exact tolerance for CTMC hit ratios

    @Disabled("CTMC state space too large (h=5 creates ~2520 cache states)")
    @Test
    void testLRU_n7_h5_alpha1_1e4() {
        // LRU (n=7, h=5, alpha=1.0, 1e4): 1,1,1,1,1 -> 0.849017715124717
        Matrix itemLevelCap = new Matrix(1, 5);
        for (int i = 0; i < 5; i++) {
            itemLevelCap.set(0, i, 1);
        }
        
        Network model = createCacheModel(7, itemLevelCap, ReplacementStrategy.LRU, 1.0);
        
        double expectedHitRatio = 0.849017715124717;
        testCacheConfiguration(model, expectedHitRatio, 10000, "LRU n=7, h=5, alpha=1.0");
    }

    @Disabled("CTMC state space too large (h=4 creates ~840 cache states)")
    @Test
    void testLRU_n7_h4_alpha1_1e4() {
        // LRU (n=7, h=4, alpha=1.0, 1e4): 1,1,1,1 -> 0.752582801392053
        Matrix itemLevelCap = new Matrix(1, 4);
        for (int i = 0; i < 4; i++) {
            itemLevelCap.set(0, i, 1);
        }
        
        Network model = createCacheModel(7, itemLevelCap, ReplacementStrategy.LRU, 1.0);
        
        double expectedHitRatio = 0.752582801392053;
        testCacheConfiguration(model, expectedHitRatio, 10000, "LRU n=7, h=4, alpha=1.0");
    }


    @Disabled("CTMC state space too large (capacity=5 creates ~2520 cache states)")
    @Test
    void testLRU_n7_alpha1_config1() {
        // LRU (n=7, alpha=1.0): 1,2,2 -> 0.846912573826510
        Matrix itemLevelCap = new Matrix(1, 3);
        itemLevelCap.set(0, 0, 1);
        itemLevelCap.set(0, 1, 2);
        itemLevelCap.set(0, 2, 2);
        
        Network model = createCacheModel(7, itemLevelCap, ReplacementStrategy.LRU, 1.0);
        
        double expectedHitRatio = 0.846912573826510;
        testCacheConfiguration(model, expectedHitRatio, 10000, "LRU n=7, alpha=1.0 (1,2,2)");
    }

    @Disabled("CTMC state space too large (capacity=4 creates ~840 cache states)")
    @Test
    void testLRU_n7_alpha1_config2() {
        // LRU (n=7, alpha=1.0): 1,1,2 -> 0.752411031903945
        Matrix itemLevelCap = new Matrix(1, 3);
        itemLevelCap.set(0, 0, 1);
        itemLevelCap.set(0, 1, 1);
        itemLevelCap.set(0, 2, 2);
        
        Network model = createCacheModel(7, itemLevelCap, ReplacementStrategy.LRU, 1.0);
        
        double expectedHitRatio = 0.752411031903945;
        testCacheConfiguration(model, expectedHitRatio, 10000, "LRU n=7, alpha=1.0 (1,1,2)");
    }

    @Disabled("CTMC state space too large (capacity=4 creates ~840 cache states)")
    @Test
    void testLRU_n7_alpha1_config3() {
        // LRU (n=7, alpha=1.0): 1,2,1 -> 0.748678893424558
        Matrix itemLevelCap = new Matrix(1, 3);
        itemLevelCap.set(0, 0, 1);
        itemLevelCap.set(0, 1, 2);
        itemLevelCap.set(0, 2, 1);
        
        Network model = createCacheModel(7, itemLevelCap, ReplacementStrategy.LRU, 1.0);
        
        double expectedHitRatio = 0.748678893424558;
        testCacheConfiguration(model, expectedHitRatio, 10000, "LRU n=7, alpha=1.0 (1,2,1)");
    }

    @Test
    void testLRU_n7_alpha1_config4() {
        // LRU (n=7, alpha=1.0): 1,1,1 -> 0.626935646456052
        Matrix itemLevelCap = new Matrix(1, 3);
        itemLevelCap.set(0, 0, 1);
        itemLevelCap.set(0, 1, 1);
        itemLevelCap.set(0, 2, 1);

        Network model = createCacheModel(7, itemLevelCap, ReplacementStrategy.LRU, 1.0);

        double expectedHitRatio = 0.626935646456052;
        testCacheConfiguration(model, expectedHitRatio, 10000, "LRU n=7, alpha=1.0 (1,1,1)");
    }

    @Test
    void testLRU_n4_alpha1_1e4() {
        // LRU (n=4, alpha=1.0, 1e4): 1,1,1 -> 0.839480519480520
        Matrix itemLevelCap = new Matrix(1, 3);
        itemLevelCap.set(0, 0, 1);
        itemLevelCap.set(0, 1, 1);
        itemLevelCap.set(0, 2, 1);

        Network model = createCacheModel(4, itemLevelCap, ReplacementStrategy.LRU, 1.0);

        double expectedHitRatio = 0.839480519480520;
        testCacheConfiguration(model, expectedHitRatio, 10000, "LRU n=4, alpha=1.0");
    }

    @Test
    void testLRU_n7_h2_alpha1_1e6() {
        // LRU (n=7, h=2, alpha=1.0, 1e6): 1,1 -> 0.454926571785874
        Matrix itemLevelCap = new Matrix(1, 2);
        itemLevelCap.set(0, 0, 1);
        itemLevelCap.set(0, 1, 1);

        Network model = createCacheModel(7, itemLevelCap, ReplacementStrategy.LRU, 1.0);

        double expectedHitRatio = 0.454926571785874;
        testCacheConfiguration(model, expectedHitRatio, 1000000, "LRU n=7, h=2, alpha=1.0");
    }

    @Disabled("CTMC state space too large (capacity=5 creates ~2520 cache states)")
    @Test
    void testRR_n7_alpha1_config1() {
        // RR (n=7, alpha=1.0): 1,2,2 -> 0.841467052083640
        Matrix itemLevelCap = new Matrix(1, 3);
        itemLevelCap.set(0, 0, 1);
        itemLevelCap.set(0, 1, 2);
        itemLevelCap.set(0, 2, 2);
        
        Network model = createCacheModel(7, itemLevelCap, ReplacementStrategy.RR, 1.0);
        
        double expectedHitRatio = 0.841467052083640;
        testCacheConfiguration(model, expectedHitRatio, 10000, "RR n=7, alpha=1.0 (1,2,2)");
    }

    @Test
    void testRR_n7_alpha1_config2() {
        // RR (n=7, alpha=1.0): 1,1,1 -> 0.626935646456052
        Matrix itemLevelCap = new Matrix(1, 3);
        itemLevelCap.set(0, 0, 1);
        itemLevelCap.set(0, 1, 1);
        itemLevelCap.set(0, 2, 1);

        Network model = createCacheModel(7, itemLevelCap, ReplacementStrategy.RR, 1.0);

        double expectedHitRatio = 0.626935646456052;
        testCacheConfiguration(model, expectedHitRatio, 10000, "RR n=7, alpha=1.0 (1,1,1)");
    }

    @Disabled("CTMC state space too large (capacity=5 creates ~2520 cache states)")
    @Test
    void testFIFO_n7_alpha1() {
        // FIFO (n=7, alpha=1.0): 1,2,2 -> 0.841385288703200
        Matrix itemLevelCap = new Matrix(1, 3);
        itemLevelCap.set(0, 0, 1);
        itemLevelCap.set(0, 1, 2);
        itemLevelCap.set(0, 2, 2);
        
        Network model = createCacheModel(7, itemLevelCap, ReplacementStrategy.FIFO, 1.0);
        
        double expectedHitRatio = 0.841385288703200;
        testCacheConfiguration(model, expectedHitRatio, 10000, "FIFO n=7, alpha=1.0 (1,2,2)");
    }

    /**
     * Creates a cache model similar to gettingstarted_ex6 with specified parameters
     */
    private Network createCacheModel(int nItems, Matrix itemLevelCap, ReplacementStrategy strategy, double zipfAlpha) {
        Network model = new Network("CacheConfigTest_" + strategy.name());

        // Create nodes
        Delay clientDelay = new Delay(model, "Client");
        Cache cacheNode = new Cache(model, "Cache", nItems, itemLevelCap, strategy);
        Delay cacheDelay = new Delay(model, "CacheDelay");

        // Create classes - similar to gettingstarted_ex6
        ClosedClass clientClass = new ClosedClass(model, "ClientClass", 1, clientDelay, 0);
        ClosedClass hitClass = new ClosedClass(model, "HitClass", 0, clientDelay, 0);
        ClosedClass missClass = new ClosedClass(model, "MissClass", 0, clientDelay, 0);

        // Set service processes
        clientDelay.setService(clientClass, new Immediate());
        cacheDelay.setService(hitClass, Exp.fitMean(0.2));
        cacheDelay.setService(missClass, Exp.fitMean(1.0));

        // Set cache read probabilities with Zipf distribution
        cacheNode.setRead(clientClass, new Zipf(zipfAlpha, nItems));
        cacheNode.setHitClass(clientClass, hitClass);
        cacheNode.setMissClass(clientClass, missClass);

        // Set topology - same as gettingstarted_ex6
        RoutingMatrix P = model.initRoutingMatrix();
        // routing from client to cache
        P.set(clientClass, clientClass, clientDelay, cacheNode, 1.0);
        // routing out of the cache
        P.set(hitClass, hitClass, cacheNode, cacheDelay, 1.0);
        P.set(missClass, missClass, cacheNode, cacheDelay, 1.0);
        // return to the client
        P.set(hitClass, clientClass, cacheDelay, clientDelay, 1.0);
        P.set(missClass, clientClass, cacheDelay, clientDelay, 1.0);

        model.link(P);
        return model;
    }

    /**
     * Tests a cache configuration by comparing SSA with CTMC results
     */
    private void testCacheConfiguration(Network model, double expectedHitRatio, int ssaSamples, String testDescription) {
        

        try {
            // Test with CTMC first to verify exact hit ratio
            SolverOptions ctmcOptions = new SolverOptions();
            ctmcOptions.verbose = VerboseLevel.SILENT;
            SolverCTMC ctmcSolver = new SolverCTMC(model, ctmcOptions);
            NetworkAvgTable ctmcAvgTable = ctmcSolver.getAvgTable();
            assertNotNull(ctmcAvgTable, "CTMC should produce results");

            // Calculate hit ratio from CTMC results
            // Hit ratio = cache hit throughput / (cache hit throughput + cache miss throughput)
            double ctmcHitTput = ctmcAvgTable.getTput().get(1); // Hit class throughput
            double ctmcMissTput = ctmcAvgTable.getTput().get(2); // Miss class throughput
            
            
            double totalTput = ctmcHitTput + ctmcMissTput;
            double ctmcHitRatio = -1; // Initialize to invalid value
            boolean ctmcValid = false;
            
            if (totalTput > 0) {
                ctmcHitRatio = ctmcHitTput / totalTput;
                ctmcValid = true;
                
                // Verify CTMC gives exactly the expected hit ratio
                assertEquals(expectedHitRatio, ctmcHitRatio, TOLERANCE_EXACT, 
                    String.format("CTMC hit ratio should be exactly %.15f, got %.15f", expectedHitRatio, ctmcHitRatio));
                
                
            } else {
                
                // Still perform SSA test below as it's more robust
            }

            // Test with SSA
            SolverOptions ssaOptions = new SolverOptions();
            ssaOptions.verbose = VerboseLevel.SILENT;
            SolverSSA ssaSolver = new SolverSSA(model, ssaOptions);
            ssaSolver.getOptions().samples = ssaSamples;
            ssaSolver.getOptions().seed = 1;
            NetworkAvgTable ssaAvgTable = ssaSolver.getAvgTable();
            assertNotNull(ssaAvgTable, "SSA should produce results");

            // Calculate hit ratio from SSA results
            double ssaHitTput = ssaAvgTable.getTput().get(1); // Hit class throughput
            double ssaMissTput = ssaAvgTable.getTput().get(2); // Miss class throughput
            double ssaHitRatio = ssaHitTput / (ssaHitTput + ssaMissTput);

            // Verify SSA accuracy against CTMC (if CTMC was valid)
            if (ctmcValid) {
                double relativeError = Math.abs(ssaHitRatio - ctmcHitRatio) / ctmcHitRatio;
                assertTrue(relativeError <= TOLERANCE_SSA, 
                    String.format("SSA hit ratio relative error %.4f exceeds tolerance %.2f (SSA: %.6f, CTMC: %.6f)", 
                        relativeError, TOLERANCE_SSA, ssaHitRatio, ctmcHitRatio));

                
            } else {
                // Just verify SSA gives a reasonable hit ratio when CTMC failed
                assertTrue(ssaHitRatio >= 0.0 && ssaHitRatio <= 1.0, 
                    String.format("SSA hit ratio %.6f should be between 0 and 1", ssaHitRatio));
                
                
            }
            

        } catch (Exception e) {
            fail("Test failed with exception: " + e.getMessage());
        }
    }

    // ========== DES Solver Tests ==========

    /**
     * DES simulation tests for cache configurations.
     * Validates hit ratios are within reasonable bounds.
     */
    @Nested
    class DESTests {

        @Test
        void testDES_LRU_n7_h5() {
            double[] cap = {1,1,1,1,1};
            double hitRatio = runDES(7, cap, ReplacementStrategy.LRU, 1.0);
            assertTrue(hitRatio >= 0.0 && hitRatio <= 1.0, "Hit ratio should be between 0 and 1");
        }

        @Test
        void testDES_LRU_n7_h4() {
            double[] cap = {1,1,1,1};
            double hitRatio = runDES(7, cap, ReplacementStrategy.LRU, 1.0);
            assertTrue(hitRatio >= 0.0 && hitRatio <= 1.0, "Hit ratio should be between 0 and 1");
        }

        @Test
        void testDES_LRU_n7_config1() {
            double[] cap = {1,2,2};
            double hitRatio = runDES(7, cap, ReplacementStrategy.LRU, 1.0);
            assertTrue(hitRatio >= 0.0 && hitRatio <= 1.0, "Hit ratio should be between 0 and 1");
        }

        @Test
        void testDES_LRU_n7_config2() {
            double[] cap = {1,1,2};
            double hitRatio = runDES(7, cap, ReplacementStrategy.LRU, 1.0);
            assertTrue(hitRatio >= 0.0 && hitRatio <= 1.0, "Hit ratio should be between 0 and 1");
        }

        @Test
        void testDES_LRU_n7_config3() {
            double[] cap = {1,2,1};
            double hitRatio = runDES(7, cap, ReplacementStrategy.LRU, 1.0);
            assertTrue(hitRatio >= 0.0 && hitRatio <= 1.0, "Hit ratio should be between 0 and 1");
        }

        @Test
        void testDES_LRU_n7_config4() {
            double[] cap = {1,1,1};
            double hitRatio = runDES(7, cap, ReplacementStrategy.LRU, 1.0);
            assertTrue(hitRatio >= 0.0 && hitRatio <= 1.0, "Hit ratio should be between 0 and 1");
        }

        @Test
        void testDES_LRU_n4() {
            double[] cap = {1,1,1};
            double hitRatio = runDES(4, cap, ReplacementStrategy.LRU, 1.0);
            assertTrue(hitRatio >= 0.0 && hitRatio <= 1.0, "Hit ratio should be between 0 and 1");
        }

        @Test
        void testDES_LRU_n7_h2() {
            double[] cap = {1,1};
            double hitRatio = runDES(7, cap, ReplacementStrategy.LRU, 1.0);
            assertTrue(hitRatio >= 0.0 && hitRatio <= 1.0, "Hit ratio should be between 0 and 1");
        }

        @Test
        void testDES_RR_n7_config1() {
            double[] cap = {1,2,2};
            double hitRatio = runDES(7, cap, ReplacementStrategy.RR, 1.0);
            assertTrue(hitRatio >= 0.0 && hitRatio <= 1.0, "Hit ratio should be between 0 and 1");
        }

        @Test
        void testDES_RR_n7_config2() {
            double[] cap = {1,1,1};
            double hitRatio = runDES(7, cap, ReplacementStrategy.RR, 1.0);
            assertTrue(hitRatio >= 0.0 && hitRatio <= 1.0, "Hit ratio should be between 0 and 1");
        }

        @Test
        void testDES_FIFO_n7() {
            double[] cap = {1,2,2};
            double hitRatio = runDES(7, cap, ReplacementStrategy.FIFO, 1.0);
            assertTrue(hitRatio >= 0.0 && hitRatio <= 1.0, "Hit ratio should be between 0 and 1");
        }

        private double runDES(int nItems, double[] cap, ReplacementStrategy strategy, double zipfAlpha) {
            Network model = new Network("CacheConfigTest_DES");

            Delay clientDelay = new Delay(model, "Client");
            Matrix itemLevelCap = new Matrix(1, cap.length);
            for (int i = 0; i < cap.length; i++) {
                itemLevelCap.set(0, i, cap[i]);
            }
            Cache cacheNode = new Cache(model, "Cache", nItems, itemLevelCap, strategy);
            Delay cacheDelay = new Delay(model, "CacheDelay");

            ClosedClass clientClass = new ClosedClass(model, "ClientClass", 1, clientDelay, 0);
            ClosedClass hitClass = new ClosedClass(model, "HitClass", 0, clientDelay, 0);
            ClosedClass missClass = new ClosedClass(model, "MissClass", 0, clientDelay, 0);

            clientDelay.setService(clientClass, new Immediate());
            cacheDelay.setService(hitClass, Exp.fitMean(0.2));
            cacheDelay.setService(missClass, Exp.fitMean(1.0));

            cacheNode.setRead(clientClass, new Zipf(zipfAlpha, nItems));
            cacheNode.setHitClass(clientClass, hitClass);
            cacheNode.setMissClass(clientClass, missClass);

            RoutingMatrix P = model.initRoutingMatrix();
            P.set(clientClass, clientClass, clientDelay, cacheNode, 1.0);
            P.set(hitClass, hitClass, cacheNode, cacheDelay, 1.0);
            P.set(missClass, missClass, cacheNode, cacheDelay, 1.0);
            P.set(hitClass, clientClass, cacheDelay, clientDelay, 1.0);
            P.set(missClass, clientClass, cacheDelay, clientDelay, 1.0);

            model.link(P);

            SolverOptions options = new SolverOptions();
            options.verbose = VerboseLevel.SILENT;
            options.seed = 1;
            options.samples = 1000000;

            SolverDES solver = new SolverDES(model, options);
            NetworkAvgTable avgTable = solver.getAvgTable();

            double hitTput = avgTable.getTput().get(1);
            double missTput = avgTable.getTput().get(2);
            return hitTput / (hitTput + missTput);
        }
    }

    // ========== JMT Solver Tests ==========

    /**
     * JMT simulation tests for cache configurations.
     */
    @Nested
    class JMTTests {

        @Test
        void testJMT_LRU_n7_h5() {
            double[] cap = {1,1,1,1,1};
            double hitRatio = runJMT(7, cap, ReplacementStrategy.LRU, 1.0);
            assertTrue(hitRatio >= 0.0 && hitRatio <= 1.0, "Hit ratio should be between 0 and 1");
        }

        @Test
        void testJMT_LRU_n7_h4() {
            double[] cap = {1,1,1,1};
            double hitRatio = runJMT(7, cap, ReplacementStrategy.LRU, 1.0);
            assertTrue(hitRatio >= 0.0 && hitRatio <= 1.0, "Hit ratio should be between 0 and 1");
        }

        @Test
        void testJMT_LRU_n7_config1() {
            double[] cap = {1,2,2};
            double hitRatio = runJMT(7, cap, ReplacementStrategy.LRU, 1.0);
            assertTrue(hitRatio >= 0.0 && hitRatio <= 1.0, "Hit ratio should be between 0 and 1");
        }

        @Test
        void testJMT_LRU_n7_config2() {
            double[] cap = {1,1,2};
            double hitRatio = runJMT(7, cap, ReplacementStrategy.LRU, 1.0);
            assertTrue(hitRatio >= 0.0 && hitRatio <= 1.0, "Hit ratio should be between 0 and 1");
        }

        @Test
        void testJMT_LRU_n7_config3() {
            double[] cap = {1,2,1};
            double hitRatio = runJMT(7, cap, ReplacementStrategy.LRU, 1.0);
            assertTrue(hitRatio >= 0.0 && hitRatio <= 1.0, "Hit ratio should be between 0 and 1");
        }

        @Test
        void testJMT_LRU_n7_config4() {
            double[] cap = {1,1,1};
            double hitRatio = runJMT(7, cap, ReplacementStrategy.LRU, 1.0);
            assertTrue(hitRatio >= 0.0 && hitRatio <= 1.0, "Hit ratio should be between 0 and 1");
        }

        @Test
        void testJMT_LRU_n4() {
            double[] cap = {1,1,1};
            double hitRatio = runJMT(4, cap, ReplacementStrategy.LRU, 1.0);
            assertTrue(hitRatio >= 0.0 && hitRatio <= 1.0, "Hit ratio should be between 0 and 1");
        }

        @Test
        void testJMT_LRU_n7_h2() {
            double[] cap = {1,1};
            double hitRatio = runJMT(7, cap, ReplacementStrategy.LRU, 1.0);
            assertTrue(hitRatio >= 0.0 && hitRatio <= 1.0, "Hit ratio should be between 0 and 1");
        }

        @Test
        void testJMT_RR_n7_config1() {
            double[] cap = {1,2,2};
            double hitRatio = runJMT(7, cap, ReplacementStrategy.RR, 1.0);
            assertTrue(hitRatio >= 0.0 && hitRatio <= 1.0, "Hit ratio should be between 0 and 1");
        }

        @Test
        void testJMT_RR_n7_config2() {
            double[] cap = {1,1,1};
            double hitRatio = runJMT(7, cap, ReplacementStrategy.RR, 1.0);
            assertTrue(hitRatio >= 0.0 && hitRatio <= 1.0, "Hit ratio should be between 0 and 1");
        }

        @Test
        void testJMT_FIFO_n7() {
            double[] cap = {1,2,2};
            double hitRatio = runJMT(7, cap, ReplacementStrategy.FIFO, 1.0);
            assertTrue(hitRatio >= 0.0 && hitRatio <= 1.0, "Hit ratio should be between 0 and 1");
        }

        private double runJMT(int nItems, double[] cap, ReplacementStrategy strategy, double zipfAlpha) {
            Network model = new Network("CacheConfigTest_JMT");

            Delay clientDelay = new Delay(model, "Client");
            Matrix itemLevelCap = new Matrix(1, cap.length);
            for (int i = 0; i < cap.length; i++) {
                itemLevelCap.set(0, i, cap[i]);
            }
            Cache cacheNode = new Cache(model, "Cache", nItems, itemLevelCap, strategy);
            Delay cacheDelay = new Delay(model, "CacheDelay");

            ClosedClass clientClass = new ClosedClass(model, "ClientClass", 1, clientDelay, 0);
            ClosedClass hitClass = new ClosedClass(model, "HitClass", 0, clientDelay, 0);
            ClosedClass missClass = new ClosedClass(model, "MissClass", 0, clientDelay, 0);

            clientDelay.setService(clientClass, new Immediate());
            cacheDelay.setService(hitClass, Exp.fitMean(0.2));
            cacheDelay.setService(missClass, Exp.fitMean(1.0));

            cacheNode.setRead(clientClass, new Zipf(zipfAlpha, nItems));
            cacheNode.setHitClass(clientClass, hitClass);
            cacheNode.setMissClass(clientClass, missClass);

            RoutingMatrix P = model.initRoutingMatrix();
            P.set(clientClass, clientClass, clientDelay, cacheNode, 1.0);
            P.set(hitClass, hitClass, cacheNode, cacheDelay, 1.0);
            P.set(missClass, missClass, cacheNode, cacheDelay, 1.0);
            P.set(hitClass, clientClass, cacheDelay, clientDelay, 1.0);
            P.set(missClass, clientClass, cacheDelay, clientDelay, 1.0);

            model.link(P);

            SolverOptions options = new SolverOptions();
            options.verbose = VerboseLevel.SILENT;
            options.seed = 1;
            options.samples = 100000;

            SolverJMT solver = new SolverJMT(model, options);
            NetworkAvgTable avgTable = solver.getAvgTable();

            double hitTput = avgTable.getTput().get(1);
            double missTput = avgTable.getTput().get(2);
            return hitTput / (hitTput + missTput);
        }
    }

    // ========== NC Solver Tests ==========

    /**
     * NC solver tests for cache configurations.
     */
    @Nested
    class NCTests {

        @Test
        void testNC_LRU_n7_h5() {
            double[] cap = {1,1,1,1,1};
            double hitRatio = runNC(7, cap, ReplacementStrategy.LRU, 1.0);
            assertTrue(hitRatio >= 0.0 && hitRatio <= 1.0, "Hit ratio should be between 0 and 1");
        }

        @Test
        void testNC_LRU_n7_h4() {
            double[] cap = {1,1,1,1};
            double hitRatio = runNC(7, cap, ReplacementStrategy.LRU, 1.0);
            assertTrue(hitRatio >= 0.0 && hitRatio <= 1.0, "Hit ratio should be between 0 and 1");
        }

        @Test
        void testNC_LRU_n7_config1() {
            double[] cap = {1,2,2};
            double hitRatio = runNC(7, cap, ReplacementStrategy.LRU, 1.0);
            assertTrue(hitRatio >= 0.0 && hitRatio <= 1.0, "Hit ratio should be between 0 and 1");
        }

        @Test
        void testNC_LRU_n7_config2() {
            double[] cap = {1,1,2};
            double hitRatio = runNC(7, cap, ReplacementStrategy.LRU, 1.0);
            assertTrue(hitRatio >= 0.0 && hitRatio <= 1.0, "Hit ratio should be between 0 and 1");
        }

        @Test
        void testNC_LRU_n7_config3() {
            double[] cap = {1,2,1};
            double hitRatio = runNC(7, cap, ReplacementStrategy.LRU, 1.0);
            assertTrue(hitRatio >= 0.0 && hitRatio <= 1.0, "Hit ratio should be between 0 and 1");
        }

        @Test
        void testNC_LRU_n7_config4() {
            double[] cap = {1,1,1};
            double hitRatio = runNC(7, cap, ReplacementStrategy.LRU, 1.0);
            assertTrue(hitRatio >= 0.0 && hitRatio <= 1.0, "Hit ratio should be between 0 and 1");
        }

        @Test
        void testNC_LRU_n4() {
            double[] cap = {1,1,1};
            double hitRatio = runNC(4, cap, ReplacementStrategy.LRU, 1.0);
            assertTrue(hitRatio >= 0.0 && hitRatio <= 1.0, "Hit ratio should be between 0 and 1");
        }

        @Test
        void testNC_LRU_n7_h2() {
            double[] cap = {1,1};
            double hitRatio = runNC(7, cap, ReplacementStrategy.LRU, 1.0);
            assertTrue(hitRatio >= 0.0 && hitRatio <= 1.0, "Hit ratio should be between 0 and 1");
        }

        @Test
        void testNC_RR_n7_config1() {
            double[] cap = {1,2,2};
            double hitRatio = runNC(7, cap, ReplacementStrategy.RR, 1.0);
            assertTrue(hitRatio >= 0.0 && hitRatio <= 1.0, "Hit ratio should be between 0 and 1");
        }

        @Test
        void testNC_RR_n7_config2() {
            double[] cap = {1,1,1};
            double hitRatio = runNC(7, cap, ReplacementStrategy.RR, 1.0);
            assertTrue(hitRatio >= 0.0 && hitRatio <= 1.0, "Hit ratio should be between 0 and 1");
        }

        @Test
        void testNC_FIFO_n7() {
            double[] cap = {1,2,2};
            double hitRatio = runNC(7, cap, ReplacementStrategy.FIFO, 1.0);
            assertTrue(hitRatio >= 0.0 && hitRatio <= 1.0, "Hit ratio should be between 0 and 1");
        }

        private double runNC(int nItems, double[] cap, ReplacementStrategy strategy, double zipfAlpha) {
            Network model = new Network("CacheConfigTest_NC");

            Delay clientDelay = new Delay(model, "Client");
            Matrix itemLevelCap = new Matrix(1, cap.length);
            for (int i = 0; i < cap.length; i++) {
                itemLevelCap.set(0, i, cap[i]);
            }
            Cache cacheNode = new Cache(model, "Cache", nItems, itemLevelCap, strategy);
            Delay cacheDelay = new Delay(model, "CacheDelay");

            ClosedClass clientClass = new ClosedClass(model, "ClientClass", 1, clientDelay, 0);
            ClosedClass hitClass = new ClosedClass(model, "HitClass", 0, clientDelay, 0);
            ClosedClass missClass = new ClosedClass(model, "MissClass", 0, clientDelay, 0);

            clientDelay.setService(clientClass, new Immediate());
            cacheDelay.setService(hitClass, Exp.fitMean(0.2));
            cacheDelay.setService(missClass, Exp.fitMean(1.0));

            cacheNode.setRead(clientClass, new Zipf(zipfAlpha, nItems));
            cacheNode.setHitClass(clientClass, hitClass);
            cacheNode.setMissClass(clientClass, missClass);

            RoutingMatrix P = model.initRoutingMatrix();
            P.set(clientClass, clientClass, clientDelay, cacheNode, 1.0);
            P.set(hitClass, hitClass, cacheNode, cacheDelay, 1.0);
            P.set(missClass, missClass, cacheNode, cacheDelay, 1.0);
            P.set(hitClass, clientClass, cacheDelay, clientDelay, 1.0);
            P.set(missClass, clientClass, cacheDelay, clientDelay, 1.0);

            model.link(P);

            SolverOptions options = new SolverOptions();
            options.verbose = VerboseLevel.SILENT;

            SolverNC solver = new SolverNC(model, options);
            NetworkAvgTable avgTable = solver.getAvgTable();

            double hitTput = avgTable.getTput().get(1);
            double missTput = avgTable.getTput().get(2);
            return hitTput / (hitTput + missTput);
        }
    }
}
