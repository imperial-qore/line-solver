package jline.lang.nodes;

import java.util.ArrayList;
import java.util.List;

import jline.lang.ClosedClass;
import jline.lang.JobClass;
import jline.lang.OpenClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import static jline.TestTools.*;
import jline.lang.constant.ReplacementStrategy;
import jline.lang.processes.Exp;
import jline.lang.processes.Immediate;
import jline.lang.processes.DiscreteSampler;
import jline.lang.processes.MarkedMAP;
import jline.lang.processes.Zipf;
import jline.lang.nodes.Source;
import jline.lang.nodes.Sink;
import jline.solvers.mva.SolverMVA;
import jline.lang.constant.SolverType;
import jline.util.matrix.MatrixCell;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.ldes.SolverLDES;
import jline.solvers.wrappers.jmt.SolverJMT;
import jline.solvers.nc.SolverNC;
import jline.solvers.ssa.SolverSSA;
import jline.solvers.NetworkAvgTable;
import jline.solvers.SolverOptions;
import jline.VerboseLevel;
import jline.util.matrix.Matrix;


import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

/**
 * Unit tests for specific cache configurations with exact CTMC hit ratios.
 * Tests verify SSA accuracy against CTMC values and exact CTMC hit ratios.
 * Based on gettingstarted_ex6 with closed-loop arrival patterns.
 */
class CacheListTest {

    private static final double TOLERANCE_EXACT = ZERO_TOL; // Exact tolerance for CTMC hit ratios

    @Test
    @Tag("slow") // ~255s
    void testLRU_n7_h4_alpha1_1e4() {
        // LRU (n=7, h=4, alpha=1.0, 1e4): 1,1,1,1 -> 0.752582801392053
        Matrix itemLevelCap = new Matrix(1, 4);
        for (int i = 0; i < 4; i++) {
            itemLevelCap.set(0, i, 1);
        }
        
        Network model = createCacheModel(7, itemLevelCap, ReplacementStrategy.LRU, 1.0);
        
        double expectedHitRatio = 0.752582801392053;
        testCacheConfiguration(model, expectedHitRatio, 10000, "LRU n=7, h=4, alpha=1.0", COARSE_TOL);
    }


    @Test
    @Tag("slow") // ~261s
    void testLRU_n7_alpha1_config2() {
        // LRU (n=7, alpha=1.0): 1,1,2 -> 0.752411031903945
        Matrix itemLevelCap = new Matrix(1, 3);
        itemLevelCap.set(0, 0, 1);
        itemLevelCap.set(0, 1, 1);
        itemLevelCap.set(0, 2, 2);
        
        Network model = createCacheModel(7, itemLevelCap, ReplacementStrategy.LRU, 1.0);
        
        double expectedHitRatio = 0.752411031903945;
        testCacheConfiguration(model, expectedHitRatio, 10000, "LRU n=7, alpha=1.0 (1,1,2)", COARSE_TOL);
    }

    @Test
    @Tag("slow") // ~267s
    void testLRU_n7_alpha1_config3() {
        // LRU (n=7, alpha=1.0): 1,2,1 -> 0.748678893424558
        Matrix itemLevelCap = new Matrix(1, 3);
        itemLevelCap.set(0, 0, 1);
        itemLevelCap.set(0, 1, 2);
        itemLevelCap.set(0, 2, 1);
        
        Network model = createCacheModel(7, itemLevelCap, ReplacementStrategy.LRU, 1.0);
        
        double expectedHitRatio = 0.748678893424558;
        testCacheConfiguration(model, expectedHitRatio, 10000, "LRU n=7, alpha=1.0 (1,2,1)", COARSE_TOL);
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
    void testRR_n7_alpha1_config2() {
        // RR (n=7, alpha=1.0): 1,1,1 -> 0.626935646456052
        Matrix itemLevelCap = new Matrix(1, 3);
        itemLevelCap.set(0, 0, 1);
        itemLevelCap.set(0, 1, 1);
        itemLevelCap.set(0, 2, 1);

        Network model = createCacheModel(7, itemLevelCap, ReplacementStrategy.RR, 1.0);

        double expectedHitRatio = 0.626935646456052;
        // RR/CLIMB/QLRU cache hit sequences are strongly autocorrelated (the cache
        // contents are sticky), so the SSA hit-ratio estimator has a small effective
        // sample size. At 1e4 the seed=1 realization sat at rel err ~0.066, just over
        // the 0.05 gate; 5e5 converges it to <0.01 (RUN-17). Do not lower this.
        testCacheConfiguration(model, expectedHitRatio, 500000, "RR n=7, alpha=1.0 (1,1,1)");
    }

    @Test
    void testHLRU_n6_alpha12_caps21() {
        // h-LRU / LRU(m) (n=6, alpha=1.2): 2,1 -> 0.700849020321130
        // (MATLAB SolverCTMC reference, 2026-07-12)
        Matrix itemLevelCap = new Matrix(1, 2);
        itemLevelCap.set(0, 0, 2);
        itemLevelCap.set(0, 1, 1);

        Network model = createCacheModel(6, itemLevelCap, ReplacementStrategy.HLRU, 1.2);

        double expectedHitRatio = 0.700849020321130;
        testCacheConfiguration(model, expectedHitRatio, 10000, "HLRU n=6, alpha=1.2 (2,1)");
    }

    @Test
    void testCLIMB_n5_alpha12_cap2() {
        // CLIMB (transposition) (n=5, alpha=1.2): 2 -> 0.597633165363664
        // (MATLAB SolverCTMC reference, 2026-07-12)
        Matrix itemLevelCap = new Matrix(1, 1);
        itemLevelCap.set(0, 0, 2);

        Network model = createCacheModel(5, itemLevelCap, ReplacementStrategy.CLIMB, 1.2);

        double expectedHitRatio = 0.597633165363664;
        // 5e5 samples: autocorrelated hit sequence needs a large N to converge the
        // SSA hit ratio within the 0.05 gate (RUN-17). Do not lower this.
        testCacheConfiguration(model, expectedHitRatio, 500000, "CLIMB n=5, alpha=1.2 (2)");
    }

    @Test
    void testQLRU_n5_alpha12_cap2_q05() {
        // q-LRU (n=5, alpha=1.2, q=0.5): 2 -> 0.575143053074336
        // (MATLAB SolverCTMC reference, 2026-07-12)
        Matrix itemLevelCap = new Matrix(1, 1);
        itemLevelCap.set(0, 0, 2);

        Network model = createCacheModel(5, itemLevelCap, ReplacementStrategy.QLRU, 1.2);
        // set the q-LRU admission probability on the cache node
        for (jline.lang.nodes.Node nd : model.getNodes()) {
            if (nd instanceof Cache) {
                ((Cache) nd).setAdmissionProb(0.5);
            }
        }

        double expectedHitRatio = 0.575143053074336;
        // 5e5 samples: autocorrelated hit sequence needs a large N to converge the
        // SSA hit ratio within the 0.05 gate (RUN-17). Do not lower this.
        testCacheConfiguration(model, expectedHitRatio, 500000, "QLRU n=5, alpha=1.2 (2), q=0.5");
    }

    @Test
    void testLRU_markedMAP_MVA() {
        // Marked MMAP source driving per-item classes of an LRU [2,1] cache:
        // phase 1 prefers items 1-2, phase 2 prefers items 3-4 (non-IRM).
        // CTMC exact 0.8108225223 and MVA (cache_ttl_lrum_map dispatch)
        // 0.8400817539; MATLAB/python references, 2026-07-12.
        int n = 4;
        Matrix caps = new Matrix(1, 2);
        caps.set(0, 0, 2);
        caps.set(0, 1, 1);
        Matrix D0 = new Matrix("[-2.1,0.1;0.2,-0.7]");
        Matrix D1 = new Matrix("[2.0,0.0;0.0,0.5]");
        double[] pfast = {0.4, 0.3, 0.2, 0.1};
        double[] pslow = {0.1, 0.2, 0.3, 0.4};
        MatrixCell cell = new MatrixCell(n + 2);
        cell.set(0, D0);
        cell.set(1, D1);
        for (int k = 0; k < n; k++) {
            Matrix Dk = new Matrix(2, 2);
            Dk.set(0, 0, D1.get(0, 0) * pfast[k]);
            Dk.set(1, 1, D1.get(1, 1) * pslow[k]);
            cell.set(2 + k, Dk);
        }
        MarkedMAP mmap = new MarkedMAP(cell, n);

        Network model = new Network("markedLRU");
        Source source = new Source(model, "Source");
        Cache cacheNode = new Cache(model, "Cache", n, caps, ReplacementStrategy.LRU);
        Sink sink = new Sink(model, "Sink");
        List<JobClass> readers = new ArrayList<JobClass>();
        for (int k = 0; k < n; k++) {
            readers.add(new OpenClass(model, "Item" + (k + 1), 0));
        }
        OpenClass hitClass = new OpenClass(model, "HitClass", 0);
        OpenClass missClass = new OpenClass(model, "MissClass", 0);
        source.setMarkedArrival(mmap, readers);
        RoutingMatrix P = model.initRoutingMatrix();
        for (int k = 0; k < n; k++) {
            Matrix onehot = new Matrix(n, 1);
            onehot.set(k, 0, 1.0);
            Matrix items = new Matrix(n, 1);
            for (int i = 0; i < n; i++) {
                items.set(i, 0, i + 1);
            }
            cacheNode.setRead(readers.get(k), new DiscreteSampler(onehot, items));
            cacheNode.setHitClass(readers.get(k), hitClass);
            cacheNode.setMissClass(readers.get(k), missClass);
            P.set(readers.get(k), readers.get(k), source, cacheNode, 1.0);
        }
        P.set(hitClass, hitClass, cacheNode, sink, 1.0);
        P.set(missClass, missClass, cacheNode, sink, 1.0);
        model.link(P);

        SolverOptions mvaOpt = new SolverOptions(SolverType.MVA);
        mvaOpt.verbose = VerboseLevel.SILENT;
        SolverMVA mva = new SolverMVA(model, mvaOpt);
        jline.solvers.NetworkAvgNodeTable t = mva.getAvgNodeTable();
        double hit = 0.0, miss = 0.0;
        for (int i = 0; i < t.getNodeNames().size(); i++) {
            if (t.getNodeNames().get(i).equals("Cache")) {
                if (t.getClassNames().get(i).equals("HitClass")) {
                    hit = t.getTput().get(i);
                } else if (t.getClassNames().get(i).equals("MissClass")) {
                    miss = t.getTput().get(i);
                }
            }
        }
        assertEquals(0.8400817539, hit / (hit + miss), 1e-6,
                "MVA marked-MAP LRU hit ratio (cache_ttl_lrum_map)");
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
        testCacheConfiguration(model, expectedHitRatio, ssaSamples, testDescription, TOLERANCE_EXACT);
    }

    private void testCacheConfiguration(Network model, double expectedHitRatio, int ssaSamples, String testDescription, double ctmcTolerance) {


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

                // Verify CTMC gives the expected hit ratio within tolerance
                assertEquals(expectedHitRatio, ctmcHitRatio, ctmcTolerance,
                    String.format("CTMC hit ratio should be %.15f, got %.15f", expectedHitRatio, ctmcHitRatio));
                
                
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
                assertTrue(relativeError <= LOOSE_COARSE_TOL,
                    String.format("SSA hit ratio relative error %.4f exceeds tolerance %.2f (SSA: %.6f, CTMC: %.6f)", 
                        relativeError, LOOSE_COARSE_TOL, ssaHitRatio, ctmcHitRatio));

                
            } else {
                // Just verify SSA gives a reasonable hit ratio when CTMC failed
                assertTrue(ssaHitRatio >= 0.0 && ssaHitRatio <= 1.0, 
                    String.format("SSA hit ratio %.6f should be between 0 and 1", ssaHitRatio));
                
                
            }
            

        } catch (Exception e) {
            fail("Test failed with exception: " + e.getMessage());
        }
    }

    // ========== LDES Solver Tests ==========

    /**
     * LDES simulation tests for cache configurations.
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

            SolverLDES solver = new SolverLDES(model, options);
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

        // NC supports only RR/FIFO replacement (the cache_prob_erec/cache_miss_spm
        // algorithms are exact for the RR==FIFO family). LRU is out of contract and
        // rejected at the featSet gate, so runNC must throw. LRU hit ratios are
        // covered by the DES/JMT/CTMC groups above.
        @Test
        void testNC_LRU_n7_h5() {
            double[] cap = {1,1,1,1,1};
            assertThrows(RuntimeException.class, () -> runNC(7, cap, ReplacementStrategy.LRU, 1.0),
                    "NC does not support LRU replacement (RR/FIFO only)");
        }

        @Test
        void testNC_LRU_n7_h4() {
            double[] cap = {1,1,1,1};
            assertThrows(RuntimeException.class, () -> runNC(7, cap, ReplacementStrategy.LRU, 1.0),
                    "NC does not support LRU replacement (RR/FIFO only)");
        }

        @Test
        void testNC_LRU_n7_config1() {
            double[] cap = {1,2,2};
            assertThrows(RuntimeException.class, () -> runNC(7, cap, ReplacementStrategy.LRU, 1.0),
                    "NC does not support LRU replacement (RR/FIFO only)");
        }

        @Test
        void testNC_LRU_n7_config2() {
            double[] cap = {1,1,2};
            assertThrows(RuntimeException.class, () -> runNC(7, cap, ReplacementStrategy.LRU, 1.0),
                    "NC does not support LRU replacement (RR/FIFO only)");
        }

        @Test
        void testNC_LRU_n7_config3() {
            double[] cap = {1,2,1};
            assertThrows(RuntimeException.class, () -> runNC(7, cap, ReplacementStrategy.LRU, 1.0),
                    "NC does not support LRU replacement (RR/FIFO only)");
        }

        @Test
        void testNC_LRU_n7_config4() {
            double[] cap = {1,1,1};
            assertThrows(RuntimeException.class, () -> runNC(7, cap, ReplacementStrategy.LRU, 1.0),
                    "NC does not support LRU replacement (RR/FIFO only)");
        }

        @Test
        void testNC_LRU_n7_h2() {
            double[] cap = {1,1};
            assertThrows(RuntimeException.class, () -> runNC(7, cap, ReplacementStrategy.LRU, 1.0),
                    "NC does not support LRU replacement (RR/FIFO only)");
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
