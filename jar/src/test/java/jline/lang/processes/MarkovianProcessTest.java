/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.VerboseLevel;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.solvers.NetworkAvgTable;
import jline.solvers.SolverOptions;
import jline.solvers.des.DESOptions;
import jline.solvers.des.SolverDES;
import jline.solvers.mva.MVAOptions;
import jline.solvers.mva.SolverMVA;
import jline.util.Maths;
import jline.util.RandomManager;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;

import java.util.List;
import java.util.Random;

import static jline.lib.butools.ph.CheckRAPRepresentationKt.checkRAPRepresentation;
import static org.junit.jupiter.api.Assertions.*;

/**
 * Comprehensive test suite for Markovian arrival processes.
 *
 * Merged from:
 * - BMAPTest: Batch Markovian Arrival Process tests
 * - RAPDistributionTest: Rational Arrival Process tests
 *
 * Test organization:
 * - BMAP Tests: Batch Markovian Arrival Process functionality
 * - RAP Tests: Rational Arrival Process functionality
 */
public class MarkovianProcessTest {

    private static final double TOLERANCE = 1e-6;
    private static final double LOOSE_TOLERANCE = 1e-3;
    private static final double SAMPLING_TOLERANCE = 0.1; // 10% tolerance for empirical samples

    @BeforeAll
    public static void setUp() {
        // Ensure MATLAB-compatible random number generation with fixed seed
        Maths.setRandomNumbersMatlab(true);
        RandomManager.setMasterSeed(23000);
    }

    // ==================== BMAP Tests ====================

    @Nested
    class BMAPTests {

        @Test
        public void testBMAPConstruction() {
            // Create simple 2-phase BMAP with batch sizes 1, 2, 3
            Matrix D0 = new Matrix(new double[][]{
                    {-2.0, 1.0},
                    {1.0, -2.0}
            });

            Matrix D1 = new Matrix(new double[][]{
                    {0.3, 0.2},
                    {0.2, 0.3}
            });

            Matrix D2 = new Matrix(new double[][]{
                    {0.2, 0.1},
                    {0.1, 0.2}
            });

            Matrix D3 = new Matrix(new double[][]{
                    {0.125, 0.075},
                    {0.075, 0.125}
            });

            BMAP bmap = new BMAP(D0, D1, D2, D3);

            assertNotNull(bmap);
            assertEquals(3, bmap.getMaxBatchSize());
            assertEquals(2, bmap.getNumberOfPhases());
        }

        @Test
        public void testBMAPFromMAPWithPMF() {
            // Create base MAP for inter-batch arrivals
            Matrix D0 = new Matrix(new double[][]{
                    {-2.0, 1.0},
                    {1.0, -2.0}
            });

            Matrix D1 = new Matrix(new double[][]{
                    {0.5, 0.5},
                    {0.5, 0.5}
            });

            // Batch size distribution
            int[] batchSizes = {1, 2, 3};
            double[] pmf = {0.3, 0.5, 0.2};

            BMAP bmap = BMAP.fromMAPWithBatchPMF(D0, D1, batchSizes, pmf);

            assertNotNull(bmap);
            assertEquals(3, bmap.getMaxBatchSize());

            // Verify mean batch size
            double meanBatchSize = 1 * 0.3 + 2 * 0.5 + 3 * 0.2;
            assertEquals(1.9, meanBatchSize, TOLERANCE);
        }

        @Test
        public void testBMAPBatchStatistics() {
            // Create BMAP with known batch size distribution
            Matrix D0 = new Matrix(new double[][]{
                    {-1.0, 0.0},
                    {0.0, -1.0}
            });

            Matrix D1 = new Matrix(new double[][]{
                    {0.3, 0.1},
                    {0.1, 0.3}
            });

            Matrix D2 = new Matrix(new double[][]{
                    {0.4, 0.2},
                    {0.2, 0.4}
            });

            BMAP bmap = new BMAP(D0, D1, D2);

            // Verify basic statistics
            assertNotNull(bmap);
            assertEquals(2, bmap.getMaxBatchSize());
        }

        @Test
        public void testBMAPSampling() {
            // Create simple BMAP and verify sampling works
            Matrix D0 = new Matrix(new double[][]{
                    {-2.0, 1.0},
                    {1.0, -2.0}
            });

            Matrix D1 = new Matrix(new double[][]{
                    {0.5, 0.5},
                    {0.5, 0.5}
            });

            BMAP bmap = new BMAP(D0, D1);

            // Sample from the BMAP
            RandomManager.setMasterSeed(42);
            Random random = RandomManager.getThreadRandomAsRandom();

            double[] samples = bmap.sample(1000, random);

            // Verify samples are non-negative
            for (double sample : samples) {
                assertTrue(sample >= 0, "BMAP samples should be non-negative");
            }

            // Compute empirical mean
            double empiricalMean = 0;
            for (double sample : samples) {
                empiricalMean += sample;
            }
            empiricalMean /= samples.length;

            // Verify empirical mean is close to theoretical mean
            double theoreticalMean = bmap.getMean();
            double relativeError = Math.abs(empiricalMean - theoreticalMean) / theoreticalMean;
            assertTrue(relativeError < SAMPLING_TOLERANCE,
                    "BMAP empirical mean should match theoretical mean within tolerance");
        }

        @Test
        public void testBMAPInQueueingNetwork() {
            // Test BMAP as arrival process in a queueing network
            Network model = new Network("BMAP/M/1");
            Source source = new Source(model, "Source");
            Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
            Sink sink = new Sink(model, "Sink");

            OpenClass jobClass = new OpenClass(model, "Class1");

            // Create BMAP arrival with mean inter-arrival time = 2.0, batch size = 1
            Matrix D0 = new Matrix(new double[][]{
                    {-0.5, 0.0},
                    {0.0, -0.5}
            });
            Matrix D1 = new Matrix(new double[][]{
                    {0.25, 0.25},
                    {0.25, 0.25}
            });
            BMAP bmap = new BMAP(D0, D1);

            source.setArrival(jobClass, bmap);
            queue.setService(jobClass, new Exp(1.0));

            RoutingMatrix P = model.initRoutingMatrix();
            P.set(jobClass, jobClass, source, queue, 1.0);
            P.set(jobClass, jobClass, queue, sink, 1.0);
            model.link(P);

            // Solve with DES
            DESOptions desOptions = new DESOptions();
            desOptions.samples = 100000;
            desOptions.seed = 23000;
            desOptions.verbose = VerboseLevel.SILENT;

            SolverDES desSolver = new SolverDES(model, desOptions);
            NetworkAvgTable desResult = desSolver.getAvgTable();

            assertNotNull(desResult);

            // Verify throughput matches arrival rate
            List<Double> tput = desResult.getTput();
            double queueTput = tput.get(1); // Queue throughput
            double expectedRate = 1.0 / bmap.getMean();
            double relativeError = Math.abs(queueTput - expectedRate) / expectedRate;
            assertTrue(relativeError < 0.1,
                    "BMAP throughput should match expected arrival rate");
        }

        @Test
        public void testBMAPValidation() {
            // Test that invalid BMAP matrices are rejected
            Matrix D0 = new Matrix(new double[][]{
                    {-2.0, 1.0},
                    {1.0, -2.0}
            });

            Matrix D1 = new Matrix(new double[][]{
                    {0.5, 0.5},
                    {0.5, 0.5}
            });

            Matrix D2_invalid = new Matrix(new double[][]{
                    {2.0, 1.0},  // Row sum > 1, should be rejected
                    {1.0, 2.0}
            });

            assertThrows(IllegalArgumentException.class, () -> {
                new BMAP(D0, D1, D2_invalid);
            });
        }
    }

    // ==================== RAP Tests ====================

    @Nested
    class RAPTests {

        @Test
        public void testRAPConstructionValid() {
            // Create valid RAP distribution (order 2)
            Matrix H0 = new Matrix(new double[][]{
                    {-2.0, 1.0},
                    {0.5, -1.5}
            });
            Matrix H1 = new Matrix(new double[][]{
                    {0.5, 0.5},
                    {0.5, 0.5}
            });

            RAP rap = new RAP(H0, H1);

            assertNotNull(rap);
            assertEquals(2, rap.getNumberOfPhases());
        }

        @Test
        public void testRAPConstructionInvalid() {
            // Test with H0 + H1 not forming valid generator (row sums != 0)
            Matrix H0 = new Matrix(new double[][]{
                    {-2.0, 1.0},
                    {0.5, -1.5}
            });
            Matrix H1 = new Matrix(new double[][]{
                    {0.5, 0.3}, // Row sum = 0.8, not 1.0
                    {0.5, 0.8}  // Row sum = 1.3, not 1.0
            });

            assertThrows(IllegalArgumentException.class, () -> {
                new RAP(H0, H1);
            });
        }

        @Test
        public void testRAPFromPoisson() {
            // Poisson process is simplest RAP (order 1)
            double rate = 2.0;
            RAP rap = RAP.fromPoisson(rate);

            assertNotNull(rap);
            assertEquals(1, rap.getNumberOfPhases());

            // Verify it's a valid RAP
            assertTrue(checkRAPRepresentation(rap.getH0(), rap.getH1(), 1e-14));

            // Mean should be 1/rate
            double mean = rap.getMean();
            assertEquals(1.0 / rate, mean, TOLERANCE);

            // SCV should be 1 (for exponential/Poisson)
            double scv = rap.getSCV();
            assertEquals(1.0, scv, TOLERANCE);
        }

        @Test
        public void testRAPFromErlang() {
            // Erlang-3 renewal process with rate 1.0
            int k = 3;
            double rate = 1.0;
            RAP rap = RAP.fromErlang(k, rate);

            assertNotNull(rap);
            assertEquals(k, rap.getNumberOfPhases());

            // Verify it's a valid RAP
            assertTrue(checkRAPRepresentation(rap.getH0(), rap.getH1(), 1e-14));

            // Mean should be k/rate
            double mean = rap.getMean();
            assertEquals((double) k / rate, mean, TOLERANCE);

            // SCV for Erlang-k should be 1/k
            double scv = rap.getSCV();
            assertEquals(1.0 / k, scv, LOOSE_TOLERANCE);
        }

        @Test
        public void testRAPFromHyperExp() {
            // Create Hyper-exponential-2 with specified parameters
            double p = 0.7;
            double lambda1 = 2.0;
            double lambda2 = 0.5;

            RAP rap = RAP.fromHyperExp(p, lambda1, lambda2);

            assertNotNull(rap);
            assertEquals(2, rap.getNumberOfPhases());

            // Verify it's a valid RAP
            assertTrue(checkRAPRepresentation(rap.getH0(), rap.getH1(), 1e-14));

            // Mean of HyperExp-2: p/lambda1 + (1-p)/lambda2
            double expectedMean = p / lambda1 + (1 - p) / lambda2;
            double mean = rap.getMean();
            assertEquals(expectedMean, mean, TOLERANCE);

            // SCV of HyperExp-2 should be > 1 (more variable than exponential)
            double scv = rap.getSCV();
            assertTrue(scv > 1.0, "HyperExp should have SCV > 1");
        }

        @Test
        public void testRAPSampling() {
            // Create RAP and verify sampling works
            Matrix H0 = new Matrix(new double[][]{
                    {-2.0, 1.0},
                    {0.5, -1.5}
            });
            Matrix H1 = new Matrix(new double[][]{
                    {0.5, 0.5},
                    {0.5, 0.5}
            });

            RAP rap = new RAP(H0, H1);

            // Sample from the RAP
            RandomManager.setMasterSeed(42);
            Random random = RandomManager.getThreadRandomAsRandom();

            double[] samples = rap.sample(1000, random);

            // Verify samples are non-negative
            for (double sample : samples) {
                assertTrue(sample >= 0, "RAP samples should be non-negative");
            }

            // Compute empirical mean
            double empiricalMean = 0;
            for (double sample : samples) {
                empiricalMean += sample;
            }
            empiricalMean /= samples.length;

            // Verify empirical mean is close to theoretical mean
            double theoreticalMean = rap.getMean();
            double relativeError = Math.abs(empiricalMean - theoreticalMean) / theoreticalMean;
            assertTrue(relativeError < SAMPLING_TOLERANCE,
                    "RAP empirical mean should match theoretical mean within tolerance");
        }

        @Test
        public void testRAPInQueueingNetwork() {
            // Test RAP as arrival process in a queueing network
            Network model = new Network("RAP/M/1");
            Source source = new Source(model, "Source");
            Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
            Sink sink = new Sink(model, "Sink");

            OpenClass jobClass = new OpenClass(model, "Class1");

            // Create RAP arrival (Erlang-2 with rate 1.0)
            RAP rap = RAP.fromErlang(2, 1.0);
            source.setArrival(jobClass, rap);
            queue.setService(jobClass, new Exp(2.0));

            RoutingMatrix P = model.initRoutingMatrix();
            P.set(jobClass, jobClass, source, queue, 1.0);
            P.set(jobClass, jobClass, queue, sink, 1.0);
            model.link(P);

            // Solve with DES
            DESOptions desOptions = new DESOptions();
            desOptions.samples = 100000;
            desOptions.seed = 23000;
            desOptions.verbose = VerboseLevel.SILENT;

            SolverDES desSolver = new SolverDES(model, desOptions);
            NetworkAvgTable desResult = desSolver.getAvgTable();

            assertNotNull(desResult);

            // Verify throughput matches arrival rate
            List<Double> tput = desResult.getTput();
            double queueTput = tput.get(1); // Queue throughput
            double expectedRate = 1.0 / rap.getMean();
            double relativeError = Math.abs(queueTput - expectedRate) / expectedRate;
            assertTrue(relativeError < 0.1,
                    "RAP throughput should match expected arrival rate");
        }

        @Test
        public void testRAPMomentMatching() {
            // Test that fitMeanAndSCV creates valid RAP with correct moments
            double targetMean = 2.0;
            double targetSCV = 0.5; // Less variable than exponential

            RAP rap = RAP.fitMeanAndSCV(targetMean, targetSCV);

            assertNotNull(rap);

            // Verify it's a valid RAP
            assertTrue(checkRAPRepresentation(rap.getH0(), rap.getH1(), 1e-14));

            // Verify mean matches
            double mean = rap.getMean();
            assertEquals(targetMean, mean, LOOSE_TOLERANCE);

            // Verify SCV matches
            double scv = rap.getSCV();
            assertEquals(targetSCV, scv, LOOSE_TOLERANCE);
        }

        @Test
        public void testRAPValidation() {
            // Test validation of RAP representation
            Matrix H0_valid = new Matrix(new double[][]{
                    {-2.0, 1.0},
                    {0.5, -1.5}
            });
            Matrix H1_valid = new Matrix(new double[][]{
                    {0.5, 0.5},
                    {0.5, 0.5}
            });

            // Valid RAP should pass validation
            assertTrue(checkRAPRepresentation(H0_valid, H1_valid, 1e-14));

            // Invalid RAP (H0 + H1 row sums != 0) should fail validation
            Matrix H0_invalid = new Matrix(new double[][]{
                    {-2.0, 1.0},
                    {0.5, -1.5}
            });
            Matrix H1_invalid = new Matrix(new double[][]{
                    {0.5, 0.3},
                    {0.5, 0.8}
            });

            assertFalse(checkRAPRepresentation(H0_invalid, H1_invalid, 1e-14));
        }
    }
}
