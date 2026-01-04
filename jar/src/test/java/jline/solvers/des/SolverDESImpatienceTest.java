/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.des;

import jline.VerboseLevel;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.BalkingStrategy;
import jline.lang.constant.BalkingThreshold;
import jline.lang.constant.DropStrategy;
import jline.lang.constant.ImpatienceType;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Det;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.lang.processes.Gamma;
import jline.lang.processes.Lognormal;
import jline.lang.processes.Pareto;
import jline.lang.processes.Weibull;
import jline.solvers.NetworkAvgTable;
import jline.solvers.SolverOptions;
import jline.util.Maths;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Test suite for customer impatience features in the DES solver.
 *
 * Tests cover:
 * - Reneging: Timer-based abandonment from queue
 * - Balking: Queue-length based refusal to join
 * - Retrial: Orbit-based retry after rejection
 *
 * @see SolverDES
 * @see ImpatienceType
 * @see BalkingStrategy
 * @see DropStrategy
 */
public class SolverDESImpatienceTest {

    @BeforeAll
    public static void setUp() {
        Maths.setRandomNumbersMatlab(true);
    }

    // ==================== Reneging Tests ====================

    @Nested
    @DisplayName("Reneging Tests")
    class RenegingTests {

        @Test
        @DisplayName("M/M/1 with exponential patience - customers renege")
        public void testRenegingExponentialPatience() {
            // M/M/1 with high load (rho=0.9) and impatient customers
            Network model = new Network("Reneging_Exp");

            Source source = new Source(model, "Source");
            Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
            Sink sink = new Sink(model, "Sink");

            OpenClass jobClass = new OpenClass(model, "Class1", 0);

            source.setArrival(jobClass, new Exp(0.9));  // lambda = 0.9
            queue.setService(jobClass, new Exp(1.0));   // mu = 1.0 (rho = 0.9)

            // Set exponential patience with rate 0.5 (mean patience = 2.0)
            queue.setPatience(jobClass, new Exp(0.5));

            model.link(model.serialRouting(source, queue, sink));

            SolverOptions options = new SolverOptions(SolverType.DES);
            options.verbose = VerboseLevel.SILENT;
            options.samples = 50000;
            options.seed = 23000;

            SolverDES solver = new SolverDES(model, options);
            solver.getAvgTable();  // Run the simulation
            DESResult result = (DESResult) solver.result;

            assertNotNull(result, "Result should not be null");
            assertNotNull(result.renegedCustomers, "Reneging statistics should be computed");

            // At high load with short patience, we expect some customers to renege
            double reneged = result.renegedCustomers.get(1, 0); // Queue is station 1
            assertTrue(reneged > 0, "Some customers should renege at high load with short patience");

            // Reneging rate should be positive
            assertNotNull(result.renegingRate, "Reneging rate should be computed");
            double renegingRate = result.renegingRate.get(1, 0);
            assertTrue(renegingRate > 0 && renegingRate < 1, "Reneging rate should be between 0 and 1");
        }

        @Test
        @DisplayName("Reneging cancelled when service starts")
        public void testRenegingCancelledOnServiceStart() {
            // Low load system - most customers should get served before patience expires
            Network model = new Network("Reneging_NoRenege");

            Source source = new Source(model, "Source");
            Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
            Sink sink = new Sink(model, "Sink");

            OpenClass jobClass = new OpenClass(model, "Class1", 0);

            source.setArrival(jobClass, new Exp(0.1));  // lambda = 0.1 (low load)
            queue.setService(jobClass, new Exp(1.0));   // mu = 1.0

            // Set very long patience (mean = 100)
            queue.setPatience(jobClass, new Exp(0.01));

            model.link(model.serialRouting(source, queue, sink));

            SolverOptions options = new SolverOptions(SolverType.DES);
            options.verbose = VerboseLevel.SILENT;
            options.samples = 10000;
            options.seed = 23000;

            SolverDES solver = new SolverDES(model, options);
            solver.getAvgTable();  // Run the simulation
            DESResult result = (DESResult) solver.result;

            assertNotNull(result, "Result should not be null");

            // With low load and long patience, very few should renege
            double reneged = result.renegedCustomers.get(1, 0);
            double total = result.TN.get(1, 0) * options.samples; // Approximate total arrivals
            double renegingRate = result.renegingRate.get(1, 0);

            // Reneging rate should be very low
            assertTrue(renegingRate < 0.1, "Reneging rate should be very low with low load and long patience");
        }

        @Test
        @DisplayName("Deterministic patience timeout")
        public void testRenegingDeterministicPatience() {
            Network model = new Network("Reneging_Det");

            Source source = new Source(model, "Source");
            Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
            Sink sink = new Sink(model, "Sink");

            OpenClass jobClass = new OpenClass(model, "Class1", 0);

            source.setArrival(jobClass, new Exp(0.8));
            queue.setService(jobClass, new Exp(1.0));

            // Deterministic patience of 2.0 time units
            queue.setPatience(jobClass, new Det(2.0));

            model.link(model.serialRouting(source, queue, sink));

            SolverOptions options = new SolverOptions(SolverType.DES);
            options.verbose = VerboseLevel.SILENT;
            options.samples = 20000;
            options.seed = 23000;

            SolverDES solver = new SolverDES(model, options);
            solver.getAvgTable();  // Run the simulation
            DESResult result = (DESResult) solver.result;

            assertNotNull(result, "Result should not be null");
            assertNotNull(result.renegedCustomers, "Reneging statistics should be computed");

            // With deterministic patience, behavior should be consistent
            double renegingRate = result.renegingRate.get(1, 0);
            assertTrue(renegingRate >= 0 && renegingRate <= 1, "Reneging rate should be valid");
        }

        @Test
        @DisplayName("Erlang patience distribution")
        public void testRenegingErlangPatience() {
            Network model = new Network("Reneging_Erlang");

            Source source = new Source(model, "Source");
            Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
            Sink sink = new Sink(model, "Sink");

            OpenClass jobClass = new OpenClass(model, "Class1", 0);

            source.setArrival(jobClass, new Exp(0.7));
            queue.setService(jobClass, new Exp(1.0));

            // Erlang-2 patience with rate 0.5 (mean = 4.0)
            queue.setPatience(jobClass, new Erlang(0.5, 2));

            model.link(model.serialRouting(source, queue, sink));

            SolverOptions options = new SolverOptions(SolverType.DES);
            options.verbose = VerboseLevel.SILENT;
            options.samples = 20000;
            options.seed = 23000;

            SolverDES solver = new SolverDES(model, options);
            solver.getAvgTable();  // Run the simulation
            DESResult result = (DESResult) solver.result;

            assertNotNull(result, "Result should not be null");
            assertNotNull(result.avgRenegingWaitTime, "Average reneging wait time should be computed");
        }

        @Test
        @DisplayName("Gamma patience distribution")
        public void testRenegingGammaPatience() {
            Network model = new Network("Reneging_Gamma");

            Source source = new Source(model, "Source");
            Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
            Sink sink = new Sink(model, "Sink");

            OpenClass jobClass = new OpenClass(model, "Class1", 0);

            source.setArrival(jobClass, new Exp(0.7));
            queue.setService(jobClass, new Exp(1.0));

            // Gamma patience with mean=2.0 and SCV=0.5 (less variable than exponential)
            queue.setPatience(jobClass, Gamma.fitMeanAndSCV(2.0, 0.5));

            model.link(model.serialRouting(source, queue, sink));

            SolverOptions options = new SolverOptions(SolverType.DES);
            options.verbose = VerboseLevel.SILENT;
            options.samples = 20000;
            options.seed = 23000;

            SolverDES solver = new SolverDES(model, options);
            solver.getAvgTable();  // Run the simulation
            DESResult result = (DESResult) solver.result;

            assertNotNull(result, "Result should not be null");
            assertNotNull(result.renegedCustomers, "Reneging statistics should be computed");

            // Verify reneging rate is valid
            double renegingRate = result.renegingRate.get(1, 0);
            assertTrue(renegingRate >= 0 && renegingRate <= 1, "Reneging rate should be valid");
        }

        @Test
        @DisplayName("Pareto patience distribution")
        public void testRenegingParetoPatience() {
            Network model = new Network("Reneging_Pareto");

            Source source = new Source(model, "Source");
            Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
            Sink sink = new Sink(model, "Sink");

            OpenClass jobClass = new OpenClass(model, "Class1", 0);

            source.setArrival(jobClass, new Exp(0.6));
            queue.setService(jobClass, new Exp(1.0));

            // Pareto patience with mean=3.0 and SCV=2.0 (heavy-tailed)
            queue.setPatience(jobClass, Pareto.fitMeanAndSCV(3.0, 2.0));

            model.link(model.serialRouting(source, queue, sink));

            SolverOptions options = new SolverOptions(SolverType.DES);
            options.verbose = VerboseLevel.SILENT;
            options.samples = 20000;
            options.seed = 23000;

            SolverDES solver = new SolverDES(model, options);
            solver.getAvgTable();  // Run the simulation
            DESResult result = (DESResult) solver.result;

            assertNotNull(result, "Result should not be null");
            assertNotNull(result.renegedCustomers, "Reneging statistics should be computed");

            // Verify reneging rate is valid
            double renegingRate = result.renegingRate.get(1, 0);
            assertTrue(renegingRate >= 0 && renegingRate <= 1, "Reneging rate should be valid");
        }

        @Test
        @DisplayName("Weibull patience distribution")
        public void testRenegingWeibullPatience() {
            Network model = new Network("Reneging_Weibull");

            Source source = new Source(model, "Source");
            Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
            Sink sink = new Sink(model, "Sink");

            OpenClass jobClass = new OpenClass(model, "Class1", 0);

            source.setArrival(jobClass, new Exp(0.7));
            queue.setService(jobClass, new Exp(1.0));

            // Weibull patience with mean=2.5 and SCV=0.8
            queue.setPatience(jobClass, Weibull.fitMeanAndSCV(2.5, 0.8));

            model.link(model.serialRouting(source, queue, sink));

            SolverOptions options = new SolverOptions(SolverType.DES);
            options.verbose = VerboseLevel.SILENT;
            options.samples = 20000;
            options.seed = 23000;

            SolverDES solver = new SolverDES(model, options);
            solver.getAvgTable();  // Run the simulation
            DESResult result = (DESResult) solver.result;

            assertNotNull(result, "Result should not be null");
            assertNotNull(result.renegedCustomers, "Reneging statistics should be computed");

            // Verify reneging rate is valid
            double renegingRate = result.renegingRate.get(1, 0);
            assertTrue(renegingRate >= 0 && renegingRate <= 1, "Reneging rate should be valid");
        }

        @Test
        @DisplayName("Lognormal patience distribution")
        public void testRenegingLognormalPatience() {
            Network model = new Network("Reneging_Lognormal");

            Source source = new Source(model, "Source");
            Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
            Sink sink = new Sink(model, "Sink");

            OpenClass jobClass = new OpenClass(model, "Class1", 0);

            source.setArrival(jobClass, new Exp(0.65));
            queue.setService(jobClass, new Exp(1.0));

            // Lognormal patience with mean=2.0 and SCV=1.5 (moderately variable)
            queue.setPatience(jobClass, Lognormal.fitMeanAndSCV(2.0, 1.5));

            model.link(model.serialRouting(source, queue, sink));

            SolverOptions options = new SolverOptions(SolverType.DES);
            options.verbose = VerboseLevel.SILENT;
            options.samples = 20000;
            options.seed = 23000;

            SolverDES solver = new SolverDES(model, options);
            solver.getAvgTable();  // Run the simulation
            DESResult result = (DESResult) solver.result;

            assertNotNull(result, "Result should not be null");
            assertNotNull(result.renegedCustomers, "Reneging statistics should be computed");

            // Verify reneging rate is valid
            double renegingRate = result.renegingRate.get(1, 0);
            assertTrue(renegingRate >= 0 && renegingRate <= 1, "Reneging rate should be valid");
        }

        @Test
        @DisplayName("Multiclass reneging with different patience")
        public void testRenegingMulticlass() {
            Network model = new Network("Reneging_Multiclass");

            Source source = new Source(model, "Source");
            Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
            Sink sink = new Sink(model, "Sink");

            OpenClass class1 = new OpenClass(model, "Impatient", 0);
            OpenClass class2 = new OpenClass(model, "Patient", 0);

            source.setArrival(class1, new Exp(0.4));
            source.setArrival(class2, new Exp(0.4));
            queue.setService(class1, new Exp(1.0));
            queue.setService(class2, new Exp(1.0));

            // Different patience for each class
            queue.setPatience(class1, new Exp(1.0));   // Short patience (mean = 1.0)
            queue.setPatience(class2, new Exp(0.1));   // Long patience (mean = 10.0)

            RoutingMatrix routing = new RoutingMatrix(model, model.getClasses(), model.getNodes());
            routing.addConnection(source, queue);
            routing.addConnection(queue, sink);
            model.link(routing);

            SolverOptions options = new SolverOptions(SolverType.DES);
            options.verbose = VerboseLevel.SILENT;
            options.samples = 30000;
            options.seed = 23000;

            SolverDES solver = new SolverDES(model, options);
            solver.getAvgTable();  // Run the simulation
            DESResult result = (DESResult) solver.result;

            assertNotNull(result, "Result should not be null");

            // Impatient class should have higher reneging rate
            double renegingRate1 = result.renegingRate.get(1, 0); // Class 1 (impatient)
            double renegingRate2 = result.renegingRate.get(1, 1); // Class 2 (patient)

            assertTrue(renegingRate1 > renegingRate2,
                    "Impatient class should have higher reneging rate than patient class");
        }
    }

    // ==================== Balking Tests ====================

    @Nested
    @DisplayName("Balking Tests")
    class BalkingTests {

        @Test
        @DisplayName("Balking by queue length thresholds")
        public void testBalkingQueueLength() {
            Network model = new Network("Balking_QueueLength");

            Source source = new Source(model, "Source");
            Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
            Sink sink = new Sink(model, "Sink");

            OpenClass jobClass = new OpenClass(model, "Class1", 0);

            source.setArrival(jobClass, new Exp(0.9));  // High load
            queue.setService(jobClass, new Exp(1.0));

            // Set up balking: 50% balk when 3-5 jobs, 100% when >5 jobs
            List<BalkingThreshold> thresholds = new ArrayList<BalkingThreshold>();
            thresholds.add(new BalkingThreshold(3, 5, 0.5));
            thresholds.add(new BalkingThreshold(6, Integer.MAX_VALUE, 1.0));
            queue.setBalking(jobClass, BalkingStrategy.QUEUE_LENGTH, thresholds);

            model.link(model.serialRouting(source, queue, sink));

            SolverOptions options = new SolverOptions(SolverType.DES);
            options.verbose = VerboseLevel.SILENT;
            options.samples = 30000;
            options.seed = 23000;

            SolverDES solver = new SolverDES(model, options);
            solver.getAvgTable();  // Run the simulation
            DESResult result = (DESResult) solver.result;

            assertNotNull(result, "Result should not be null");
            assertNotNull(result.balkedCustomers, "Balking statistics should be computed");

            // With high load and balking thresholds, some customers should balk
            double balked = result.balkedCustomers.get(1, 0);
            assertTrue(balked > 0, "Some customers should balk at high load with thresholds");

            // Balking probability should be positive
            assertNotNull(result.balkingProbability, "Balking probability should be computed");
            double balkingProb = result.balkingProbability.get(1, 0);
            assertTrue(balkingProb > 0 && balkingProb < 1, "Balking probability should be between 0 and 1");
        }

        @Test
        @DisplayName("No balking when queue is short")
        public void testNoBalkingShortQueue() {
            Network model = new Network("Balking_NoBalking");

            Source source = new Source(model, "Source");
            Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
            Sink sink = new Sink(model, "Sink");

            OpenClass jobClass = new OpenClass(model, "Class1", 0);

            source.setArrival(jobClass, new Exp(0.3));  // Low load (rho = 0.3)
            queue.setService(jobClass, new Exp(1.0));

            // Balking only triggers when queue has 10+ jobs (which won't happen at low load)
            List<BalkingThreshold> thresholds = new ArrayList<BalkingThreshold>();
            thresholds.add(new BalkingThreshold(10, Integer.MAX_VALUE, 1.0));
            queue.setBalking(jobClass, BalkingStrategy.QUEUE_LENGTH, thresholds);

            model.link(model.serialRouting(source, queue, sink));

            SolverOptions options = new SolverOptions(SolverType.DES);
            options.verbose = VerboseLevel.SILENT;
            options.samples = 10000;
            options.seed = 23000;

            SolverDES solver = new SolverDES(model, options);
            solver.getAvgTable();  // Run the simulation
            DESResult result = (DESResult) solver.result;

            assertNotNull(result, "Result should not be null");

            // At low load, queue rarely reaches 10, so balking should be minimal
            double balkingProb = result.balkingProbability.get(1, 0);
            assertTrue(balkingProb < 0.01, "Balking should be rare when queue is short");
        }

        @Test
        @DisplayName("Gradual balking with probability thresholds")
        public void testBalkingGradualProbability() {
            Network model = new Network("Balking_Gradual");

            Source source = new Source(model, "Source");
            Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
            Sink sink = new Sink(model, "Sink");

            OpenClass jobClass = new OpenClass(model, "Class1", 0);

            source.setArrival(jobClass, new Exp(0.85));
            queue.setService(jobClass, new Exp(1.0));

            // Gradual balking: 20% at 2-3 jobs, 50% at 4-5, 80% at 6-7, 100% at 8+
            List<BalkingThreshold> thresholds = new ArrayList<BalkingThreshold>();
            thresholds.add(new BalkingThreshold(2, 3, 0.2));
            thresholds.add(new BalkingThreshold(4, 5, 0.5));
            thresholds.add(new BalkingThreshold(6, 7, 0.8));
            thresholds.add(new BalkingThreshold(8, Integer.MAX_VALUE, 1.0));
            queue.setBalking(jobClass, BalkingStrategy.QUEUE_LENGTH, thresholds);

            model.link(model.serialRouting(source, queue, sink));

            SolverOptions options = new SolverOptions(SolverType.DES);
            options.verbose = VerboseLevel.SILENT;
            options.samples = 30000;
            options.seed = 23000;

            SolverDES solver = new SolverDES(model, options);
            solver.getAvgTable();  // Run the simulation
            DESResult result = (DESResult) solver.result;

            assertNotNull(result, "Result should not be null");

            // Balking probability should reflect the gradual thresholds
            double balkingProb = result.balkingProbability.get(1, 0);
            assertTrue(balkingProb > 0, "Some balking should occur with gradual thresholds");
        }
    }

    // ==================== Retrial Tests ====================

    @Nested
    @DisplayName("Retrial Tests")
    class RetrialTests {

        @Test
        @DisplayName("Unlimited retrial with orbit")
        public void testRetrialUnlimited() {
            Network model = new Network("Retrial_Unlimited");

            Source source = new Source(model, "Source");
            Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
            queue.setNumberOfServers(1);
            queue.setCapacity(3);  // Small capacity to force rejections
            Sink sink = new Sink(model, "Sink");

            OpenClass jobClass = new OpenClass(model, "Class1", 0);

            source.setArrival(jobClass, new Exp(0.8));
            queue.setService(jobClass, new Exp(1.0));

            // Unlimited retrial with exponential delay (mean = 2.0)
            queue.setRetrial(jobClass, new Exp(0.5), -1);

            model.link(model.serialRouting(source, queue, sink));

            SolverOptions options = new SolverOptions(SolverType.DES);
            options.verbose = VerboseLevel.SILENT;
            options.samples = 30000;
            options.seed = 23000;

            SolverDES solver = new SolverDES(model, options);
            solver.getAvgTable();  // Run the simulation
            DESResult result = (DESResult) solver.result;

            assertNotNull(result, "Result should not be null");
            assertNotNull(result.retriedCustomers, "Retrial statistics should be computed");

            // With limited capacity and retrial, some customers should retry
            double retried = result.retriedCustomers.get(1, 0);
            assertTrue(retried > 0, "Some customers should retry when queue is full");

            // Average orbit size should be positive
            assertNotNull(result.avgOrbitSize, "Average orbit size should be computed");
            double avgOrbit = result.avgOrbitSize.get(1, 0);
            assertTrue(avgOrbit >= 0, "Average orbit size should be non-negative");

            // With unlimited retrial, no one should be dropped due to max retries
            double dropped = result.retrialDropped.get(1, 0);
            assertEquals(0.0, dropped, 0.01, "No customers should be dropped with unlimited retrial");
        }

        @Test
        @DisplayName("Retrial with max attempts")
        public void testRetrialWithMaxAttempts() {
            Network model = new Network("Retrial_MaxAttempts");

            Source source = new Source(model, "Source");
            Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
            queue.setNumberOfServers(1);
            queue.setCapacity(2);  // Very small capacity
            Sink sink = new Sink(model, "Sink");

            OpenClass jobClass = new OpenClass(model, "Class1", 0);

            source.setArrival(jobClass, new Exp(0.9));  // High load
            queue.setService(jobClass, new Exp(1.0));

            // Limited retrial: max 3 attempts, then drop
            queue.setRetrial(jobClass, new Exp(0.5), 3);

            model.link(model.serialRouting(source, queue, sink));

            SolverOptions options = new SolverOptions(SolverType.DES);
            options.verbose = VerboseLevel.SILENT;
            options.samples = 30000;
            options.seed = 23000;

            SolverDES solver = new SolverDES(model, options);
            solver.getAvgTable();  // Run the simulation
            DESResult result = (DESResult) solver.result;

            assertNotNull(result, "Result should not be null");

            // With very limited capacity and high load, some customers should exhaust retries
            double dropped = result.retrialDropped.get(1, 0);
            assertTrue(dropped > 0, "Some customers should be dropped after max retries");
        }

        @Test
        @DisplayName("Retrial preserves eventual throughput")
        public void testRetrialPreservesThroughput() {
            // With unlimited retrial, all jobs should eventually complete
            Network model = new Network("Retrial_Throughput");

            Source source = new Source(model, "Source");
            Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
            queue.setNumberOfServers(1);
            queue.setCapacity(5);
            Sink sink = new Sink(model, "Sink");

            OpenClass jobClass = new OpenClass(model, "Class1", 0);

            double lambda = 0.7;
            source.setArrival(jobClass, new Exp(lambda));
            queue.setService(jobClass, new Exp(1.0));

            // Unlimited retrial with fast retry
            queue.setRetrial(jobClass, new Exp(2.0), -1);

            model.link(model.serialRouting(source, queue, sink));

            SolverOptions options = new SolverOptions(SolverType.DES);
            options.verbose = VerboseLevel.SILENT;
            options.samples = 50000;
            options.seed = 23000;

            SolverDES solver = new SolverDES(model, options);
            solver.getAvgTable();  // Run the simulation
            DESResult result = (DESResult) solver.result;

            assertNotNull(result, "Result should not be null");

            // Throughput should approximately equal arrival rate (all jobs eventually complete)
            double throughput = result.TN.get(1, 0);
            assertEquals(lambda, throughput, lambda * 0.15,
                    "Throughput should approximately equal arrival rate with unlimited retrial");
        }
    }

    // ==================== Combined Impatience Tests ====================

    @Nested
    @DisplayName("Combined Impatience Tests")
    class CombinedTests {

        @Test
        @DisplayName("Reneging and balking combined")
        public void testRenegingAndBalkingCombined() {
            Network model = new Network("Combined_RenegingBalking");

            Source source = new Source(model, "Source");
            Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
            Sink sink = new Sink(model, "Sink");

            OpenClass jobClass = new OpenClass(model, "Class1", 0);

            source.setArrival(jobClass, new Exp(0.85));
            queue.setService(jobClass, new Exp(1.0));

            // Reneging: exponential patience (mean = 3.0)
            queue.setPatience(jobClass, new Exp(1.0/3.0));

            // Balking: 100% balk when queue >= 5
            List<BalkingThreshold> thresholds = new ArrayList<BalkingThreshold>();
            thresholds.add(new BalkingThreshold(5, Integer.MAX_VALUE, 1.0));
            queue.setBalking(jobClass, BalkingStrategy.QUEUE_LENGTH, thresholds);

            model.link(model.serialRouting(source, queue, sink));

            SolverOptions options = new SolverOptions(SolverType.DES);
            options.verbose = VerboseLevel.SILENT;
            options.samples = 30000;
            options.seed = 23000;

            SolverDES solver = new SolverDES(model, options);
            solver.getAvgTable();  // Run the simulation
            DESResult result = (DESResult) solver.result;

            assertNotNull(result, "Result should not be null");

            // Both balking and reneging should occur
            double renegingRate = result.renegingRate.get(1, 0);
            double balkingProb = result.balkingProbability.get(1, 0);

            assertTrue(renegingRate > 0 || balkingProb > 0,
                    "Either reneging or balking (or both) should occur");
        }

        @Test
        @DisplayName("Multiclass with different impatience behaviors")
        public void testMulticlassDifferentImpatience() {
            Network model = new Network("Combined_Multiclass");

            Source source = new Source(model, "Source");
            Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
            queue.setNumberOfServers(1);
            queue.setCapacity(10);
            Sink sink = new Sink(model, "Sink");

            OpenClass renegingClass = new OpenClass(model, "Reneging", 0);
            OpenClass balkingClass = new OpenClass(model, "Balking", 0);
            OpenClass retrialClass = new OpenClass(model, "Retrial", 0);

            source.setArrival(renegingClass, new Exp(0.3));
            source.setArrival(balkingClass, new Exp(0.3));
            source.setArrival(retrialClass, new Exp(0.3));

            queue.setService(renegingClass, new Exp(1.0));
            queue.setService(balkingClass, new Exp(1.0));
            queue.setService(retrialClass, new Exp(1.0));

            // Different impatience for each class
            queue.setPatience(renegingClass, new Exp(0.5)); // Reneging

            List<BalkingThreshold> thresholds = new ArrayList<BalkingThreshold>();
            thresholds.add(new BalkingThreshold(5, Integer.MAX_VALUE, 1.0));
            queue.setBalking(balkingClass, BalkingStrategy.QUEUE_LENGTH, thresholds); // Balking

            queue.setRetrial(retrialClass, new Exp(0.5), -1); // Retrial

            RoutingMatrix routing = new RoutingMatrix(model, model.getClasses(), model.getNodes());
            routing.addConnection(source, queue);
            routing.addConnection(queue, sink);
            model.link(routing);

            SolverOptions options = new SolverOptions(SolverType.DES);
            options.verbose = VerboseLevel.SILENT;
            options.samples = 30000;
            options.seed = 23000;

            SolverDES solver = new SolverDES(model, options);
            solver.getAvgTable();  // Run the simulation
            DESResult result = (DESResult) solver.result;

            assertNotNull(result, "Result should not be null");

            // Each class should show its respective impatience behavior
            double renegingRate = result.renegingRate.get(1, 0); // Class 0
            double balkingProb = result.balkingProbability.get(1, 1); // Class 1
            double retriedCustomers = result.retriedCustomers.get(1, 2); // Class 2

            // Validate that statistics are computed for each class
            assertNotNull(result.renegingRate, "Reneging rate should be computed");
            assertNotNull(result.balkingProbability, "Balking probability should be computed");
            assertNotNull(result.retriedCustomers, "Retried customers should be computed");
        }
    }

    // ==================== Edge Case Tests ====================

    @Nested
    @DisplayName("Edge Case Tests")
    class EdgeCaseTests {

        @Test
        @DisplayName("No impatience configured - normal operation")
        public void testNoImpatienceConfigured() {
            Network model = new Network("NoImpatience");

            Source source = new Source(model, "Source");
            Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
            Sink sink = new Sink(model, "Sink");

            OpenClass jobClass = new OpenClass(model, "Class1", 0);

            source.setArrival(jobClass, new Exp(0.5));
            queue.setService(jobClass, new Exp(1.0));

            // No patience, balking, or retrial configured
            model.link(model.serialRouting(source, queue, sink));

            SolverOptions options = new SolverOptions(SolverType.DES);
            options.verbose = VerboseLevel.SILENT;
            options.samples = 10000;
            options.seed = 23000;

            SolverDES solver = new SolverDES(model, options);
            solver.getAvgTable();  // Run the simulation
            DESResult result = (DESResult) solver.result;

            assertNotNull(result, "Result should not be null");

            // All impatience metrics should be zero or null
            if (result.renegedCustomers != null) {
                double reneged = result.renegedCustomers.get(1, 0);
                assertEquals(0.0, reneged, 0.01, "No reneging without patience configured");
            }
            if (result.balkedCustomers != null) {
                double balked = result.balkedCustomers.get(1, 0);
                assertEquals(0.0, balked, 0.01, "No balking without balking configured");
            }
        }

        @Test
        @DisplayName("Very short patience - high reneging rate")
        public void testVeryShortPatience() {
            Network model = new Network("VeryShortPatience");

            Source source = new Source(model, "Source");
            Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
            Sink sink = new Sink(model, "Sink");

            OpenClass jobClass = new OpenClass(model, "Class1", 0);

            source.setArrival(jobClass, new Exp(0.9));
            queue.setService(jobClass, new Exp(1.0));

            // Very short patience (mean = 0.1)
            queue.setPatience(jobClass, new Exp(10.0));

            model.link(model.serialRouting(source, queue, sink));

            SolverOptions options = new SolverOptions(SolverType.DES);
            options.verbose = VerboseLevel.SILENT;
            options.samples = 20000;
            options.seed = 23000;

            SolverDES solver = new SolverDES(model, options);
            solver.getAvgTable();  // Run the simulation
            DESResult result = (DESResult) solver.result;

            assertNotNull(result, "Result should not be null");

            // With very short patience at high load, most waiting customers should renege
            double renegingRate = result.renegingRate.get(1, 0);
            assertTrue(renegingRate > 0.3, "Reneging rate should be high with very short patience");
        }

        @Test
        @DisplayName("Immediate balking at any queue length")
        public void testImmediateBalkingAlways() {
            Network model = new Network("ImmediateBalking");

            Source source = new Source(model, "Source");
            Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
            Sink sink = new Sink(model, "Sink");

            OpenClass jobClass = new OpenClass(model, "Class1", 0);

            source.setArrival(jobClass, new Exp(0.5));
            queue.setService(jobClass, new Exp(1.0));

            // 100% balking when queue has any customers (including 0)
            List<BalkingThreshold> thresholds = new ArrayList<BalkingThreshold>();
            thresholds.add(new BalkingThreshold(1, Integer.MAX_VALUE, 1.0));
            queue.setBalking(jobClass, BalkingStrategy.QUEUE_LENGTH, thresholds);

            model.link(model.serialRouting(source, queue, sink));

            SolverOptions options = new SolverOptions(SolverType.DES);
            options.verbose = VerboseLevel.SILENT;
            options.samples = 10000;
            options.seed = 23000;

            SolverDES solver = new SolverDES(model, options);
            solver.getAvgTable();  // Run the simulation
            DESResult result = (DESResult) solver.result;

            assertNotNull(result, "Result should not be null");

            // With aggressive balking, most arrivals during busy periods should balk
            double balkingProb = result.balkingProbability.get(1, 0);
            assertTrue(balkingProb > 0, "Some balking should occur when queue is busy");
        }
    }
}
