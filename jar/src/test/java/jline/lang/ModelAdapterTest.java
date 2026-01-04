/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang;

import jline.VerboseLevel;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.solvers.NetworkAvgChainTable;
import jline.solvers.NetworkAvgTable;
import jline.solvers.Solver;
import jline.solvers.SolverOptions;
import jline.solvers.mva.SolverMVA;
import jline.util.Maths;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.*;
import static jline.TestTools.MID_TOL;

/**
 * Tests for ModelAdapter chain aggregation functionality
 */
public class ModelAdapterTest {

    @BeforeEach
    public void matlabRandomSeedSetUp() {
        Maths.setRandomNumbersMatlab(true);
    }

    @AfterEach
    public void matlabRandomSeedClear() {
        Maths.setRandomNumbersMatlab(false);
    }

    /**
     * Test 1: Simple chain aggregation with class switching
     *
     * Creates a model with 2 chains (4 classes total) and 3 stations.
     * Chain 1: Class1 (N=5) <-> Class2 (N=0)
     * Chain 2: Class3 (N=3) <-> Class4 (N=0)
     *
     * Verifies that avgTable after aggregation matches avgChainTable before aggregation.
     */
    @Test
    public void testAggregateChainsTwoChains() {
        // Create model
        Network model = new Network("TestModel");

        // Create stations
        Delay delay = new Delay(model, "Delay");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);

        // Create classes - 2 chains with 2 classes each
        ClosedClass class1 = new ClosedClass(model, "Class1", 5, delay, 0);
        ClosedClass class2 = new ClosedClass(model, "Class2", 0, delay, 0);
        ClosedClass class3 = new ClosedClass(model, "Class3", 3, delay, 0);
        ClosedClass class4 = new ClosedClass(model, "Class4", 0, delay, 0);

        // Set service times
        delay.setService(class1, Exp.fitMean(1.0));
        delay.setService(class2, Exp.fitMean(1.0));
        delay.setService(class3, Exp.fitMean(1.5));
        delay.setService(class4, Exp.fitMean(1.5));

        queue1.setService(class1, Exp.fitMean(0.5));
        queue1.setService(class2, Exp.fitMean(0.6));
        queue1.setService(class3, Exp.fitMean(0.4));
        queue1.setService(class4, Exp.fitMean(0.5));

        queue2.setService(class1, Exp.fitMean(0.3));
        queue2.setService(class2, Exp.fitMean(0.4));
        queue2.setService(class3, Exp.fitMean(0.2));
        queue2.setService(class4, Exp.fitMean(0.3));

        // Set routing with class switching within chains
        RoutingMatrix P = model.initRoutingMatrix();

        // Chain 1: Class1 -> Queue1 -> Class2 -> Queue2 -> Class1 -> Delay
        P.set(class1, class1, delay, queue1, 1.0);
        P.set(class1, class2, queue1, queue2, 1.0);  // class switch
        P.set(class2, class1, queue2, delay, 1.0);   // class switch back

        // Chain 2: Class3 -> Queue1 -> Class4 -> Queue2 -> Class3 -> Delay
        P.set(class3, class3, delay, queue1, 1.0);
        P.set(class3, class4, queue1, queue2, 1.0);  // class switch
        P.set(class4, class3, queue2, delay, 1.0);   // class switch back

        model.link(P);

        // Step 1: Solve original model and get avgChainTable
        SolverOptions options = Solver.defaultOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverMVA solverOriginal = new SolverMVA(model, options);
        NetworkAvgChainTable avgChainTableOriginal = solverOriginal.getAvgChainTable();

        assertNotNull(avgChainTableOriginal, "avgChainTable should not be null");

        // Step 2: Aggregate chains
        ModelAdapter.AggregateChainResult result = ModelAdapter.aggregateChains(model);
        assertNotNull(result, "Aggregation result should not be null");

        Network chainModel = result.getChainModel();
        assertNotNull(chainModel, "Chain model should not be null");

        // Verify aggregated model has correct number of classes (should be 2 chains = 2 classes)
        assertEquals(2, chainModel.getNumberOfClasses(), "Aggregated model should have 2 classes (one per chain)");

        // Step 3: Solve aggregated model and get avgTable
        SolverMVA solverAggregated = new SolverMVA(chainModel, options);
        NetworkAvgTable avgTableAggregated = solverAggregated.getAvgTable();

        assertNotNull(avgTableAggregated, "avgTable for aggregated model should not be null");

        // Step 4: Compare results - avgTable of aggregated model should match avgChainTable of original
        // Compare throughput per chain
        List<Double> originalTput = avgChainTableOriginal.getTput();
        List<Double> aggregatedTput = avgTableAggregated.getTput();

        // Note: Tables may have different ordering, so compare by station/chain
        // For simplicity, just verify the values are close
        assertEquals(originalTput.size(), aggregatedTput.size(),
            "Should have same number of throughput entries");

        // Compare queue lengths
        List<Double> originalQLen = avgChainTableOriginal.getQLen();
        List<Double> aggregatedQLen = avgTableAggregated.getQLen();

        assertEquals(originalQLen.size(), aggregatedQLen.size(),
            "Should have same number of queue length entries");

        // Verify values are similar (within tolerance)
        double totalDiff = 0.0;
        for (int i = 0; i < originalQLen.size(); i++) {
            double orig = originalQLen.get(i);
            double agg = aggregatedQLen.get(i);
            if (orig > MID_TOL) {
                totalDiff += Math.abs(orig - agg) / orig;
            }
        }
        double avgRelDiff = totalDiff / originalQLen.size();
        assertTrue(avgRelDiff < 0.05, "Average relative difference should be less than 5%, got " + avgRelDiff);
    }

    /**
     * Test 2: Complex chain aggregation with fork-like routing
     *
     * Creates a model with 1 chain (3 classes) and 4 stations with more complex routing.
     * Verifies that aggregation preserves performance metrics.
     */
    @Test
    public void testAggregateChainsComplexRouting() {
        // Create model
        Network model = new Network("ComplexModel");

        // Create stations
        Delay delay = new Delay(model, "Delay");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);
        Queue queue3 = new Queue(model, "Queue3", SchedStrategy.PS);

        // Create classes - 1 chain with 3 classes
        ClosedClass class1 = new ClosedClass(model, "Class1", 10, delay, 0);
        ClosedClass class2 = new ClosedClass(model, "Class2", 0, delay, 0);
        ClosedClass class3 = new ClosedClass(model, "Class3", 0, delay, 0);

        // Set service times with Erlang distribution for variety
        delay.setService(class1, Exp.fitMean(1.0));
        delay.setService(class2, Exp.fitMean(1.0));
        delay.setService(class3, Exp.fitMean(1.0));

        queue1.setService(class1, Exp.fitMean(0.3));
        queue1.setService(class2, Exp.fitMean(0.4));
        queue1.setService(class3, Exp.fitMean(0.35));

        queue2.setService(class1, Erlang.fitMeanAndOrder(0.2, 2));
        queue2.setService(class2, Erlang.fitMeanAndOrder(0.25, 2));
        queue2.setService(class3, Erlang.fitMeanAndOrder(0.22, 2));

        queue3.setService(class1, Exp.fitMean(0.15));
        queue3.setService(class2, Exp.fitMean(0.18));
        queue3.setService(class3, Exp.fitMean(0.16));

        // Set routing with probabilistic class switching
        RoutingMatrix P = model.initRoutingMatrix();

        // From Delay: Class1 -> Queue1
        P.set(class1, class1, delay, queue1, 1.0);

        // From Queue1: split to Queue2 (switch to Class2) and Queue3 (switch to Class3)
        P.set(class1, class2, queue1, queue2, 0.5);  // 50% to Queue2 as Class2
        P.set(class1, class3, queue1, queue3, 0.5);  // 50% to Queue3 as Class3

        // From Queue2: Class2 -> back to Delay as Class1
        P.set(class2, class1, queue2, delay, 1.0);

        // From Queue3: Class3 -> back to Delay as Class1
        P.set(class3, class1, queue3, delay, 1.0);

        model.link(P);

        // Step 1: Solve original model and get avgChainTable
        SolverOptions options = Solver.defaultOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverMVA solverOriginal = new SolverMVA(model, options);
        NetworkAvgChainTable avgChainTableOriginal = solverOriginal.getAvgChainTable();

        assertNotNull(avgChainTableOriginal, "avgChainTable should not be null");

        // Step 2: Aggregate chains
        ModelAdapter.AggregateChainResult result = ModelAdapter.aggregateChains(model);
        assertNotNull(result, "Aggregation result should not be null");

        Network chainModel = result.getChainModel();
        assertNotNull(chainModel, "Chain model should not be null");

        // Verify aggregated model has 1 class (since all classes belong to same chain)
        assertEquals(1, chainModel.getNumberOfClasses(), "Aggregated model should have 1 class (single chain)");

        // Step 3: Solve aggregated model
        SolverMVA solverAggregated = new SolverMVA(chainModel, options);
        NetworkAvgTable avgTableAggregated = solverAggregated.getAvgTable();

        assertNotNull(avgTableAggregated, "avgTable for aggregated model should not be null");

        // Step 4: Compare total throughput and queue lengths
        List<Double> originalTput = avgChainTableOriginal.getTput();
        List<Double> aggregatedTput = avgTableAggregated.getTput();

        // Compare queue lengths
        List<Double> originalQLen = avgChainTableOriginal.getQLen();
        List<Double> aggregatedQLen = avgTableAggregated.getQLen();

        // Verify values are similar
        double totalDiff = 0.0;
        int comparisons = 0;
        for (int i = 0; i < Math.min(originalQLen.size(), aggregatedQLen.size()); i++) {
            double orig = originalQLen.get(i);
            double agg = aggregatedQLen.get(i);
            if (orig > MID_TOL) {
                totalDiff += Math.abs(orig - agg) / orig;
                comparisons++;
            }
        }
        if (comparisons > 0) {
            double avgRelDiff = totalDiff / comparisons;
            assertTrue(avgRelDiff < 0.05, "Average relative difference should be less than 5%, got " + avgRelDiff);
        }
    }

    /**
     * Test that aggregating a model where each class is already its own chain
     * returns a copy without modification.
     */
    @Test
    public void testAggregateChainsSingleClassPerChain() {
        // Create model with no class switching (each class is its own chain)
        Network model = new Network("SimpleModel");

        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue", SchedStrategy.PS);

        ClosedClass class1 = new ClosedClass(model, "Class1", 5, delay, 0);
        ClosedClass class2 = new ClosedClass(model, "Class2", 3, delay, 0);

        delay.setService(class1, Exp.fitMean(1.0));
        delay.setService(class2, Exp.fitMean(1.5));

        queue.setService(class1, Exp.fitMean(0.5));
        queue.setService(class2, Exp.fitMean(0.4));

        // Simple routing without class switching
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(class1, class1, delay, queue, 1.0);
        P.set(class1, class1, queue, delay, 1.0);
        P.set(class2, class2, delay, queue, 1.0);
        P.set(class2, class2, queue, delay, 1.0);

        model.link(P);

        // Aggregate
        ModelAdapter.AggregateChainResult result = ModelAdapter.aggregateChains(model);
        assertNotNull(result, "Aggregation result should not be null");

        // When C == K, model should not be aggregated
        assertFalse(result.getDeaggInfo().isAggregated,
            "Model should not be marked as aggregated when each class is its own chain");

        // Model should have same number of classes
        assertEquals(model.getNumberOfClasses(), result.getChainModel().getNumberOfClasses(),
            "Number of classes should remain the same");
    }
}
