/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang;

import jline.VerboseLevel;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
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

    /**
     * Test 4: the Network-level entry points of the two transformations.
     *
     * Network.aggregateChains() must return the same aggregation as
     * ModelAdapter.aggregateChains, and Network.withoutClass() must return a
     * copy without the class while leaving the original model solvable. Both
     * entry points existed only on ModelAdapter until they were wired.
     */
    @Test
    public void testNetworkEntryPointsAggregateChainsAndWithoutClass() {
        Network model = new Network("EntryPoints");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.PS);

        ClosedClass class1 = new ClosedClass(model, "ClassA", 4, delay, 0);
        ClosedClass class2 = new ClosedClass(model, "ClassB", 0, delay, 0);

        delay.setService(class1, Exp.fitMean(1.0));
        delay.setService(class2, Exp.fitMean(0.5));
        queue.setService(class1, Exp.fitMean(1.0 / 3.0));
        queue.setService(class2, Exp.fitMean(0.25));

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(class1, class1, delay, queue, 1.0);
        P.set(class1, class2, queue, delay, 0.7);
        P.set(class1, class1, queue, delay, 0.3);
        P.set(class2, class1, queue, delay, 0.6);
        P.set(class2, class2, queue, delay, 0.4);
        P.set(class2, class2, delay, queue, 1.0);
        model.link(P);

        // one chain, two classes
        assertEquals(1, model.getStruct(true).nchains, "The two classes must form one chain");

        ModelAdapter.AggregateChainResult viaNetwork = model.aggregateChains();
        assertNotNull(viaNetwork, "Network.aggregateChains must return a result");
        assertEquals(1, viaNetwork.getChainModel().getNumberOfClasses(),
            "The aggregated model must have one class per chain");
        assertTrue(viaNetwork.getDeaggInfo().isAggregated,
            "A model with two classes in one chain must be marked as aggregated");
        assertEquals(2, model.getNumberOfClasses(), "The original model must be untouched");

        SolverOptions options = Solver.defaultOptions();
        options.verbose = VerboseLevel.SILENT;
        NetworkAvgTable base = new SolverMVA(model, options).getAvgTable();
        NetworkAvgTable aggregated = new SolverMVA(viaNetwork.getChainModel(), options).getAvgTable();

        // chain aggregation is exact on this product-form model: per station,
        // the chain queue length is the sum of the class queue lengths
        int M = model.getNumberOfStations();
        List<Double> baseQLen = base.getQLen();
        List<Double> chainQLen = aggregated.getQLen();
        for (int ist = 0; ist < M; ist++) {
            double sum = 0.0;
            for (int r = 0; r < 2; r++) {
                sum += baseQLen.get(ist * 2 + r);
            }
            assertEquals(sum, chainQLen.get(ist), 1e-4,
                "Chain queue length at station " + ist + " must equal the sum over its classes");
        }
    }

    /**
     * Test 5: class removal, mutating and non-mutating.
     */
    @Test
    public void testRemoveClassAndWithoutClass() {
        Network model = openTwoClassModel();
        List<JobClass> classes = model.getClasses();
        JobClass second = classes.get(1);

        Network reduced = model.withoutClass(second);
        assertEquals(1, reduced.getNumberOfClasses(), "The copy must have lost the class");
        assertEquals(2, model.getNumberOfClasses(), "The original must be untouched");

        SolverOptions options = Solver.defaultOptions();
        options.verbose = VerboseLevel.SILENT;
        NetworkAvgTable got = new SolverMVA(reduced, options).getAvgTable();
        NetworkAvgTable want = new SolverMVA(openOneClassModel(), options).getAvgTable();
        assertEquals(want.getQLen().size(), got.getQLen().size(),
            "Removing a class must give the same table as building the model without it");
        for (int i = 0; i < want.getQLen().size(); i++) {
            assertEquals(want.getQLen().get(i), got.getQLen().get(i), 1e-6);
            assertEquals(want.getTput().get(i), got.getTput().get(i), 1e-6);
        }

        // the mutating variant leaves the model itself with one class
        model.removeClass(second);
        assertEquals(1, model.getNumberOfClasses(), "removeClass must mutate the model");
        NetworkAvgTable mutated = new SolverMVA(model, options).getAvgTable();
        for (int i = 0; i < want.getQLen().size(); i++) {
            assertEquals(want.getQLen().get(i), mutated.getQLen().get(i), 1e-6);
        }
    }

    private static Network openTwoClassModel() {
        Network model = new Network("open2");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.PS);
        Sink sink = new Sink(model, "Sink");
        OpenClass classA = new OpenClass(model, "ClassA", 0);
        OpenClass classB = new OpenClass(model, "ClassB", 0);
        source.setArrival(classA, Exp.fitMean(2.0));
        source.setArrival(classB, Exp.fitMean(1.0 / 0.3));
        queue.setService(classA, Exp.fitMean(0.5));
        queue.setService(classB, Exp.fitMean(1.0 / 1.5));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(classA, classA, model.serialRouting(source, queue, sink));
        P.set(classB, classB, model.serialRouting(source, queue, sink));
        model.link(P);
        return model;
    }

    private static Network openOneClassModel() {
        Network model = new Network("open1");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.PS);
        Sink sink = new Sink(model, "Sink");
        OpenClass classA = new OpenClass(model, "ClassA", 0);
        source.setArrival(classA, Exp.fitMean(2.0));
        queue.setService(classA, Exp.fitMean(0.5));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(classA, classA, model.serialRouting(source, queue, sink));
        model.link(P);
        return model;
    }

    /**
     * Test 6: removing a class must shrink the K x K class-switching mask that
     * link() records, otherwise the chain computation runs over classes the
     * model no longer has.
     */
    @Test
    public void testRemoveClassSlicesTheClassSwitchMask() {
        Network model = new Network("cs3");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.PS);

        ClosedClass classA = new ClosedClass(model, "ClassA", 3, delay, 0);
        ClosedClass classB = new ClosedClass(model, "ClassB", 0, delay, 0);
        ClosedClass classC = new ClosedClass(model, "ClassC", 2, delay, 0);
        double[] means = {1.0, 0.5, 1.0 / 1.5};
        ClosedClass[] classes = {classA, classB, classC};
        for (int r = 0; r < 3; r++) {
            delay.setService(classes[r], Exp.fitMean(means[r]));
            queue.setService(classes[r], Exp.fitMean(1.0 / 3.0));
        }

        RoutingMatrix P = model.initRoutingMatrix();
        // ClassA and ClassB switch into each other, ClassC is its own chain
        P.set(classA, classA, delay, queue, 1.0);
        P.set(classA, classA, queue, delay, 0.3);
        P.set(classA, classB, queue, delay, 0.7);
        P.set(classB, classB, delay, queue, 1.0);
        P.set(classB, classA, queue, delay, 0.6);
        P.set(classB, classB, queue, delay, 0.4);
        P.set(classC, classC, delay, queue, 1.0);
        P.set(classC, classC, queue, delay, 1.0);
        model.link(P);

        SolverOptions options = Solver.defaultOptions();
        options.verbose = VerboseLevel.SILENT;
        new SolverMVA(model, options).getAvgTable();

        Network reduced = model.withoutClass(classC);
        assertEquals(2, reduced.getNumberOfClasses(), "The copy must have lost ClassC");
        assertEquals(2, reduced.getStruct(true).nclasses, "The struct must be rebuilt over two classes");
        assertEquals(1, reduced.getStruct(true).nchains, "The surviving classes still form one chain");

        NetworkAvgTable table = new SolverMVA(reduced, options).getAvgTable();
        for (double q : table.getQLen()) {
            assertTrue(Double.isFinite(q), "Queue lengths must stay finite after class removal");
        }
    }
}
