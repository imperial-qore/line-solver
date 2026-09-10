/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.fes;

import jline.VerboseLevel;
import jline.lang.*;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.*;
import jline.lang.processes.Exp;
import jline.solvers.NetworkAvgTable;
import jline.solvers.Solver;
import jline.solvers.SolverOptions;
import jline.solvers.mva.SolverMVA;
import jline.util.Maths;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;
import static jline.TestTools.COARSE_TOL;

/**
 * Tests for FES (Flow-Equivalent Server) aggregation via ModelAdapter.aggregateFES().
 *
 * Mirrors the MATLAB tests test_fes_single_class.m and test_fes_multi_class.m.
 */
public class FESAggregatorTest {

    @BeforeEach
    public void matlabRandomSeedSetUp() {
        Maths.setRandomNumbersMatlab(true);
    }

    @AfterEach
    public void matlabRandomSeedClear() {
        Maths.setRandomNumbersMatlab(false);
    }

    /**
     * Test 1: Single-class FES aggregation on a 4-station tandem network.
     *
     * Delay -> Q1 -> Q2 -> Q3 -> Delay (closed, N=5)
     * Service times: Delay=5.0, Q1=1.5, Q2=1.0, Q3=0.8
     * Aggregate Q1+Q2 into FES, solve, compare throughput at Delay.
     * Single-class FES should be exact.
     */
    @Test
    public void testFESSingleClass() {
        // Build original 4-station tandem
        Network model = new Network("FES_SingleClass");

        Delay delay = new Delay(model, "Delay");
        Queue q1 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue q2 = new Queue(model, "Queue2", SchedStrategy.PS);
        Queue q3 = new Queue(model, "Queue3", SchedStrategy.PS);

        ClosedClass class1 = new ClosedClass(model, "Class1", 5, delay, 0);

        delay.setService(class1, Exp.fitMean(5.0));
        q1.setService(class1, Exp.fitMean(1.5));
        q2.setService(class1, Exp.fitMean(1.0));
        q3.setService(class1, Exp.fitMean(0.8));

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(class1, class1, delay, q1, 1.0);
        P.set(class1, class1, q1, q2, 1.0);
        P.set(class1, class1, q2, q3, 1.0);
        P.set(class1, class1, q3, delay, 1.0);
        model.link(P);

        // Solve original model
        SolverOptions options = Solver.defaultOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverMVA solverOrig = new SolverMVA(model, options);
        NetworkAvgTable origTable = solverOrig.getAvgTable();
        assertNotNull(origTable, "Original avgTable should not be null");

        // Get throughput at Delay (station index 0, class index 0)
        List<Double> origTput = origTable.getTput();
        double origDelayTput = origTput.get(0); // Delay, Class1

        // Aggregate Q1+Q2 into FES
        List<Station> stationSubset = Arrays.asList(q1, q2);
        FESResult fesResult = ModelAdapter.aggregateFES(model, stationSubset);
        assertNotNull(fesResult, "FES result should not be null");

        Network fesModel = fesResult.getFesModel();
        assertNotNull(fesModel, "FES model should not be null");

        // Solve FES model
        SolverMVA solverFES = new SolverMVA(fesModel, options);
        NetworkAvgTable fesTable = solverFES.getAvgTable();
        assertNotNull(fesTable, "FES avgTable should not be null");

        // Find Delay throughput in FES model (Delay is first complement station)
        List<Double> fesTput = fesTable.getTput();
        double fesDelayTput = fesTput.get(0); // Delay, Class1

        // Single-class FES should be exact: assert within 1%
        double relError = Math.abs(origDelayTput - fesDelayTput) / origDelayTput;
        assertTrue(relError < COARSE_TOL,
            "Single-class FES throughput should match within 1%: original=" + origDelayTput
                + ", FES=" + fesDelayTput + ", relError=" + relError);

        // Verify deaggInfo fields
        FESDeaggInfo deaggInfo = fesResult.getDeaggInfo();
        assertNotNull(deaggInfo, "DeaggInfo should not be null");
        assertEquals(2, deaggInfo.subsetIndices.length,
            "subsetIndices should have 2 entries (Q1, Q2)");
        assertEquals(2, deaggInfo.complementIndices.length,
            "complementIndices should have 2 entries (Delay, Q3)");
        assertTrue(deaggInfo.fesNodeIdx >= 0,
            "fesNodeIdx should be non-negative");
    }

    /**
     * Test 2: Multi-class FES aggregation on a 4-station tandem network.
     *
     * Delay -> Q1 -> Q2 -> Q3 -> Delay (closed, N1=3, N2=2)
     * Service times per class: Delay(5.0,4.0), Q1(1.5,2.0), Q2(1.0,1.2), Q3(0.8,1.0)
     * Aggregate Q1+Q2 into FES; verify against the exact NC solver.
     * (Approximate MVA is inaccurate for the state-dependent FES rate.)
     */
    @Test
    public void testFESMultiClass() {
        // Build original 4-station tandem
        Network model = new Network("FES_MultiClass");

        Delay delay = new Delay(model, "Delay");
        Queue q1 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue q2 = new Queue(model, "Queue2", SchedStrategy.PS);
        Queue q3 = new Queue(model, "Queue3", SchedStrategy.PS);

        ClosedClass class1 = new ClosedClass(model, "Class1", 3, delay, 0);
        ClosedClass class2 = new ClosedClass(model, "Class2", 2, delay, 0);

        delay.setService(class1, Exp.fitMean(5.0));
        delay.setService(class2, Exp.fitMean(4.0));

        q1.setService(class1, Exp.fitMean(1.5));
        q1.setService(class2, Exp.fitMean(2.0));

        q2.setService(class1, Exp.fitMean(1.0));
        q2.setService(class2, Exp.fitMean(1.2));

        q3.setService(class1, Exp.fitMean(0.8));
        q3.setService(class2, Exp.fitMean(1.0));

        RoutingMatrix P = model.initRoutingMatrix();
        // Class1 tandem routing
        P.set(class1, class1, delay, q1, 1.0);
        P.set(class1, class1, q1, q2, 1.0);
        P.set(class1, class1, q2, q3, 1.0);
        P.set(class1, class1, q3, delay, 1.0);
        // Class2 tandem routing
        P.set(class2, class2, delay, q1, 1.0);
        P.set(class2, class2, q1, q2, 1.0);
        P.set(class2, class2, q2, q3, 1.0);
        P.set(class2, class2, q3, delay, 1.0);
        model.link(P);

        // Solve original model
        SolverOptions options = Solver.defaultOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverMVA solverOrig = new SolverMVA(model, options);
        NetworkAvgTable origTable = solverOrig.getAvgTable();
        assertNotNull(origTable, "Original avgTable should not be null");

        List<Double> origTput = origTable.getTput();
        // Delay throughputs: index 0 = (Delay, Class1), index 1 = (Delay, Class2)
        double origTput1 = origTput.get(0);
        double origTput2 = origTput.get(1);

        // Aggregate Q1+Q2 into FES
        List<Station> stationSubset = Arrays.asList(q1, q2);
        FESResult fesResult = ModelAdapter.aggregateFES(model, stationSubset);
        assertNotNull(fesResult, "FES result should not be null");

        Network fesModel = fesResult.getFesModel();
        assertNotNull(fesModel, "FES model should not be null");

        // Solve FES model with NC (exact). The FES composite rate is state
        // dependent, so approximate MVA (linearizer) is not accurate here (~40%
        // off, in both LINE-Java and MATLAB); the exact normalizing-constant
        // solver reproduces the original throughput to machine precision.
        SolverOptions ncOpts = Solver.defaultOptions();
        ncOpts.verbose = VerboseLevel.SILENT;
        ncOpts.method = "exact";
        NetworkAvgTable fesTable = new jline.solvers.nc.SolverNC(fesModel, ncOpts).getAvgTable();
        assertNotNull(fesTable, "FES avgTable should not be null");

        List<Double> fesTput = fesTable.getTput();
        double fesTput1 = fesTput.get(0);
        double fesTput2 = fesTput.get(1);

        double relError1 = Math.abs(origTput1 - fesTput1) / origTput1;
        double relError2 = Math.abs(origTput2 - fesTput2) / origTput2;

        assertTrue(relError1 < COARSE_TOL,
            "Multi-class FES throughput (Class1) should match (NC exact): original=" + origTput1
                + ", FES=" + fesTput1 + ", relError=" + relError1);
        assertTrue(relError2 < COARSE_TOL,
            "Multi-class FES throughput (Class2) should match (NC exact): original=" + origTput2
                + ", FES=" + fesTput2 + ", relError=" + relError2);

        // Verify deaggInfo fields
        FESDeaggInfo deaggInfo = fesResult.getDeaggInfo();
        assertNotNull(deaggInfo, "DeaggInfo should not be null");
        assertEquals(2, deaggInfo.subsetIndices.length,
            "subsetIndices should have 2 entries (Q1, Q2)");
        assertEquals(2, deaggInfo.complementIndices.length,
            "complementIndices should have 2 entries (Delay, Q3)");
        assertTrue(deaggInfo.fesNodeIdx >= 0,
            "fesNodeIdx should be non-negative");
    }

    /**
     * Test 3: Validation tests for FES aggregation error handling.
     *
     * - Aggregating all stations should throw IllegalArgumentException.
     * - Open class model should throw IllegalArgumentException.
     */
    @Test
    public void testFESValidation() {
        // --- Test: aggregating all stations should fail ---
        Network model = new Network("FES_AllStations");

        Delay delay = new Delay(model, "Delay");
        Queue q1 = new Queue(model, "Queue1", SchedStrategy.PS);

        ClosedClass class1 = new ClosedClass(model, "Class1", 5, delay, 0);

        delay.setService(class1, Exp.fitMean(1.0));
        q1.setService(class1, Exp.fitMean(0.5));

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(class1, class1, delay, q1, 1.0);
        P.set(class1, class1, q1, delay, 1.0);
        model.link(P);

        List<Station> allStations = Arrays.asList(delay, q1);
        assertThrows(IllegalArgumentException.class,
            () -> ModelAdapter.aggregateFES(model, allStations),
            "Aggregating all stations should throw IllegalArgumentException");

        // --- Test: open class model should fail ---
        Network openModel = new Network("FES_Open");

        Source source = new Source(openModel, "Source");
        Queue queue = new Queue(openModel, "Queue", SchedStrategy.PS);
        Sink sink = new Sink(openModel, "Sink");

        OpenClass openClass = new OpenClass(openModel, "OpenClass");

        source.setArrival(openClass, Exp.fitMean(1.0));
        queue.setService(openClass, Exp.fitMean(0.5));

        RoutingMatrix Popen = openModel.initRoutingMatrix();
        Popen.set(openClass, openClass, source, queue, 1.0);
        Popen.set(openClass, openClass, queue, sink, 1.0);
        openModel.link(Popen);

        List<Station> openSubset = Arrays.asList(queue);
        assertThrows(IllegalArgumentException.class,
            () -> ModelAdapter.aggregateFES(openModel, openSubset),
            "FES aggregation on open model should throw IllegalArgumentException");
    }
}
