/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.bench.mqn;

import jline.util.matrix.Matrix;
import jline.bench.BenchmarkSolvers;
import jline.bench.BenchmarkUtils;
import jline.lang.ClosedClass;
import jline.lang.OpenClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.processes.Exp;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.solvers.SolverResult;
import jline.solvers.Solver;
import jline.solvers.NetworkSolver;
import jline.solvers.jmt.SolverJMT;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Template for MQN (Mixed Queueing Network) benchmark implementations.
 * Aligned with MATLAB benchmarks in line-test.git/bench/bench_MQN_*.
 *
 * Model structure:
 * - Open class: Source -> Delay -> Queue1 -> Queue2 -> Sink
 * - Closed class: Delay -> Queue1 -> Queue2 -> Delay (circular)
 */
public class BenchMQNTemplate {

    /**
     * Run a standard MQN benchmark (1 open + 1 closed class)
     *
     * @param benchmarkName  Name of the benchmark
     * @param iteration      Iteration number
     * @param openArrival    Arrival rate for open class
     * @param closedPop      Population for closed class
     * @param sched          Scheduling strategy
     * @return Map containing error metrics
     */
    public static Map<String, Object> runBenchmark(String benchmarkName, int iteration,
                                                   double openArrival, int closedPop,
                                                   SchedStrategy sched) {
        Network model = new Network(benchmarkName + "_" + iteration);

        // Create nodes
        Source source = new Source(model, "Source");
        Delay delay = new Delay(model, "Delay");
        Queue queue1 = new Queue(model, "Queue1", sched);
        Queue queue2 = new Queue(model, "Queue2", sched);
        Sink sink = new Sink(model, "Sink");

        // Create classes
        OpenClass openClass = new OpenClass(model, "OpenClass", 0);
        ClosedClass closedClass = new ClosedClass(model, "ClosedClass", closedPop, delay);

        // Set arrival for open class
        source.setArrival(openClass, new Exp(openArrival));

        // Set service times (aligned with MATLAB benchmark)
        delay.setService(openClass, new Exp(2.0));
        delay.setService(closedClass, new Exp(1.0));

        queue1.setService(openClass, new Exp(1.0));
        queue1.setService(closedClass, new Exp(1.2));

        queue2.setService(openClass, new Exp(1.0));
        queue2.setService(closedClass, new Exp(1.0));

        // Set routing
        RoutingMatrix P = model.initRoutingMatrix();

        // Open class routing: Source -> Delay -> Queue1 -> Queue2 -> Sink
        P.set(openClass, source, delay, 1.0);
        P.set(openClass, delay, queue1, 1.0);
        P.set(openClass, queue1, queue2, 1.0);
        P.set(openClass, queue2, sink, 1.0);

        // Closed class routing: Delay -> Queue1 -> Queue2 -> Delay
        P.set(closedClass, delay, queue1, 1.0);
        P.set(closedClass, queue1, queue2, 1.0);
        P.set(closedClass, queue2, delay, 1.0);

        model.link(P);

        return runSolversAndCollectErrors(benchmarkName, model);
    }

    /**
     * Run MQN benchmark with higher closed population
     */
    public static Map<String, Object> runHighPopBenchmark(String benchmarkName, int iteration,
                                                          double openArrival, int closedPop,
                                                          SchedStrategy sched) {
        return runBenchmark(benchmarkName, iteration, openArrival, closedPop, sched);
    }

    /**
     * Run a tandem MQN benchmark (longer chain)
     */
    public static Map<String, Object> runTandemBenchmark(String benchmarkName, int iteration,
                                                          double openArrival, int closedPop,
                                                          SchedStrategy sched) {
        Network model = new Network(benchmarkName + "_" + iteration);

        // Create nodes with longer chain
        Source source = new Source(model, "Source");
        Delay delay = new Delay(model, "Delay");
        Queue queue1 = new Queue(model, "Queue1", sched);
        Queue queue2 = new Queue(model, "Queue2", sched);
        Queue queue3 = new Queue(model, "Queue3", sched);
        Sink sink = new Sink(model, "Sink");

        OpenClass openClass = new OpenClass(model, "OpenClass", 0);
        ClosedClass closedClass = new ClosedClass(model, "ClosedClass", closedPop, delay);

        source.setArrival(openClass, new Exp(openArrival));

        delay.setService(openClass, new Exp(2.0));
        delay.setService(closedClass, new Exp(1.0));

        queue1.setService(openClass, new Exp(1.5));
        queue1.setService(closedClass, new Exp(1.2));

        queue2.setService(openClass, new Exp(1.0));
        queue2.setService(closedClass, new Exp(1.0));

        queue3.setService(openClass, new Exp(1.2));
        queue3.setService(closedClass, new Exp(1.1));

        RoutingMatrix P = model.initRoutingMatrix();

        // Open class: Source -> Delay -> Queue1 -> Queue2 -> Queue3 -> Sink
        P.set(openClass, source, delay, 1.0);
        P.set(openClass, delay, queue1, 1.0);
        P.set(openClass, queue1, queue2, 1.0);
        P.set(openClass, queue2, queue3, 1.0);
        P.set(openClass, queue3, sink, 1.0);

        // Closed class: Delay -> Queue1 -> Queue2 -> Queue3 -> Delay
        P.set(closedClass, delay, queue1, 1.0);
        P.set(closedClass, queue1, queue2, 1.0);
        P.set(closedClass, queue2, queue3, 1.0);
        P.set(closedClass, queue3, delay, 1.0);

        model.link(P);

        return runSolversAndCollectErrors(benchmarkName, model);
    }

    /**
     * Run all solvers and collect error metrics
     */
    private static Map<String, Object> runSolversAndCollectErrors(String benchmarkName, Network model) {
        SolverJMT simSolver = BenchmarkSolvers.getSimulationSolver(model);
        SolverResult simResult = simSolver.getAvg();

        List<Solver> solvers = BenchmarkSolvers.getBenchmarkSolvers(model);

        Map<String, Double> errQ = new HashMap<>();
        Map<String, Double> errU = new HashMap<>();
        Map<String, Double> errR = new HashMap<>();
        Map<String, Double> errT = new HashMap<>();

        for (Solver solver : solvers) {
            String solverName = solver.getClass().getSimpleName();

            try {
                SolverResult result = ((NetworkSolver) solver).getAvg();

                double errorQ = BenchmarkUtils.maxErrorOnSum(result.QN, simResult.QN);
                double errorU = BenchmarkUtils.utilizationError(result.UN, simResult.UN, model);
                double errorR = BenchmarkUtils.mape(result.RN, simResult.RN);
                double errorT = BenchmarkUtils.mape(result.TN, simResult.TN);

                errQ.put(solverName, errorQ);
                errU.put(solverName, errorU);
                errR.put(solverName, errorR);
                errT.put(solverName, errorT);

            } catch (Exception e) {
                errQ.put(solverName, Double.MAX_VALUE);
                errU.put(solverName, Double.MAX_VALUE);
                errR.put(solverName, Double.MAX_VALUE);
                errT.put(solverName, Double.MAX_VALUE);
            }
        }

        Map<String, Object> results = new HashMap<>();
        results.put("ERRQ", errQ);
        results.put("ERRU", errU);
        results.put("ERRR", errR);
        results.put("ERRT", errT);
        results.put("model", model);

        MQNResultsFormatter.addGlobalResult(benchmarkName, results);

        return results;
    }
}
