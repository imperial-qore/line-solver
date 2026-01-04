/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.bench.oqn;

import jline.util.matrix.Matrix;
import jline.bench.BenchmarkSolvers;
import jline.bench.BenchmarkUtils;
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
 * Template for OQN (Open Queueing Network) benchmark implementations.
 * Aligned with MATLAB benchmarks in line-test.git/bench/bench_OQN_*.
 *
 * Model structure: Source -> Delay -> Queue1 -> Queue2 -> Sink
 */
public class BenchOQNTemplate {

    /**
     * Run a single-class OQN benchmark
     *
     * @param benchmarkName Name of the benchmark
     * @param iteration     Iteration number for reproducibility
     * @param arrivalRate   Arrival rate (lambda) - controls load level
     * @param sched         Scheduling strategy (PS or FCFS)
     * @return Map containing error metrics
     */
    public static Map<String, Object> runBenchmark(String benchmarkName, int iteration,
                                                   double arrivalRate, SchedStrategy sched) {
        // Create model: Source -> Delay -> Queue1 -> Queue2 -> Sink
        Network model = new Network(benchmarkName + "_" + iteration);

        // Create nodes
        Source source = new Source(model, "Source");
        Delay delay = new Delay(model, "Delay");
        Queue queue1 = new Queue(model, "Queue1", sched);
        Queue queue2 = new Queue(model, "Queue2", sched);
        Sink sink = new Sink(model, "Sink");

        // Create open class
        OpenClass openClass = new OpenClass(model, "Class1", 0);

        // Set arrival process at source
        source.setArrival(openClass, new Exp(arrivalRate));

        // Set service times (aligned with MATLAB benchmark)
        delay.setService(openClass, new Exp(2.0));
        queue1.setService(openClass, new Exp(1.0));
        queue2.setService(openClass, new Exp(1.0));

        // Set routing: Source -> Delay -> Queue1 -> Queue2 -> Sink
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(openClass, source, delay, 1.0);
        P.set(openClass, delay, queue1, 1.0);
        P.set(openClass, queue1, queue2, 1.0);
        P.set(openClass, queue2, sink, 1.0);
        model.link(P);

        return runSolversAndCollectErrors(benchmarkName, model);
    }

    /**
     * Run a multiclass OQN benchmark
     *
     * @param benchmarkName Name of the benchmark
     * @param iteration     Iteration number
     * @param arrivalRate1  Arrival rate for class 1
     * @param arrivalRate2  Arrival rate for class 2
     * @param sched         Scheduling strategy
     * @return Map containing error metrics
     */
    public static Map<String, Object> runMulticlassBenchmark(String benchmarkName, int iteration,
                                                              double arrivalRate1, double arrivalRate2,
                                                              SchedStrategy sched) {
        Network model = new Network(benchmarkName + "_" + iteration);

        // Create nodes
        Source source = new Source(model, "Source");
        Delay delay = new Delay(model, "Delay");
        Queue queue1 = new Queue(model, "Queue1", sched);
        Queue queue2 = new Queue(model, "Queue2", sched);
        Sink sink = new Sink(model, "Sink");

        // Create two open classes
        OpenClass class1 = new OpenClass(model, "Class1", 0);
        OpenClass class2 = new OpenClass(model, "Class2", 0);

        // Set arrivals
        source.setArrival(class1, new Exp(arrivalRate1));
        source.setArrival(class2, new Exp(arrivalRate2));

        // Set service times
        delay.setService(class1, new Exp(2.0));
        delay.setService(class2, new Exp(2.5));
        queue1.setService(class1, new Exp(1.0));
        queue1.setService(class2, new Exp(1.2));
        queue2.setService(class1, new Exp(1.0));
        queue2.setService(class2, new Exp(1.0));

        // Set routing
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(class1, source, delay, 1.0);
        P.set(class1, delay, queue1, 1.0);
        P.set(class1, queue1, queue2, 1.0);
        P.set(class1, queue2, sink, 1.0);

        P.set(class2, source, delay, 1.0);
        P.set(class2, delay, queue1, 1.0);
        P.set(class2, queue1, queue2, 1.0);
        P.set(class2, queue2, sink, 1.0);
        model.link(P);

        return runSolversAndCollectErrors(benchmarkName, model);
    }

    /**
     * Run a tandem OQN benchmark (longer chain)
     *
     * @param benchmarkName Name of the benchmark
     * @param iteration     Iteration number
     * @param arrivalRate   Arrival rate
     * @param sched         Scheduling strategy
     * @return Map containing error metrics
     */
    public static Map<String, Object> runTandemBenchmark(String benchmarkName, int iteration,
                                                          double arrivalRate, SchedStrategy sched) {
        Network model = new Network(benchmarkName + "_" + iteration);

        // Create nodes: Source -> Delay -> Queue1 -> Queue2 -> Queue3 -> Sink
        Source source = new Source(model, "Source");
        Delay delay = new Delay(model, "Delay");
        Queue queue1 = new Queue(model, "Queue1", sched);
        Queue queue2 = new Queue(model, "Queue2", sched);
        Queue queue3 = new Queue(model, "Queue3", sched);
        Sink sink = new Sink(model, "Sink");

        OpenClass openClass = new OpenClass(model, "Class1", 0);

        source.setArrival(openClass, new Exp(arrivalRate));
        delay.setService(openClass, new Exp(2.0));
        queue1.setService(openClass, new Exp(1.5));
        queue2.setService(openClass, new Exp(1.0));
        queue3.setService(openClass, new Exp(1.2));

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(openClass, source, delay, 1.0);
        P.set(openClass, delay, queue1, 1.0);
        P.set(openClass, queue1, queue2, 1.0);
        P.set(openClass, queue2, queue3, 1.0);
        P.set(openClass, queue3, sink, 1.0);
        model.link(P);

        return runSolversAndCollectErrors(benchmarkName, model);
    }

    /**
     * Run all solvers and collect error metrics
     */
    private static Map<String, Object> runSolversAndCollectErrors(String benchmarkName, Network model) {
        // Get simulation results for comparison
        SolverJMT simSolver = BenchmarkSolvers.getSimulationSolver(model);
        SolverResult simResult = simSolver.getAvg();

        // Test all available solvers
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

        OQNResultsFormatter.addGlobalResult(benchmarkName, results);

        return results;
    }
}
