/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.bench.fj;

import jline.util.matrix.Matrix;
import jline.bench.BenchmarkSolvers;
import jline.bench.BenchmarkUtils;
import jline.lang.OpenClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.processes.Exp;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Fork;
import jline.lang.nodes.Join;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Source;
import jline.lang.nodes.Sink;
import jline.solvers.SolverResult;
import jline.solvers.Solver;
import jline.solvers.NetworkSolver;
import jline.solvers.jmt.SolverJMT;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Benchmark for Fork-Join Open Networks
 * Tests both homogeneous and heterogeneous configurations
 */
public class BenchFJ_Open {
    
    public static class BenchmarkResult {
        public String variant;
        public int iteration;
        public Map<String, Double> errQ = new HashMap<>();
        public Map<String, Double> errU = new HashMap<>();
        public Map<String, Double> errR = new HashMap<>();
        public Map<String, Double> errT = new HashMap<>();
    }
    
    /**
     * Run homogeneous FCFS benchmarks
     */
    public static void runHomogeneousFCFS() {
        System.out.println("\n=== Running FJ Open Homogeneous FCFS Benchmarks ===");
        
        // Run variants 1 and 2
        runBenchmark("homogeneous_FCFS", 1, SchedStrategy.FCFS, true, 2);
        runBenchmark("homogeneous_FCFS", 2, SchedStrategy.FCFS, true, 2);
    }
    
    /**
     * Run homogeneous PS benchmarks
     */
    public static void runHomogeneousPS() {
        System.out.println("\n=== Running FJ Open Homogeneous PS Benchmarks ===");
        
        // Run variants 1-4
        runBenchmark("homogeneous_PS", 1, SchedStrategy.PS, true, 2);
        runBenchmark("homogeneous_PS", 2, SchedStrategy.PS, true, 3);
        runBenchmark("homogeneous_PS", 3, SchedStrategy.PS, true, 2);
        runBenchmark("homogeneous_PS", 4, SchedStrategy.PS, true, 3);
    }
    
    /**
     * Run heterogeneous FCFS benchmarks
     */
    public static void runHeterogeneousFCFS() {
        System.out.println("\n=== Running FJ Open Heterogeneous FCFS Benchmarks ===");
        
        // Run variants 1-3
        runBenchmark("heterogeneous_FCFS", 1, SchedStrategy.FCFS, false, 2);
        runBenchmark("heterogeneous_FCFS", 2, SchedStrategy.FCFS, false, 3);
        runBenchmark("heterogeneous_FCFS", 3, SchedStrategy.FCFS, false, 2);
    }
    
    /**
     * Run heterogeneous PS benchmarks
     */
    public static void runHeterogeneousPS() {
        System.out.println("\n=== Running FJ Open Heterogeneous PS Benchmarks ===");
        
        // Run variants 1-5
        runBenchmark("heterogeneous_PS", 1, SchedStrategy.PS, false, 2);
        runBenchmark("heterogeneous_PS", 2, SchedStrategy.PS, false, 3);
        runBenchmark("heterogeneous_PS", 3, SchedStrategy.PS, false, 2);
        runBenchmark("heterogeneous_PS", 4, SchedStrategy.PS, false, 2);
        runBenchmark("heterogeneous_PS", 5, SchedStrategy.PS, false, 3);
    }
    
    /**
     * Core benchmark implementation
     */
    private static BenchmarkResult runBenchmark(String variant, int iteration, 
                                               SchedStrategy sched, boolean homogeneous, int numQueues) {
        String benchName = "bench_FJ_o_" + variant + "_" + iteration;
        System.out.println("\nRunning " + benchName);
        
        BenchmarkResult result = new BenchmarkResult();
        result.variant = variant;
        result.iteration = iteration;
        
        // Create model
        Network model = new Network(benchName);
        
        // Create stations
        Source source = new Source(model, "source");
        Delay delay = new Delay(model, "delay");
        Fork fork = new Fork(model, "fork");
        Join join = new Join(model, "join", fork);
        Sink sink = new Sink(model, "sink");
        
        // Create queue stations
        List<Queue> queues = new ArrayList<>();
        for (int i = 1; i <= numQueues; i++) {
            queues.add(new Queue(model, "queue" + i, sched));
        }
        
        // Create job classes - open classes with arrival rates
        OpenClass class1 = new OpenClass(model, "class1");
        OpenClass class2 = new OpenClass(model, "class2");
        
        // Set arrival rates
        source.setArrival(class1, new Exp(0.3));
        source.setArrival(class2, new Exp(0.2));
        
        // Set service distributions at delay
        delay.setService(class1, new Exp(0.5));
        delay.setService(class2, new Exp(0.8));
        
        // Set service distributions at queues
        if (homogeneous) {
            // All queues have same rates for each class
            for (Queue queue : queues) {
                queue.setService(class1, new Exp(1.5));
                queue.setService(class2, new Exp(2.0));
            }
        } else {
            // Heterogeneous - different rates at each queue
            double[] ratesClass1 = {1.5, 1.1, 2.5}; // Extend if needed
            double[] ratesClass2 = {2.0, 1.8, 2.2}; // Extend if needed
            
            for (int i = 0; i < queues.size(); i++) {
                queues.get(i).setService(class1, new Exp(ratesClass1[i % ratesClass1.length]));
                queues.get(i).setService(class2, new Exp(ratesClass2[i % ratesClass2.length]));
            }
        }
        
        // Set routing
        RoutingMatrix P = model.initRoutingMatrix();
        
        // Source -> Delay
        P.set(class1, source, delay, 1.0);
        P.set(class2, source, delay, 1.0);
        
        // Delay -> Fork
        P.set(class1, delay, fork, 1.0);
        P.set(class2, delay, fork, 1.0);
        
        // Fork -> Queues
        for (Queue queue : queues) {
            P.set(class1, fork, queue, 1.0);
            P.set(class2, fork, queue, 1.0);
        }
        
        // Queues -> Join
        for (Queue queue : queues) {
            P.set(class1, queue, join, 1.0);
            P.set(class2, queue, join, 1.0);
        }
        
        // Join -> Sink
        P.set(class1, join, sink, 1.0);
        P.set(class2, join, sink, 1.0);
        
        model.link(P);
        
        // Get simulation results for comparison
        SolverJMT simSolver = BenchmarkSolvers.getForkJoinSimulationSolver(model);
        SolverResult simResult = simSolver.getAvg();
        
        // Test fork-join solvers
        List<Solver> solvers = BenchmarkSolvers.getForkJoinBenchmarkSolvers(model);
        
        for (Solver solver : solvers) {
            String solverName = solver.getName();
            if (solver.getOptions() != null && solver.getOptions().config != null) {
                solverName += "_" + solver.getOptions().config.fork_join;
            }
            
            System.out.println("  Testing " + solverName);
            
            try {
                SolverResult solverResult = ((NetworkSolver) solver).getAvg();
                
                // Calculate errors
                double errorQ = BenchmarkUtils.maxErrorOnSum(solverResult.QN, simResult.QN);
                double errorU = BenchmarkUtils.utilizationError(solverResult.UN, simResult.UN, model);
                double errorR = BenchmarkUtils.mape(solverResult.RN, simResult.RN);
                double errorT = BenchmarkUtils.mape(solverResult.TN, simResult.TN);
                
                result.errQ.put(solverName, errorQ);
                result.errU.put(solverName, errorU);
                result.errR.put(solverName, errorR);
                result.errT.put(solverName, errorT);
                
                System.out.printf("    Errors - Q: %.6f, U: %.6f, R: %.6f, T: %.6f%n", 
                                errorQ, errorU, errorR, errorT);
                
            } catch (Exception e) {
                System.out.println("    ERROR: " + e.getMessage());
                result.errQ.put(solverName, Double.MAX_VALUE);
                result.errU.put(solverName, Double.MAX_VALUE);
                result.errR.put(solverName, Double.MAX_VALUE);
                result.errT.put(solverName, Double.MAX_VALUE);
            }
        }
        
        return result;
    }
    
    /**
     * Run all open FJ benchmarks
     */
    public static void runAll() {
        runHomogeneousFCFS();
        runHomogeneousPS();
        runHeterogeneousFCFS();
        runHeterogeneousPS();
    }
    
    public static void main(String[] args) {
        jline.GlobalConstants.Verbose = jline.VerboseLevel.SILENT;
        runAll();
    }
}