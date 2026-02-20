/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.bench.fj;

import jline.util.matrix.Matrix;
import jline.bench.BenchmarkSolvers;
import jline.bench.BenchmarkUtils;
import jline.lang.ClosedClass;
import jline.lang.OpenClass;
import jline.lang.JobClass;
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

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Benchmark for Fork-Join Nested Networks
 * Tests nested fork-join structures with closed, open, and mixed configurations
 */
public class BenchFJ_Nested {
    
    public static class BenchmarkResult {
        public String networkType;
        public String schedType;
        public int iteration;
        public Map<String, Double> errQ = new HashMap<>();
        public Map<String, Double> errU = new HashMap<>();
        public Map<String, Double> errR = new HashMap<>();
        public Map<String, Double> errT = new HashMap<>();
    }
    
    /**
     * Run nested closed FCFS benchmarks
     */
    public static void runNestedClosedFCFS() {
        System.out.println("\n=== Running FJ Nested Closed FCFS Benchmarks ===");
        for (int i = 1; i <= 4; i++) {
            runBenchmark("c", SchedStrategy.FCFS, i);
        }
    }
    
    /**
     * Run nested closed PS benchmarks
     */
    public static void runNestedClosedPS() {
        System.out.println("\n=== Running FJ Nested Closed PS Benchmarks ===");
        for (int i = 1; i <= 4; i++) {
            runBenchmark("c", SchedStrategy.PS, i);
        }
    }
    
    /**
     * Run nested open FCFS benchmarks
     */
    public static void runNestedOpenFCFS() {
        System.out.println("\n=== Running FJ Nested Open FCFS Benchmarks ===");
        for (int i = 1; i <= 4; i++) {
            runBenchmark("o", SchedStrategy.FCFS, i);
        }
    }
    
    /**
     * Run nested open PS benchmarks
     */
    public static void runNestedOpenPS() {
        System.out.println("\n=== Running FJ Nested Open PS Benchmarks ===");
        for (int i = 1; i <= 4; i++) {
            runBenchmark("o", SchedStrategy.PS, i);
        }
    }
    
    /**
     * Run nested mixed FCFS benchmarks
     */
    public static void runNestedMixedFCFS() {
        System.out.println("\n=== Running FJ Nested Mixed FCFS Benchmarks ===");
        for (int i = 1; i <= 4; i++) {
            runBenchmark("m", SchedStrategy.FCFS, i);
        }
    }
    
    /**
     * Run nested mixed PS benchmarks
     */
    public static void runNestedMixedPS() {
        System.out.println("\n=== Running FJ Nested Mixed PS Benchmarks ===");
        for (int i = 1; i <= 4; i++) {
            runBenchmark("m", SchedStrategy.PS, i);
        }
    }
    
    /**
     * Core benchmark implementation for nested fork-join networks
     */
    private static BenchmarkResult runBenchmark(String networkType, SchedStrategy sched, int iteration) {
        String schedName = sched == SchedStrategy.FCFS ? "FCFS" : "PS";
        String benchName = "bench_FJ_nested_" + networkType + "_" + schedName + "_" + iteration;
        System.out.println("\nRunning " + benchName);
        
        BenchmarkResult result = new BenchmarkResult();
        result.networkType = networkType;
        result.schedType = schedName;
        result.iteration = iteration;
        
        // Create model
        Network model = new Network(benchName);
        
        // Create stations based on network type
        Source source = null;
        Sink sink = null;
        Delay delay = new Delay(model, "delay");
        
        if (!networkType.equals("c")) { // Open or mixed networks need source/sink
            source = new Source(model, "source");
            sink = new Sink(model, "sink");
        }
        
        // Create main fork-join structure
        Fork fork = new Fork(model, "fork");
        Join join = new Join(model, "join", fork);
        
        // Create nested fork-join structure
        Fork fork_nested = new Fork(model, "fork_nested");
        Join join_nested = new Join(model, "join_nested", fork_nested);
        
        // Create queue stations
        Queue queue1 = new Queue(model, "queue1", sched);
        Queue queue2 = new Queue(model, "queue2", sched);
        Queue queue3 = new Queue(model, "queue3", sched);
        
        // Create job classes based on network type
        JobClass class1 = null;
        JobClass class2 = null;
        
        if (networkType.equals("c")) { // Closed
            class1 = new ClosedClass(model, "class1", 10, delay);
            class2 = new ClosedClass(model, "class2", 10, delay);
        } else if (networkType.equals("o")) { // Open
            class1 = new OpenClass(model, "class1");
            class2 = new OpenClass(model, "class2");
            source.setArrival((OpenClass)class1, new Exp(0.3));
            source.setArrival((OpenClass)class2, new Exp(0.2));
        } else { // Mixed
            class1 = new ClosedClass(model, "class1", 10, delay);
            class2 = new OpenClass(model, "class2");
            source.setArrival((OpenClass)class2, new Exp(0.3));
        }
        
        // Set service distributions at delay
        delay.setService(class1, new Exp(0.5));
        delay.setService(class2, new Exp(0.8));
        
        // Set service distributions at queues
        // Different rates for different iterations
        double[][] rates = {
            {1.5, 1.2, 1.8}, // iteration 1
            {2.0, 1.5, 2.5}, // iteration 2
            {1.0, 2.0, 1.5}, // iteration 3
            {1.8, 1.3, 2.2}  // iteration 4
        };
        
        int idx = (iteration - 1) % rates.length;
        queue1.setService(class1, new Exp(rates[idx][0]));
        queue1.setService(class2, new Exp(rates[idx][0] * 1.2));
        
        queue2.setService(class1, new Exp(rates[idx][1]));
        queue2.setService(class2, new Exp(rates[idx][1] * 1.1));
        
        queue3.setService(class1, new Exp(rates[idx][2]));
        queue3.setService(class2, new Exp(rates[idx][2] * 0.9));
        
        // Set routing
        RoutingMatrix P = model.initRoutingMatrix();
        
        // Set routing based on network type
        if (networkType.equals("c")) { // Closed
            // Delay -> Fork
            P.set(class1, delay, fork, 1.0);
            P.set(class2, delay, fork, 1.0);
            
            // Main fork splits to queue1 and nested fork
            P.set(class1, fork, queue1, 1.0);
            P.set(class1, fork, fork_nested, 1.0);
            P.set(class2, fork, queue1, 1.0);
            P.set(class2, fork, fork_nested, 1.0);
            
            // Nested fork to queue2 and queue3
            P.set(class1, fork_nested, queue2, 1.0);
            P.set(class1, fork_nested, queue3, 1.0);
            P.set(class2, fork_nested, queue2, 1.0);
            P.set(class2, fork_nested, queue3, 1.0);
            
            // Queues to joins
            P.set(class1, queue1, join, 1.0);
            P.set(class2, queue1, join, 1.0);
            P.set(class1, queue2, join_nested, 1.0);
            P.set(class2, queue2, join_nested, 1.0);
            P.set(class1, queue3, join_nested, 1.0);
            P.set(class2, queue3, join_nested, 1.0);
            
            // Nested join to main join
            P.set(class1, join_nested, join, 1.0);
            P.set(class2, join_nested, join, 1.0);
            
            // Main join back to delay
            P.set(class1, join, delay, 1.0);
            P.set(class2, join, delay, 1.0);
            
        } else if (networkType.equals("o")) { // Open
            // Source -> Delay
            P.set(class1, source, delay, 1.0);
            P.set(class2, source, delay, 1.0);
            
            // Delay -> Fork
            P.set(class1, delay, fork, 1.0);
            P.set(class2, delay, fork, 1.0);
            
            // Same fork-join routing as closed
            P.set(class1, fork, queue1, 1.0);
            P.set(class1, fork, fork_nested, 1.0);
            P.set(class2, fork, queue1, 1.0);
            P.set(class2, fork, fork_nested, 1.0);
            
            P.set(class1, fork_nested, queue2, 1.0);
            P.set(class1, fork_nested, queue3, 1.0);
            P.set(class2, fork_nested, queue2, 1.0);
            P.set(class2, fork_nested, queue3, 1.0);
            
            P.set(class1, queue1, join, 1.0);
            P.set(class2, queue1, join, 1.0);
            P.set(class1, queue2, join_nested, 1.0);
            P.set(class2, queue2, join_nested, 1.0);
            P.set(class1, queue3, join_nested, 1.0);
            P.set(class2, queue3, join_nested, 1.0);
            
            P.set(class1, join_nested, join, 1.0);
            P.set(class2, join_nested, join, 1.0);
            
            // Main join to sink
            P.set(class1, join, sink, 1.0);
            P.set(class2, join, sink, 1.0);
            
        } else { // Mixed
            // Closed class routing (delay -> fork -> ... -> join -> delay)
            P.set(class1, delay, fork, 1.0);
            P.set(class1, fork, queue1, 1.0);
            P.set(class1, fork, fork_nested, 1.0);
            P.set(class1, fork_nested, queue2, 1.0);
            P.set(class1, fork_nested, queue3, 1.0);
            P.set(class1, queue1, join, 1.0);
            P.set(class1, queue2, join_nested, 1.0);
            P.set(class1, queue3, join_nested, 1.0);
            P.set(class1, join_nested, join, 1.0);
            P.set(class1, join, delay, 1.0);
            
            // Open class routing (source -> delay -> fork -> ... -> join -> sink)
            P.set(class2, source, delay, 1.0);
            P.set(class2, delay, fork, 1.0);
            P.set(class2, fork, queue1, 1.0);
            P.set(class2, fork, fork_nested, 1.0);
            P.set(class2, fork_nested, queue2, 1.0);
            P.set(class2, fork_nested, queue3, 1.0);
            P.set(class2, queue1, join, 1.0);
            P.set(class2, queue2, join_nested, 1.0);
            P.set(class2, queue3, join_nested, 1.0);
            P.set(class2, join_nested, join, 1.0);
            P.set(class2, join, sink, 1.0);
        }
        
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
     * Run all nested FJ benchmarks
     */
    public static void runAll() {
        runNestedClosedFCFS();
        runNestedClosedPS();
        runNestedOpenFCFS();
        runNestedOpenPS();
        runNestedMixedFCFS();
        runNestedMixedPS();
    }
    
    public static void main(String[] args) {
        jline.GlobalConstants.Verbose = jline.VerboseLevel.SILENT;
        runAll();
    }
}