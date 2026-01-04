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

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Benchmark for Fork-Join Mixed Networks (combination of open and closed classes)
 * Tests both homogeneous and heterogeneous configurations
 */
public class BenchFJ_Mixed {
    
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
        System.out.println("\n=== Running FJ Mixed Homogeneous FCFS Benchmarks ===");
        
        // Run variants 1 and 2
        runBenchmark("homogeneous_FCFS", 1, SchedStrategy.FCFS, true, 2);
        runBenchmark("homogeneous_FCFS", 2, SchedStrategy.FCFS, true, 2);
    }
    
    /**
     * Run homogeneous PS benchmarks
     */
    public static void runHomogeneousPS() {
        System.out.println("\n=== Running FJ Mixed Homogeneous PS Benchmarks ===");
        
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
        System.out.println("\n=== Running FJ Mixed Heterogeneous FCFS Benchmarks ===");
        
        // Run variants 1-3
        runBenchmark("heterogeneous_FCFS", 1, SchedStrategy.FCFS, false, 2);
        runBenchmark("heterogeneous_FCFS", 2, SchedStrategy.FCFS, false, 3);
        runBenchmark("heterogeneous_FCFS", 3, SchedStrategy.FCFS, false, 2);
    }
    
    /**
     * Run heterogeneous PS benchmarks
     */
    public static void runHeterogeneousPS() {
        System.out.println("\n=== Running FJ Mixed Heterogeneous PS Benchmarks ===");
        
        // Run variants 1-5
        runBenchmark("heterogeneous_PS", 1, SchedStrategy.PS, false, 2);
        runBenchmark("heterogeneous_PS", 2, SchedStrategy.PS, false, 3);
        runBenchmark("heterogeneous_PS", 3, SchedStrategy.PS, false, 2);
        runBenchmark("heterogeneous_PS", 4, SchedStrategy.PS, false, 2);
        runBenchmark("heterogeneous_PS", 5, SchedStrategy.PS, false, 3);
    }
    
    /**
     * Core benchmark implementation for mixed networks
     */
    private static BenchmarkResult runBenchmark(String variant, int iteration, 
                                               SchedStrategy sched, boolean homogeneous, int numQueues) {
        String benchName = "bench_FJ_m_" + variant + "_" + iteration;
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
        
        // Create job classes - one closed and one open class
        ClosedClass closedClass = new ClosedClass(model, "closedClass", 10, delay);
        OpenClass openClass = new OpenClass(model, "openClass");
        
        // Set arrival rate for open class
        source.setArrival(openClass, new Exp(0.3));
        
        // Set service distributions at delay
        delay.setService(closedClass, new Exp(0.5));
        delay.setService(openClass, new Exp(0.8));
        
        // Set service distributions at queues
        if (homogeneous) {
            // All queues have same rates for each class
            for (Queue queue : queues) {
                queue.setService(closedClass, new Exp(1.5));
                queue.setService(openClass, new Exp(2.0));
            }
        } else {
            // Heterogeneous - different rates at each queue
            double[] ratesClosed = {1.5, 1.1, 2.5}; // Extend if needed
            double[] ratesOpen = {2.0, 1.8, 2.2}; // Extend if needed
            
            for (int i = 0; i < queues.size(); i++) {
                queues.get(i).setService(closedClass, new Exp(ratesClosed[i % ratesClosed.length]));
                queues.get(i).setService(openClass, new Exp(ratesOpen[i % ratesOpen.length]));
            }
        }
        
        // Set routing
        RoutingMatrix P = model.initRoutingMatrix();
        
        // Open class routing: Source -> Delay -> Fork -> Queues -> Join -> Sink
        P.set(openClass, source, delay, 1.0);
        P.set(openClass, delay, fork, 1.0);
        for (Queue queue : queues) {
            P.set(openClass, fork, queue, 1.0);
        }
        for (Queue queue : queues) {
            P.set(openClass, queue, join, 1.0);
        }
        P.set(openClass, join, sink, 1.0);
        
        // Closed class routing: Delay -> Fork -> Queues -> Join -> Delay (loop)
        P.set(closedClass, delay, fork, 1.0);
        for (Queue queue : queues) {
            P.set(closedClass, fork, queue, 1.0);
        }
        for (Queue queue : queues) {
            P.set(closedClass, queue, join, 1.0);
        }
        P.set(closedClass, join, delay, 1.0); // Loop back for closed class
        
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
     * Run all mixed FJ benchmarks
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