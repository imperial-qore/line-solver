/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.bench.cqn;

import java.util.Map;
import jline.bench.cqn.CQNResultsFormatter;

/**
 * Benchmark for Closed Queueing Networks with Processor Sharing (PS) scheduling
 * Light load configuration
 */
public class BenchCQN_PS {
    
    public static class BenchmarkResult {
        public double errorQ;
        public double errorU;
        public double errorR;
        public double errorT;
        
        public BenchmarkResult(double errorQ, double errorU, double errorR, double errorT) {
            this.errorQ = errorQ;
            this.errorU = errorU;
            this.errorR = errorR;
            this.errorT = errorT;
        }
    }
    
    public static void runLightLoad() {
        runLightLoad(1);
    }
    
    public static void runLightLoad(int iteration) {
        System.out.println("=== Running CQN PS Light Load Benchmark (iteration " + iteration + ") ===");
        
        // Standard variant
        BenchCQNTemplate.runBenchmark("CQN_PS_lightload", iteration, 5, 5, 1, false, false);
        
        // High CV variant
        BenchCQNTemplate.runBenchmark("CQN_PS_hicv_lightload", iteration, 5, 5, 1, true, false);
        
        // Multi-server variant
        BenchCQNTemplate.runBenchmark("CQN_PS_multiserver_lightload", iteration, 5, 5, 2, false, false);
        
        // Multi-server high CV variant
        BenchCQNTemplate.runBenchmark("CQN_PS_multiserver_hicv_lightload", iteration, 5, 5, 2, true, false);
        
        // Random mapping variant
        BenchCQNTemplate.runBenchmark("CQN_PS_rm_lightload", iteration, 5, 5, 1, false, true);
        
        // Random mapping high CV variant
        BenchCQNTemplate.runBenchmark("CQN_PS_rm_hicv_lightload", iteration, 5, 5, 1, true, true);
        
        // Random mapping multi-server variant
        BenchCQNTemplate.runBenchmark("CQN_PS_rm_multiserver_lightload", iteration, 5, 5, 2, false, true);
        
        // Random mapping multi-server high CV variant
        BenchCQNTemplate.runBenchmark("CQN_PS_rm_multiserver_hicv_lightload", iteration, 5, 5, 2, true, true);
    }
    
    public static void runMidLoad() {
        runMidLoad(1);
    }
    
    public static void runMidLoad(int iteration) {
        System.out.println("=== Running CQN PS Mid Load Benchmark (iteration " + iteration + ") ===");
        
        // Standard variant
        BenchCQNTemplate.runBenchmark("CQN_PS_midload", iteration, 15, 15, 1, false, false);
        
        // High CV variant
        BenchCQNTemplate.runBenchmark("CQN_PS_hicv_midload", iteration, 15, 15, 1, true, false);
        
        // Multi-server variant
        BenchCQNTemplate.runBenchmark("CQN_PS_multiserver_midload", iteration, 15, 15, 2, false, false);
        
        // Multi-server high CV variant
        BenchCQNTemplate.runBenchmark("CQN_PS_multiserver_hicv_midload", iteration, 15, 15, 2, true, false);
        
        // Random mapping variant
        BenchCQNTemplate.runBenchmark("CQN_PS_rm_midload", iteration, 15, 15, 1, false, true);
        
        // Random mapping high CV variant
        BenchCQNTemplate.runBenchmark("CQN_PS_rm_hicv_midload", iteration, 15, 15, 1, true, true);
        
        // Random mapping multi-server variant
        BenchCQNTemplate.runBenchmark("CQN_PS_rm_multiserver_midload", iteration, 15, 15, 2, false, true);
        
        // Random mapping multi-server high CV variant
        BenchCQNTemplate.runBenchmark("CQN_PS_rm_multiserver_hicv_midload", iteration, 15, 15, 2, true, true);
    }
    
    public static void runHighLoad() {
        runHighLoad(1);
    }
    
    public static void runHighLoad(int iteration) {
        System.out.println("=== Running CQN PS High Load Benchmark (iteration " + iteration + ") ===");
        
        // Standard variant
        BenchCQNTemplate.runBenchmark("CQN_PS_highload", iteration, 30, 30, 1, false, false);
        
        // High CV variant
        BenchCQNTemplate.runBenchmark("CQN_PS_hicv_highload", iteration, 30, 30, 1, true, false);
        
        // Multi-server variant
        BenchCQNTemplate.runBenchmark("CQN_PS_multiserver_highload", iteration, 30, 30, 2, false, false);
        
        // Multi-server high CV variant
        BenchCQNTemplate.runBenchmark("CQN_PS_multiserver_hicv_highload", iteration, 30, 30, 2, true, false);
        
        // Random mapping variant
        BenchCQNTemplate.runBenchmark("CQN_PS_rm_highload", iteration, 30, 30, 1, false, true);
        
        // Random mapping high CV variant
        BenchCQNTemplate.runBenchmark("CQN_PS_rm_hicv_highload", iteration, 30, 30, 1, true, true);
        
        // Random mapping multi-server variant
        BenchCQNTemplate.runBenchmark("CQN_PS_rm_multiserver_highload", iteration, 30, 30, 2, false, true);
        
        // Random mapping multi-server high CV variant
        BenchCQNTemplate.runBenchmark("CQN_PS_rm_multiserver_hicv_highload", iteration, 30, 30, 2, true, true);
    }
    
    /**
     * Run only the repairmen (RM) benchmarks with formatted output
     */
    public static void runRMBenchmarks() {
        System.out.println("=== Running CQN PS Repairmen (RM) Benchmarks ===");
        
        // Run RM-specific benchmarks for all load levels
        runRMLightLoad();
        runRMMidLoad();
        runRMHighLoad();
        
        // Print all accumulated results
        CQNResultsFormatter.printAndClearGlobalResults();
    }
    
    /**
     * Run only RM light load benchmarks
     */
    public static void runRMLightLoad() {
        // Only run RM variants
        BenchCQNTemplate.runBenchmark("bench_CQN_PS_rm_lightload", 1, 5, 5, 1, false, true);
        BenchCQNTemplate.runBenchmark("bench_CQN_PS_rm_hicv_lightload", 1, 5, 5, 1, true, true);
        BenchCQNTemplate.runBenchmark("bench_CQN_PS_rm_multiserver_lightload", 1, 5, 5, 2, false, true);
        BenchCQNTemplate.runBenchmark("bench_CQN_PS_rm_multiserver_hicv_lightload", 1, 5, 5, 2, true, true);
    }
    
    /**
     * Run only RM mid load benchmarks
     */
    public static void runRMMidLoad() {
        // Only run RM variants
        BenchCQNTemplate.runBenchmark("bench_CQN_PS_rm_midload", 1, 15, 15, 1, false, true);
        BenchCQNTemplate.runBenchmark("bench_CQN_PS_rm_hicv_midload", 1, 15, 15, 1, true, true);
        BenchCQNTemplate.runBenchmark("bench_CQN_PS_rm_multiserver_midload", 1, 15, 15, 2, false, true);
        BenchCQNTemplate.runBenchmark("bench_CQN_PS_rm_multiserver_hicv_midload", 1, 15, 15, 2, true, true);
    }
    
    /**
     * Run only RM high load benchmarks
     */
    public static void runRMHighLoad() {
        // Only run RM variants
        BenchCQNTemplate.runBenchmark("bench_CQN_PS_rm_highload", 1, 30, 30, 1, false, true);
        BenchCQNTemplate.runBenchmark("bench_CQN_PS_rm_hicv_highload", 1, 30, 30, 1, true, true);
        BenchCQNTemplate.runBenchmark("bench_CQN_PS_rm_multiserver_highload", 1, 30, 30, 2, false, true);
        BenchCQNTemplate.runBenchmark("bench_CQN_PS_rm_multiserver_hicv_highload", 1, 30, 30, 2, true, true);
    }
    
    public static void main(String[] args) {
        jline.GlobalConstants.Verbose = jline.VerboseLevel.SILENT;
        // Check if we should run only RM benchmarks
        if (args.length > 0 && "rm".equals(args[0])) {
            runRMBenchmarks();
        } else {
            // Run all benchmarks
            runLightLoad();
            runMidLoad();
            runHighLoad();
            
            // Print all accumulated results
            CQNResultsFormatter.printAndClearGlobalResults();
        }
    }
}