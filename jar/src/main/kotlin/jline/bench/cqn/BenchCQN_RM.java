/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.bench.cqn;

import jline.bench.cqn.CQNResultsFormatter;

/**
 * Dedicated benchmark runner for CQN Repairmen (RM) models
 * Provides MATLAB-style formatted output
 */
public class BenchCQN_RM {
    
    /**
     * Run all repairmen benchmarks for both FCFS and PS scheduling
     */
    public static void runAll() {
        System.out.println("=== Running All CQN Repairmen (RM) Benchmarks ===\n");
        
        // Run FCFS RM benchmarks
        runFCFSBenchmarks();
        
        // Run PS RM benchmarks  
        runPSBenchmarks();
        
        // Print all accumulated results in MATLAB format
        CQNResultsFormatter.printAndClearGlobalResults();
    }
    
    /**
     * Run only PS repairmen benchmarks
     */
    public static void runPSBenchmarks() {
        System.out.println("--- Running PS Repairmen Benchmarks ---");
        
        // Light load benchmarks
        BenchCQNTemplate.runBenchmark("bench_CQN_PS_rm_lightload", 1, 5, 5, 1, false, true);
        BenchCQNTemplate.runBenchmark("bench_CQN_PS_rm_hicv_lightload", 1, 5, 5, 1, true, true);
        BenchCQNTemplate.runBenchmark("bench_CQN_PS_rm_multiserver_lightload", 1, 5, 5, 2, false, true);
        BenchCQNTemplate.runBenchmark("bench_CQN_PS_rm_multiserver_hicv_lightload", 1, 5, 5, 2, true, true);
        
        // Mid load benchmarks
        BenchCQNTemplate.runBenchmark("bench_CQN_PS_rm_midload", 1, 15, 15, 1, false, true);
        BenchCQNTemplate.runBenchmark("bench_CQN_PS_rm_hicv_midload", 1, 15, 15, 1, true, true);
        BenchCQNTemplate.runBenchmark("bench_CQN_PS_rm_multiserver_midload", 1, 15, 15, 2, false, true);
        BenchCQNTemplate.runBenchmark("bench_CQN_PS_rm_multiserver_hicv_midload", 1, 15, 15, 2, true, true);
        
        // High load benchmarks
        BenchCQNTemplate.runBenchmark("bench_CQN_PS_rm_highload", 1, 30, 30, 1, false, true);
        BenchCQNTemplate.runBenchmark("bench_CQN_PS_rm_hicv_highload", 1, 30, 30, 1, true, true);
        BenchCQNTemplate.runBenchmark("bench_CQN_PS_rm_multiserver_highload", 1, 30, 30, 2, false, true);
        BenchCQNTemplate.runBenchmark("bench_CQN_PS_rm_multiserver_hicv_highload", 1, 30, 30, 2, true, true);
    }
    
    /**
     * Run only FCFS repairmen benchmarks
     */
    public static void runFCFSBenchmarks() {
        System.out.println("--- Running FCFS Repairmen Benchmarks ---");
        
        // Light load benchmarks
        BenchCQNTemplate.runBenchmark("bench_CQN_FCFS_rm_lightload", 1, 5, 5, 1, false, true);
        BenchCQNTemplate.runBenchmark("bench_CQN_FCFS_rm_hicv_lightload", 1, 5, 5, 1, true, true);
        BenchCQNTemplate.runBenchmark("bench_CQN_FCFS_rm_multiserver_lightload", 1, 5, 5, 2, false, true);
        BenchCQNTemplate.runBenchmark("bench_CQN_FCFS_rm_multiserver_hicv_lightload", 1, 5, 5, 2, true, true);
        
        // Mid load benchmarks
        BenchCQNTemplate.runBenchmark("bench_CQN_FCFS_rm_midload", 1, 15, 15, 1, false, true);
        BenchCQNTemplate.runBenchmark("bench_CQN_FCFS_rm_hicv_midload", 1, 15, 15, 1, true, true);
        BenchCQNTemplate.runBenchmark("bench_CQN_FCFS_rm_multiserver_midload", 1, 15, 15, 2, false, true);
        BenchCQNTemplate.runBenchmark("bench_CQN_FCFS_rm_multiserver_hicv_midload", 1, 15, 15, 2, true, true);
        
        // High load benchmarks
        BenchCQNTemplate.runBenchmark("bench_CQN_FCFS_rm_highload", 1, 30, 30, 1, false, true);
        BenchCQNTemplate.runBenchmark("bench_CQN_FCFS_rm_hicv_highload", 1, 30, 30, 1, true, true);
        BenchCQNTemplate.runBenchmark("bench_CQN_FCFS_rm_multiserver_highload", 1, 30, 30, 2, false, true);
        BenchCQNTemplate.runBenchmark("bench_CQN_FCFS_rm_multiserver_hicv_highload", 1, 30, 30, 2, true, true);
    }
    
    /**
     * Run PS repairmen benchmarks only, with immediate formatted output
     */
    public static void runPSOnly() {
        System.out.println("=== Running CQN PS Repairmen (RM) Benchmarks ===\n");
        runPSBenchmarks();
        CQNResultsFormatter.printAndClearGlobalResults();
    }
    
    /**
     * Run FCFS repairmen benchmarks only, with immediate formatted output
     */
    public static void runFCFSOnly() {
        System.out.println("=== Running CQN FCFS Repairmen (RM) Benchmarks ===\n");
        runFCFSBenchmarks();
        CQNResultsFormatter.printAndClearGlobalResults();
    }
    
    /**
     * Main method with command line options
     */
    public static void main(String[] args) {
        jline.GlobalConstants.Verbose = jline.VerboseLevel.SILENT;
        if (args.length > 0) {
            switch (args[0].toLowerCase()) {
                case "ps":
                    runPSOnly();
                    break;
                case "fcfs":
                    runFCFSOnly();
                    break;
                case "all":
                default:
                    runAll();
                    break;
            }
        } else {
            runAll();
        }
    }
}