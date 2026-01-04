/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.bench.fj;

import jline.lang.constant.SchedStrategy;
import java.util.Map;

/**
 * Benchmark for Fork-Join networks with PS scheduling
 * Comprehensive benchmark suite with all configurations
 */
public class BenchFJ_PS {
    
    // Single class variants
    public static void runScLightLoad(int iteration) {
        System.out.println("=== Running FJ PS SC Light Load (iteration " + iteration + ") ===");
        BenchFJ_Template.runBenchmark("bench_FJ_PS_sc_lightload", iteration, 
                                     SchedStrategy.PS, BenchFJ_Template.FJConfig.SC, false, 1);
    }
    
    public static void runScMidLoad(int iteration) {
        System.out.println("=== Running FJ PS SC Mid Load (iteration " + iteration + ") ===");
        BenchFJ_Template.runBenchmark("bench_FJ_PS_sc_midload", iteration, 
                                     SchedStrategy.PS, BenchFJ_Template.FJConfig.SC, false, 3);
    }
    
    public static void runScHighLoad(int iteration) {
        System.out.println("=== Running FJ PS SC High Load (iteration " + iteration + ") ===");
        BenchFJ_Template.runBenchmark("bench_FJ_PS_sc_highload", iteration, 
                                     SchedStrategy.PS, BenchFJ_Template.FJConfig.SC, false, 5);
    }
    
    // Single class high CV variants
    public static void runScHicvLightLoad(int iteration) {
        System.out.println("=== Running FJ PS SC HiCV Light Load (iteration " + iteration + ") ===");
        BenchFJ_Template.runBenchmark("bench_FJ_PS_sc_hicv_lightload", iteration, 
                                     SchedStrategy.PS, BenchFJ_Template.FJConfig.SC, true, 1);
    }
    
    public static void runScHicvMidLoad(int iteration) {
        System.out.println("=== Running FJ PS SC HiCV Mid Load (iteration " + iteration + ") ===");
        BenchFJ_Template.runBenchmark("bench_FJ_PS_sc_hicv_midload", iteration, 
                                     SchedStrategy.PS, BenchFJ_Template.FJConfig.SC, true, 3);
    }
    
    public static void runScHicvHighLoad(int iteration) {
        System.out.println("=== Running FJ PS SC HiCV High Load (iteration " + iteration + ") ===");
        BenchFJ_Template.runBenchmark("bench_FJ_PS_sc_hicv_highload", iteration, 
                                     SchedStrategy.PS, BenchFJ_Template.FJConfig.SC, true, 5);
    }
    
    // Multi-class variants
    public static void runMcLightLoad(int iteration) {
        System.out.println("=== Running FJ PS MC Light Load (iteration " + iteration + ") ===");
        BenchFJ_Template.runBenchmark("bench_FJ_PS_mc_lightload", iteration, 
                                     SchedStrategy.PS, BenchFJ_Template.FJConfig.MC, false, 1);
    }
    
    public static void runMcMidLoad(int iteration) {
        System.out.println("=== Running FJ PS MC Mid Load (iteration " + iteration + ") ===");
        BenchFJ_Template.runBenchmark("bench_FJ_PS_mc_midload", iteration, 
                                     SchedStrategy.PS, BenchFJ_Template.FJConfig.MC, false, 3);
    }
    
    public static void runMcHighLoad(int iteration) {
        System.out.println("=== Running FJ PS MC High Load (iteration " + iteration + ") ===");
        BenchFJ_Template.runBenchmark("bench_FJ_PS_mc_highload", iteration, 
                                     SchedStrategy.PS, BenchFJ_Template.FJConfig.MC, false, 5);
    }
    
    // Multi-class high CV variants
    public static void runMcHicvLightLoad(int iteration) {
        System.out.println("=== Running FJ PS MC HiCV Light Load (iteration " + iteration + ") ===");
        BenchFJ_Template.runBenchmark("bench_FJ_PS_mc_hicv_lightload", iteration, 
                                     SchedStrategy.PS, BenchFJ_Template.FJConfig.MC, true, 1);
    }
    
    public static void runMcHicvMidLoad(int iteration) {
        System.out.println("=== Running FJ PS MC HiCV Mid Load (iteration " + iteration + ") ===");
        BenchFJ_Template.runBenchmark("bench_FJ_PS_mc_hicv_midload", iteration, 
                                     SchedStrategy.PS, BenchFJ_Template.FJConfig.MC, true, 3);
    }
    
    public static void runMcHicvHighLoad(int iteration) {
        System.out.println("=== Running FJ PS MC HiCV High Load (iteration " + iteration + ") ===");
        BenchFJ_Template.runBenchmark("bench_FJ_PS_mc_hicv_highload", iteration, 
                                     SchedStrategy.PS, BenchFJ_Template.FJConfig.MC, true, 5);
    }
    
    // MC2 variants
    public static void runMc2LightLoad(int iteration) {
        System.out.println("=== Running FJ PS MC2 Light Load (iteration " + iteration + ") ===");
        BenchFJ_Template.runBenchmark("bench_FJ_PS_mc2_lightload", iteration, 
                                     SchedStrategy.PS, BenchFJ_Template.FJConfig.MC2, false, 1);
    }
    
    public static void runMc2MidLoad(int iteration) {
        System.out.println("=== Running FJ PS MC2 Mid Load (iteration " + iteration + ") ===");
        BenchFJ_Template.runBenchmark("bench_FJ_PS_mc2_midload", iteration, 
                                     SchedStrategy.PS, BenchFJ_Template.FJConfig.MC2, false, 3);
    }
    
    public static void runMc2HighLoad(int iteration) {
        System.out.println("=== Running FJ PS MC2 High Load (iteration " + iteration + ") ===");
        BenchFJ_Template.runBenchmark("bench_FJ_PS_mc2_highload", iteration, 
                                     SchedStrategy.PS, BenchFJ_Template.FJConfig.MC2, false, 5);
    }
    
    // MC2 high CV variants
    public static void runMc2HicvLightLoad(int iteration) {
        System.out.println("=== Running FJ PS MC2 HiCV Light Load (iteration " + iteration + ") ===");
        BenchFJ_Template.runBenchmark("bench_FJ_PS_mc2_hicv_lightload", iteration, 
                                     SchedStrategy.PS, BenchFJ_Template.FJConfig.MC2, true, 1);
    }
    
    public static void runMc2HicvMidLoad(int iteration) {
        System.out.println("=== Running FJ PS MC2 HiCV Mid Load (iteration " + iteration + ") ===");
        BenchFJ_Template.runBenchmark("bench_FJ_PS_mc2_hicv_midload", iteration, 
                                     SchedStrategy.PS, BenchFJ_Template.FJConfig.MC2, true, 3);
    }
    
    public static void runMc2HicvHighLoad(int iteration) {
        System.out.println("=== Running FJ PS MC2 HiCV High Load (iteration " + iteration + ") ===");
        BenchFJ_Template.runBenchmark("bench_FJ_PS_mc2_hicv_highload", iteration, 
                                     SchedStrategy.PS, BenchFJ_Template.FJConfig.MC2, true, 5);
    }
    
    // MC3 variants
    public static void runMc3LightLoad(int iteration) {
        System.out.println("=== Running FJ PS MC3 Light Load (iteration " + iteration + ") ===");
        BenchFJ_Template.runBenchmark("bench_FJ_PS_mc3_lightload", iteration, 
                                     SchedStrategy.PS, BenchFJ_Template.FJConfig.MC3, false, 1);
    }
    
    public static void runMc3MidLoad(int iteration) {
        System.out.println("=== Running FJ PS MC3 Mid Load (iteration " + iteration + ") ===");
        BenchFJ_Template.runBenchmark("bench_FJ_PS_mc3_midload", iteration, 
                                     SchedStrategy.PS, BenchFJ_Template.FJConfig.MC3, false, 3);
    }
    
    public static void runMc3HighLoad(int iteration) {
        System.out.println("=== Running FJ PS MC3 High Load (iteration " + iteration + ") ===");
        BenchFJ_Template.runBenchmark("bench_FJ_PS_mc3_highload", iteration, 
                                     SchedStrategy.PS, BenchFJ_Template.FJConfig.MC3, false, 5);
    }
    
    public static void runMc3HicvLightLoad(int iteration) {
        System.out.println("=== Running FJ PS MC3 High CV Light Load (iteration " + iteration + ") ===");
        BenchFJ_Template.runBenchmark("bench_FJ_PS_mc3_hicv_lightload", iteration, 
                                     SchedStrategy.PS, BenchFJ_Template.FJConfig.MC3, true, 1);
    }
    
    public static void runMc3HicvMidLoad(int iteration) {
        System.out.println("=== Running FJ PS MC3 High CV Mid Load (iteration " + iteration + ") ===");
        BenchFJ_Template.runBenchmark("bench_FJ_PS_mc3_hicv_midload", iteration, 
                                     SchedStrategy.PS, BenchFJ_Template.FJConfig.MC3, true, 3);
    }
    
    public static void runMc3HicvHighLoad(int iteration) {
        System.out.println("=== Running FJ PS MC3 High CV High Load (iteration " + iteration + ") ===");
        BenchFJ_Template.runBenchmark("bench_FJ_PS_mc3_hicv_highload", iteration, 
                                     SchedStrategy.PS, BenchFJ_Template.FJConfig.MC3, true, 5);
    }
    
    // MC4 variants
    public static void runMc4LightLoad(int iteration) {
        System.out.println("=== Running FJ PS MC4 Light Load (iteration " + iteration + ") ===");
        BenchFJ_Template.runBenchmark("bench_FJ_PS_mc4_lightload", iteration, 
                                     SchedStrategy.PS, BenchFJ_Template.FJConfig.MC4, false, 1);
    }
    
    public static void runMc4MidLoad(int iteration) {
        System.out.println("=== Running FJ PS MC4 Mid Load (iteration " + iteration + ") ===");
        BenchFJ_Template.runBenchmark("bench_FJ_PS_mc4_midload", iteration, 
                                     SchedStrategy.PS, BenchFJ_Template.FJConfig.MC4, false, 3);
    }
    
    public static void runMc4HighLoad(int iteration) {
        System.out.println("=== Running FJ PS MC4 High Load (iteration " + iteration + ") ===");
        BenchFJ_Template.runBenchmark("bench_FJ_PS_mc4_highload", iteration, 
                                     SchedStrategy.PS, BenchFJ_Template.FJConfig.MC4, false, 5);
    }
    
    public static void runMc4HicvLightLoad(int iteration) {
        System.out.println("=== Running FJ PS MC4 High CV Light Load (iteration " + iteration + ") ===");
        BenchFJ_Template.runBenchmark("bench_FJ_PS_mc4_hicv_lightload", iteration, 
                                     SchedStrategy.PS, BenchFJ_Template.FJConfig.MC4, true, 1);
    }
    
    public static void runMc4HicvMidLoad(int iteration) {
        System.out.println("=== Running FJ PS MC4 High CV Mid Load (iteration " + iteration + ") ===");
        BenchFJ_Template.runBenchmark("bench_FJ_PS_mc4_hicv_midload", iteration, 
                                     SchedStrategy.PS, BenchFJ_Template.FJConfig.MC4, true, 3);
    }
    
    public static void runMc4HicvHighLoad(int iteration) {
        System.out.println("=== Running FJ PS MC4 High CV High Load (iteration " + iteration + ") ===");
        BenchFJ_Template.runBenchmark("bench_FJ_PS_mc4_hicv_highload", iteration, 
                                     SchedStrategy.PS, BenchFJ_Template.FJConfig.MC4, true, 5);
    }
    
    /**
     * Run all PS benchmarks
     */
    public static void runAll() {
        for (int i = 1; i <= 10; i++) {
            // Single class
            runScLightLoad(i);
            runScMidLoad(i);
            runScHighLoad(i);
            runScHicvLightLoad(i);
            runScHicvMidLoad(i);
            runScHicvHighLoad(i);
            
            // Multi-class
            runMcLightLoad(i);
            runMcMidLoad(i);
            runMcHighLoad(i);
            runMcHicvLightLoad(i);
            runMcHicvMidLoad(i);
            runMcHicvHighLoad(i);
            
            // MC2
            runMc2LightLoad(i);
            runMc2MidLoad(i);
            runMc2HighLoad(i);
            runMc2HicvLightLoad(i);
            runMc2HicvMidLoad(i);
            runMc2HicvHighLoad(i);
            
            // MC3
            runMc3LightLoad(i);
            runMc3MidLoad(i);
            runMc3HighLoad(i);
            runMc3HicvLightLoad(i);
            runMc3HicvMidLoad(i);
            runMc3HicvHighLoad(i);
            
            // MC4
            runMc4LightLoad(i);
            runMc4MidLoad(i);
            runMc4HighLoad(i);
            runMc4HicvLightLoad(i);
            runMc4HicvMidLoad(i);
            runMc4HicvHighLoad(i);
        }
    }
    
    public static void main(String[] args) {
        jline.GlobalConstants.Verbose = jline.VerboseLevel.SILENT;
        runAll();
    }
}