/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.bench.oqn;

import jline.lang.constant.SchedStrategy;

/**
 * Benchmark for Open Queueing Networks with Processor Sharing (PS) scheduling.
 * Aligned with MATLAB benchmarks: bench_OQN_PS_*.m
 */
public class BenchOQN_PS {

    // Load levels based on MATLAB benchmarks (arrival rates for ~30%, ~60%, ~90% utilization)
    private static final double LIGHT_LOAD = 0.3;
    private static final double MID_LOAD = 0.6;
    private static final double HIGH_LOAD = 0.9;

    public static void runLightLoad() {
        runLightLoad(1);
    }

    public static void runLightLoad(int iteration) {
        BenchOQNTemplate.runBenchmark("bench_OQN_PS_lightload", iteration, LIGHT_LOAD, SchedStrategy.PS);
    }

    public static void runMidLoad() {
        runMidLoad(1);
    }

    public static void runMidLoad(int iteration) {
        BenchOQNTemplate.runBenchmark("bench_OQN_PS_midload", iteration, MID_LOAD, SchedStrategy.PS);
    }

    public static void runHighLoad() {
        runHighLoad(1);
    }

    public static void runHighLoad(int iteration) {
        BenchOQNTemplate.runBenchmark("bench_OQN_PS_highload", iteration, HIGH_LOAD, SchedStrategy.PS);
    }

    public static void runMulticlassLightLoad() {
        runMulticlassLightLoad(1);
    }

    public static void runMulticlassLightLoad(int iteration) {
        BenchOQNTemplate.runMulticlassBenchmark("bench_OQN_PS_multiclass_lightload", iteration,
                0.15, 0.15, SchedStrategy.PS);
    }

    public static void runMulticlassMidLoad() {
        runMulticlassMidLoad(1);
    }

    public static void runMulticlassMidLoad(int iteration) {
        BenchOQNTemplate.runMulticlassBenchmark("bench_OQN_PS_multiclass_midload", iteration,
                0.3, 0.3, SchedStrategy.PS);
    }

    public static void runMulticlassHighLoad() {
        runMulticlassHighLoad(1);
    }

    public static void runMulticlassHighLoad(int iteration) {
        BenchOQNTemplate.runMulticlassBenchmark("bench_OQN_PS_multiclass_highload", iteration,
                0.45, 0.45, SchedStrategy.PS);
    }

    public static void runTandemLightLoad() {
        runTandemLightLoad(1);
    }

    public static void runTandemLightLoad(int iteration) {
        BenchOQNTemplate.runTandemBenchmark("bench_OQN_PS_tandem_lightload", iteration,
                LIGHT_LOAD, SchedStrategy.PS);
    }

    public static void runTandemMidLoad() {
        runTandemMidLoad(1);
    }

    public static void runTandemMidLoad(int iteration) {
        BenchOQNTemplate.runTandemBenchmark("bench_OQN_PS_tandem_midload", iteration,
                MID_LOAD, SchedStrategy.PS);
    }

    public static void runTandemHighLoad() {
        runTandemHighLoad(1);
    }

    public static void runTandemHighLoad(int iteration) {
        BenchOQNTemplate.runTandemBenchmark("bench_OQN_PS_tandem_highload", iteration,
                HIGH_LOAD, SchedStrategy.PS);
    }

    /**
     * Run all PS benchmarks
     */
    public static void runAll() {
        for (int it = 1; it <= 10; it++) {
            runLightLoad(it);
            runMidLoad(it);
            runHighLoad(it);
            runMulticlassLightLoad(it);
            runMulticlassMidLoad(it);
            runMulticlassHighLoad(it);
            runTandemLightLoad(it);
            runTandemMidLoad(it);
            runTandemHighLoad(it);
        }
        OQNResultsFormatter.printAndClearGlobalResults();
    }

    public static void main(String[] args) {
        jline.GlobalConstants.Verbose = jline.VerboseLevel.SILENT;
        runAll();
    }
}
