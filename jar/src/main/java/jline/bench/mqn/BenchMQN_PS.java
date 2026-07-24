/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.bench.mqn;

import jline.lang.constant.SchedStrategy;

/**
 * Benchmark for Mixed Queueing Networks with Processor Sharing (PS) scheduling.
 * Aligned with MATLAB benchmarks: bench_MQN_PS_*.m
 */
public class BenchMQN_PS {

    // Load levels (open arrival rates)
    private static final double LIGHT_LOAD = 0.15;
    private static final double MID_LOAD = 0.3;
    private static final double HIGH_LOAD = 0.45;

    // Closed class populations
    private static final int LOW_POP = 2;
    private static final int HIGH_POP = 5;

    public static void runLightLoad() {
        runLightLoad(1);
    }

    public static void runLightLoad(int iteration) {
        BenchMQNTemplate.runBenchmark("bench_MQN_PS_lightload", iteration, LIGHT_LOAD, LOW_POP, SchedStrategy.PS);
    }

    public static void runMidLoad() {
        runMidLoad(1);
    }

    public static void runMidLoad(int iteration) {
        BenchMQNTemplate.runBenchmark("bench_MQN_PS_midload", iteration, MID_LOAD, LOW_POP, SchedStrategy.PS);
    }

    public static void runHighLoad() {
        runHighLoad(1);
    }

    public static void runHighLoad(int iteration) {
        BenchMQNTemplate.runBenchmark("bench_MQN_PS_highload", iteration, HIGH_LOAD, LOW_POP, SchedStrategy.PS);
    }

    public static void runHipopLightLoad() {
        runHipopLightLoad(1);
    }

    public static void runHipopLightLoad(int iteration) {
        BenchMQNTemplate.runHighPopBenchmark("bench_MQN_PS_hipop_lightload", iteration,
                LIGHT_LOAD, HIGH_POP, SchedStrategy.PS);
    }

    public static void runHipopMidLoad() {
        runHipopMidLoad(1);
    }

    public static void runHipopMidLoad(int iteration) {
        BenchMQNTemplate.runHighPopBenchmark("bench_MQN_PS_hipop_midload", iteration,
                MID_LOAD, HIGH_POP, SchedStrategy.PS);
    }

    public static void runHipopHighLoad() {
        runHipopHighLoad(1);
    }

    public static void runHipopHighLoad(int iteration) {
        BenchMQNTemplate.runHighPopBenchmark("bench_MQN_PS_hipop_highload", iteration,
                HIGH_LOAD, HIGH_POP, SchedStrategy.PS);
    }

    public static void runTandemLightLoad() {
        runTandemLightLoad(1);
    }

    public static void runTandemLightLoad(int iteration) {
        BenchMQNTemplate.runTandemBenchmark("bench_MQN_PS_tandem_lightload", iteration,
                LIGHT_LOAD, LOW_POP, SchedStrategy.PS);
    }

    public static void runTandemMidLoad() {
        runTandemMidLoad(1);
    }

    public static void runTandemMidLoad(int iteration) {
        BenchMQNTemplate.runTandemBenchmark("bench_MQN_PS_tandem_midload", iteration,
                MID_LOAD, LOW_POP, SchedStrategy.PS);
    }

    public static void runTandemHighLoad() {
        runTandemHighLoad(1);
    }

    public static void runTandemHighLoad(int iteration) {
        BenchMQNTemplate.runTandemBenchmark("bench_MQN_PS_tandem_highload", iteration,
                HIGH_LOAD, LOW_POP, SchedStrategy.PS);
    }

    /**
     * Run all PS benchmarks
     */
    public static void runAll() {
        for (int it = 1; it <= 10; it++) {
            runLightLoad(it);
            runMidLoad(it);
            runHighLoad(it);
            runHipopLightLoad(it);
            runHipopMidLoad(it);
            runHipopHighLoad(it);
            runTandemLightLoad(it);
            runTandemMidLoad(it);
            runTandemHighLoad(it);
        }
        MQNResultsFormatter.printAndClearGlobalResults();
    }

    public static void main(String[] args) {
        jline.GlobalConstants.Verbose = jline.VerboseLevel.SILENT;
        runAll();
    }
}
