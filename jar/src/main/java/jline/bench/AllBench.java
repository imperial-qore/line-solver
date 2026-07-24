/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.bench;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.bench.cqn.BenchCQN_PS;
import jline.bench.cqn.BenchCQN_PS1;
import jline.bench.cqn.BenchCQN_PS2;
import jline.bench.cqn.BenchCQN_PS3;
import jline.bench.cqn.BenchCQN_PS4;
import jline.bench.cqn.BenchCQN_RM;
import jline.bench.cqn.BenchCQN_FCFS;
import jline.bench.fj.BenchFJ_FCFS;
import jline.bench.fj.BenchFJ_PS;
import jline.bench.fj.BenchFJ_Closed;
import jline.bench.fj.BenchFJ_Open;
import jline.bench.fj.BenchFJ_Mixed;
import jline.bench.fj.BenchFJ_Nested;
import jline.bench.lqn.SimpleLQNBenchmark;
import jline.bench.lqn.BenchLQN_Custom;
import jline.bench.lqn.BenchLQN_Default;
import jline.bench.lqn.BenchLQN_Fluid;
import jline.bench.lqn.BenchLQN_LQNS;
import jline.bench.lqn.BenchLQN_MVA;
import jline.bench.lqn.BenchLQN_NC;
import jline.bench.lqn.BenchLQN_SRVN;
import jline.bench.oqn.BenchOQN_PS;
import jline.bench.oqn.BenchOQN_FCFS;
import jline.bench.oqn.OQNResultsFormatter;
import jline.bench.mqn.BenchMQN_PS;
import jline.bench.mqn.BenchMQN_FCFS;
import jline.bench.mqn.MQNResultsFormatter;

/**
 * Main benchmark runner - runs all benchmark suites
 */
public class AllBench {
    
    public static void main(String[] args) {
        // Set global verbosity to silent to suppress all solver warnings
        GlobalConstants.Verbose = VerboseLevel.SILENT;

        // Also configure Java logging to suppress warnings
        java.util.logging.LogManager.getLogManager().reset();
        java.util.logging.Logger globalLogger = java.util.logging.Logger.getLogger("");
        globalLogger.setLevel(java.util.logging.Level.OFF);

        System.out.println("=== LINE Java Benchmarks ===");
        System.out.println("Starting benchmark suite...\n");

        // Run OQN benchmarks
        runOQNBenchmarks();

        // Run MQN benchmarks
        runMQNBenchmarks();

        // Run CQN benchmarks
        runCQNBenchmarks();

        // Run Fork-Join benchmarks
        runFJBenchmarks();

        // Run LQN benchmarks
        runLQNBenchmarks();

        // Run Random model benchmarks - temporarily disabled due to compilation issues
        // runRandomBenchmarks();

        System.out.println("\n=== All Benchmarks Complete ===");
    }
    
    private static void runOQNBenchmarks() {
        System.out.println("\n--- Running OQN Benchmarks ---");

        // OQN PS benchmarks
        System.out.println("\n>> OQN PS Benchmarks <<");
        BenchOQN_PS.runAll();

        // OQN FCFS benchmarks
        System.out.println("\n>> OQN FCFS Benchmarks <<");
        BenchOQN_FCFS.runAll();

        OQNResultsFormatter.printAndClearGlobalResults();
    }

    private static void runMQNBenchmarks() {
        System.out.println("\n--- Running MQN Benchmarks ---");

        // MQN PS benchmarks
        System.out.println("\n>> MQN PS Benchmarks <<");
        BenchMQN_PS.runAll();

        // MQN FCFS benchmarks
        System.out.println("\n>> MQN FCFS Benchmarks <<");
        BenchMQN_FCFS.runAll();

        MQNResultsFormatter.printAndClearGlobalResults();
    }

    private static void runCQNBenchmarks() {
        System.out.println("\n--- Running CQN Benchmarks ---");
        
        // PS benchmarks
        System.out.println("\n>> CQN PS Benchmarks <<");
        for (int i = 1; i <= 10; i++) {
            BenchCQN_PS.runLightLoad(i);
            BenchCQN_PS.runMidLoad(i);
            BenchCQN_PS.runHighLoad(i);
        }
        
        // PS1 benchmarks
        System.out.println("\n>> CQN PS1 Benchmarks <<");
        for (int i = 1; i <= 10; i++) {
            BenchCQN_PS1.runLightLoad(i);
            BenchCQN_PS1.runMidLoad(i);
            BenchCQN_PS1.runHighLoad(i);
        }
        
        // PS2 benchmarks
        System.out.println("\n>> CQN PS2 Benchmarks <<");
        for (int i = 1; i <= 10; i++) {
            BenchCQN_PS2.runLightLoad(i);
            BenchCQN_PS2.runMidLoad(i);
            BenchCQN_PS2.runHighLoad(i);
        }
        
        // PS3 benchmarks
        System.out.println("\n>> CQN PS3 Benchmarks <<");
        for (int i = 1; i <= 10; i++) {
            BenchCQN_PS3.runLightLoad(i);
            BenchCQN_PS3.runMidLoad(i);
            BenchCQN_PS3.runHighLoad(i);
        }
        
        // PS4 benchmarks
        System.out.println("\n>> CQN PS4 Benchmarks <<");
        for (int i = 1; i <= 10; i++) {
            BenchCQN_PS4.runLightLoad(i);
            BenchCQN_PS4.runMidLoad(i);
            BenchCQN_PS4.runHighLoad(i);
        }
        
        // FCFS benchmarks
        System.out.println("\n>> CQN FCFS Benchmarks <<");
        for (int i = 1; i <= 10; i++) {
            BenchCQN_FCFS.runLightLoad(i);
            BenchCQN_FCFS.runMidLoad(i);
            BenchCQN_FCFS.runHighLoad(i);
        }
        
        // Repairmen (RM) benchmarks with formatted output
        System.out.println("\n>> CQN Repairmen (RM) Benchmarks <<");
        BenchCQN_RM.runAll();
    }
    
    private static void runFJBenchmarks() {
        System.out.println("\n--- Running Fork-Join Benchmarks ---");
        
        // FCFS Fork-Join benchmarks
        System.out.println("\n>> FJ FCFS Benchmarks <<");
        BenchFJ_FCFS.runAll();
        
        // PS Fork-Join benchmarks
        System.out.println("\n>> FJ PS Benchmarks <<");
        BenchFJ_PS.runAll();
        
        // Closed Fork-Join benchmarks
        System.out.println("\n>> FJ Closed Benchmarks <<");
        BenchFJ_Closed.runAll();
        
        // Open Fork-Join benchmarks
        System.out.println("\n>> FJ Open Benchmarks <<");
        BenchFJ_Open.runAll();
        
        // Mixed Fork-Join benchmarks
        System.out.println("\n>> FJ Mixed Benchmarks <<");
        BenchFJ_Mixed.runAll();
        
        // Nested Fork-Join benchmarks
        System.out.println("\n>> FJ Nested Benchmarks <<");
        BenchFJ_Nested.runAll();
    }
    
    private static void runLQNBenchmarks() {
        System.out.println("\n--- Running LQN Benchmarks ---");
        
        // Run all comprehensive LQN solver benchmarks
        System.out.println("\n>> LQN Basic Benchmarks <<");
        SimpleLQNBenchmark.runAll();
        
        System.out.println("\n>> LQN MVA Solver Benchmarks <<");
        BenchLQN_MVA.runAll();
        
        System.out.println("\n>> LQN NC Solver Benchmarks <<");
        BenchLQN_NC.runAll();
        
        System.out.println("\n>> LQN Fluid Solver Benchmarks <<");
        BenchLQN_Fluid.runAll();
        
        System.out.println("\n>> LQN Default Solver Benchmarks <<");
        BenchLQN_Default.runAll();
        
        System.out.println("\n>> LQN Custom Solver Benchmarks <<");
        BenchLQN_Custom.runAll();
        
        System.out.println("\n>> LQN LQNS Solver Benchmarks <<");
        BenchLQN_LQNS.runAll();
        
        System.out.println("\n>> LQN SRVN Solver Benchmarks <<");
        BenchLQN_SRVN.runAll();
    }
    
    private static void runRandomBenchmarks() {
        System.out.println("\n--- Running Random Model Benchmarks ---");
        System.out.println(">> Random benchmarks temporarily disabled due to compilation issues <<");
        // System.out.println("\n>> Random Model Generation and Testing <<");
        // BenchRand.runAll();
    }
}