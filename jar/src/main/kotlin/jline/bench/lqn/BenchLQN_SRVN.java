/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.bench.lqn;

import jline.lang.layered.LayeredNetwork;
import jline.lang.Network;
import jline.lang.constant.SolverType;
import jline.solvers.ln.SolverLN;
import jline.solvers.ln.SolverFactory;
import jline.solvers.NetworkSolver;
import jline.solvers.SolverOptions;
import jline.solvers.AvgTable;
import jline.VerboseLevel;
import jline.solvers.jmt.SolverJMT;

import java.io.File;
import java.io.FilenameFilter;
import java.util.*;

/**
 * LQN benchmark using SRVN solver (Stochastic Rendezvous Networks)
 * Note: Uses JMT as a proxy since native SRVN is external
 */
public class BenchLQN_SRVN {
    
    /**
     * SRVN-style solver factory for LQN benchmarks
     */
    public static class SRVNSolverFactory implements SolverFactory {
        @Override
        public NetworkSolver at(Network model) {
            SolverOptions options = new SolverOptions(SolverType.JMT);
            options.verbose = VerboseLevel.SILENT;
            options.samples = 100000; // SRVN-style high precision
            options.seed = 12345; // Fixed seed for reproducibility
            return new SolverJMT(model, options);
        }
    }
    
    /**
     * Run a single LQN model with SRVN solver
     */
    public static boolean runModel(String modelPath) {
        try {
            System.out.println("  Loading model: " + new File(modelPath).getName());
            
            // Load the LQN model
            LayeredNetwork lqn = LayeredNetwork.parseXML(modelPath);
            
            // Create solver options
            SolverOptions options = new SolverOptions(SolverType.LN);
            options.verbose = VerboseLevel.SILENT;
            options.iter_max = 100;
            options.tol = 1e-6;
            
            // Create solver with SRVN-style configuration
            SolverLN solver = new SolverLN(lqn, new SRVNSolverFactory(), options);
            
            // Get results
            AvgTable avgTable = solver.getAvgTable();
            
            if (avgTable != null && solver.result.iter < options.iter_max) {
                System.out.println("    SUCCESS - SRVN-style solver converged in " + solver.result.iter + " iterations");
                return true;
            } else {
                System.out.println("    FAILED - SRVN-style solver did not converge");
                return false;
            }
            
        } catch (Exception e) {
            System.out.println("    ERROR: " + e.getMessage());
            return false;
        }
    }
    
    /**
     * Run benchmarks for a specific configuration with SRVN solver
     */
    public static void runConfiguration(int clients, int layers, int tasks, int processors, String modelsPath) {
        System.out.printf("\nSRVN Solver - Configuration: C%d_L%d_T%d_P%d\n", clients, layers, tasks, processors);
        
        // Find models matching this configuration
        String pattern = String.format("model_C%d_L%d_T%d_P%d_", clients, layers, tasks, processors);
        File modelsDir = new File(modelsPath);
        
        FilenameFilter filter = (dir, name) -> name.startsWith(pattern) && name.endsWith(".lqnx");
        File[] modelFiles = modelsDir.listFiles(filter);
        
        if (modelFiles == null || modelFiles.length == 0) {
            System.out.println("  No models found");
            return;
        }
        
        int successful = 0;
        int total = Math.min(modelFiles.length, 3); // Test first 3 models (SRVN can be slow)
        
        for (int i = 0; i < total; i++) {
            if (runModel(modelFiles[i].getAbsolutePath())) {
                successful++;
            }
        }
        
        System.out.printf("  SRVN Results: %d/%d successful\n", successful, total);
    }
    
    /**
     * Run all LQN benchmarks with SRVN solver
     */
    public static void runAll() {
        String modelsPath = "/data/line-test.git/bench/bench_LQN/models";
        
        if (!new File(modelsPath).exists()) {
            System.out.println("LQN models directory not found: " + modelsPath);
            return;
        }
        
        System.out.println("=== Running LQN SRVN Solver Benchmarks ===");
        System.out.println("Note: Using JMT as SRVN-style solver proxy");
        
        // Test smaller set of configurations due to SRVN complexity
        int[] layers = {1, 2};
        int[] tasks = {1, 2};
        int[] processors = {1, 2};
        
        for (int l : layers) {
            for (int t : tasks) {
                for (int p : processors) {
                    runConfiguration(1, l, t, p, modelsPath);
                }
            }
        }
        
        System.out.println("\n=== LQN SRVN Benchmarks Complete ===");
    }
    
    /**
     * Main method for testing
     */
    public static void main(String[] args) {
        jline.GlobalConstants.Verbose = jline.VerboseLevel.SILENT;
        runAll();
    }
}