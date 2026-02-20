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
import jline.solvers.mva.SolverMVA;
import jline.solvers.nc.SolverNC;
import jline.solvers.fluid.SolverFluid;
import jline.solvers.jmt.SolverJMT;

import java.io.File;
import java.io.FilenameFilter;
import java.util.*;

/**
 * Simplified LQN benchmark implementation
 */
public class SimpleLQNBenchmark {
    
    /**
     * Custom solver factory for LQN benchmarks
     */
    public static class BenchmarkSolverFactory implements SolverFactory {
        private final SolverType solverType;
        
        public BenchmarkSolverFactory(SolverType solverType) {
            this.solverType = solverType;
        }
        
        @Override
        public NetworkSolver at(Network model) {
            SolverOptions options = new SolverOptions(solverType);
            options.verbose = VerboseLevel.SILENT;
            
            if (solverType == SolverType.MVA) {
                return new SolverMVA(model, options);
            } else if (solverType == SolverType.NC) {
                return new SolverNC(model, options);
            } else if (solverType == SolverType.FLUID) {
                return new SolverFluid(model, options);
            } else if (solverType == SolverType.JMT) {
                options.samples = 100000;
                return new SolverJMT(model, options);
            } else {
                // Default to MVA for other solver types
                System.out.println("    Warning: Solver type " + solverType + " not explicitly supported, using MVA");
                return new SolverMVA(model, options);
            }
        }
    }
    
    /**
     * Run a single LQN model with a specific solver
     */
    public static boolean runModel(String modelPath, SolverType solverType) {
        try {
            System.out.println("  Loading model: " + new File(modelPath).getName());
            
            // Load the LQN model
            LayeredNetwork lqn = LayeredNetwork.parseXML(modelPath);
            
            // Create solver options
            SolverOptions options = new SolverOptions(SolverType.LN);
            options.verbose = VerboseLevel.SILENT;
            options.iter_max = 100;
            options.tol = 1e-6;
            
            // Create solver
            SolverLN solver = new SolverLN(lqn, new BenchmarkSolverFactory(solverType), options);
            
            // Get results
            AvgTable avgTable = solver.getAvgTable();
            
            if (avgTable != null && solver.result.iter < options.iter_max) {
                System.out.println("    SUCCESS - Converged in " + solver.result.iter + " iterations");
                return true;
            } else {
                System.out.println("    FAILED - Did not converge");
                return false;
            }
            
        } catch (Exception e) {
            System.out.println("    ERROR: " + e.getMessage());
            return false;
        }
    }
    
    /**
     * Run benchmarks for a specific configuration
     */
    public static void runConfiguration(int clients, int layers, int tasks, int processors, 
                                      String modelsPath) {
        System.out.printf("\nConfiguration: C%d_L%d_T%d_P%d\n", clients, layers, tasks, processors);
        
        // Find models matching this configuration
        String pattern = String.format("model_C%d_L%d_T%d_P%d_", clients, layers, tasks, processors);
        File modelsDir = new File(modelsPath);
        
        FilenameFilter filter = (dir, name) -> name.startsWith(pattern) && name.endsWith(".lqnx");
        File[] modelFiles = modelsDir.listFiles(filter);
        
        if (modelFiles == null || modelFiles.length == 0) {
            System.out.println("  No models found");
            return;
        }
        
        // Test with different solvers
        SolverType[] solvers = {SolverType.MVA, SolverType.NC, SolverType.FLUID};
        
        for (SolverType solverType : solvers) {
            System.out.println("\n  Testing with " + solverType + ":");
            
            int successful = 0;
            int total = Math.min(modelFiles.length, 3); // Test first 3 models
            
            for (int i = 0; i < total; i++) {
                if (runModel(modelFiles[i].getAbsolutePath(), solverType)) {
                    successful++;
                }
            }
            
            System.out.printf("  Results: %d/%d successful\n", successful, total);
        }
    }
    
    /**
     * Run all LQN benchmarks
     */
    public static void runAll() {
        String modelsPath = "/data/line-test.git/bench/bench_LQN/models";
        
        if (!new File(modelsPath).exists()) {
            System.out.println("LQN models directory not found: " + modelsPath);
            return;
        }
        
        System.out.println("=== Running LQN Benchmarks ===");
        
        // Test a subset of configurations
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
        
        System.out.println("\n=== LQN Benchmarks Complete ===");
    }
    
    /**
     * Main method for testing
     */
    public static void main(String[] args) {
        jline.GlobalConstants.Verbose = jline.VerboseLevel.SILENT;
        runAll();
    }
}