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

import java.io.File;
import java.io.FilenameFilter;
import java.util.*;

/**
 * LQN benchmark using Custom solver configurations with various parameters
 */
public class BenchLQN_Custom {
    
    /**
     * Custom solver factory with configurable parameters
     */
    public static class CustomSolverFactory implements SolverFactory {
        private final SolverType solverType;
        private final double tolerance;
        private final int maxIterations;
        
        public CustomSolverFactory(SolverType solverType, double tolerance, int maxIterations) {
            this.solverType = solverType;
            this.tolerance = tolerance;
            this.maxIterations = maxIterations;
        }
        
        @Override
        public NetworkSolver at(Network model) {
            SolverOptions options = new SolverOptions(solverType);
            options.verbose = VerboseLevel.SILENT;
            options.tol = tolerance;
            options.iter_max = maxIterations;
            
            if (solverType == SolverType.MVA) {
                return new SolverMVA(model, options);
            } else if (solverType == SolverType.NC) {
                return new SolverNC(model, options);
            } else if (solverType == SolverType.FLUID) {
                return new SolverFluid(model, options);
            } else {
                // Default to MVA
                return new SolverMVA(model, options);
            }
        }
    }
    
    /**
     * Run a single LQN model with Custom solver configuration
     */
    public static boolean runModel(String modelPath, SolverType solverType, double tolerance, int maxIterations) {
        try {
            System.out.println("  Loading model: " + new File(modelPath).getName());
            
            // Load the LQN model
            LayeredNetwork lqn = LayeredNetwork.parseXML(modelPath);
            
            // Create solver options
            SolverOptions options = new SolverOptions(SolverType.LN);
            options.verbose = VerboseLevel.SILENT;
            options.iter_max = maxIterations;
            options.tol = tolerance;
            
            // Create solver with Custom configuration
            SolverLN solver = new SolverLN(lqn, new CustomSolverFactory(solverType, tolerance, maxIterations), options);
            
            // Get results
            AvgTable avgTable = solver.getAvgTable();
            
            if (avgTable != null && solver.result.iter < options.iter_max) {
                System.out.println("    SUCCESS - Custom " + solverType + " solver converged in " + solver.result.iter + " iterations");
                return true;
            } else {
                System.out.println("    FAILED - Custom " + solverType + " solver did not converge");
                return false;
            }
            
        } catch (Exception e) {
            System.out.println("    ERROR: " + e.getMessage());
            return false;
        }
    }
    
    /**
     * Run benchmarks for a specific configuration with Custom solver
     */
    public static void runConfiguration(int clients, int layers, int tasks, int processors, String modelsPath) {
        System.out.printf("\nCustom Solver - Configuration: C%d_L%d_T%d_P%d\n", clients, layers, tasks, processors);
        
        // Find models matching this configuration
        String pattern = String.format("model_C%d_L%d_T%d_P%d_", clients, layers, tasks, processors);
        File modelsDir = new File(modelsPath);
        
        FilenameFilter filter = (dir, name) -> name.startsWith(pattern) && name.endsWith(".lqnx");
        File[] modelFiles = modelsDir.listFiles(filter);
        
        if (modelFiles == null || modelFiles.length == 0) {
            System.out.println("  No models found");
            return;
        }
        
        // Test different custom configurations
        Object[][] configs = {
            {SolverType.MVA, 1e-6, 50},
            {SolverType.MVA, 1e-8, 200},
            {SolverType.NC, 1e-6, 50},
            {SolverType.FLUID, 1e-5, 100}
        };
        
        for (Object[] config : configs) {
            SolverType solverType = (SolverType) config[0];
            double tolerance = (Double) config[1];
            int maxIter = (Integer) config[2];
            
            System.out.printf("  Testing Custom %s (tol=%.1e, maxIter=%d):\n", solverType, tolerance, maxIter);
            
            int successful = 0;
            int total = Math.min(modelFiles.length, 3); // Test first 3 models per config
            
            for (int i = 0; i < total; i++) {
                if (runModel(modelFiles[i].getAbsolutePath(), solverType, tolerance, maxIter)) {
                    successful++;
                }
            }
            
            System.out.printf("    Results: %d/%d successful\n", successful, total);
        }
    }
    
    /**
     * Run all LQN benchmarks with Custom solver
     */
    public static void runAll() {
        String modelsPath = "/data/line-test.git/bench/bench_LQN/models";
        
        if (!new File(modelsPath).exists()) {
            System.out.println("LQN models directory not found: " + modelsPath);
            return;
        }
        
        System.out.println("=== Running LQN Custom Solver Benchmarks ===");
        
        // Test subset of configurations with various custom settings
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
        
        System.out.println("\n=== LQN Custom Benchmarks Complete ===");
    }
    
    /**
     * Main method for testing
     */
    public static void main(String[] args) {
        jline.GlobalConstants.Verbose = jline.VerboseLevel.SILENT;
        runAll();
    }
}