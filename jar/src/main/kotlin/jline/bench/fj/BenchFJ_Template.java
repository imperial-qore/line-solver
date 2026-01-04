/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.bench.fj;

import jline.util.matrix.Matrix;
import jline.bench.BenchmarkSolvers;
import jline.bench.BenchmarkUtils;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.processes.Exp;
import jline.lang.processes.Cox2;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Fork;
import jline.lang.nodes.Join;
import jline.lang.nodes.Queue;
import jline.solvers.SolverResult;
import jline.solvers.Solver;
import jline.solvers.NetworkSolver;
import jline.solvers.jmt.SolverJMT;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Template for Fork-Join benchmark implementations
 * Supports all configurations: sc, mc, mc2, mc3, mc4 with standard and high CV
 */
public class BenchFJ_Template {
    
    public enum FJConfig {
        SC,   // Single class
        MC,   // Multi-class (2 classes)
        MC2,  // Multi-class with extra station after join
        MC3,  // Multi-class with extra station in sequence on one branch
        MC4   // Multi-class with more complex topology
    }
    
    /**
     * Core benchmark implementation
     */
    public static Map<String, Object> runBenchmark(String benchmarkName, int iteration, 
                                                   SchedStrategy sched, FJConfig config, 
                                                   boolean highCV, int population) {
        
        // Create model
        Network model = new Network(benchmarkName + "_" + iteration);
        
        // Create base stations
        Delay delay1 = new Delay(model, "Delay1");
        Fork fork = new Fork(model, "Fork");
        Join join = new Join(model, "Join", fork);
        Queue queue1 = new Queue(model, "Queue1", sched);
        Queue queue2 = new Queue(model, "Queue2", sched);
        Queue queue3 = null;
        Queue queue4 = null;
        Delay delay2 = null;
        
        // Create additional stations based on configuration
        if (config == FJConfig.MC2) {
            // Extra station after join
            queue3 = new Queue(model, "Queue3", sched);
        } else if (config == FJConfig.MC3) {
            // Extra station in sequence with queue2
            queue3 = new Queue(model, "Queue3", sched);
        } else if (config == FJConfig.MC4) {
            // MC4: Two extra queues and one extra delay
            // Topology: Delay1 -> Fork -> (Queue1, Queue2, Queue4)
            //           Queue4 -> Delay2 -> Join
            //           Join -> Queue3 -> Delay1
            queue3 = new Queue(model, "Queue3", sched);
            queue4 = new Queue(model, "Queue4", sched);
            delay2 = new Delay(model, "Delay2");
        }
        
        // Create job classes
        List<ClosedClass> classes = new ArrayList<>();
        int numClasses = (config == FJConfig.SC) ? 1 : 2;
        
        for (int i = 1; i <= numClasses; i++) {
            classes.add(new ClosedClass(model, "class" + i, population, delay1));
        }
        
        // Generate service rates
        int numStations = 3; // base stations
        if (config == FJConfig.MC2 || config == FJConfig.MC3) numStations = 4;
        if (config == FJConfig.MC4) numStations = 6; // MC4 has 6 stations total
        
        Matrix rate = BenchmarkUtils.randGallery(numStations, numClasses, iteration);
        
        // Set service distributions
        for (int c = 0; c < numClasses; c++) {
            ClosedClass jobClass = classes.get(c);
            
            // In MATLAB multi-class benchmarks, all classes use the same service rate (column 0)
            // Only single-class uses different rates per class
            int rateCol = (config == FJConfig.SC) ? c : 0;
            
            if (highCV) {
                // High coefficient of variation (SCV = 10)
                delay1.setService(jobClass, Cox2.fitMeanAndSCV(1.0 / rate.get(0, rateCol), 10.0));
                queue1.setService(jobClass, Cox2.fitMeanAndSCV(1.0 / rate.get(1, rateCol), 10.0));
                queue2.setService(jobClass, Cox2.fitMeanAndSCV(1.0 / rate.get(2, rateCol), 10.0));
                
                if (queue3 != null) {
                    queue3.setService(jobClass, Cox2.fitMeanAndSCV(1.0 / rate.get(3, rateCol), 10.0));
                }
                if (queue4 != null) {
                    queue4.setService(jobClass, Cox2.fitMeanAndSCV(1.0 / rate.get(4, rateCol), 10.0));
                }
                if (delay2 != null) {
                    delay2.setService(jobClass, Cox2.fitMeanAndSCV(1.0 / rate.get(5, rateCol), 10.0));
                }
            } else {
                // Exponential distribution
                delay1.setService(jobClass, new Exp(rate.get(0, rateCol)));
                queue1.setService(jobClass, new Exp(rate.get(1, rateCol)));
                queue2.setService(jobClass, new Exp(rate.get(2, rateCol)));
                
                if (queue3 != null) {
                    queue3.setService(jobClass, new Exp(rate.get(3, rateCol)));
                }
                if (queue4 != null) {
                    queue4.setService(jobClass, new Exp(rate.get(4, rateCol)));
                }
                if (delay2 != null) {
                    delay2.setService(jobClass, new Exp(rate.get(5, rateCol)));
                }
            }
        }
        
        // Set routing based on configuration
        RoutingMatrix P = model.initRoutingMatrix();
        
        for (ClosedClass jobClass : classes) {
            // Common routing: delay -> fork
            P.set(jobClass, delay1, fork, 1.0);
            
            // Fork to queues
            P.set(jobClass, fork, queue1, 1.0);
            
            if (config == FJConfig.SC || config == FJConfig.MC || config == FJConfig.MC2) {
                // Standard fork to both queues
                P.set(jobClass, fork, queue2, 1.0);
                
                // Queues to join
                P.set(jobClass, queue1, join, 1.0);
                P.set(jobClass, queue2, join, 1.0);
                
                if (config == FJConfig.MC2) {
                    // MC2: Class 1 goes through queue3, Class 2 bypasses it
                    int classIndex = classes.indexOf(jobClass);
                    if (classIndex == 0) {
                        // Class 1: Join to queue3, then back to delay
                        P.set(jobClass, join, queue3, 1.0);
                        P.set(jobClass, queue3, delay1, 1.0);
                    } else {
                        // Class 2: Bypass queue3, go directly to delay
                        P.set(jobClass, join, delay1, 1.0);
                    }
                } else {
                    // Join back to delay
                    P.set(jobClass, join, delay1, 1.0);
                }
                
            } else if (config == FJConfig.MC3) {
                // Fork to queue2
                P.set(jobClass, fork, queue2, 1.0);
                
                // Queue2 to queue3 (sequential)
                P.set(jobClass, queue2, queue3, 1.0);
                
                // Both branches to join
                P.set(jobClass, queue1, join, 1.0);
                P.set(jobClass, queue3, join, 1.0);
                
                // Join back to delay
                P.set(jobClass, join, delay1, 1.0);
                
            } else if (config == FJConfig.MC4) {
                // MC4: Fork to THREE branches
                P.set(jobClass, fork, queue2, 1.0);  // Branch 2
                P.set(jobClass, fork, queue4, 1.0);  // Branch 3
                
                // Branch 1: queue1 -> join
                P.set(jobClass, queue1, join, 1.0);
                
                // Branch 2: queue2 -> join
                P.set(jobClass, queue2, join, 1.0);
                
                // Branch 3: queue4 -> delay2 -> join
                P.set(jobClass, queue4, delay2, 1.0);
                P.set(jobClass, delay2, join, 1.0);
                
                // Join -> queue3 -> delay1
                P.set(jobClass, join, queue3, 1.0);
                P.set(jobClass, queue3, delay1, 1.0);
            }
        }
        
        model.link(P);
        
        // Get simulation results for comparison
        SolverJMT simSolver = BenchmarkSolvers.getForkJoinSimulationSolver(model);
        SolverResult simResult = simSolver.getAvg();
        
        // Test fork-join solvers
        List<Solver> solvers = BenchmarkSolvers.getForkJoinBenchmarkSolvers(model);
        
        Map<String, Double> errQ = new HashMap<>();
        Map<String, Double> errU = new HashMap<>();
        Map<String, Double> errR = new HashMap<>();
        Map<String, Double> errT = new HashMap<>();
        
        for (Solver solver : solvers) {
            String solverName = solver.getName();
            if (solver.getOptions() != null && solver.getOptions().config != null) {
                solverName += "_" + solver.getOptions().config.fork_join;
            }
            
            try {
                SolverResult result = ((NetworkSolver) solver).getAvg();
                
                // Calculate errors
                double errorQ = BenchmarkUtils.maxErrorOnSum(result.QN, simResult.QN);
                double errorU = BenchmarkUtils.utilizationError(result.UN, simResult.UN, model);
                double errorR = BenchmarkUtils.mape(result.RN, simResult.RN);
                double errorT = BenchmarkUtils.mape(result.TN, simResult.TN);
                
                errQ.put(solverName, errorQ);
                errU.put(solverName, errorU);
                errR.put(solverName, errorR);
                errT.put(solverName, errorT);
                
            } catch (Exception e) {
                // Solver failed - record as max error
                errQ.put(solverName, Double.MAX_VALUE);
                errU.put(solverName, Double.MAX_VALUE);
                errR.put(solverName, Double.MAX_VALUE);
                errT.put(solverName, Double.MAX_VALUE);
            }
        }
        
        // Store results
        Map<String, Object> results = new HashMap<>();
        results.put("ERRQ", errQ);
        results.put("ERRU", errU);
        results.put("ERRR", errR);
        results.put("ERRT", errT);
        results.put("model", model);
        
        return results;
    }
}