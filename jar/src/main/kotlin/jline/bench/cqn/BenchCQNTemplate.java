/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.bench.cqn;

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
import jline.lang.nodes.Queue;
import jline.solvers.SolverResult;
import jline.solvers.Solver;
import jline.solvers.NetworkSolver;
import jline.solvers.jmt.SolverJMT;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import jline.bench.cqn.CQNResultsFormatter;

/**
 * Template for CQN benchmark implementations
 */
public class BenchCQNTemplate {
    
    /**
     * Core benchmark implementation template
     */
    public static Map<String, Object> runBenchmark(String benchmarkName, int iteration, int N1, int N2, 
                                                   int numServers, boolean highCV, boolean randomMapping) {
        // Default to PS scheduling if not specified in benchmark name
        SchedStrategy sched = benchmarkName.contains("FCFS") ? SchedStrategy.FCFS : SchedStrategy.PS;
        return runBenchmark(benchmarkName, iteration, N1, N2, numServers, highCV, randomMapping, sched);
    }
    
    /**
     * Core benchmark implementation template with explicit scheduling strategy
     */
    public static Map<String, Object> runBenchmark(String benchmarkName, int iteration, int N1, int N2, 
                                                   int numServers, boolean highCV, boolean randomMapping,
                                                   SchedStrategy sched) {
        
        // Create model
        Network model = new Network(benchmarkName + "_" + iteration);
        
        // Create stations
        Delay station1 = new Delay(model, "Delay");
        Queue station2, station3 = null;
        
        if (randomMapping) {
            // Random mapping: only 2 stations total
            station2 = new Queue(model, "Queue1", sched);
        } else {
            // Standard: 3 stations total  
            station2 = new Queue(model, "Queue1", sched);
            station3 = new Queue(model, "Queue2", sched);
        }
        
        // Create job classes
        ClosedClass jobclass1 = new ClosedClass(model, "ClassA", N1, station1);
        ClosedClass jobclass2 = new ClosedClass(model, "ClassB", N2, station1);
        
        // Generate service rates
        int numStations = randomMapping ? 2 : 3;
        Matrix rate = BenchmarkUtils.randGallery(numStations, 2, iteration);
        
        // Set service distributions
        if (highCV) {
            // High coefficient of variation (SCV = 10)
            station1.setService(jobclass1, Cox2.fitMeanAndSCV(1.0 / rate.get(0, 0), 10.0));
            station1.setService(jobclass2, Cox2.fitMeanAndSCV(1.0 / rate.get(0, 1), 10.0));
            
            station2.setService(jobclass1, Cox2.fitMeanAndSCV(1.0 / rate.get(1, 0), 10.0));
            station2.setService(jobclass2, Cox2.fitMeanAndSCV(1.0 / rate.get(1, 1), 10.0));
            
            if (station3 != null) {
                station3.setService(jobclass1, Cox2.fitMeanAndSCV(1.0 / rate.get(2, 0), 10.0));
                station3.setService(jobclass2, Cox2.fitMeanAndSCV(1.0 / rate.get(2, 1), 10.0));
            }
        } else {
            // Exponential distribution
            station1.setService(jobclass1, new Exp(rate.get(0, 0)));
            station1.setService(jobclass2, new Exp(rate.get(0, 1)));
            
            station2.setService(jobclass1, new Exp(rate.get(1, 0)));
            station2.setService(jobclass2, new Exp(rate.get(1, 1)));
            
            if (station3 != null) {
                station3.setService(jobclass1, new Exp(rate.get(2, 0)));
                station3.setService(jobclass2, new Exp(rate.get(2, 1)));
            }
        }
        
        // Set number of servers
        station2.setNumberOfServers(numServers);
        if (station3 != null) {
            station3.setNumberOfServers(numServers);
        }
        
        // Set routing
        RoutingMatrix P = model.initRoutingMatrix();
        
        if (randomMapping) {
            // Circular routing for 2-station random mapping
            P.set(jobclass1, station1, station2, 1.0);
            P.set(jobclass1, station2, station1, 1.0);
            
            P.set(jobclass2, station1, station2, 1.0);
            P.set(jobclass2, station2, station1, 1.0);
        } else {
            // Standard 3-station routing
            // Class 1 routing: Delay -> 60% Queue1, 40% Queue2 -> back to Delay
            P.set(jobclass1, station1, station2, 0.6);
            P.set(jobclass1, station1, station3, 0.4);
            P.set(jobclass1, station2, station1, 1.0);
            P.set(jobclass1, station3, station1, 1.0);
            
            // Class 2 routing: Delay -> 100% Queue1 -> back to Delay
            P.set(jobclass2, station1, station2, 1.0);
            P.set(jobclass2, station2, station1, 1.0);
            P.set(jobclass2, station3, station1, 1.0);
        }
        
        model.link(P);
        
        // Get simulation results for comparison
        SolverJMT simSolver = BenchmarkSolvers.getSimulationSolver(model);
        SolverResult simResult = simSolver.getAvg();
        
        // Test all available solvers
        List<Solver> solvers = BenchmarkSolvers.getBenchmarkSolvers(model);
        
        Map<String, Double> errQ = new HashMap<>();
        Map<String, Double> errU = new HashMap<>();
        Map<String, Double> errR = new HashMap<>();
        Map<String, Double> errT = new HashMap<>();
        
        for (Solver solver : solvers) {
            String solverName = solver.getClass().getSimpleName();
            
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
        
        // Format and add to global results accumulator
        CQNResultsFormatter.addGlobalResult(benchmarkName, results);
        
        return results;
    }
}