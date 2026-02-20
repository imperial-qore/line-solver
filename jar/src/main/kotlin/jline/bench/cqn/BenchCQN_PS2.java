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

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * CQN PS2 Benchmarks
 */
public class BenchCQN_PS2 {
    
    /**
     * CQN PS2 Light Load benchmark
     */
    public static Map<String, Object> runLightLoad(int iteration) {
        return runBenchmark(iteration, 5, 7, 1, false, false);
    }
    
    /**
     * CQN PS2 Mid Load benchmark
     */
    public static Map<String, Object> runMidLoad(int iteration) {
        return runBenchmark(iteration, 5, 7, 1, false, false);
    }
    
    /**
     * CQN PS2 High Load benchmark
     */
    public static Map<String, Object> runHighLoad(int iteration) {
        return runBenchmark(iteration, 5, 7, 1, false, false);
    }
    
    /**
     * CQN PS2 High CV Light Load benchmark
     */
    public static Map<String, Object> runHicvLightLoad(int iteration) {
        return runBenchmark(iteration, 5, 7, 1, true, false);
    }
    
    /**
     * CQN PS2 High CV Mid Load benchmark
     */
    public static Map<String, Object> runHicvMidLoad(int iteration) {
        return runBenchmark(iteration, 5, 7, 1, true, false);
    }
    
    /**
     * CQN PS2 High CV High Load benchmark
     */
    public static Map<String, Object> runHicvHighLoad(int iteration) {
        return runBenchmark(iteration, 5, 7, 1, true, false);
    }
    
    /**
     * CQN PS2 Multiserver Light Load benchmark
     */
    public static Map<String, Object> runMultiserverLightLoad(int iteration) {
        return runBenchmark(iteration, 5, 7, 10, false, false);
    }
    
    /**
     * CQN PS2 Multiserver Mid Load benchmark
     */
    public static Map<String, Object> runMultiserverMidLoad(int iteration) {
        return runBenchmark(iteration, 5, 7, 10, false, false);
    }
    
    /**
     * CQN PS2 Multiserver High Load benchmark
     */
    public static Map<String, Object> runMultiserverHighLoad(int iteration) {
        return runBenchmark(iteration, 5, 7, 10, false, false);
    }
    
    /**
     * CQN PS2 Random Mapping Light Load benchmark
     */
    public static Map<String, Object> runRmLightLoad(int iteration) {
        // Note: For PS2 rm configurations, N1=2, N2=2
        return runBenchmark(iteration, 2, 2, 1, false, true);
    }
    
    /**
     * CQN PS2 Random Mapping Mid Load benchmark
     */
    public static Map<String, Object> runRmMidLoad(int iteration) {
        return runBenchmark(iteration, 2, 2, 1, false, true);
    }
    
    /**
     * CQN PS2 Random Mapping High Load benchmark
     */
    public static Map<String, Object> runRmHighLoad(int iteration) {
        return runBenchmark(iteration, 2, 2, 1, false, true);
    }
    
    /**
     * Core benchmark implementation for CQN PS2
     */
    private static Map<String, Object> runBenchmark(int iteration, int N1, int N2, 
                                                   int numServers, boolean highCV, boolean randomMapping) {
        return BenchCQNTemplate.runBenchmark("CQN_PS2", iteration, N1, N2, numServers, highCV, randomMapping);
    }
}