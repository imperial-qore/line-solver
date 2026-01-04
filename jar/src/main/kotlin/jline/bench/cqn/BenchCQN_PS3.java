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
 * CQN PS3 Benchmarks
 */
public class BenchCQN_PS3 {
    
    /**
     * CQN PS3 Light Load benchmark
     */
    public static Map<String, Object> runLightLoad(int iteration) {
        return runBenchmark(iteration, 3, 9, 1, false, false);
    }
    
    /**
     * CQN PS3 Mid Load benchmark
     */
    public static Map<String, Object> runMidLoad(int iteration) {
        return runBenchmark(iteration, 3, 9, 1, false, false);
    }
    
    /**
     * CQN PS3 High Load benchmark
     */
    public static Map<String, Object> runHighLoad(int iteration) {
        return runBenchmark(iteration, 3, 9, 1, false, false);
    }
    
    /**
     * CQN PS3 High CV Light Load benchmark
     */
    public static Map<String, Object> runHicvLightLoad(int iteration) {
        return runBenchmark(iteration, 3, 9, 1, true, false);
    }
    
    /**
     * CQN PS3 High CV Mid Load benchmark
     */
    public static Map<String, Object> runHicvMidLoad(int iteration) {
        return runBenchmark(iteration, 3, 9, 1, true, false);
    }
    
    /**
     * CQN PS3 High CV High Load benchmark
     */
    public static Map<String, Object> runHicvHighLoad(int iteration) {
        return runBenchmark(iteration, 3, 9, 1, true, false);
    }
    
    /**
     * CQN PS3 Multiserver Light Load benchmark
     */
    public static Map<String, Object> runMultiserverLightLoad(int iteration) {
        return runBenchmark(iteration, 3, 9, 10, false, false);
    }
    
    /**
     * CQN PS3 Multiserver Mid Load benchmark
     */
    public static Map<String, Object> runMultiserverMidLoad(int iteration) {
        return runBenchmark(iteration, 3, 9, 10, false, false);
    }
    
    /**
     * CQN PS3 Multiserver High Load benchmark
     */
    public static Map<String, Object> runMultiserverHighLoad(int iteration) {
        return runBenchmark(iteration, 3, 9, 10, false, false);
    }
    
    /**
     * CQN PS3 Random Mapping Light Load benchmark
     */
    public static Map<String, Object> runRmLightLoad(int iteration) {
        return runBenchmark(iteration, 3, 9, 1, false, true);
    }
    
    /**
     * CQN PS3 Random Mapping Mid Load benchmark
     */
    public static Map<String, Object> runRmMidLoad(int iteration) {
        return runBenchmark(iteration, 3, 9, 1, false, true);
    }
    
    /**
     * CQN PS3 Random Mapping High Load benchmark
     */
    public static Map<String, Object> runRmHighLoad(int iteration) {
        return runBenchmark(iteration, 3, 9, 1, false, true);
    }
    
    /**
     * Core benchmark implementation for CQN PS3
     */
    private static Map<String, Object> runBenchmark(int iteration, int N1, int N2, 
                                                   int numServers, boolean highCV, boolean randomMapping) {
        return BenchCQNTemplate.runBenchmark("CQN_PS3", iteration, N1, N2, numServers, highCV, randomMapping);
    }
}