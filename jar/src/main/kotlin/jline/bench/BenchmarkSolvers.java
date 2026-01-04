/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.bench;

import jline.lang.Network;
import jline.VerboseLevel;
import jline.solvers.Solver;
import jline.solvers.SolverOptions;
import jline.solvers.auto.SolverAUTO;
import jline.solvers.fluid.SolverFluid;
import jline.solvers.jmt.SolverJMT;
import jline.solvers.mam.SolverMAM;
import jline.solvers.mva.SolverMVA;
import jline.solvers.nc.SolverNC;
import jline.solvers.qns.SolverQNS;

import java.util.ArrayList;
import java.util.List;

/**
 * Benchmark solver initialization utility
 */
public class BenchmarkSolvers {
    
    public static List<Solver> getBenchmarkSolvers(Network model) {
        List<Solver> solvers = new ArrayList<>();

        SolverOptions options = new SolverOptions();
        options.verbose = VerboseLevel.SILENT;
        options.seed = 23000;

        solvers.add(new SolverFluid(model, options));
        solvers.add(new SolverMVA(model, options));
        solvers.add(new SolverNC(model, options));
        solvers.add(new SolverAUTO(model, options));
        solvers.add(new SolverQNS(model, options));
        solvers.add(new SolverMAM(model, options));

        return solvers;
    }
    
    public static SolverJMT getSimulationSolver(Network model) {
        SolverOptions options = new SolverOptions();
        options.verbose = VerboseLevel.SILENT;
        options.seed = 23000;
        options.samples = 1000000; // 1e6
        options.method = "jmt.comom";
        
        return new SolverJMT(model, options);
    }
    
    /**
     * Get Fork-Join specific benchmark solvers
     * Tests different fork-join approximation methods
     */
    public static List<Solver> getForkJoinBenchmarkSolvers(Network model) {
        List<Solver> solvers = new ArrayList<>();
        
        // MVA with heidelberger-trivedi method
        SolverOptions mvaOptions1 = new SolverOptions();
        mvaOptions1.verbose = VerboseLevel.SILENT;
        mvaOptions1.seed = 23000;
        mvaOptions1.config.fork_join = "heidelberger-trivedi";
        solvers.add(new SolverMVA(model, mvaOptions1));
        
        // MVA with mmt method
        SolverOptions mvaOptions2 = new SolverOptions();
        mvaOptions2.verbose = VerboseLevel.SILENT;
        mvaOptions2.seed = 23000;
        mvaOptions2.config.fork_join = "mmt";
        solvers.add(new SolverMVA(model, mvaOptions2));
        
        return solvers;
    }
    
    /**
     * Get Fork-Join simulation solver with appropriate settings
     */
    public static SolverJMT getForkJoinSimulationSolver(Network model) {
        SolverOptions options = new SolverOptions();
        options.verbose = VerboseLevel.SILENT;
        options.seed = 23000;
        options.samples = 1000000; // 1e6
        options.method = "jmt.comom";
        
        return new SolverJMT(model, options);
    }
}