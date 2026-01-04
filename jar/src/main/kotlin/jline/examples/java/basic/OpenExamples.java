/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java.basic;

import jline.VerboseLevel;
import jline.lang.Network;
import jline.solvers.ctmc.CTMC;
import jline.solvers.fluid.FLD;
import jline.solvers.jmt.JMT;
import jline.solvers.mam.MAM;
import jline.solvers.mva.MVA;
import jline.solvers.nc.NC;
import jline.solvers.ssa.SSA;
import java.util.Scanner;

/**
 * Open queueing network examples mirroring the Kotlin notebooks in openQN.
 * <p>
 * This class contains Java implementations that mirror the Kotlin notebook examples
 * found in jar/src/main/kotlin/jline/examples/kotlin/basic/openQN/. Each method 
 * demonstrates a specific open queueing network concept using models from the basic package.
 * <p>
 * The examples cover:
 * - Basic open networks with multiple solver comparisons
 * - Class switching and routing patterns
 * - Multi-class systems with complex topologies
 * - Trace-driven service and empirical distributions
 * - One-line network specifications
 * - Virtual sinks and probabilistic routing
 */
public class OpenExamples {

    private static final Scanner scanner = new Scanner(System.in);

    private static void pauseForUser() {
        // Skip pause if running in non-interactive mode (e.g., Maven exec)
        if (System.console() == null) {
            System.out.println("\n[Running in non-interactive mode, continuing...]");
            return;
        }
        System.out.println("\nPress Enter to continue to next example...");
        try {
            scanner.nextLine();
        } catch (Exception e) {
            // Ignore scanner errors in case of pipe or redirection
        }
    }

    /**
     * Basic open queueing network (oqn_basic.ipynb).
     * <p>
     * Demonstrates a basic open queueing network with hyperexponential service
     * time distribution and comparison of multiple solvers.
     * <p>
     * Features:
     * - Open network structure: jobs arrive from external source and depart through sink
     * - Hyperexponential service distribution at delay station
     * - Multiple solver comparison (CTMC, JMT, SSA, Fluid, MVA, NC)
     * - Performance analysis with steady-state metrics
     * 
     * @throws Exception if any solver fails
     */
    public static void oqn_basic() throws Exception {
        Network model = OpenModel.oqn_basic();
        
        // Create array of different solvers to compare
        CTMC solverCTMC = new CTMC(model, "keep", false, "cutoff", 10);
        JMT solverJMT = new JMT(model, "seed", 23000, "verbose", VerboseLevel.STD, "keep", false);
        SSA solverSSA = new SSA(model, "seed", 23000, "verbose", VerboseLevel.STD, "samples", 10000);
        FLD solverFluid = new FLD(model);
        MVA solverMVA = new MVA(model);
        NC solverNC = new NC(model);
        
        // Execute each solver and display results
        try {
            solverCTMC.getAvgTable().print();
        } catch (Exception e) {
            System.out.println("Solver failed: " + e.getMessage());
        }
        
        try {
            solverJMT.getAvgTable().print();
        } catch (Exception e) {
            System.out.println("Solver failed: " + e.getMessage());
        }
        
        try {
            solverSSA.getAvgTable().print();
        } catch (Exception e) {
            System.out.println("Solver failed: " + e.getMessage());
        }
        
        try {
            solverFluid.getAvgTable().print();
        } catch (Exception e) {
            System.out.println("Solver failed: " + e.getMessage());
        }
        
        try {
            solverMVA.getAvgTable().print();
        } catch (Exception e) {
            System.out.println("Solver failed: " + e.getMessage());
        }
        
        try {
            solverNC.getAvgTable().print();
        } catch (Exception e) {
            System.out.println("Solver failed: " + e.getMessage());
        }
        pauseForUser();
    }
    
    /**
     * Three-class open network with class switching (oqn_cs_routing.ipynb).
     * <p>
     * Demonstrates class switching behavior in open networks where jobs can
     * transform from one class to another during processing.
     * <p>
     * Features:
     * - Three open classes: Class A and B arrive, Class C created via switching
     * - Class switching from A→C and B→C at ClassSwitch node
     * - Two PS queues with different service rates per class
     * - Multiple solver comparison (CTMC, Fluid, MVA, MAM, NC, JMT, SSA)
     * 
     * @throws Exception if solver fails
     */
    public static void oqn_cs_routing() throws Exception {
        Network model = OpenModel.oqn_cs_routing();
        
        // Create multiple solvers as in notebook
        CTMC solverCTMC = new CTMC(model, "keep", true, "verbose", 1, "cutoff", new int[][]{{1,1,0}, {3,3,0}, {0,0,3}});
        FLD solverFluid = new FLD(model, "keep", true, "verbose", 1);
        MVA solverMVA = new MVA(model, "keep", true, "verbose", 1);
        MAM solverMAM = new MAM(model, "keep", true, "verbose", 1);
        NC solverNC = new NC(model, "keep", true, "verbose", 1);
        JMT solverJMT = new JMT(model, "keep", true, "verbose", 1, "seed", 23000, "samples", 100000);
        SSA solverSSA = new SSA(model, "keep", true, "verbose", 1, "seed", 23000, "samples", 100000);
        
        Object[] solvers = {solverCTMC, solverFluid, solverMVA, solverMAM, solverNC, solverJMT, solverSSA};
        
        // Execute all solvers and collect results
        for (int i = 0; i < solvers.length; i++) {
            try {
                if (solvers[i] instanceof CTMC) {
                    ((CTMC) solvers[i]).getAvgTable().print();
                } else if (solvers[i] instanceof FLD) {
                    ((FLD) solvers[i]).getAvgTable().print();
                } else if (solvers[i] instanceof MVA) {
                    ((MVA) solvers[i]).getAvgTable().print();
                } else if (solvers[i] instanceof MAM) {
                    ((MAM) solvers[i]).getAvgTable().print();
                } else if (solvers[i] instanceof NC) {
                    ((NC) solvers[i]).getAvgTable().print();
                } else if (solvers[i] instanceof JMT) {
                    ((JMT) solvers[i]).getAvgTable().print();
                } else if (solvers[i] instanceof SSA) {
                    ((SSA) solvers[i]).getAvgTable().print();
                }
            } catch (Exception e) {
                System.out.println("Solver failed: " + e.getMessage());
            }
        }
        pauseForUser();
    }
    
    /**
     * Complex multi-class open network with four queues (oqn_fourqueues.ipynb).
     * <p>
     * Demonstrates a web service architecture with multiple storage systems
     * and complex feedback routing patterns.
     * <p>
     * Features:
     * - Three open classes with different arrival rates and priorities
     * - Four queues: WebServer (FCFS), Storage1 (FCFS), Storage2 (PS), Storage3 (FCFS)
     * - Complex feedback routing with 25% probability splits
     * - Multiple solver comparison (CTMC, MVA, MAM, JMT) - matches MATLAB coverage
     * 
     * @throws Exception if solver fails
     */
    public static void oqn_fourqueues() throws Exception {
        Network model = OpenModel.oqn_fourqueues();
        
        // Create multiple solvers as in MATLAB version (some commented out in MATLAB)
        CTMC solverCTMC = new CTMC(model, "seed", 23000, "cutoff", 1);
        // FLD solverFluid = new FLD(model, "seed", 23000); // Commented out in MATLAB
        MVA solverMVA = new MVA(model, "seed", 23000);
        MAM solverMAM = new MAM(model, "seed", 23000);
        JMT solverJMT = new JMT(model, "seed", 23000, "samples", 1000000);
        // SSA solverSSA = new SSA(model, "seed", 23000); // Commented out in MATLAB
        // NC solverNC = new NC(model, "seed", 23000); // Commented out in MATLAB
        
        Object[] solvers = {solverCTMC, solverMVA, solverMAM, solverJMT};
        String[] solverNames = {"CTMC", "MVA", "MAM", "JMT"};
        
        // Execute all solvers and collect results
        for (int i = 0; i < solvers.length; i++) {
            try {
                if (solvers[i] instanceof CTMC) {
                    ((CTMC) solvers[i]).getAvgTable().print();
                } else if (solvers[i] instanceof MVA) {
                    ((MVA) solvers[i]).getAvgTable().print();
                } else if (solvers[i] instanceof MAM) {
                    ((MAM) solvers[i]).getAvgTable().print();
                } else if (solvers[i] instanceof JMT) {
                    ((JMT) solvers[i]).getAvgTable().print();
                }
            } catch (Exception e) {
                System.out.println("Solver failed: " + e.getMessage());
            }
        }
        pauseForUser();
    }
    
    /**
     * One-line tandem PS network specification (oqn_oneline.ipynb).
     * <p>
     * Demonstrates compact network specification using matrix-based constructors
     * for processor sharing queues with delays.
     * <p>
     * Features:
     * - Matrix-based constructor for PS queues with delay
     * - Lambda matrix defines arrival rates for multiple classes
     * - D matrix defines service demands at stations
     * - Z matrix defines service times at delay stations
     * 
     * @throws Exception if solver fails
     */
    public static void oqn_oneline() throws Exception {
        Network model = OpenModel.oqn_oneline();
        
        // Create solver as in MATLAB version
        MVA solver = new MVA(model);
        solver.getAvgTable().print();
        
        pauseForUser();
    }
    
    /**
     * Open network with trace-driven service (oqn_trace_driven.ipynb).
     * <p>
     * Demonstrates empirical service time distributions driven by trace files,
     * useful for modeling real workload patterns.
     * <p>
     * Features:
     * - Single open class with exponential arrivals
     * - Queue service driven by trace file (example_trace.txt)
     * - Empirical service time distributions from real data
     * - Simple Source → Queue → Sink topology
     * 
     * @throws Exception if solver fails
     */
    public static void oqn_trace_driven() throws Exception {
        Network model = OpenModel.oqn_trace_driven();
        JMT solver = new JMT(model, "seed", 23000);
        
        solver.getAvgTable().print();
        pauseForUser();
    }
    
    /**
     * Open network with virtual sinks (oqn_vsinks.ipynb).
     * <p>
     * Demonstrates probabilistic routing to multiple exit points using
     * virtual sinks for different departure streams.
     * <p>
     * Features:
     * - Two open classes with different routing patterns
     * - Class1: 60% to VSink1, 40% to VSink2
     * - Class2: 10% to VSink1, 90% to VSink2
     * - Router nodes as intermediate destinations
     * - Multiple exit points from the network
     * 
     * @throws Exception if solver fails
     */
    public static void oqn_vsinks() throws Exception {
        Network model = OpenModel.oqn_vsinks();
        
        // Create solvers as in MATLAB version
        MVA solverMVA = new MVA(model);
        MAM solverMAM = new MAM(model);
        NC solverNC = new NC(model);
        
        // MVA solver results
        solverMVA.getAvgTable().print();
        solverMVA.getAvgNodeTable().print();
        
        // MAM solver results
        solverMAM.getAvgTable().print();
        solverMAM.getAvgNodeTable().print();
        
        // NC solver results
        solverNC.getAvgTable().print();
        solverNC.getAvgNodeTable().print();
        
        pauseForUser();
    }

    /**
     * Main method demonstrating all open queueing network examples.
     * <p>
     * Executes all example methods to showcase the different open network
     * concepts and solution approaches available in LINE.
     * 
     * @param args command line arguments (not used)
     * @throws Exception if any example fails
     */
    public static void main(String[] args) throws Exception {
        System.out.println("\n=== Running example: oqn_basic ===");
        try {
            oqn_basic();
        } catch (Exception e) {
            System.err.println("oqn_basic failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: oqn_cs_routing ===");
        try {
            oqn_cs_routing();
        } catch (Exception e) {
            System.err.println("oqn_cs_routing failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: oqn_fourqueues ===");
        try {
            oqn_fourqueues();
        } catch (Exception e) {
            System.err.println("oqn_fourqueues failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: oqn_oneline ===");
        try {
            oqn_oneline();
        } catch (Exception e) {
            System.err.println("oqn_oneline failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: oqn_trace_driven ===");
        try {
            oqn_trace_driven();
        } catch (Exception e) {
            System.err.println("oqn_trace_driven failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: oqn_vsinks ===");
        try {
            oqn_vsinks();
        } catch (Exception e) {
            System.err.println("oqn_vsinks failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        scanner.close();
    }
}