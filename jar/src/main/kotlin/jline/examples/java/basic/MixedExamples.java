/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java.basic;

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
 * Mixed queueing network examples mirroring the Kotlin notebooks in mixedQN.
 * <p>
 * This class contains Java implementations that mirror the Kotlin notebook examples
 * found in jar/src/main/kotlin/jline/examples/kotlin/basic/mixedQN/. Each method 
 * demonstrates a specific mixed queueing network concept using models from the basic package.
 * <p>
 * The examples cover:
 * - Basic mixed networks combining open and closed classes
 * - Multi-server systems with different scheduling strategies
 * - Performance comparison between FCFS and PS scheduling
 * - Large population systems and scalability analysis
 */
public class MixedExamples {

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
     * Basic mixed open/closed network (mqn_basic.ipynb).
     * <p>
     * Demonstrates a fundamental mixed network combining open arrivals
     * with a circulating closed population, comparing multiple solvers.
     * <p>
     * Features:
     * - One closed class (2 jobs) and one open class
     * - Delay node with Erlang service for closed, HyperExp for open
     * - PS queue shared by both classes with different service distributions
     * - Multiple solver comparison (CTMC, JMT, SSA, Fluid, MVA, NC, MAM)
     * 
     * @throws Exception if solver fails
     */
    public static void mqn_basic() throws Exception {
        Network model = MixedModel.mqn_basic();
        
        // Create multiple solvers for comparison as in the notebook
        CTMC solverCTMC = new CTMC(model, "cutoff", 3);
        JMT solverJMT = new JMT(model, "seed", 23000, "verbose", true, "keep", false);
        SSA solverSSA = new SSA(model, "seed", 23000, "verbose", false, "samples", 10000);
        FLD solverFluid = new FLD(model);
        MVA solverMVA = new MVA(model);
        NC solverNC = new NC(model);
        MAM solverMAM = new MAM(model);
        
        // Solve and display results for each solver
        Object[] solvers = {solverCTMC, solverJMT, solverSSA, solverFluid, solverMVA, solverNC, solverMAM};
        
        for (Object solverObj : solvers) {
            try {
                if (solverObj instanceof CTMC) {
                    CTMC solver = (CTMC) solverObj;
                    solver.getAvgTable().print();
                } else if (solverObj instanceof JMT) {
                    JMT solver = (JMT) solverObj;
                    solver.getAvgTable().print();
                } else if (solverObj instanceof SSA) {
                    SSA solver = (SSA) solverObj;
                    solver.getAvgTable().print();
                } else if (solverObj instanceof FLD) {
                    FLD solver = (FLD) solverObj;
                    solver.getAvgTable().print();
                } else if (solverObj instanceof MVA) {
                    MVA solver = (MVA) solverObj;
                    solver.getAvgTable().print();
                } else if (solverObj instanceof NC) {
                    NC solver = (NC) solverObj;
                    solver.getAvgTable().print();
                } else if (solverObj instanceof MAM) {
                    MAM solver = (MAM) solverObj;
                    solver.getAvgTable().print();
                }
            } catch (Exception e) {
                System.out.println("Error with solver: " + e.getMessage());
            }
        }
        pauseForUser();
    }
    
    /**
     * Mixed network with FCFS multi-server queues (mqn_multiserver_fcfs.ipynb).
     * <p>
     * Demonstrates mixed networks with multiple servers using FCFS scheduling,
     * showing capacity effects on mixed traffic.
     * <p>
     * Features:
     * - 5 multi-server FCFS queues with different capacities (1-5 servers)
     * - Closed class (3 jobs) and open class with separate routing
     * - Multiple solver comparison (CTMC, JMT, SSA, MVA, MAM)
     * - Different service rates for each class
     * 
     * @throws Exception if solver fails
     */
    public static void mqn_multiserver_fcfs() throws Exception {
        Network model = MixedModel.mqn_multiserver_fcfs();
        
        // Create solvers with appropriate configurations as in notebook
        try {
            CTMC solverCTMC = new CTMC(model, "cutoff", 3);
            solverCTMC.getAvgTable().print();
        } catch (Exception e) {
            System.out.println("CTMC solver not available: " + e.getMessage());
        }
        
        try {
            JMT solverJMT = new JMT(model, "samples", 100000, "seed", 23000);
            solverJMT.getAvgTable().print();
        } catch (Exception e) {
            System.out.println("JMT solver not available: " + e.getMessage());
        }
        
        try {
            SSA solverSSA = new SSA(model, "cutoff", 3, "seed", 23000);
            solverSSA.getAvgTable().print();
        } catch (Exception e) {
            System.out.println("SSA solver not available: " + e.getMessage());
        }
        
        try {
            MVA solverMVA = new MVA(model);
            solverMVA.getAvgTable().print();
        } catch (Exception e) {
            System.out.println("MVA solver not available: " + e.getMessage());
        }
        
        try {
            MAM solverMAM = new MAM(model, "cutoff", 3, "seed", 23000, "keep", false, "verbose", true);
            solverMAM.getAvgTable().print();
        } catch (Exception e) {
            System.out.println("MAM solver not available: " + e.getMessage());
        }
        pauseForUser();
    }
    
    /**
     * Mixed network with PS multi-server queues (mqn_multiserver_ps.ipynb).
     * <p>
     * Demonstrates mixed networks with processor sharing scheduling
     * across multiple servers (1-5 servers per queue).
     * <p>
     * Features:
     * - 5 multi-server PS queues with different capacities
     * - Closed class (3 jobs) and open class with separate routing
     * - Different service rates optimized for server counts
     * - Multiple solver comparison (CTMC, JMT, MVA, MAM)
     * 
     * @throws Exception if solver fails
     */
    public static void mqn_multiserver_ps() throws Exception {
        Network model = MixedModel.mqn_multiserver_ps();
        
        // Create solvers with appropriate configurations as in notebook
        try {
            CTMC solverCTMC = new CTMC(model, "cutoff", 3);
            solverCTMC.getAvgTable().print();
        } catch (Exception e) {
            System.out.println("CTMC solver not available: " + e.getMessage());
        }
        
        try {
            JMT solverJMT = new JMT(model, "samples", 100000, "seed", 23000);
            solverJMT.getAvgTable().print();
        } catch (Exception e) {
            System.out.println("JMT solver not available: " + e.getMessage());
        }
        
        try {
            MVA solverMVA = new MVA(model, "method", "exact");
            solverMVA.getAvgTable().print();
        } catch (Exception e) {
            System.out.println("MVA solver not available: " + e.getMessage());
        }
        
        try {
            MAM solverMAM = new MAM(model, "seed", 23000, "keep", false, "verbose", true);
            solverMAM.getAvgTable().print();
        } catch (Exception e) {
            System.out.println("MAM solver not available: " + e.getMessage());
        }
        pauseForUser();
    }
    
    /**
     * Mixed FCFS network with large closed population (mqn_singleserver_fcfs.ipynb).
     * <p>
     * Demonstrates scalability with large closed populations and
     * high-variability open arrivals under FCFS scheduling.
     * <p>
     * Features:
     * - Large closed population (100 jobs) with 4 FCFS queues
     * - High-variability APH arrivals (mean=3, SCV=64) for open class
     * - CTMC/SSA avoided due to infinite state space with high population
     * - Multiple solver comparison (JMT, MVA, NC, MAM)
     * 
     * @throws Exception if solver fails
     */
    public static void mqn_singleserver_fcfs() throws Exception {
        Network model = MixedModel.mqn_singleserver_fcfs();
        
        // Create solvers as in notebook - CTMC/SSA omitted due to high population
        try {
            JMT solverJMT = new JMT(model, "cutoff", 3, "keep", false, "verbose", true, "seed", 23000, "samples", 20000);
            solverJMT.getAvgTable().print();
        } catch (Exception e) {
            System.out.println("JMT solver not available: " + e.getMessage());
        }
        
        try {
            MVA solverMVA = new MVA(model, "method", "lin", "keep", false, "verbose", true, "seed", 23000);
            solverMVA.getAvgTable().print();
        } catch (Exception e) {
            System.out.println("MVA solver not available: " + e.getMessage());
        }
        
        try {
            NC solverNC = new NC(model, "keep", false, "verbose", true, "seed", 23000);
            solverNC.getAvgTable().print();
        } catch (Exception e) {
            System.out.println("NC solver not available: " + e.getMessage());
        }
        
        try {
            MAM solverMAM = new MAM(model, "keep", false, "verbose", true, "seed", 23000);
            solverMAM.getAvgTable().print();
        } catch (Exception e) {
            System.out.println("MAM solver not available: " + e.getMessage());
        }
        pauseForUser();
    }
    
    /**
     * Mixed PS network with simplified routing (mqn_singleserver_ps.ipynb).
     * <p>
     * Demonstrates processor sharing in mixed networks with
     * large closed populations and simplified open routing.
     * <p>
     * Features:
     * - Large closed population (100 jobs) with 4 PS queues
     * - Exponential arrivals (rate=0.3) for open class
     * - PS scheduling provides fairness under high congestion
     * - Multiple solver comparison (JMT, SSA, Fluid, MVA, NC, MAM)
     * 
     * @throws Exception if solver fails
     */
    public static void mqn_singleserver_ps() throws Exception {
        Network model = MixedModel.mqn_singleserver_ps();
        
        // Create solvers as in notebook - comprehensive comparison
        try {
            JMT solverJMT = new JMT(model, "cutoff", 3, "keep", false, "verbose", true, "seed", 23000, "samples", 100000);
            solverJMT.getAvgTable().print();
        } catch (Exception e) {
            System.out.println("JMT solver not available: " + e.getMessage());
        }
        
        try {
            SSA solverSSA = new SSA(model, "seed", 23000, "verbose", false, "samples", 100000);
            solverSSA.getAvgTable().print();
        } catch (Exception e) {
            System.out.println("SSA solver not available: " + e.getMessage());
        }
        
        try {
            FLD solverFluid = new FLD(model, "cutoff", 3, "keep", false, "verbose", true, "seed", 23000, "samples", 100000);
            solverFluid.getAvgTable().print();
        } catch (Exception e) {
            System.out.println("Fluid solver not available: " + e.getMessage());
        }
        
        try {
            MVA solverMVA = new MVA(model, "method", "lin");
            solverMVA.getAvgTable().print();
        } catch (Exception e) {
            System.out.println("MVA solver not available: " + e.getMessage());
        }
        
        try {
            NC solverNC = new NC(model, "cutoff", 3, "keep", false, "verbose", true, "seed", 23000, "samples", 100000);
            solverNC.getAvgTable().print();
        } catch (Exception e) {
            System.out.println("NC solver not available: " + e.getMessage());
        }
        
        try {
            MAM solverMAM = new MAM(model, "cutoff", 3, "keep", false, "verbose", true, "seed", 23000, "samples", 100000);
            solverMAM.getAvgTable().print();
        } catch (Exception e) {
            System.out.println("MAM solver not available: " + e.getMessage());
        }
        pauseForUser();
    }

    /**
     * Main method demonstrating all mixed queueing network examples.
     * <p>
     * Executes all example methods to showcase the different mixed network
     * concepts and scheduling comparisons available in LINE.
     * 
     * @param args command line arguments (not used)
     * @throws Exception if any example fails
     */
    public static void main(String[] args) throws Exception {
        System.out.println("\n=== Running example: mqn_basic ===");
        try {
            mqn_basic();
        } catch (Exception e) {
            System.err.println("mqn_basic failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: mqn_multiserver_ps ===");
        try {
            mqn_multiserver_ps();
        } catch (Exception e) {
            System.err.println("mqn_multiserver_ps failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: mqn_multiserver_fcfs ===");
        try {
            mqn_multiserver_fcfs();
        } catch (Exception e) {
            System.err.println("mqn_multiserver_fcfs failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: mqn_singleserver_ps ===");
        try {
            mqn_singleserver_ps();
        } catch (Exception e) {
            System.err.println("mqn_singleserver_ps failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        scanner.close();
    }
}
