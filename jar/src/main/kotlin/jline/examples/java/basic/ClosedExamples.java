/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java.basic;

import jline.lang.Network;
import jline.lang.constant.SchedStrategy;
import jline.solvers.ctmc.CTMC;
import jline.solvers.fluid.FLD;
import jline.solvers.jmt.JMT;
import jline.solvers.mam.MAM;
import jline.solvers.mva.MVA;
import jline.solvers.nc.NC;
import jline.solvers.ssa.SSA;

import java.util.Scanner;

/**
 * Closed queueing network examples mirroring the Kotlin notebooks in closedQN.
 * <p>
 * This class contains Java implementations that mirror the Kotlin notebook examples
 * found in jar/src/main/kotlin/jline/examples/kotlin/basic/closedQN/. Each method 
 * demonstrates a specific closed queueing network concept using models from the basic package.
 * <p>
 * The examples cover:
 * - BCMP theorem validation with different scheduling strategies
 * - Multi-server systems and repairman models
 * - Various service distributions (Erlang, Hyperexponential, MMPP)
 * - Discriminatory processor sharing (DPS) scheduling
 * - Multi-class systems with complex routing
 * - One-line network specifications
 */
public class ClosedExamples {
    
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
     * BCMP theorem validation example (cqn_bcmp_theorem.ipynb).
     * <p>
     * Demonstrates the BCMP theorem by comparing different scheduling strategies
     * (PS, FCFS, LCFS-PR) and showing that certain performance measures remain invariant.
     * <p>
     * Features:
     * - Multiple job classes (2 classes, 2 jobs each) in closed network
     * - Comparison of PS, FCFS, and LCFS-PR scheduling strategies
     * - Validation of product-form solution properties
     * - Uses CTMC solver for exact analysis
     * 
     * @throws Exception if solver fails
     */
    public static void cqn_bcmp_theorem() throws Exception {
        // Test PS scheduling model
        Network modelPS = ClosedModel.cqn_bcmp_theorem_ps();
        try {
            CTMC solverPS = new CTMC(modelPS);
            solverPS.getAvgTable().print();
        } catch (Exception e) {
            System.out.println("Error solving PS model: " + e.getMessage());
        }
        
        // Test FCFS scheduling model
        Network modelFCFS = ClosedModel.cqn_bcmp_theorem_fcfs();
        try {
            CTMC solverFCFS = new CTMC(modelFCFS);
            solverFCFS.getAvgTable().print();
        } catch (Exception e) {
            System.out.println("Error solving FCFS model: " + e.getMessage());
        }
        
        // Test LCFS-PR scheduling model
        Network modelLCFSPR = ClosedModel.cqn_bcmp_theorem_lcfspr();
        CTMC solverLCFSPR = new CTMC(modelLCFSPR);
        solverLCFSPR.getAvgTable().print();
        
        pauseForUser();
    }
    
    /**
     * MMPP2 service process example (cqn_mmpp2_service.ipynb).
     * <p>
     * Demonstrates the use of Markov Modulated Poisson Process (MMPP)
     * service distributions in closed networks.
     * <p>
     * Features:
     * - MMPP2 service process with state-dependent rates
     * - Bursty service patterns with correlation
     * - Closed network with complex service behavior
     * - Performance analysis under correlated service
     * 
     * @throws Exception if solver fails
     */
    public static void cqn_mmpp2_service() throws Exception {
        Network model = ClosedModel.cqn_mmpp2_service();
        JMT solver = new JMT(model, "seed", 23000, "verbose", 1, "samples", 20000);
        
        solver.getAvgTable().print();
        
        pauseForUser();
    }
    
    /**
     * Multi-server closed network (cqn_multiserver.ipynb).
     * <p>
     * Demonstrates closed networks with multi-server stations,
     * showing capacity effects on performance.
     * <p>
     * Features:
     * - Multiple servers at queueing stations
     * - Capacity constraints in closed networks
     * - Server utilization and queueing effects
     * - Performance comparison with single-server systems
     * 
     * @throws Exception if solver fails
     */
    public static void cqn_multiserver() throws Exception {
        Network model = ClosedModel.cqn_multiserver();
        MVA solver = new MVA(model);
        
        solver.getAvgTable().print();
        
        pauseForUser();
    }
    
    /**
     * One-line closed network specification (cqn_oneline.ipynb).
     * <p>
     * Demonstrates compact specification of closed networks using
     * matrix-based constructors.
     * <p>
     * Features:
     * - Matrix-based closed network constructor
     * - N vector defines job populations per class
     * - S matrix defines service times at stations
     * - Compact notation for standard topologies
     * 
     * @throws Exception if solver fails
     */
    public static void cqn_oneline() throws Exception {
        Network model = ClosedModel.cqn_oneline();
        MVA solver = new MVA(model, "method", "exact");
        
        solver.getAvgTable().print();
        
        pauseForUser();
    }
    
    /**
     * Repairman model (cqn_repairmen.ipynb).
     * <p>
     * Demonstrates the classic repairman model where machines
     * break down and compete for repair resources.
     * <p>
     * Features:
     * - Machine population circulating between operational and repair
     * - Limited repair capacity creates resource contention
     * - Delay station represents operational (working) time
     * - Queue represents repair facility with limited servers
     * 
     * @throws Exception if solver fails
     */
    public static void cqn_repairmen() throws Exception {
        Network model = ClosedModel.cqn_repairmen();
        MVA solver = new MVA(model);
        
        solver.getAvgTable().print();
        
        pauseForUser();
    }
    
    /**
     * Multi-class repairman model (cqn_repairmen_multi.ipynb).
     * <p>
     * Extends the repairman model to multiple machine types
     * with different failure and repair characteristics.
     * <p>
     * Features:
     * - Multiple machine classes with different characteristics
     * - Class-dependent failure rates and repair times
     * - Shared repair facility for all machine types
     * - Resource contention across machine classes
     * 
     * @throws Exception if solver fails
     */
    public static void cqn_repairmen_multi() throws Exception {
        Network model = ClosedModel.cqn_repairmen_multi();
        MVA solver = new MVA(model);
        
        solver.getAvgTable().print();
        
        pauseForUser();
    }
    
    /**
     * Discriminatory processor sharing (cqn_scheduling_dps.ipynb).
     * <p>
     * Demonstrates discriminatory processor sharing (DPS) scheduling
     * where different job classes receive different service shares.
     * <p>
     * Features:
     * - DPS scheduling with class-dependent weights
     * - Proportional service allocation based on class priorities
     * - Performance differentiation across job classes
     * - Weighted fair queueing behavior
     * 
     * @throws Exception if solver fails
     */
    public static void cqn_scheduling_dps() throws Exception {
        Network model = ClosedModel.cqn_scheduling_dps();
        JMT solver = new JMT(model, "seed", 23000);
        
        solver.getAvgTable().print();
        
        pauseForUser();
    }
    
    /**
     * Three-class hyperexponential system (cqn_threeclass_hyperl.ipynb).
     * <p>
     * Demonstrates a three-class closed network with hyperexponential
     * service distributions, showing high variability effects.
     * <p>
     * Features:
     * - Three job classes with different populations
     * - Hyperexponential service distributions (high variability)
     * - Complex routing patterns across multiple stations
     * - Performance analysis under high service variability
     * 
     * @throws Exception if solver fails
     */
    public static void cqn_threeclass_hyperl() throws Exception {
        Network model = ClosedModel.cqn_threeclass_hyperl();
        MVA solver = new MVA(model);
        
        solver.getAvgTable().print();
        
        pauseForUser();
    }
    
    /**
     * Two-class Erlang system (cqn_twoclass_erl.ipynb).
     * <p>
     * Demonstrates a two-class closed network with Erlang service
     * distributions, showing low variability effects.
     * <p>
     * Features:
     * - Two job classes with different characteristics
     * - Erlang service distributions (low variability)
     * - Comparison with exponential service systems
     * - Impact of service variability on performance
     * 
     * @throws Exception if solver fails
     */
    public static void cqn_twoclass_erl() throws Exception {
        Network model = ClosedModel.cqn_twoclass_erl();
        MVA solver = new MVA(model);
        
        solver.getAvgTable().print();
        
        pauseForUser();
    }
    
    /**
     * Two-class hyperexponential system (cqn_twoclass_hyperl.ipynb).
     * <p>
     * Demonstrates a two-class closed network with hyperexponential
     * service distributions for both classes.
     * <p>
     * Features:
     * - Two job classes with hyperexponential service
     * - High service variability for all classes
     * - Performance comparison between classes
     * - Variability effects on class interactions
     * 
     * @throws Exception if solver fails
     */
    public static void cqn_twoclass_hyperl() throws Exception {
        Network model = ClosedModel.cqn_twoclass_hyperl();
        MVA solver = new MVA(model);
        
        solver.getAvgTable().print();
        
        pauseForUser();
    }
    
    /**
     * Two queues with reducible routing (cqn_twoqueues.ipynb).
     * <p>
     * Demonstrates a closed network with reducible routing matrix,
     * containing probabilistic routing patterns.
     * <p>
     * Features:
     * - Reducible routing matrix with probabilistic transitions
     * - Self-loops and probabilistic routing decisions
     * - Delay station with two queue alternatives
     * - Analysis of routing probability effects
     * 
     * @throws Exception if solver fails
     */
    public static void cqn_twoqueues() throws Exception {
        Network model = ClosedModel.cqn_twoqueues();
        MVA solver = new MVA(model);
        
        solver.getAvgTable().print();
        
        pauseForUser();
    }
    
    /**
     * Two queues with multi-class system (cqn_twoqueues_multi.ipynb).
     * <p>
     * Extends the two-queue model to multiple classes with
     * different routing patterns per class.
     * <p>
     * Features:
     * - Multiple job classes with class-specific routing
     * - Different probabilistic patterns per class
     * - Class interactions in shared queueing stations
     * - Performance differentiation across classes
     * 
     * @throws Exception if solver fails
     */
    public static void cqn_twoqueues_multi() throws Exception {
        Network model = ClosedModel.cqn_twoqueues_multi();
        MVA solver = new MVA(model);
        
        solver.getAvgTable().print();
        
        pauseForUser();
    }

    /**
     * Main method demonstrating all closed queueing network examples.
     * <p>
     * Executes all example methods to showcase the different closed network
     * concepts and solution approaches available in LINE.
     * 
     * @param args command line arguments (not used)
     * @throws Exception if any example fails
     */
    public static void main(String[] args) throws Exception {
        System.out.println("\n=== Running example: cqn_bcmp_theorem ===");
        try {
            cqn_bcmp_theorem();
        } catch (Exception e) {
            System.err.println("cqn_bcmp_theorem failed: " + e.getMessage());
            e.printStackTrace();
        }

        System.out.println("\n=== Running example: cqn_mmpp2_service ===");
        try {
            cqn_mmpp2_service();
        } catch (Exception e) {
            System.err.println("cqn_mmpp2_service failed: " + e.getMessage());
            e.printStackTrace();
        }

        System.out.println("\n=== Running example: cqn_multiserver ===");
        try {
            cqn_multiserver();
        } catch (Exception e) {
            System.err.println("cqn_multiserver failed: " + e.getMessage());
            e.printStackTrace();
        }

        System.out.println("\n=== Running example: cqn_oneline ===");
        try {
            cqn_oneline();
        } catch (Exception e) {
            System.err.println("cqn_oneline failed: " + e.getMessage());
            e.printStackTrace();
        }

        System.out.println("\n=== Running example: cqn_repairmen ===");
        try {
            cqn_repairmen();
        } catch (Exception e) {
            System.err.println("cqn_repairmen failed: " + e.getMessage());
            e.printStackTrace();
        }

        System.out.println("\n=== Running example: cqn_repairmen_multi ===");
        try {
            cqn_repairmen_multi();
        } catch (Exception e) {
            System.err.println("cqn_repairmen_multi failed: " + e.getMessage());
            e.printStackTrace();
        }

        System.out.println("\n=== Running example: cqn_scheduling_dps ===");
        try {
            cqn_scheduling_dps();
        } catch (Exception e) {
            System.err.println("cqn_scheduling_dps failed: " + e.getMessage());
            e.printStackTrace();
        }

        System.out.println("\n=== Running example: cqn_threeclass_hyperl ===");
        try {
            cqn_threeclass_hyperl();
        } catch (Exception e) {
            System.err.println("cqn_threeclass_hyperl failed: " + e.getMessage());
            e.printStackTrace();
        }

        System.out.println("\n=== Running example: cqn_twoclass_erl ===");
        try {
            cqn_twoclass_erl();
        } catch (Exception e) {
            System.err.println("cqn_twoclass_erl failed: " + e.getMessage());
            e.printStackTrace();
        }

        System.out.println("\n=== Running example: cqn_twoclass_hyperl ===");
        try {
            cqn_twoclass_hyperl();
        } catch (Exception e) {
            System.err.println("cqn_twoclass_hyperl failed: " + e.getMessage());
            e.printStackTrace();
        }

        System.out.println("\n=== Running example: cqn_twoqueues ===");
        try {
            cqn_twoqueues();
        } catch (Exception e) {
            System.err.println("cqn_twoqueues failed: " + e.getMessage());
            e.printStackTrace();
        }

        System.out.println("\n=== Running example: cqn_twoqueues_multi ===");
        try {
            cqn_twoqueues_multi();
        } catch (Exception e) {
            System.err.println("cqn_twoqueues_multi failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        scanner.close();
    }
}