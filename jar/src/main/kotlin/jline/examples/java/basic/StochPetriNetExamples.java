/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java.basic;

import jline.lang.Network;
import jline.solvers.ctmc.CTMC;
import jline.solvers.jmt.JMT;
import java.util.Scanner;

/**
 * Stochastic Petri net examples mirroring the Kotlin notebooks in stochPetriNet.
 * <p>
 * This class contains Java implementations that mirror the Kotlin notebook examples
 * found in jar/src/main/kotlin/jline/examples/kotlin/basic/stochPetriNet/. Each method 
 * demonstrates a specific Petri net concept using models from the basic package.
 * <p>
 * The examples cover:
 * - Basic open and closed Petri net structures
 * - Multiple firing modes and batch processing
 * - Inhibiting conditions and complex token flows
 * - Various stochastic distributions in transitions
 * - Multi-class token systems
 */
public class StochPetriNetExamples {

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
     * Basic closed stochastic Petri net (spn_basic_closed.ipynb).
     * <p>
     * Demonstrates a simple closed Petri net with cyclic token flow
     * between two places using multiple solvers.
     * <p>
     * Features:
     * - 3 tokens initially in Place1, cycling between Place1 and Place2
     * - Exponential firing rates for transitions T1 and T2
     * - Multiple solver comparison (CTMC and JMT)
     * - Token conservation analysis
     */
    public static void spn_basic_closed() throws Exception {
        Network model = StochPetriNetModel.spn_basic_closed();
        
        // CTMC solver (exact for small Petri nets)
//        try {
//            CTMC solverCtmc = new CTMC(model);
//            solverCtmc.getAvgTable().print();
//        } catch (Exception e) {
//            System.out.println("CTMC solver error: " + e.getMessage());
//        }

        // JMT solver (simulation)
        try {
            JMT solverJmt = new JMT(model, "seed", 23000);
            solverJmt.getAvgTable().print();
        } catch (Exception e) {
            System.out.println("JMT solver not available: " + e.getMessage());
        }
        pauseForUser();
    }
    
    /**
     * Basic open stochastic Petri net (spn_basic_open.ipynb).
     * <p>
     * Shows fundamental open Petri net structure with
     * source, place, transition, and sink.
     */
    public static void spn_basic_open() throws Exception {
        Network model = StochPetriNetModel.spn_basic_open();
        JMT solver = new JMT(model, "seed", 23000);
        
        solver.getAvgTable().print();
        pauseForUser();
    }
    
    /**
     * Batch processing Petri net (spn_twomodes.ipynb).
     * <p>
     * Demonstrates batch token processing where transitions
     * require and produce multiple tokens.
     */
    public static void spn_twomodes() throws Exception {
        Network model = StochPetriNetModel.spn_twomodes();
        JMT solver = new JMT(model, "seed", 23000);
        
        solver.getAvgTable().print();
        pauseForUser();
    }
    
    /**
     * Competing transitions Petri net (spn_fourmodes.ipynb).
     * <p>
     * Shows resource competition between transitions
     * requiring different numbers of tokens.
     */
    public static void spn_fourmodes() throws Exception {
        Network model = StochPetriNetModel.spn_fourmodes();
        JMT solver = new JMT(model, "seed", 23000);
        
        solver.getAvgTable().print();
        pauseForUser();
    }
    
    /**
     * Inhibiting transitions Petri net (spn_inhibiting.ipynb).
     * <p>
     * Demonstrates inhibitor arcs that prevent firing
     * when tokens are present in certain places.
     */
    public static void spn_inhibiting() throws Exception {
        Network model = StochPetriNetModel.spn_inhibiting();
        JMT solver = new JMT(model, "seed", 23000);
        
        solver.getAvgTable().print();
        pauseForUser();
    }
    
    /**
     * Closed Petri net with two places (spn_closed_twoplaces.ipynb).
     * <p>
     * Simple closed system demonstrating token circulation
     * between two places with different service rates.
     */
    public static void spn_closed_twoplaces() throws Exception {
        Network model = StochPetriNetModel.spn_closed_twoplaces();
        
        // CTMC solver for exact solution
        try {
            CTMC solverCtmc = new CTMC(model);
            solverCtmc.getAvgTable().print();
        } catch (Exception e) {
            System.out.println("CTMC solver error: " + e.getMessage());
        }
        
        // JMT solver for simulation
        try {
            JMT solverJmt = new JMT(model, "seed", 23000);
            solverJmt.getAvgTable().print();
        } catch (Exception e) {
            System.out.println("JMT solver not available: " + e.getMessage());
        }
        pauseForUser();
    }
    
    /**
     * Closed Petri net with four places (spn_closed_fourplaces.ipynb).
     * <p>
     * More complex closed system with tokens cycling through
     * four different places using multiple transitions.
     */
    public static void spn_closed_fourplaces() throws Exception {
        Network model = StochPetriNetModel.spn_closed_fourplaces();
        
        // CTMC solver for exact solution
        try {
            CTMC solverCtmc = new CTMC(model);
            solverCtmc.getAvgTable().print();
        } catch (Exception e) {
            System.out.println("CTMC solver error: " + e.getMessage());
        }
        
        // JMT solver for simulation
        try {
            JMT solverJmt = new JMT(model, "seed", 23000);
            solverJmt.getAvgTable().print();
        } catch (Exception e) {
            System.out.println("JMT solver not available: " + e.getMessage());
        }
        pauseForUser();
    }
    
    /**
     * Open Petri net with seven places (spn_open_sevenplaces.ipynb).
     * <p>
     * Complex open system demonstrating token flow through
     * multiple places with varied routing probabilities.
     */
    public static void spn_open_sevenplaces() throws Exception {
        Network model = StochPetriNetModel.spn_open_sevenplaces();
        JMT solver = new JMT(model, "seed", 23000);
        
        solver.getAvgTable().print();
        pauseForUser();
    }

    /**
     * Main method demonstrating selected Petri net examples.
     */
    public static void main(String[] args) throws Exception {
        System.out.println("\n=== Running example: spn_basic_closed ===");
        try {
            spn_basic_closed();
        } catch (Exception e) {
            System.err.println("spn_basic_closed failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: spn_basic_open ===");
        try {
            spn_basic_open();
        } catch (Exception e) {
            System.err.println("spn_basic_open failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: spn_twomodes ===");
        try {
            spn_twomodes();
        } catch (Exception e) {
            System.err.println("spn_twomodes failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: spn_fourmodes ===");
        try {
            spn_fourmodes();
        } catch (Exception e) {
            System.err.println("spn_fourmodes failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: spn_inhibiting ===");
        try {
            spn_inhibiting();
        } catch (Exception e) {
            System.err.println("spn_inhibiting failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: spn_closed_twoplaces ===");
        try {
            spn_closed_twoplaces();
        } catch (Exception e) {
            System.err.println("spn_closed_twoplaces failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: spn_closed_fourplaces ===");
        try {
            spn_closed_fourplaces();
        } catch (Exception e) {
            System.err.println("spn_closed_fourplaces failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: spn_open_sevenplaces ===");
        try {
            spn_open_sevenplaces();
        } catch (Exception e) {
            System.err.println("spn_open_sevenplaces failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        scanner.close();
    }
}