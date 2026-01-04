/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java.basic;

import jline.lang.Network;
import jline.solvers.ctmc.CTMC;
import jline.solvers.jmt.JMT;
import jline.solvers.mva.MVA;
import jline.solvers.ssa.SSA;
import java.util.Scanner;

/**
 * Priority queueing examples mirroring the Kotlin notebooks in prioModel.
 * <p>
 * This class contains Java implementations that mirror the Kotlin notebook examples
 * found in jar/src/main/kotlin/jline/examples/kotlin/basic/prioModel/. Each method 
 * demonstrates a specific priority queueing concept using models from the basic package.
 * <p>
 * The examples cover:
 * - Head-of-line (HOL) priority scheduling in open and closed networks
 * - Processor sharing with priority (PS-PRIO) disciplines
 * - Performance comparison across different priority classes
 * - Priority effects with identical service requirements
 */
public class PrioExamples {

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
     * Closed priority network with HOL scheduling (prio_hol_closed.ipynb).
     * <p>
     * Demonstrates head-of-line priority scheduling in closed networks
     * with multiple job classes having different priorities.
     */
    public static void prio_hol_closed() throws Exception {
        Network model = PrioModel.prio_hol_closed();
        
        // Create solvers as in MATLAB version (lines 82, 84)
        MVA solverMVA = new MVA(model, "seed", 23000, "cutoff", 1, "samples", 10000);
        JMT solverJMT = new JMT(model, "seed", 23000, "cutoff", 1, "samples", 10000);
        
        solverMVA.getAvgTable().print();
        solverJMT.getAvgTable().print();
        
        pauseForUser();
    }
    
    /**
     * Open priority network with HOL scheduling (prio_hol_open.ipynb).
     * <p>
     * Shows head-of-line priority effects in open networks
     * with different arrival rates per priority class.
     */
    public static void prio_hol_open() throws Exception {
        Network model = PrioModel.prio_hol_open();
        
        // Create solvers as in MATLAB version (lines 69, 71, 73, 74)
        CTMC solverCTMC = new CTMC(model, "seed", 23000, "cutoff", 1, "samples", 10000);
        MVA solverMVA = new MVA(model, "seed", 23000, "cutoff", 1, "samples", 10000);
        JMT solverJMT = new JMT(model, "seed", 23000, "cutoff", 1, "samples", 10000);
        SSA solverSSA = new SSA(model, "seed", 23000, "cutoff", 1, "samples", 10000);
        
        solverCTMC.getAvgTable().print();
        solverMVA.getAvgTable().print();
        solverJMT.getAvgTable().print();
        solverSSA.getAvgTable().print();
        
        pauseForUser();
    }
    
    /**
     * Priority network with GPSPRIO scheduling (prio_identical.ipynb).
     * <p>
     * Demonstrates GPS (Generalized Processor Sharing) with priorities
     * using weighted fair sharing across multiple job classes.
     */
    public static void prio_identical() throws Exception {
        Network model = PrioModel.prio_identical();
        
        // Create solvers as in MATLAB version (lines 35, 36)
        CTMC solverCTMC = new CTMC(model);
        JMT solverJMT = new JMT(model, "seed", 23000, "verbose", true, "samples", 30000);
        
        solverCTMC.getAvgTable().print();
        solverJMT.getAvgTable().print();
        
        pauseForUser();
    }
    
    /**
     * Processor sharing with priority (prio_psprio.ipynb).
     * <p>
     * Shows processor sharing combined with priority-based
     * resource allocation across job classes.
     */
    public static void prio_psprio() throws Exception {
        Network model = PrioModel.prio_psprio();
        
        // Create solvers as in MATLAB version (lines 27, 28, 29)
        CTMC solverCTMC = new CTMC(model);
        JMT solverJMT = new JMT(model, "seed", 23000, "verbose", true, "samples", 5000);
        SSA solverSSA = new SSA(model, "seed", 23000, "verbose", true, "samples", 5000);
        
        solverCTMC.getAvgTable().print();
        solverJMT.getAvgTable().print();
        solverSSA.getAvgTable().print();
        
        pauseForUser();
    }

    /**
     * Main method demonstrating all priority queueing examples.
     */
    public static void main(String[] args) throws Exception {
        try {
            prio_hol_open();
        } catch (Exception e) {
            System.err.println("prio_hol_open failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        try {
            prio_hol_closed();
        } catch (Exception e) {
            System.err.println("prio_hol_closed failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        try {
            prio_identical();
        } catch (Exception e) {
            System.err.println("prio_identical failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        try {
            prio_psprio();
        } catch (Exception e) {
            System.err.println("prio_psprio failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        scanner.close();
    }
}