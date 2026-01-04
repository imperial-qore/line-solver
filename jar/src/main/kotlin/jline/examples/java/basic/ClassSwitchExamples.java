/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java.basic;

import jline.lang.Network;
import jline.solvers.mva.MVA;
import java.util.Scanner;

/**
 * Class switching examples mirroring the Kotlin notebooks in classSwitching.
 */
public class ClassSwitchExamples {

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

    public static void cs_implicit() throws Exception {
        Network model = ClassSwitchingModel.cs_implicit();
        MVA solver = new MVA(model);
        
        solver.getAvgChainTable().print();
        pauseForUser();
    }
    
    public static void cs_multi_diamond() throws Exception {
        Network model = ClassSwitchingModel.cs_multi_diamond();
        MVA solver = new MVA(model);
        
        solver.getAvgChainTable().print();
        pauseForUser();
    }
    
    public static void cs_single_diamond() throws Exception {
        Network model = ClassSwitchingModel.cs_single_diamond();
        MVA solver = new MVA(model);
        
        solver.getAvgChainTable().print();
        pauseForUser();
    }
    
    public static void cs_transient_class() throws Exception {
        Network model = ClassSwitchingModel.cs_transient_class();
        MVA solver = new MVA(model);
        
        solver.getAvgChainTable().print();
        pauseForUser();
    }

    public static void main(String[] args) throws Exception {
        System.out.println("\n=== Running example: cs_implicit ===");
        try {
            cs_implicit();
        } catch (Exception e) {
            System.err.println("cs_implicit failed: " + e.getMessage());
            e.printStackTrace();
        }
        System.out.println("\n=== Running example: cs_single_diamond ===");
        try {
            cs_single_diamond();
        } catch (Exception e) {
            System.err.println("cs_single_diamond failed: " + e.getMessage());
            e.printStackTrace();
        }
        System.out.println("\n=== Running example: cs_multi_diamond ===");
        try {
            cs_multi_diamond();
        } catch (Exception e) {
            System.err.println("cs_multi_diamond failed: " + e.getMessage());
            e.printStackTrace();
        }
        System.out.println("\n=== Running example: cs_transient_class ===");
        try {
            cs_transient_class();
        } catch (Exception e) {
            System.err.println("cs_transient_class failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        scanner.close();
    }
}
