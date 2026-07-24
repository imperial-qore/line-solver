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
        model.printRoutingMatrix();

        jline.lang.NetworkStruct sn = model.getStruct(true);
        System.out.println("\n=== rt (station routing matrix) ===");
        sn.rt.print();

        // Debug SCC
        int K = sn.nclasses;
        int M = sn.nstateful;
        jline.util.matrix.Matrix inchain = sn.inchain.get(0);
        java.util.List<Integer> cols = new java.util.ArrayList<>();
        for (int ist = 0; ist < M; ist++) {
            for (int ik = 0; ik < inchain.length(); ik++) {
                cols.add(ist * K + (int) inchain.get(ik));
            }
        }
        jline.util.matrix.Matrix Pchain = new jline.util.matrix.Matrix(cols.size(), cols.size());
        for (int row = 0; row < cols.size(); row++) {
            for (int col = 0; col < cols.size(); col++) {
                Pchain.set(row, col, sn.rt.get(cols.get(row), cols.get(col)));
            }
        }
        jline.util.graph.DirectedGraph graph = new jline.util.graph.DirectedGraph(Pchain);
        jline.util.graph.DirectedGraph.SCCResult sccResult = graph.stronglyconncomp();
        System.out.println("\n=== Java SCC ===");
        for (int i = 0; i < sccResult.I.length; i++) System.out.print(sccResult.I[i] + " ");
        System.out.println("\n=== Java isrec ===");
        for (boolean b : sccResult.recurrent) System.out.print((b ? "1" : "0") + " ");
        System.out.println();

        System.out.println("\n=== visits ===");
        for (int c = 0; c < sn.nchains; c++) {
            System.out.println("Chain " + c + ":");
            if (sn.visits.get(c) != null) {
                sn.visits.get(c).print();
            }
        }

        MVA solver = new MVA(model);

        System.out.println("\n=== Per-Class Table ===");
        solver.getAvgTable().print();
        System.out.println("\n=== Chain Table ===");
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
