/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java.basic;

import jline.VerboseLevel;
import jline.lang.constant.SolverType;
import jline.lang.layered.LayeredNetwork;
import jline.solvers.AvgTable;
import jline.solvers.ln.LN;
import jline.solvers.wrappers.lqns.LQNS;
import jline.solvers.mva.MVA;
import jline.solvers.nc.NC;
import jline.solvers.SolverOptions;
import java.util.Scanner;

/**
 * Layered network examples mirroring the example notebooks in layeredModel.
 * <p>
 * This class contains Java implementations that mirror the example notebooks
 * found in jar/src/main/java/jline/examples/java/basic/layeredModel/. Each method
 * demonstrates a specific layered network concept using models from the basic package.
 * <p>
 * The examples cover:
 * - Basic layered network structures with processors and tasks
 * - Activity precedence patterns and synchronous calls
 * - Multi-solver approaches for layered networks
 * - Enterprise application modeling patterns
 * - BPMN-style workflow representations
 */
public class LayeredExamples {

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
     * Basic layered network (lqn_basic.ipynb).
     * <p>
     * Demonstrates fundamental layered network concepts with
     * processors, tasks, and activity precedence.
     */
    /**
     * A processor whose servers are not interchangeable (lqn_server_pools).
     * <p>
     * Prints the fully-compatible pool, which is the neutral declaration, and
     * then the compatibility graph, under which neither task reaches more than
     * two of the three servers.
     */
    public static void lqn_server_pools() throws Exception {
        SolverOptions opt = LN.defaultOptions();
        // a compatibility declaration is a station rate law, which only the
        // class-switching layerings can carry
        opt.method = "srvn.cs";

        System.out.println("--- homogeneous pool on P1 ---");
        new LN(LayeredModel.lqn_server_pools(false), opt).getAvgTable().print();

        System.out.println("--- compatibility pool on P1 ---");
        new LN(LayeredModel.lqn_server_pools(true), opt).getAvgTable().print();

        pauseForUser();
    }

    public static void lqn_basic() throws Exception {
        LayeredNetwork model = LayeredModel.lqn_basic();

        // The reference solves this with the layered solver on its default layers.
        new LN(model).getAvgTable().print();

        pauseForUser();
    }
    
    /**
     * Serial layered network (lqn_serial.ipynb).
     * <p>
     * Shows serial activity precedence in layered networks
     * with synchronous call patterns.
     */
    public static void lqn_serial() throws Exception {
        LayeredNetwork model = LayeredModel.lqn_serial();
        
        // Create LQNS solver as in MATLAB (lines 32-34)
        SolverOptions options = LQNS.defaultOptions();
        options.keep = true;
        LQNS solver = new LQNS(model, options);
        AvgTable avgTable = solver.getAvgTable();
        avgTable.print();
        
        // Get raw avg tables as in MATLAB (lines 37-39)
        AvgTable[] rawTables = solver.getRawAvgTables();
        rawTables[0].print();
        if (rawTables.length > 1) {
            rawTables[1].print();
        }
        
        // LN as in MATLAB (line 41)
        LN solverLN = new LN(model, SolverType.MVA);
        AvgTable avgTableLN = solverLN.getAvgTable();
        avgTableLN.print();
        
        pauseForUser();
    }
    
    /**
     * Multi-solver layered network (lqn_multi_solvers.ipynb).
     * <p>
     * Demonstrates different solver approaches for layered networks
     * with infinite server capacity.
     */
    public static void lqn_multi_solvers() throws Exception {
        LayeredNetwork model = LayeredModel.lqn_multi_solvers();
        
        // LQNS solver as in MATLAB (lines 27-29)
        SolverOptions lqnsOptions = LQNS.defaultOptions();
        lqnsOptions.keep = true;
        lqnsOptions.verbose = VerboseLevel.STD;
        LQNS lqnsSolver = new LQNS(model, lqnsOptions);
        AvgTable avgTableLQNS = lqnsSolver.getAvgTable();
        avgTableLQNS.print();
        
        // LN with MVA solver in each layer (lines 37-39)
        SolverOptions lnOptions = LN.defaultOptions();
        lnOptions.verbose = VerboseLevel.SILENT;
        lnOptions.seed = 2300;
        SolverOptions mvaOptions = MVA.defaultOptions();
        mvaOptions.verbose = VerboseLevel.SILENT;
        LN solverLN_MVA = new LN(model, (subModel) -> new MVA(subModel, mvaOptions), lnOptions);
        AvgTable avgTableLN_MVA = solverLN_MVA.getAvgTable();
        avgTableLN_MVA.print();
        
        // LN with NC solver in each layer (lines 47-49)
        SolverOptions ncOptions = NC.defaultOptions();
        ncOptions.verbose = VerboseLevel.SILENT;
        LN solverLN_NC = new LN(model, (subModel) -> new NC(subModel, ncOptions), lnOptions);
        AvgTable avgTableLN_NC = solverLN_NC.getAvgTable();
        avgTableLN_NC.print();
        
        pauseForUser();
    }
    
    /**
     * Two-task layered network (lqn_twotasks.ipynb).
     * <p>
     * Shows interaction between multiple tasks with
     * synchronous call patterns.
     */
    public static void lqn_twotasks() throws Exception {
        LayeredNetwork model = LayeredModel.lqn_twotasks();

        try {
            new LQNS(model).getAvgTable().print();
        } catch (Exception e) {
            System.out.println("LQNS failed: " + e.getMessage());
        }
        // NC LAYERS, NOT THE DEFAULT: MVA layers and NC layers are different fixed
        // points, so the layer solver the reference pins is part of the golden.
        new LN(model, (subModel) -> new NC(subModel, "verbose", false)).getAvgTable().print();

        pauseForUser();
    }
    
    /**
     * BPMN-style layered network (lqn_bpmn.ipynb).
     * <p>
     * Shows fork-join patterns with OrFork, AndFork,
     * OrJoin, and AndJoin activity precedence.
     */
    public static void lqn_bpmn() throws Exception {
        LayeredNetwork model = LayeredModel.lqn_bpmn();
        
        LQNS solver = new LQNS(model);
        AvgTable avgTable = solver.getAvgTable();
        avgTable.print();
        
        pauseForUser();
    }
    
    /**
     * Layered network with a setup / delay-off task (lqn_setup.ipynb).
     * <p>
     * The servers of task F2 switch off when idle and pay a setup
     * time when a request reactivates them.
     */
    public static void lqn_setup() throws Exception {
        LayeredNetwork model = LayeredModel.lqn_setup();

        // The reference drives MVA layers, which is what its golden holds.
        new LN(model, SolverType.MVA).getAvgTable().print();

        pauseForUser();
    }
    
    /**
     * Workflow layered network (lqn_workflows.ipynb).
     * <p>
     * Demonstrates loop and fork-join precedence patterns
     * with nested synchronous calls.
     */
    public static void lqn_workflows() throws Exception {
        LayeredNetwork model = LayeredModel.lqn_workflows();

        try {
            new LQNS(model).getAvgTable().print();
        } catch (Exception e) {
            System.out.println("LQNS failed: " + e.getMessage());
        }
        new LN(model).getAvgTable().print();

        pauseForUser();
    }
    
    /**
     * OFBiz-style layered network (lqn_ofbiz.ipynb).
     * <p>
     * Enterprise application model inspired by Apache OFBiz
     * with database and application layers.
     */
    public static void lqn_ofbiz() throws Exception {
        LayeredNetwork model = LayeredModel.lqn_ofbiz();

        try {
            new LQNS(model).getAvgTable().print();
        } catch (Exception e) {
            System.out.println("LQNS failed: " + e.getMessage());
        }
        // NC layers: the golden keys this row LN(NC).
        new LN(model, (subModel) -> new NC(subModel, "verbose", false)).getAvgTable().print();

        pauseForUser();
    }

    /**
     * Sock Shop microservice layered network (lqn_sockshop).
     * <p>
     * Demonstrates a multi-tier microservice architecture with
     * processor replication, fan-in/fan-out, and PS scheduling.
     */
    public static void lqn_sockshop() throws Exception {
        LayeredNetwork model = LayeredModel.lqn_sockshop();

        LN solverLN = new LN(model, SolverType.MVA);
        AvgTable avgTable = solverLN.getAvgTable();
        avgTable.print();

        pauseForUser();
    }

    /**
     * Main method demonstrating selected layered network examples.
     */
    public static void main(String[] args) throws Exception {
        try {
            lqn_basic();
        } catch (Exception e) {
            System.err.println("lqn_basic failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        try {
            lqn_serial();
        } catch (Exception e) {
            System.err.println("lqn_serial failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        try {
            lqn_multi_solvers();
        } catch (Exception e) {
            System.err.println("lqn_multi_solvers failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        try {
            lqn_twotasks();
        } catch (Exception e) {
            System.err.println("lqn_twotasks failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        try {
            lqn_bpmn();
        } catch (Exception e) {
            System.err.println("lqn_bpmn failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        try {
            lqn_setup();
        } catch (Exception e) {
            System.err.println("lqn_setup failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        try {
            lqn_workflows();
        } catch (Exception e) {
            System.err.println("lqn_workflows failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        try {
            lqn_ofbiz();
        } catch (Exception e) {
            System.err.println("lqn_ofbiz failed: " + e.getMessage());
            e.printStackTrace();
        }

        try {
            lqn_sockshop();
        } catch (Exception e) {
            System.err.println("lqn_sockshop failed: " + e.getMessage());
            e.printStackTrace();
        }

        scanner.close();
    }
}
