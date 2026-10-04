/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java.basic;

import jline.VerboseLevel;
import jline.lang.Network;
import jline.lang.processes.MMAPt;
import jline.solvers.SolverResult;
import jline.util.matrix.Matrix;
import jline.solvers.ctmc.CTMC;
import jline.solvers.fluid.FLD;
import jline.solvers.wrappers.jmt.JMT;
import jline.solvers.mam.MAM;
import jline.solvers.mva.MVA;
import jline.solvers.nc.NC;
import jline.solvers.ssa.SSA;
import jline.solvers.ldes.LDES;
import java.util.Scanner;

/**
 * Open queueing network examples mirroring the example notebooks in openQN.
 * <p>
 * This class contains Java implementations that mirror the example notebooks
 * found in jar/src/main/java/jline/examples/java/basic/openQN/. Each method
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
    /**
     * Non-homogeneous Poisson arrivals (oqn_nhpp.ipynb).
     *
     * <p>The golden holds the discrete-event row; the transient fluid solve
     * beside it is what the reference prints second.
     */
    public static void oqn_nhpp() throws Exception {
        Network model = OpenNHPPModel.oqn_nhpp();

        try {
            new LDES(model, "seed", 1234, "samples", 100000).getAvgTable().print();
        } catch (Exception e) {
            System.out.println("LDES failed: " + e.getMessage());
        }
        try {
            new FLD(model, "timespan", new double[]{0, 12}, "verbose", 0).getAvgTable().print();
        } catch (Exception e) {
            System.out.println("FLD failed: " + e.getMessage());
        }
        pauseForUser();
    }

    /**
     * A time-inhomogeneous MAP as the SERVICE process (oqn_mapt_service.m / .py / .cpp).
     *
     * <p>The walk is taken from the SERVICE START epoch and the modulating phase
     * carries across successive services, so the stream is correlated as well as
     * non-stationary. The constant-schedule pair beside it is the degeneracy: a
     * MAP_t whose segments are identical must reproduce the ordinary MAP.
     *
     * @throws Exception if the solver fails
     */
    public static void oqn_mapt_service() throws Exception {
        new LDES(OpenMAPtServiceModel.oqn_mapt_service(),
                "seed", 1234, "samples", 1000000).getAvgTable().print();

        System.out.println("\nconstant schedule vs ordinary MAP (must agree):");
        new LDES(OpenMAPtServiceModel.constantSchedule(),
                "seed", 1234, "samples", 400000).getAvgTable().print();
        new LDES(OpenMAPtServiceModel.homogeneousReference(),
                "seed", 1234, "samples", 400000).getAvgTable().print();
        pauseForUser();
    }

    /**
     * A MARKED, time-inhomogeneous MAP arrival stream (oqn_mmapt.m / .py / .cpp).
     *
     * <p>Both segments carry an aggregate rate of 4, so the total arrival stream
     * is statistically identical throughout; what changes is the SPLIT, 9:1
     * towards Morning by day and 1:9 towards Evening by night. The cycle average
     * is a flat 2:2, so a run showing anything else is reading the marks from the
     * segment in force rather than from the average, which is the whole point of
     * the family.
     *
     * @throws Exception if the solver fails
     */
    public static void oqn_mmapt() throws Exception {
        MMAPt arrival = OpenMMAPtModel.arrival();
        System.out.printf("Aggregate arrival rate over the cycle: %.4f%n",
                arrival.getTimeAverageRate());
        double[] lam = arrival.getTimeAverageMarkRates();
        System.out.printf("Per-mark rates over the cycle:         %.4f %.4f  "
                + "(they sum to the aggregate)%n", lam[0], lam[1]);
        System.out.println("The cycle average is a 2:2 split, so a run that shows anything "
                + "else is reading\nthe marks from the segment in force rather than from "
                + "the average.\n");

        // LDES simulates the marked schedule directly: the mark is decided by
        // WHICH block's transition fired, inside the same competing-transitions
        // draw that ends the interval, so the labelling costs no extra randomness.
        new LDES(OpenMMAPtModel.example(), "seed", 23000, "samples", 400000)
                .getAvgTable().print();
        pauseForUser();
    }

    /**
     * A BATCH, MARKED, time-inhomogeneous MAP arrival stream (oqn_bmmapt.m / .py / .cpp).
     *
     * <p>Every segment fires epochs at rate 4, so the epoch stream is identical
     * throughout and only the labels move: Premium arrives in PAIRS by day and singly
     * by night, Economy the other way round. The JOB rate therefore differs from the
     * EPOCH rate, and differs per class within a segment even though the epoch rate
     * does not, which a model reading either label off the time-averaged matrices
     * would miss entirely.
     *
     * @throws Exception if the solver fails
     */
    public static void oqn_bmmapt() throws Exception {
        new LDES(OpenBMMAPtModel.example(), "seed", 23000, "samples", 400000)
                .getAvgTable().print();
        pauseForUser();
    }

    /**
     * Transient MAP/MAP/1 (mam_transient_mapmap1.m / .py / .cpp).
     *
     * <p>WHICH ENGINE RUNS is not chosen here and is not chosen by the method
     * either: SolverMAM.getTranAvg forces "ldqbd" and the transient QBD
     * applicability test then sends a correlated MAP arrival with a correlated
     * MAP service to the Laplace-domain transient QBD. That is the same
     * dispatch the MATLAB and Python references take, so the four codebases run
     * the same algorithm on the same model rather than agreeing by coincidence.
     *
     * @throws Exception if the solver fails
     */
    public static void mam_transient_mapmap1() throws Exception {
        Network model = OpenModel.mam_transient_mapmap1();

        // The horizon is named once and printed from here, NOT from
        // SolverResult.t: nothing in jline.solvers.mam populates that field, so
        // reading it is a null dereference. InitStateTwins, the other Java
        // transient example, avoids it for the same reason. The other three
        // codebases print the curve's own last time, which is this horizon.
        final double horizon = 40.0;
        MAM solver = new MAM(model, "timespan", new double[]{0, horizon});
        solver.getTranAvg();
        SolverResult res = solver.getResults();

        // Station 1 is the Queue; the Source carries no transient curve.
        //
        // A TRANSIENT CURVE IS TWO COLUMNS, [metric, t], so the last VALUE is
        // (rows-1, 0) and not the last ELEMENT: `get(nelem-1)` walks the flat
        // row-major array and lands on the last TIME instead. This example
        // printed E[N]=U=Tput=40.00000 -- the horizon, three times -- and the
        // discrepancy had never been seen because the run before it never
        // finished (see Solver_mam_transient_qbd's capacity test). The reference
        // reads `QNt{2,1}.metric(end)`, which is this column.
        Matrix qt = res.QNt[1][0];
        Matrix ut = res.UNt[1][0];
        Matrix tt = res.TNt[1][0];
        int last = qt.getNumRows() - 1;
        System.out.println("MAP/MAP/1 transient (rho=0.6), start empty:");
        System.out.printf("  t=%5.1f  E[N]=%.5f  U=%.5f  Tput=%.5f%n",
                horizon, qt.get(last, 0), ut.get(last, 0), tt.get(last, 0));
        System.out.println("  steady-state E[N] approaches 1.5458 as t -> inf.");

        // THE PARITY GOLDEN IS READ FROM THESE THREE LINES, not from the one
        // above. getTranAvg returns a curve, and no result table carries its
        // last point, so the derived reader anchors on this label exactly as it
        // does on init_state_*'s SteadyStateQLen line. All four codebases print
        // the same three labels so the cells line up.
        System.out.printf("TranEndQLen[MAM/Queue]: %.6f%n", qt.get(last, 0));
        System.out.printf("TranEndUtil[MAM/Queue]: %.6f%n", ut.get(last, 0));
        System.out.printf("TranEndTput[MAM/Queue]: %.6f%n", tt.get(last, 0));
        pauseForUser();
    }

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
        
        try {
            new MAM(model).getAvgTable().print();
        } catch (Exception e) {
            System.out.println("Solver failed: " + e.getMessage());
        }
        
        try {
            new LDES(model, "seed", 23000, "samples", 200000).getAvgTable().print();
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
        SSA solverSSA = new SSA(model, "keep", true, "verbose", 1, "seed", 23000, "samples", 1000000);
        LDES solverLDES = new LDES(model, "keep", true, "verbose", 1, "seed", 23000, "samples", 100000);
        
        Object[] solvers = {solverCTMC, solverFluid, solverMVA, solverMAM, solverNC, solverJMT, solverSSA,
                solverLDES};
        
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
                } else if (solvers[i] instanceof LDES) {
                    ((LDES) solvers[i]).getAvgTable().print();
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

        try {
            new JMT(model, "seed", 23000).getAvgTable().print();
        } catch (Exception e) {
            System.out.println("JMT failed: " + e.getMessage());
        }
        try {
            new LDES(model, "seed", 23000).getAvgTable().print();
        } catch (Exception e) {
            System.out.println("LDES failed: " + e.getMessage());
        }
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
        System.out.println("\n=== Running example: mam_transient_mapmap1 ===");
        try {
            mam_transient_mapmap1();
        } catch (Exception e) {
            System.err.println("mam_transient_mapmap1 failed: " + e.getMessage());
            e.printStackTrace();
        }
        
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