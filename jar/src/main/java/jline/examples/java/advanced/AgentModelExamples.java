package jline.examples.java.advanced;

import jline.lang.Network;
import jline.solvers.NetworkSolver;
import jline.solvers.SolverOptions;
import jline.solvers.mam.MAM;
import jline.solvers.ctmc.CTMC;
import jline.solvers.mva.MVA;
import java.util.Scanner;

/**
 * Examples demonstrating MAM with RCAT/INAP method for agent-based analysis.
 *
 * The RCAT (Reversed Compound Agent Theorem) algorithm decomposes queueing
 * networks into interacting stochastic processes. INAP (Iterative Numerical
 * Approximation Procedure) efficiently solves the resulting fixed-point equations.
 *
 * References:
 * - Marin and Rota-Bulo', "A Mean-Field Analysis of a Class of Interactive
 *   Distributed Systems", MASCOTS 2009
 * - Harrison and Llado, "Stochastic bounds and product form solutions using
 *   RCAT", ICPE 2011
 */
public class AgentModelExamples {

    private static final Scanner scanner = new Scanner(System.in);

    private static void pauseForUser() {
        if (System.console() == null) {
            System.out.println("\n[Running in non-interactive mode, continuing...]");
            return;
        }
        System.out.println("\nPress Enter to continue to next example...");
        try {
            scanner.nextLine();
        } catch (Exception e) {
            // Ignore scanner errors
        }
    }

    /**
     * Open tandem queue (M/M/1 -> M/M/1) example.
     *
     * Demonstrates MAM with INAP method on a simple open network with two queues in series.
     * Compares results against analytical M/M/1 formulas and MVA.
     *
     * @throws Exception if the solver encounters an error
     */
    public static void ag_tandem_open() throws Exception {
        System.out.println("=== Open Tandem Queue (M/M/1 -> M/M/1) ===\n");

        Network model = AgentModel.tandemOpen();

        // Analytical solution for comparison
        double lambda = 0.5;
        double mu1 = 1.0;
        double mu2 = 1.5;
        double U1_exact = lambda / mu1;
        double U2_exact = lambda / mu2;
        double Q1_exact = U1_exact / (1 - U1_exact);
        double Q2_exact = U2_exact / (1 - U2_exact);

        System.out.println("Analytical (M/M/1):");
        System.out.printf("  Queue1: U=%.4f, Q=%.4f%n", U1_exact, Q1_exact);
        System.out.printf("  Queue2: U=%.4f, Q=%.4f%n%n", U2_exact, Q2_exact);

        // MAM with INAP method
        System.out.println("MAM (method=inap):");
        NetworkSolver solverInap = new MAM(model, "inap");
        solverInap.getAvgTable().print();

        // MAM with exact method
        System.out.println("\nMAM (method=exact):");
        NetworkSolver solverExact = new MAM(model, "exact");
        solverExact.getAvgTable().print();

        // MVA for comparison
        System.out.println("\nMVA:");
        NetworkSolver solverMVA = new MVA(model);
        solverMVA.getAvgTable().print();

        pauseForUser();
    }

    /**
     * Closed network with two PS queues.
     *
     * Demonstrates MAM with INAP method on a closed queueing network with processor-sharing
     * discipline. Compares results against MVA and CTMC.
     *
     * @throws Exception if the solver encounters an error
     */
    public static void ag_closed_network() throws Exception {
        System.out.println("=== Closed Network (2 PS Queues, 10 jobs) ===\n");

        Network model = AgentModel.closedNetwork();

        // MAM with INAP method
        System.out.println("MAM (method=inap):");
        NetworkSolver solverInap = new MAM(model, "inap");
        solverInap.getAvgTable().print();

        // MAM with exact method
        System.out.println("\nMAM (method=exact):");
        NetworkSolver solverExact = new MAM(model, "exact");
        solverExact.getAvgTable().print();

        // MVA for comparison
        System.out.println("\nMVA:");
        NetworkSolver solverMVA = new MVA(model);
        solverMVA.getAvgTable().print();

        // CTMC for exact results
        System.out.println("\nCTMC (exact):");
        NetworkSolver solverCTMC = new CTMC(model);
        solverCTMC.getAvgTable().print();

        pauseForUser();
    }

    /**
     * Multiclass closed network.
     *
     * Demonstrates MAM with INAP method on a multiclass closed network. RCAT creates
     * separate processes for each (station, class) pair.
     *
     * @throws Exception if the solver encounters an error
     */
    public static void ag_multiclass_closed() throws Exception {
        System.out.println("=== Multiclass Closed Network ===");
        System.out.println("Class 1: 5 jobs, Class 2: 3 jobs\n");

        Network model = AgentModel.multiclassClosed();

        // MAM with INAP method
        System.out.println("MAM (method=inap):");
        NetworkSolver solverInap = new MAM(model, "inap");
        solverInap.getAvgTable().print();

        // MAM with exact method
        System.out.println("\nMAM (method=exact):");
        NetworkSolver solverExact = new MAM(model, "exact");
        solverExact.getAvgTable().print();

        // MVA for comparison
        System.out.println("\nMVA:");
        NetworkSolver solverMVA = new MVA(model);
        solverMVA.getAvgTable().print();

        pauseForUser();
    }

    /**
     * Jackson network with probabilistic routing.
     *
     * Demonstrates MAM with INAP method on an open Jackson network with feedback routing.
     * RCAT models job transfers as synchronization actions between processes.
     *
     * @throws Exception if the solver encounters an error
     */
    public static void ag_jackson_network() throws Exception {
        System.out.println("=== Jackson Network (3 Queues) ===\n");

        Network model = AgentModel.jacksonNetwork();

        // MAM with INAP method
        System.out.println("MAM (method=inap):");
        NetworkSolver solverInap = new MAM(model, "inap");
        solverInap.getAvgTable().print();

        // MAM with exact method
        System.out.println("\nMAM (method=exact):");
        NetworkSolver solverExact = new MAM(model, "exact");
        solverExact.getAvgTable().print();

        // MVA for comparison
        System.out.println("\nMVA:");
        NetworkSolver solverMVA = new MVA(model);
        solverMVA.getAvgTable().print();

        pauseForUser();
    }

    /**
     * G-network (Gelenbe network) with negative customers.
     *
     * Demonstrates MAM with INAP method on a G-network where negative customers (signals)
     * remove jobs from queues. This models scenarios like job cancellations
     * or service interrupts.
     *
     * Reference: Gelenbe, E. (1991). "Product-form queueing networks with
     *            negative and positive customers", Journal of Applied Probability
     *
     * @throws Exception if the solver encounters an error
     */
    public static void ag_gnetwork() throws Exception {
        System.out.println("=== G-Network (Gelenbe Network) with Negative Customers ===");
        System.out.println("Positive arrival rate: 1.0");
        System.out.println("Negative signal rate: 0.3");
        System.out.println("Service rates: mu1=2.0, mu2=3.0\n");

        Network model = AgentModel.gNetwork();

        // MAM with INAP method
        System.out.println("MAM (method=inap):");
        NetworkSolver solverInap = new MAM(model, "inap");
        solverInap.getAvgTable().print();

        // Theoretical insight
        System.out.println("\nNote: In G-networks, negative customers reduce the effective");
        System.out.println("      load at target queues by removing jobs upon arrival.");
        System.out.println("      The utilization at Queue2 is lower than it would be");
        System.out.println("      without negative signals due to job removals.");

        pauseForUser();
    }

    /**
     * Run all agent model examples.
     *
     * @param args command line arguments (unused)
     * @throws Exception if any solver encounters an error
     */
    public static void main(String[] args) throws Exception {
        System.out.println("========================================");
        System.out.println("   MAM (RCAT/INAP) Examples");
        System.out.println("========================================\n");

        ag_tandem_open();
        ag_closed_network();
        ag_multiclass_closed();
        ag_jackson_network();
        ag_gnetwork();

        System.out.println("\nAll examples completed.");
    }
}
