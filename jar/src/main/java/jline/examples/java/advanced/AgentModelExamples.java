package jline.examples.java.advanced;

import jline.lang.Network;
import jline.solvers.NetworkSolver;
import jline.solvers.SolverOptions;
import jline.solvers.ag.AG;
import jline.solvers.ctmc.CTMC;
import jline.solvers.mva.MVA;
import java.util.Scanner;

/**
 * Examples demonstrating SolverAG, LINE's agent-based solver, and its INAP methods.
 *
 * SolverAG decomposes queueing
 * networks into interacting stochastic processes. INAP (Iterative Numerical
 * Approximation Procedure) efficiently solves the resulting fixed-point equations.
 *
 * These methods used to be reached as SolverMAM('inap'); they are SolverAG's
 * since those methods moved there, and SolverMAM no longer lists them.
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
     * Demonstrates SolverAG with the INAP method on a simple open network with two queues in series.
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

        // AG with INAP method
        System.out.println("AG (method=inap):");
        NetworkSolver solverInap = new AG(model, "inap");
        solverInap.getAvgTable().print();

        // AG with exact method
        System.out.println("\nAG (method=exact):");
        NetworkSolver solverExact = new AG(model, "exact");
        solverExact.getAvgTable().print();

        // MVA for comparison
        System.out.println("\nMVA:");
        NetworkSolver solverMVA = new MVA(model);
        solverMVA.getAvgTable().print();

        pauseForUser();
    }

    /**
     * Open tandem queue with PHASE-TYPE service.
     *
     * Queue1 sees the Poisson source directly, so it is an isolated M/PH/1 and
     * its mean queue length must be the Pollaczek-Khinchine value whatever the
     * reversed-rate iteration does. Both stations carry the same mean service
     * time and differ only in variability, which the earlier scalar
     * birth-death construction could not see: it returned the M/M/1 answer for
     * either.
     *
     * @throws Exception if the solver encounters an error
     */
    public static void ag_tandem_phasetype() throws Exception {
        System.out.println("=== Open Tandem with Phase-Type Service ===\n");

        Network model = AgentModel.tandemPhaseType();

        double lambda = 0.5, meanS = 1.0, scv1 = 0.5;
        double rho = lambda * meanS;
        double pk1 = rho + rho * rho * (1 + scv1) / (2 * (1 - rho));
        System.out.println("Queue1 is an isolated M/Er2/1, so its exact mean queue length is the");
        System.out.printf("Pollaczek-Khinchine value %.6f. The M/M/1 reading would be %.6f.%n%n",
                pk1, rho / (1 - rho));

        System.out.println("AG (method=inap):");
        new AG(model, "inap").getAvgTable().print();

        // 'inapinf' additionally drops the maxStates truncation, solving each
        // open component on its infinite state space through Neuts' rate
        // matrix R. That matters most at Queue2, whose service law has the
        // heavier tail.
        System.out.println("\nAG (method=inapinf):");
        new AG(model, "inapinf").getAvgTable().print();

        pauseForUser();
    }

    /**
     * Closed network with two PS queues.
     *
     * Demonstrates SolverAG with the INAP method on a closed queueing network with processor-sharing
     * discipline. Compares results against MVA and CTMC.
     *
     * @throws Exception if the solver encounters an error
     */
    public static void ag_closed_network() throws Exception {
        System.out.println("=== Closed Network (2 PS Queues, 10 jobs) ===\n");

        Network model = AgentModel.closedNetwork();

        // AG with INAP method
        System.out.println("AG (method=inap):");
        NetworkSolver solverInap = new AG(model, "inap");
        solverInap.getAvgTable().print();

        // AG with exact method
        System.out.println("\nAG (method=exact):");
        NetworkSolver solverExact = new AG(model, "exact");
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
     * Demonstrates SolverAG with the INAP method on a multiclass closed network. The solver creates
     * separate processes for each (station, class) pair.
     *
     * @throws Exception if the solver encounters an error
     */
    public static void ag_multiclass_closed() throws Exception {
        System.out.println("=== Multiclass Closed Network ===");
        System.out.println("Class 1: 5 jobs, Class 2: 3 jobs\n");

        Network model = AgentModel.multiclassClosed();

        // AG with INAP method
        System.out.println("AG (method=inap):");
        NetworkSolver solverInap = new AG(model, "inap");
        solverInap.getAvgTable().print();

        // AG with exact method
        System.out.println("\nAG (method=exact):");
        NetworkSolver solverExact = new AG(model, "exact");
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
     * Demonstrates SolverAG with the INAP method on an open Jackson network with feedback routing.
     * The solver models job transfers as synchronization actions between processes.
     *
     * @throws Exception if the solver encounters an error
     */
    public static void ag_jackson_network() throws Exception {
        System.out.println("=== Jackson Network (3 Queues) ===\n");

        Network model = AgentModel.jacksonNetwork();

        // AG with INAP method
        System.out.println("AG (method=inap):");
        NetworkSolver solverInap = new AG(model, "inap");
        solverInap.getAvgTable().print();

        // AG with exact method
        System.out.println("\nAG (method=exact):");
        NetworkSolver solverExact = new AG(model, "exact");
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
     * Demonstrates SolverAG with the INAP method on a G-network where negative customers (signals)
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

        // AG with INAP method
        System.out.println("AG (method=inap):");
        NetworkSolver solverInap = new AG(model, "inap");
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
        System.out.println("   AG (agent-based) Examples");
        System.out.println("========================================\n");

        ag_tandem_open();
        ag_tandem_phasetype();
        ag_closed_network();
        ag_multiclass_closed();
        ag_jackson_network();
        ag_gnetwork();

        System.out.println("\nAll examples completed.");
    }
}
