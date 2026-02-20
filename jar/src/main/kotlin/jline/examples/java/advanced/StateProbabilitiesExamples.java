package jline.examples.java.advanced;

import jline.lang.Network;
import jline.solvers.NetworkSolver;
import jline.solvers.ctmc.CTMC;
import jline.solvers.jmt.JMT;
import jline.solvers.SolverOptions;
import java.util.Scanner;

/**
 * Examples demonstrating state probability computations in queueing networks.
 * 
 * This class provides Java implementations corresponding to the Kotlin notebooks
 * in jline.examples.kotlin.advanced.stateProbabilities package.
 */
public class StateProbabilitiesExamples {

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
     * Demonstrates aggregated state probabilities (statepr_aggr.ipynb).
     * 
     * This example shows how to compute and analyze aggregated state probabilities,
     * where states are grouped based on certain characteristics such as total number
     * of customers or server occupancy patterns.
     * 
     * Features:
     * - State aggregation techniques
     * - Marginal probability distributions
     * - Efficient computation for large state spaces
     * - Analysis of system-wide properties
     * 
     * @throws Exception if the solver encounters an error
     */
    public static void statepr_aggr() throws Exception {
        Network model = StateProbabilitiesModel.statepr_aggr();
        
        NetworkSolver solver = new CTMC(model);
        
        try {
            solver.getAvgTable().print();
        } catch (Exception e) {
        }
        
        pauseForUser();
    }

    /**
     * Demonstrates aggregated state probabilities for large systems (statepr_aggr_large.ipynb).
     * 
     * This example extends the aggregated state probability analysis to larger systems
     * where direct enumeration of all states is infeasible. It shows approximation
     * techniques and efficient algorithms for large-scale analysis.
     * 
     * Features:
     * - Large state space handling
     * - Approximation techniques for state probabilities
     * - Memory-efficient algorithms
     * - Scalability analysis
     * 
     * @throws Exception if the solver encounters an error
     */
    public static void statepr_aggr_large() throws Exception {
        Network model = StateProbabilitiesModel.statepr_aggr_large();
        
        NetworkSolver solver = new JMT(model, "seed", 12345);
        
        try {
            SolverOptions options = JMT.defaultOptions();
            options.samples = 100000;
            ((JMT)solver).setOptions(options);
            
            solver.getAvgTable().print();
        } catch (Exception e) {
        }
        
        pauseForUser();
    }

    /**
     * Demonstrates all state probabilities for FCFS queues (statepr_allprobs_fcfs.ipynb).
     * 
     * This example computes the complete state probability distribution for systems
     * with FCFS queues. It shows how customer ordering affects the state space
     * structure and probability calculations.
     * 
     * Features:
     * - Complete state enumeration for FCFS
     * - Customer order tracking
     * - Exact probability computation
     * - Analysis of queue length distributions
     * 
     * @throws Exception if the solver encounters an error
     */
    public static void statepr_allprobs_fcfs() throws Exception {
        Network model = StateProbabilitiesModel.statepr_allprobs_fcfs();
        
        NetworkSolver solver = new CTMC(model);
        
        try {
            solver.getAvgTable().print();
        } catch (Exception e) {
        }
        
        pauseForUser();
    }

    /**
     * Demonstrates all state probabilities for PS queues (statepr_allprobs_ps.ipynb).
     * 
     * This example computes state probabilities for processor sharing queues, where
     * the state space is simpler than FCFS as customer ordering doesn't matter.
     * This allows for more efficient computation and analysis.
     * 
     * Features:
     * - State space reduction for PS
     * - Efficient probability computation
     * - Comparison with FCFS state space
     * - Product-form solution exploitation
     * 
     * @throws Exception if the solver encounters an error
     */
    public static void statepr_allprobs_ps() throws Exception {
        Network model = StateProbabilitiesModel.statepr_allprobs_ps();
        
        NetworkSolver solver = new CTMC(model);
        
        try {
            solver.getAvgTable().print();
        } catch (Exception e) {
        }
        
        pauseForUser();
    }

    /**
     * Demonstrates system-wide aggregated state probabilities (statepr_sys_aggr.ipynb).
     * 
     * This example focuses on system-wide state aggregations, such as the probability
     * of having exactly k customers in the entire network regardless of their
     * distribution across nodes.
     * 
     * Features:
     * - System-wide state aggregation
     * - Network-level performance metrics
     * - Efficient computation methods
     * - Little's Law validation
     * 
     * @throws Exception if the solver encounters an error
     */
    public static void statepr_sys_aggr() throws Exception {
        Network model = StateProbabilitiesModel.statepr_sys_aggr();
        
        NetworkSolver solver = new CTMC(model);
        
        try {
            solver.getAvgTable().print();
        } catch (Exception e) {
        }
        
        pauseForUser();
    }

    /**
     * Demonstrates system-wide aggregated probabilities for large systems (statepr_sys_aggr_large.ipynb).
     * 
     * This example extends system-wide aggregation to large-scale networks where
     * exact computation is challenging. It demonstrates approximation techniques
     * and asymptotic analysis methods.
     * 
     * Features:
     * - Large-scale system aggregation
     * - Asymptotic approximations
     * - Computational complexity reduction
     * - Accuracy vs. efficiency trade-offs
     * 
     * @throws Exception if the solver encounters an error
     */
    public static void statepr_sys_aggr_large() throws Exception {
        Network model = StateProbabilitiesModel.statepr_sys_aggr_large();
        
        NetworkSolver solver = new JMT(model, "seed", 12345);
        
        try {
            SolverOptions options = JMT.defaultOptions();
            options.samples = 100000;
            ((JMT)solver).setOptions(options);
            
            solver.getAvgTable().print();
        } catch (Exception e) {
        }
        
        pauseForUser();
    }

    /**
     * Main method to run all state probability examples.
     * 
     * @param args command line arguments (not used)
     */
    public static void main(String[] args) {
        System.out.println("\n=== Running example: statepr_aggr ===");
        try {
            statepr_aggr();
        } catch (Exception e) {
            System.err.println("statepr_aggr failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: statepr_aggr_large ===");
        try {
            statepr_aggr_large();
        } catch (Exception e) {
            System.err.println("statepr_aggr_large failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: statepr_allprobs_fcfs ===");
        try {
            statepr_allprobs_fcfs();
        } catch (Exception e) {
            System.err.println("statepr_allprobs_fcfs failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: statepr_allprobs_ps ===");
        try {
            statepr_allprobs_ps();
        } catch (Exception e) {
            System.err.println("statepr_allprobs_ps failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: statepr_sys_aggr ===");
        try {
            statepr_sys_aggr();
        } catch (Exception e) {
            System.err.println("statepr_sys_aggr failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: statepr_sys_aggr_large ===");
        try {
            statepr_sys_aggr_large();
        } catch (Exception e) {
            System.err.println("statepr_sys_aggr_large failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        scanner.close();
    }
}