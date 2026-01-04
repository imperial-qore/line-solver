package jline.examples.java.advanced;

import jline.lang.Network;
import jline.solvers.NetworkSolver;
import jline.solvers.jmt.JMT;
import jline.solvers.SolverOptions;
import java.util.Scanner;

/**
 * Examples demonstrating cyclic polling systems.
 * 
 * This class provides Java implementations corresponding to the Kotlin notebooks
 * in jline.examples.kotlin.advanced.cyclicPolling package.
 */
public class CyclicPollingExamples {

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
     * Demonstrates exhaustive polling with deterministic service times (polling_exhaustive_det.ipynb).
     * 
     * This example models a polling system where the server visits queues in cyclic order
     * and serves all customers at a queue before moving to the next (exhaustive service).
     * The service times are deterministic, allowing analysis of worst-case behaviors.
     * 
     * Features:
     * - Multiple queues with POLLING scheduling strategy
     * - Exhaustive service discipline
     * - Deterministic service time distributions
     * - Cyclic server movement pattern
     * 
     * @throws Exception if the solver encounters an error
     */
    public static void polling_exhaustive_det() throws Exception {
        Network model = CyclicPollingModel.polling_exhaustive_det();
        
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
     * Demonstrates exhaustive polling with exponential service times (polling_exhaustive_exp.ipynb).
     * 
     * This example is similar to the deterministic case but uses exponential service times,
     * which is more common in practice and allows for analytical approximations. The
     * exhaustive service discipline ensures fairness by completely serving each queue.
     * 
     * Features:
     * - Multiple queues with POLLING scheduling strategy
     * - Exhaustive service discipline
     * - Exponential service time distributions
     * - Analysis of mean waiting times and queue lengths
     * 
     * @throws Exception if the solver encounters an error
     */
    public static void polling_exhaustive_exp() throws Exception {
        Network model = CyclicPollingModel.polling_exhaustive_exp();
        
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
     * Demonstrates gated polling service discipline (polling_gated.ipynb).
     * 
     * In gated polling, the server only serves customers that were present when it arrived
     * at the queue. New arrivals during service must wait for the next cycle. This provides
     * better control over service times but may increase waiting times.
     * 
     * Features:
     * - Multiple queues with POLLING scheduling strategy
     * - Gated service discipline
     * - Comparison with exhaustive service
     * - Analysis of cycle times and fairness
     * 
     * @throws Exception if the solver encounters an error
     */
    public static void polling_gated() throws Exception {
        Network model = CyclicPollingModel.polling_gated();
        
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
     * Demonstrates k-limited polling service discipline (polling_klimited.ipynb).
     * 
     * In k-limited polling, the server serves at most k customers at each queue before
     * moving to the next. This provides a compromise between exhaustive and 1-limited
     * (pure cyclic) service, balancing fairness and efficiency.
     * 
     * Features:
     * - Multiple queues with POLLING scheduling strategy
     * - K-limited service discipline with configurable k
     * - Analysis of the impact of k on performance
     * - Trade-off between fairness and efficiency
     * 
     * @throws Exception if the solver encounters an error
     */
    public static void polling_klimited() throws Exception {
        Network model = CyclicPollingModel.polling_klimited();
        
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
     * Main method to run all cyclic polling examples.
     * 
     * @param args command line arguments (not used)
     */
    public static void main(String[] args) {
        System.out.println("\n=== Running example: polling_exhaustive_det ===");
        try {
            polling_exhaustive_det();
        } catch (Exception e) {
            System.err.println("polling_exhaustive_det failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: polling_exhaustive_exp ===");
        try {
            polling_exhaustive_exp();
        } catch (Exception e) {
            System.err.println("polling_exhaustive_exp failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: polling_gated ===");
        try {
            polling_gated();
        } catch (Exception e) {
            System.err.println("polling_gated failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: polling_klimited ===");
        try {
            polling_klimited();
        } catch (Exception e) {
            System.err.println("polling_klimited failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        scanner.close();
    }
}