package jline.examples.java.advanced;

import jline.lang.Network;
import jline.solvers.NetworkSolver;
import jline.solvers.ctmc.CTMC;
import jline.solvers.jmt.JMT;
import jline.solvers.SolverOptions;
import java.util.Scanner;

/**
 * Examples demonstrating initial state configurations in queueing networks.
 * 
 * This class provides Java implementations corresponding to the Kotlin notebooks
 * in jline.examples.kotlin.advanced.initState package.
 */
public class InitStateExamples {

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
     * Demonstrates initial state with FCFS queues and exponential service times (init_state_fcfs_exp.ipynb).
     * 
     * This example shows how to specify and analyze the transient behavior of a queueing
     * network starting from a specific initial state. The FCFS discipline with exponential
     * service times allows for Markovian analysis of the transient behavior.
     * 
     * Features:
     * - FCFS (First-Come-First-Served) scheduling discipline
     * - Exponential service time distributions
     * - Specific initial customer distribution
     * - Transient and steady-state analysis comparison
     * 
     * @throws Exception if the solver encounters an error
     */
    public static void init_state_fcfs_exp() throws Exception {
        Network model = InitStateModel.init_state_fcfs_exp();
        
        NetworkSolver[] solvers = new NetworkSolver[] {
            new JMT(model, "seed", 12345),
            new CTMC(model)
        };
        
        for (NetworkSolver solver : solvers) {
            try {
                
                if (solver instanceof JMT) {
                    SolverOptions options = JMT.defaultOptions();
                    options.samples = 50000;
                    ((JMT)solver).setOptions(options);
                }
                
                solver.getAvgTable().print();
            } catch (Exception e) {
            }
        }
        
        pauseForUser();
    }

    /**
     * Demonstrates initial state with FCFS queues and non-exponential service times (init_state_fcfs_nonexp.ipynb).
     * 
     * This example extends the FCFS initial state analysis to non-exponential service times.
     * Non-exponential distributions make the transient analysis more complex as the system
     * is no longer Markovian, requiring simulation or approximation techniques.
     * 
     * Features:
     * - FCFS scheduling with non-exponential service times
     * - Erlang and deterministic distributions
     * - Impact of service time variability on transient behavior
     * - Simulation-based transient analysis
     * 
     * @throws Exception if the solver encounters an error
     */
    public static void init_state_fcfs_nonexp() throws Exception {
        Network model = InitStateModel.init_state_fcfs_nonexp();
        
        NetworkSolver solver = new JMT(model, "seed", 12345);
        
        try {
            SolverOptions options = JMT.defaultOptions();
            options.samples = 50000;
            ((JMT)solver).setOptions(options);
            
            solver.getAvgTable().print();
        } catch (Exception e) {
        }
        
        pauseForUser();
    }

    /**
     * Demonstrates initial state with processor sharing (PS) discipline (init_state_ps.ipynb).
     * 
     * This example analyzes initial state configurations in processor sharing queues.
     * PS discipline provides different transient behavior compared to FCFS, as all
     * customers receive service simultaneously at reduced rates.
     * 
     * Features:
     * - Processor Sharing (PS) scheduling discipline
     * - Analysis of fairness in transient regime
     * - Comparison with FCFS initial state behavior
     * - Impact on response time distributions
     * 
     * @throws Exception if the solver encounters an error
     */
    public static void init_state_ps() throws Exception {
        Network model = InitStateModel.init_state_ps();
        
        NetworkSolver[] solvers = new NetworkSolver[] {
            new JMT(model, "seed", 12345),
            new CTMC(model)
        };
        
        for (NetworkSolver solver : solvers) {
            try {
                
                if (solver instanceof JMT) {
                    SolverOptions options = JMT.defaultOptions();
                    options.samples = 50000;
                    ((JMT)solver).setOptions(options);
                }
                
                solver.getAvgTable().print();
            } catch (Exception e) {
            }
        }
        
        pauseForUser();
    }

    /**
     * Main method to run all initial state examples.
     * 
     * @param args command line arguments (not used)
     */
    public static void main(String[] args) {
        System.out.println("\n=== Running example: init_state_fcfs_exp ===");
        try {
            init_state_fcfs_exp();
        } catch (Exception e) {
            System.err.println("init_state_fcfs_exp failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: init_state_fcfs_nonexp ===");
        try {
            init_state_fcfs_nonexp();
        } catch (Exception e) {
            System.err.println("init_state_fcfs_nonexp failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: init_state_ps ===");
        try {
            init_state_ps();
        } catch (Exception e) {
            System.err.println("init_state_ps failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        scanner.close();
    }
}