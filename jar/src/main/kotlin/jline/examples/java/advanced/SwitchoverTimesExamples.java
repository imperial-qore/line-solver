package jline.examples.java.advanced;

import jline.lang.Network;
import jline.solvers.NetworkSolver;
import jline.solvers.jmt.JMT;
import jline.solvers.SolverOptions;
import java.util.Scanner;

/**
 * Examples demonstrating switchover times in queueing systems.
 * 
 * This class provides Java implementations corresponding to the Kotlin notebooks
 * in jline.examples.kotlin.advanced.switchoverTimes package.
 */
public class SwitchoverTimesExamples {

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
     * Demonstrates basic switchover time modeling (switchover_basic.ipynb).
     * 
     * This example shows how to model systems where there is a time penalty when
     * switching between different types of work or when a server moves between
     * queues. Switchover times are common in manufacturing, computer systems,
     * and communication networks.
     * 
     * Features:
     * - Switchover time between job classes
     * - Setup time modeling
     * - Impact on system throughput and response time
     * - Optimization of switching policies
     * 
     * @throws Exception if the solver encounters an error
     */
    public static void switchover_basic() throws Exception {
        Network model = SwitchoverTimesModel.switchover_basic();
        
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
     * Main method to run all switchover time examples.
     * 
     * @param args command line arguments (not used)
     */
    public static void main(String[] args) {
        System.out.println("\n=== Running example: switchover_basic ===");
        try {
            switchover_basic();
        } catch (Exception e) {
            System.err.println("switchover_basic failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        scanner.close();
    }
}