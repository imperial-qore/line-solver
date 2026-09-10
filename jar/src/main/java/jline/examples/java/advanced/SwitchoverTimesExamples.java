package jline.examples.java.advanced;

import jline.lang.Network;
import jline.solvers.NetworkSolver;
import jline.solvers.wrappers.jmt.JMT;
import java.util.Scanner;

/**
 * Examples demonstrating switchover times in queueing systems.
 * 
 * This class provides Java implementations corresponding to the example notebooks
 * in jline.examples.java.advanced.switchoverTimes package.
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
        
        // The reference's own seed and its own run length, which is the engine
        // default: setOptions(defaultOptions()) would REPLACE the seed set here
        // with 0, and a zero seed is drawn at random, so the row would not repeat.
        NetworkSolver solver = new JMT(model, "seed", 23000, "keep", true);
        
        try {
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