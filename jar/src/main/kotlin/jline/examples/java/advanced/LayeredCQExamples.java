package jline.examples.java.advanced;

import jline.lang.Network;
import jline.lang.constant.SolverType;
import jline.lang.layered.LayeredNetwork;
import jline.solvers.NetworkSolver;
import jline.solvers.ln.LN;
import jline.solvers.SolverOptions;
import java.util.Scanner;

/**
 * Examples demonstrating layered queueing networks with contention queues (CQ).
 * 
 * This class provides Java implementations corresponding to the Kotlin notebooks
 * in jline.examples.kotlin.advanced.layeredCQ package.
 */
public class LayeredCQExamples {

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
     * Demonstrates a layered CQ model with a single host (lcq_singlehost.ipynb).
     * 
     * This example shows a basic layered queueing network where software tasks execute
     * on a single host. The model captures both software contention (at the task level)
     * and hardware contention (at the host level), typical in client-server architectures.
     * 
     * Features:
     * - Single host with multiple software tasks
     * - Nested queueing structure (software and hardware layers)
     * - Analysis of software bottlenecks vs. hardware bottlenecks
     * - Mean value analysis for layered networks
     * 
     * @throws Exception if the solver encounters an error
     */
    public static void lcq_singlehost() throws Exception {
        LayeredNetwork model = LayeredCQModel.lcq_singlehost();
        
        LN solver = new LN(model, SolverType.MVA);
        
        try {
            solver.getAvgTable().print();
        } catch (Exception e) {
        }
        
        pauseForUser();
    }

    /**
     * Demonstrates a layered CQ model with three hosts (lcq_threehosts.ipynb).
     * 
     * This example extends the single host case to a distributed system with three hosts.
     * Different software tasks are deployed on different hosts, creating a more complex
     * interaction pattern typical of multi-tier architectures.
     * 
     * Features:
     * - Three hosts with distributed software tasks
     * - Inter-host communication patterns
     * - Analysis of distributed system bottlenecks
     * - Load balancing considerations in layered networks
     * 
     * @throws Exception if the solver encounters an error
     */
    public static void lcq_threehosts() throws Exception {
        LayeredNetwork model = LayeredCQModel.lcq_threehosts();
        
        LN solver = new LN(model, SolverType.MVA);
        
        try {
            solver.getAvgTable().print();
        } catch (Exception e) {
        }
        
        pauseForUser();
    }

    /**
     * Main method to run all layered CQ examples.
     * 
     * @param args command line arguments (not used)
     */
    public static void main(String[] args) {
        System.out.println("\n=== Running example: lcq_singlehost ===");
        try {
            lcq_singlehost();
        } catch (Exception e) {
            System.err.println("lcq_singlehost failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: lcq_threehosts ===");
        try {
            lcq_threehosts();
        } catch (Exception e) {
            System.err.println("lcq_threehosts failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        scanner.close();
    }
}