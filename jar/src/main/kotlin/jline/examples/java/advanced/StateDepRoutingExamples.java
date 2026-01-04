package jline.examples.java.advanced;

import jline.lang.Network;
import jline.solvers.NetworkSolver;
import jline.solvers.jmt.JMT;
import jline.solvers.SolverOptions;
import java.util.Scanner;

/**
 * Examples demonstrating state-dependent routing in queueing networks.
 * 
 * This class provides Java implementations corresponding to the Kotlin notebooks
 * in jline.examples.kotlin.advanced.stateDepRouting package.
 */
public class StateDepRoutingExamples {

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
     * Demonstrates state-dependent routing in a closed network (sdroute_closed.ipynb).
     * 
     * This example shows how routing decisions can depend on the current state of the
     * system, such as queue lengths or server utilization. This models load balancing
     * strategies that adapt to current system conditions.
     * 
     * Features:
     * - Routing probabilities dependent on queue states
     * - Dynamic load balancing in closed networks
     * - Comparison with static routing strategies
     * - Analysis of adaptation effectiveness
     * 
     * @throws Exception if the solver encounters an error
     */
    public static void sdroute_closed() throws Exception {
        Network model = StateDepRoutingModel.sdroute_closed();
        
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
     * Demonstrates state-dependent routing in an open network (sdroute_open.ipynb).
     * 
     * This example extends state-dependent routing to open networks where customers
     * arrive from external sources. The routing decisions adapt to system load to
     * minimize response times or balance utilization.
     * 
     * Features:
     * - Adaptive routing for external arrivals
     * - Queue-length based routing decisions
     * - Overflow and alternate path routing
     * - Performance under varying loads
     * 
     * @throws Exception if the solver encounters an error
     */
    public static void sdroute_open() throws Exception {
        Network model = StateDepRoutingModel.sdroute_open();
        
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
     * Demonstrates state-dependent routing with two classes in a closed network (sdroute_twoclasses_closed.ipynb).
     * 
     * This example shows how state-dependent routing can be class-specific, allowing
     * different customer types to have different routing strategies based on their
     * priorities or service requirements.
     * 
     * Features:
     * - Class-specific routing strategies
     * - Multi-class state-dependent decisions
     * - Priority-aware load balancing
     * - Analysis of class interaction effects
     * 
     * @throws Exception if the solver encounters an error
     */
    public static void sdroute_twoclasses_closed() throws Exception {
        Network model = StateDepRoutingModel.sdroute_twoclasses_closed();
        
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
     * Main method to run all state-dependent routing examples.
     * 
     * @param args command line arguments (not used)
     */
    public static void main(String[] args) {
        System.out.println("\n=== Running example: sdroute_closed ===");
        try {
            sdroute_closed();
        } catch (Exception e) {
            System.err.println("sdroute_closed failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: sdroute_open ===");
        try {
            sdroute_open();
        } catch (Exception e) {
            System.err.println("sdroute_open failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: sdroute_twoclasses_closed ===");
        try {
            sdroute_twoclasses_closed();
        } catch (Exception e) {
            System.err.println("sdroute_twoclasses_closed failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        scanner.close();
    }
}