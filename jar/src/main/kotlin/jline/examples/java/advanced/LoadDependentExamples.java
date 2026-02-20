package jline.examples.java.advanced;

import jline.lang.Network;
import jline.solvers.NetworkSolver;
import jline.solvers.ctmc.CTMC;
import jline.solvers.jmt.JMT;
import jline.solvers.mva.MVA;
import jline.solvers.nc.NC;
import jline.solvers.SolverOptions;
import java.util.Scanner;

/**
 * Examples demonstrating load-dependent queueing behavior.
 * 
 * This class provides Java implementations corresponding to the Kotlin notebooks
 * in jline.examples.kotlin.advanced.loadDependent package.
 */
public class LoadDependentExamples {

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
     * Demonstrates load-dependent behavior with class dependence (ld_class_dependence.ipynb).
     * 
     * This example shows how service rates can depend on the mix of customer classes
     * present in the queue. This models situations where different types of customers
     * create different levels of resource contention.
     * 
     * Features:
     * - Service rates dependent on class composition
     * - Multi-class load-dependent queues
     * - Analysis of class interaction effects
     * - Performance isolation between classes
     * 
     * @throws Exception if the solver encounters an error
     */
    public static void ld_class_dependence() throws Exception {
        Network model = LoadDependentModel.ld_class_dependence();
        
        NetworkSolver[] solvers = new NetworkSolver[] {
            new CTMC(model),
            new JMT(model, "seed", 12345)
        };
        
        for (NetworkSolver solver : solvers) {
            try {
                
                if (solver instanceof JMT) {
                    SolverOptions options = JMT.defaultOptions();
                    options.samples = 100000;
                    ((JMT)solver).setOptions(options);
                }
                
                solver.getAvgTable().print();
            } catch (Exception e) {
            }
        }
        
        pauseForUser();
    }

    /**
     * Demonstrates load-dependent multi-server with FCFS discipline (ld_multiserver_fcfs.ipynb).
     * 
     * This example models a multi-server queue where the effective service rate depends
     * on the number of customers in the system. Under FCFS, this creates complex
     * interactions between waiting times and system occupancy.
     * 
     * Features:
     * - Multi-server queue with load-dependent rates
     * - FCFS scheduling discipline
     * - Modeling of server activation/deactivation
     * - Analysis of economies/diseconomies of scale
     * 
     * @throws Exception if the solver encounters an error
     */
    public static void ld_multiserver_fcfs() throws Exception {
        Network model = LoadDependentModel.ld_multiserver_fcfs();
        
        NetworkSolver[] solvers = new NetworkSolver[] {
            new CTMC(model),
            new JMT(model, "seed", 12345)
        };
        
        for (NetworkSolver solver : solvers) {
            try {
                
                if (solver instanceof JMT) {
                    SolverOptions options = JMT.defaultOptions();
                    options.samples = 100000;
                    ((JMT)solver).setOptions(options);
                }
                
                solver.getAvgTable().print();
            } catch (Exception e) {
            }
        }
        
        pauseForUser();
    }

    /**
     * Demonstrates load-dependent multi-server with processor sharing (ld_multiserver_ps.ipynb).
     * 
     * This example shows load-dependent behavior in a processor sharing environment.
     * PS discipline distributes service capacity equally among all customers, and the
     * load dependence models how this capacity changes with system load.
     * 
     * Features:
     * - Multi-server queue with processor sharing
     * - Load-dependent service capacity
     * - Fair sharing of variable capacity
     * - Comparison with FCFS load dependence
     * 
     * @throws Exception if the solver encounters an error
     */
    public static void ld_multiserver_ps() throws Exception {
        Network model = LoadDependentModel.ld_multiserver_ps();
        
        NetworkSolver[] solvers = new NetworkSolver[] {
            new NC(model, "method", "rd"),
            new MVA(model, "method", "exact"),
            new JMT(model, "seed", 12345)
        };
        
        for (NetworkSolver solver : solvers) {
            try {
                
                if (solver instanceof JMT) {
                    SolverOptions options = JMT.defaultOptions();
                    options.samples = 100000;
                    ((JMT)solver).setOptions(options);
                }
                
                solver.getAvgTable().print();
            } catch (Exception e) {
            }
        }
        
        pauseForUser();
    }

    /**
     * Demonstrates load-dependent PS with two classes (ld_multiserver_ps_twoclasses.ipynb).
     * 
     * This example extends the PS load-dependent case to multiple customer classes.
     * Each class can have different impacts on the load-dependent behavior, modeling
     * heterogeneous workloads with different resource requirements.
     * 
     * Features:
     * - Two-class processor sharing system
     * - Class-specific load dependence effects
     * - Analysis of fairness across classes
     * - Performance differentiation strategies
     * 
     * @throws Exception if the solver encounters an error
     */
    public static void ld_multiserver_ps_twoclasses() throws Exception {
        Network model = LoadDependentModel.ld_multiserver_ps_twoclasses();
        
        NetworkSolver[] solvers = new NetworkSolver[] {
            new NC(model, "method", "rd"),
            new MVA(model, "method", "exact"),
            new JMT(model, "seed", 12345)
        };
        
        for (NetworkSolver solver : solvers) {
            try {
                
                if (solver instanceof JMT) {
                    SolverOptions options = JMT.defaultOptions();
                    options.samples = 100000;
                    ((JMT)solver).setOptions(options);
                }
                
                solver.getAvgTable().print();
            } catch (Exception e) {
            }
        }
        
        pauseForUser();
    }

    /**
     * Main method to run all load-dependent examples.
     * 
     * @param args command line arguments (not used)
     */
    public static void main(String[] args) {
        System.out.println("\n=== Running example: ld_class_dependence ===");
        try {
            ld_class_dependence();
        } catch (Exception e) {
            System.err.println("ld_class_dependence failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: ld_multiserver_fcfs ===");
        try {
            ld_multiserver_fcfs();
        } catch (Exception e) {
            System.err.println("ld_multiserver_fcfs failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: ld_multiserver_ps ===");
        try {
            ld_multiserver_ps();
        } catch (Exception e) {
            System.err.println("ld_multiserver_ps failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: ld_multiserver_ps_twoclasses ===");
        try {
            ld_multiserver_ps_twoclasses();
        } catch (Exception e) {
            System.err.println("ld_multiserver_ps_twoclasses failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        scanner.close();
    }
}