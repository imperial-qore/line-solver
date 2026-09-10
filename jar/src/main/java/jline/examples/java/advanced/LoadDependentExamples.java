package jline.examples.java.advanced;

import jline.api.fes.FESResult;
import jline.lang.Network;
import jline.lang.ModelAdapter;
import jline.lang.nodes.Station;
import jline.solvers.NetworkSolver;
import jline.solvers.ctmc.CTMC;
import jline.solvers.wrappers.jmt.JMT;
import jline.solvers.mva.MVA;
import jline.solvers.nc.NC;
import jline.solvers.SolverOptions;
import java.util.Scanner;

/**
 * Examples demonstrating load-dependent queueing behavior.
 * 
 * This class provides Java implementations corresponding to the example notebooks
 * in jline.examples.java.advanced.loadDependent package.
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

        // JMT is not solved here, as in ld_class_dependence.m: the JSIM writer
        // has no representation for the class-dependence handle, so SolverJMT
        // rejects the model rather than silently solving it unscaled.
        new CTMC(model, "exact").getAvgTable().print();
        new MVA(model, "method", "qd").getAvgTable().print();

        pauseForUser();
    }

    /**
     * Demonstrates joint (non-product-form) load dependence (ld_joint_dependence.m).
     *
     * <p>The station rate reads the class-1 marginal only and is shared by every
     * class, so it is NOT the product-form beta_{i,r}(n_{i,r}) of
     * {@link #ld_class_dependence()}; the exact CTMC is the reference against
     * which the QD-AMVA approximation is read.</p>
     *
     * @throws Exception if the solver encounters an error
     */
    public static void ld_joint_dependence() throws Exception {
        Network model = LoadDependentModel.ld_joint_dependence();

        // JMT is not solved here, as in ld_joint_dependence.m: the JSIM writer
        // has no representation for the joint-dependence handle.
        new CTMC(model, "exact").getAvgTable().print();
        new MVA(model, "method", "qd").getAvgTable().print();

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

        // The reference's order on the load-dependent model: CTMC, then exact MVA,
        // then JMT at seed 23000. Its golden holds the FIRST table each solver
        // produced, so the order and the pinned method are part of the answer.
        try {
            new CTMC(model).getAvgTable().print();
        } catch (Exception e) {
            System.out.println("CTMC failed: " + e.getMessage());
        }
        try {
            new MVA(model, "method", "exact").getAvgTable().print();
        } catch (Exception e) {
            System.out.println("MVA failed: " + e.getMessage());
        }
        try {
            new JMT(model, "seed", 23000, "samples", 100000).getAvgTable().print();
        } catch (Exception e) {
            System.out.println("JMT failed: " + e.getMessage());
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
                    // The solver's OWN options, not a fresh defaultOptions(): the latter
                    // draws a RANDOM seed in its constructor (SolverOptions ->
                    // RandomManager.generateRandomSeed), so replacing the object wholesale
                    // discards the seed passed to the constructor above and makes this
                    // example irreproducible -- which is what left the ld_multiserver_ps
                    // parity row disagreeing with its golden on a Monte-Carlo margin.
                    SolverOptions options = solver.getOptions();
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
                    // The solver's OWN options, not a fresh defaultOptions(): the latter
                    // draws a RANDOM seed in its constructor (SolverOptions ->
                    // RandomManager.generateRandomSeed), so replacing the object wholesale
                    // discards the seed passed to the constructor above and makes this
                    // example irreproducible -- which is what left the ld_multiserver_ps
                    // parity row disagreeing with its golden on a Monte-Carlo margin.
                    SolverOptions options = solver.getOptions();
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
     * The subset every FES entry aggregates, named on the model it belongs to.
     *
     * <p>The stations are looked up BY NAME rather than by index: the reference
     * scripts name them (Queue1/Queue2, or Q1/Q2/Q3) and a positional lookup
     * would silently aggregate a different subset if the tandem were ever
     * reordered, which is the one error this transform cannot report.
     */
    private static java.util.List<Station> subsetOf(Network model, String[] names) {
        java.util.List<Station> subset = new java.util.ArrayList<Station>();
        for (String name : names) {
            for (Station station : model.getStations()) {
                if (station.getName().equals(name)) {
                    subset.add(station);
                    break;
                }
            }
        }
        if (subset.size() != names.length)
            throw new IllegalArgumentException("FES subset: a named station is not in the model");
        return subset;
    }

    /**
     * Flow-equivalent-server aggregation of a two-class tandem (fes_aggregation).
     *
     * <p>Queue1 and Queue2 are replaced by one limited-class-dependent station
     * whose per-class rates are the isolated subnetwork's throughputs. For a
     * closed product-form network Norton's theorem makes this EXACT, so the
     * aggregated throughput reproduces the original's rather than approximating
     * it, which is what the comparison below prints.
     *
     * @throws Exception if the solver encounters an error
     */
    public static void fes_aggregation() throws Exception {
        Network model = LoadDependentModel.fes_aggregation();
        System.out.println("MVA (original):");
        new MVA(model, "method", "exact").getAvgTable().print();

        FESResult fes = ModelAdapter.aggregateFES(
                model, subsetOf(model, new String[] {"Queue1", "Queue2"}));
        System.out.println("MVA (FES model):");
        new MVA(fes.fesModel, "method", "exact").getAvgTable().print();
        pauseForUser();
    }

    /**
     * Flow-equivalent-server aggregation of a single-class tandem (fes_single_class).
     *
     * @throws Exception if the solver encounters an error
     */
    public static void fes_single_class() throws Exception {
        Network model = LoadDependentModel.fes_single_class();
        System.out.println("MVA (original):");
        new MVA(model, "method", "exact").getAvgTable().print();

        FESResult fes = ModelAdapter.aggregateFES(
                model, subsetOf(model, new String[] {"Queue1", "Queue2"}));
        System.out.println("MVA (FES model):");
        new MVA(fes.fesModel, "method", "exact").getAvgTable().print();
        pauseForUser();
    }

    /**
     * Norton's theorem on a single-class tandem, solved by convolution
     * (ld_fes_singleclass).
     *
     * <p>All THREE queues are aggregated, leaving the think time beside one FES,
     * and the aggregate is solved with exact NC: the FES rates are Sauer's
     * chain-dependent service rates, which the convolution consumes directly.
     *
     * @throws Exception if the solver encounters an error
     */
    public static void ld_fes_singleclass() throws Exception {
        Network model = LoadDependentModel.ld_fes_singleclass();
        System.out.println("MVA (original):");
        new MVA(model, "method", "exact").getAvgTable().print();

        FESResult fes = ModelAdapter.aggregateFES(
                model, subsetOf(model, new String[] {"Q1", "Q2", "Q3"}));
        System.out.println("NC (FES model):");
        new NC(fes.fesModel, "method", "exact").getAvgTable().print();
        pauseForUser();
    }

    /**
     * The same aggregation with two classes (ld_fes_multiclass).
     *
     * @throws Exception if the solver encounters an error
     */
    public static void ld_fes_multiclass() throws Exception {
        Network model = LoadDependentModel.ld_fes_multiclass();
        System.out.println("MVA (original):");
        new MVA(model, "method", "exact").getAvgTable().print();

        FESResult fes = ModelAdapter.aggregateFES(
                model, subsetOf(model, new String[] {"Q1", "Q2", "Q3"}));
        System.out.println("NC (FES model):");
        new NC(fes.fesModel, "method", "exact").getAvgTable().print();
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

        System.out.println("\n=== Running example: fes_aggregation ===");
        try {
            fes_aggregation();
        } catch (Exception e) {
            System.err.println("fes_aggregation failed: " + e.getMessage());
            e.printStackTrace();
        }

        System.out.println("\n=== Running example: fes_single_class ===");
        try {
            fes_single_class();
        } catch (Exception e) {
            System.err.println("fes_single_class failed: " + e.getMessage());
            e.printStackTrace();
        }

        System.out.println("\n=== Running example: ld_fes_singleclass ===");
        try {
            ld_fes_singleclass();
        } catch (Exception e) {
            System.err.println("ld_fes_singleclass failed: " + e.getMessage());
            e.printStackTrace();
        }

        System.out.println("\n=== Running example: ld_fes_multiclass ===");
        try {
            ld_fes_multiclass();
        } catch (Exception e) {
            System.err.println("ld_fes_multiclass failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        scanner.close();
    }
}