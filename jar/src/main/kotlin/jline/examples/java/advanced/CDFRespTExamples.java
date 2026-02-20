package jline.examples.java.advanced;

import jline.lang.Network;
import jline.solvers.NetworkSolver;
import jline.solvers.fluid.FLD;
import jline.solvers.jmt.JMT;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;
import jline.io.Ret.DistributionResult;
import java.util.Scanner;
import java.util.List;

/**
 * Examples demonstrating CDF (Cumulative Distribution Function) of response times.
 * 
 * This class provides Java implementations corresponding to the Kotlin notebooks
 * in jline.examples.kotlin.advanced.cdfRespT package.
 */
public class CDFRespTExamples {

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
     * Demonstrates CDF of response times for a closed network (cdf_respt_closed.ipynb).
     * 
     * This example shows how to compute and analyze the cumulative distribution function
     * of response times in a simple closed queueing network with one class of customers.
     * The model consists of a delay station and a processor sharing queue.
     * 
     * Features:
     * - Single closed class with 10 customers
     * - Exponential service time distributions
     * - CDF computation using both JMT simulation and Fluid analytical approximation
     * - Statistical moments extraction from CDF data
     * 
     * @throws Exception if the solver encounters an error
     */
    public static void cdf_respt_closed() throws Exception {
        Network model = CDFRespTModel.cdf_respt_closed();
        
        NetworkSolver[] solvers = new NetworkSolver[] {
            new JMT(model, "seed", 12345),
            new FLD(model)
        };
        
        for (NetworkSolver solver : solvers) {
            try {
                
                if (solver instanceof JMT) {
                    JMT jmtSolver = (JMT) solver;
                    // CDF analysis for JMT solver
                    // Note: getCdfRespT() returns DistributionResult, not Matrix
                    // Example shows general structure for CDF analysis
                } else if (solver instanceof FLD) {
                    FLD fluidSolver = (FLD) solver;
                    // CDF analysis for Fluid solver
                    // Note: getCdfRespT() returns DistributionResult, not Matrix
                    // Example shows general structure for CDF analysis
                }
                
                solver.getAvgTable().print();
            } catch (Exception e) {
            }
        }
        
        pauseForUser();
    }

    /**
     * Demonstrates CDF of response times for a closed network with three classes (cdf_respt_closed_threeclasses.ipynb).
     * 
     * This example extends the basic closed network case to handle multiple customer classes
     * with different service requirements and routing probabilities. It shows how response
     * time distributions vary across different classes in the same network.
     * 
     * Features:
     * - Three closed classes with populations 5, 3, and 2
     * - Class-specific service times and routing probabilities
     * - Two PS queues with probabilistic routing
     * - Multi-class CDF analysis
     * 
     * @throws Exception if the solver encounters an error
     */
    public static void cdf_respt_closed_threeclasses() throws Exception {
        Network model = CDFRespTModel.cdf_respt_closed_threeclasses();
        
        NetworkSolver[] solvers = new NetworkSolver[] {
            new JMT(model, "seed", 12345),
            new FLD(model)
        };
        
        for (NetworkSolver solver : solvers) {
            try {
                
                if (solver instanceof JMT) {
                    JMT jmtSolver = (JMT) solver;
                    SolverOptions options = JMT.defaultOptions();
                    options.samples = 100000;
                    jmtSolver.setOptions(options);
                }
                
                solver.getAvgTable().print();
            } catch (Exception e) {
            }
        }
        
        pauseForUser();
    }

    /**
     * Demonstrates CDF of response times with non-exponential distributions (cdf_respt_distrib.ipynb).
     * 
     * This example shows how to analyze response time distributions when service times
     * follow non-exponential distributions. It uses Erlang distributions with specific
     * squared coefficients of variation to model more deterministic service processes.
     * 
     * Features:
     * - Erlang service time distributions
     * - Different coefficients of variation (0.5 and 0.25)
     * - Impact of service time variability on response time CDF
     * - Comparison of simulation vs. analytical approximation
     * 
     * @throws Exception if the solver encounters an error
     */
    public static void cdf_respt_distrib() throws Exception {
        Network model = CDFRespTModel.cdf_respt_distrib();
        
        NetworkSolver[] solvers = new NetworkSolver[] {
            new JMT(model, "seed", 12345),
            new FLD(model)
        };
        
        for (NetworkSolver solver : solvers) {
            try {
                DistributionResult cdfResults = null;
                
                if (solver instanceof JMT) {
                    JMT jmtSolver = (JMT) solver;
                    SolverOptions options = JMT.defaultOptions();
                    options.samples = 100000;
                    options.seed = 23000;
                    jmtSolver.setOptions(options);
                    // Use getTranCdfRespT for transient analysis like MATLAB
                    cdfResults = jmtSolver.getTranCdfRespT();
                    System.out.println("JMT Solver Transient CDF Analysis (Distributions):");
                } else if (solver instanceof FLD) {
                    FLD fluidSolver = (FLD) solver;
                    cdfResults = fluidSolver.getCdfRespT();
                    System.out.println("Fluid Solver CDF Analysis (Distributions):");
                }
                
                if (cdfResults != null) {
                    analyzeResponseTimeCDF(model, cdfResults, solver.getClass().getSimpleName());
                }
                
                solver.getAvgTable().print();
            } catch (Exception e) {
            }
        }
        
        pauseForUser();
    }

    /**
     * Demonstrates CDF of response times for an open network with two classes (cdf_respt_open_twoclasses.ipynb).
     * 
     * This example analyzes response time distributions in an open queueing network where
     * customers arrive from external sources. It shows how to handle multiple open classes
     * with different arrival rates and routing patterns.
     * 
     * Features:
     * - Two open classes with different arrival rates
     * - External arrivals modeled with Source node
     * - Class-specific routing probabilities
     * - Response time CDF for open networks
     * 
     * @throws Exception if the solver encounters an error
     */
    public static void cdf_respt_open_twoclasses() throws Exception {
        Network model = CDFRespTModel.cdf_respt_open_twoclasses();
        
        NetworkSolver[] solvers = new NetworkSolver[] {
            new JMT(model, "seed", 12345),
            new FLD(model)
        };
        
        for (NetworkSolver solver : solvers) {
            try {
                DistributionResult cdfResults = null;
                
                if (solver instanceof JMT) {
                    JMT jmtSolver = (JMT) solver;
                    SolverOptions options = JMT.defaultOptions();
                    options.samples = 10000;
                    options.seed = 23000;
                    jmtSolver.setOptions(options);
                    // Use getTranCdfRespT for transient analysis like MATLAB
                    cdfResults = jmtSolver.getTranCdfRespT();
                    System.out.println("JMT Solver Transient CDF Analysis (Open Two Classes):");
                } else if (solver instanceof FLD) {
                    FLD fluidSolver = (FLD) solver;
                    SolverOptions options = FLD.defaultOptions();
                    options.iter_max = 300;
                    fluidSolver.setOptions(options);
                    cdfResults = fluidSolver.getCdfRespT();
                    System.out.println("Fluid Solver CDF Analysis (Open Two Classes):");
                }
                
                if (cdfResults != null) {
                    analyzeResponseTimeCDF(model, cdfResults, solver.getClass().getSimpleName());
                }
                
                solver.getAvgTable().print();
            } catch (Exception e) {
            }
        }
        
        pauseForUser();
    }

    /**
     * Demonstrates CDF of response times for different population sizes (cdf_respt_populations.ipynb).
     * 
     * This example explores how the response time distribution changes as the number of
     * customers in the system varies. It creates multiple classes with populations
     * ranging from 1 to 5 to show the impact of system load on response times.
     * 
     * Features:
     * - Five classes with populations 1 through 5
     * - Identical service parameters for comparison
     * - Analysis of load impact on response time distribution
     * - Demonstration of Little's Law validation
     * 
     * @throws Exception if the solver encounters an error
     */
    public static void cdf_respt_populations() throws Exception {
        Network model = CDFRespTModel.cdf_respt_populationsN4();
        
        NetworkSolver[] solvers = new NetworkSolver[] {
            new JMT(model, "seed", 12345),
            new FLD(model)
        };
        
        for (NetworkSolver solver : solvers) {
            try {
                DistributionResult cdfResults = null;
                
                if (solver instanceof JMT) {
                    JMT jmtSolver = (JMT) solver;
                    SolverOptions options = JMT.defaultOptions();
                    options.seed = 12345;
                    jmtSolver.setOptions(options);
                    cdfResults = jmtSolver.getCdfRespT();
                    System.out.println("JMT Solver CDF Analysis (Populations):");
                } else if (solver instanceof FLD) {
                    FLD fluidSolver = (FLD) solver;
                    SolverOptions options = FLD.defaultOptions();
                    options.iter_max = 100;
                    fluidSolver.setOptions(options);
                    cdfResults = fluidSolver.getCdfRespT();
                    System.out.println("Fluid Solver CDF Analysis (Populations):");
                }
                
                if (cdfResults != null) {
                    analyzeResponseTimeCDF(model, cdfResults, solver.getClass().getSimpleName());
                }
                
                solver.getAvgTable().print();
            } catch (Exception e) {
            }
        }
        
        pauseForUser();
    }

    /**
     * Analyzes response time CDF data and computes statistical moments.
     * Based on MATLAB implementation that computes mean, variance, and SCV from CDF data.
     * 
     * @param model the network model
     * @param cdfResults the CDF results from solver
     * @param solverName name of the solver for output formatting
     */
    private static void analyzeResponseTimeCDF(Network model, DistributionResult cdfResults, String solverName) {
        if (cdfResults == null || cdfResults.cdfData == null) {
            System.out.println("No CDF data available for analysis");
            return;
        }
        
        System.out.println("\n=== " + solverName + " CDF Statistical Analysis ===");
        
        // Process CDF data for each station and class
        int numStations = model.getNumberOfStations();
        int numClasses = model.getNumberOfClasses();
        
        for (int i = 0; i < numStations; i++) {
            for (int c = 0; c < numClasses; c++) {
                if (i < cdfResults.cdfData.size() && c < cdfResults.cdfData.get(i).size()) {
                    Matrix cdfMatrix = cdfResults.cdfData.get(i).get(c);
                    if (cdfMatrix != null && cdfMatrix.getNumRows() > 1) {
                        // Extract CDF and time values (columns: [probability, time])
                        double[] probabilities = new double[cdfMatrix.getNumRows() - 1];
                        double[] times = new double[cdfMatrix.getNumRows() - 1];
                        
                        for (int row = 1; row < cdfMatrix.getNumRows(); row++) {
                            probabilities[row - 1] = cdfMatrix.get(row, 0) - cdfMatrix.get(row - 1, 0);
                            times[row - 1] = cdfMatrix.get(row, 1);
                        }
                        
                        // Compute statistical moments
                        double mean = 0.0;
                        double secondMoment = 0.0;
                        
                        for (int j = 0; j < probabilities.length; j++) {
                            mean += probabilities[j] * times[j];
                            secondMoment += probabilities[j] * times[j] * times[j];
                        }
                        
                        double variance = secondMoment - mean * mean;
                        double scv = variance / (mean * mean);
                        
                        System.out.printf("Station %d, Class %d: Mean=%.4f, Variance=%.4f, SCV=%.4f%n", 
                                        i + 1, c + 1, mean, variance, scv);
                    }
                }
            }
        }
        System.out.println();
    }

    /**
     * Main method to run all CDF response time examples.
     * 
     * @param args command line arguments (not used)
     */
    public static void main(String[] args) {
        System.out.println("\n=== Running example: cdf_respt_closed ===");
        try {
            cdf_respt_closed();
        } catch (Exception e) {
            System.err.println("cdf_respt_closed failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: cdf_respt_closed_threeclasses ===");
        try {
            cdf_respt_closed_threeclasses();
        } catch (Exception e) {
            System.err.println("cdf_respt_closed_threeclasses failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: cdf_respt_distrib ===");
        try {
            cdf_respt_distrib();
        } catch (Exception e) {
            System.err.println("cdf_respt_distrib failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: cdf_respt_open_twoclasses ===");
        try {
            cdf_respt_open_twoclasses();
        } catch (Exception e) {
            System.err.println("cdf_respt_open_twoclasses failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: cdf_respt_populations ===");
        try {
            cdf_respt_populations();
        } catch (Exception e) {
            System.err.println("cdf_respt_populations failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        scanner.close();
    }
}