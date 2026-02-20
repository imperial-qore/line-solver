/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

// Copyright (c) 2012-2026, Imperial College London
// All rights reserved.

package jline.solvers;

import static jline.io.InputOutputKt.line_warning;

import jline.lang.Ensemble;
import jline.lang.Network;
import jline.lang.constant.SolverType;
import jline.VerboseLevel;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

// Abstract class for solvers applicable to Ensemble models
public abstract class EnsembleSolver extends Solver {

    protected Network[] ensemble;

    protected NetworkSolver[] solvers;

    protected Map<Integer, Map<Integer, SolverResult>> results;
    
    protected ExecutorService threadPool;
    protected int numThreads;

    protected EnsembleSolver(String name) {
        this(null, name, EnsembleSolver.defaultOptions());
    }

    protected EnsembleSolver(Ensemble ensModel, String name, SolverOptions options) {
        super(name);
        this.options = options;
        if (ensModel != null) {
            this.ensemble = ensModel.getEnsemble().toArray(new Network[0]);
            this.solvers = new NetworkSolver[ensemble.length];
        }
        this.results = new HashMap<>();
        
        // Initialize thread pool following SSA solver pattern
        this.numThreads = (int) Math.ceil(Runtime.getRuntime().availableProcessors() / 2.0);
        this.threadPool = Executors.newFixedThreadPool(numThreads);
    }

    protected EnsembleSolver(Ensemble ensModel, String name) {
        this(ensModel, name, EnsembleSolver.defaultOptions());
    }

    // Ensemble Solver options
    public static SolverOptions defaultOptions() {
        return new SolverOptions(SolverType.ENV);
    }

    // Operations within an iteration
    protected abstract SolverResult analyze(int it, int e);

    // Convergence test at iteration it
    protected abstract boolean converged(int it);

    // Operations after iterations are completed
    protected abstract void finish();

    protected abstract AvgTable getEnsembleAvg();

    private int getIteration() {
        return results.size();
    }

    public int getNumberOfModels() {
        return ensemble.length;
    }
    
    public void setNumThreads(int numThreads) {
        this.numThreads = numThreads;
        // Recreate thread pool with new thread count
        if (threadPool != null) {
            threadPool.shutdown();
            try {
                threadPool.awaitTermination(30, TimeUnit.SECONDS);
            } catch (InterruptedException ex) {
                threadPool.shutdownNow();
            }
        }
        this.threadPool = Executors.newFixedThreadPool(numThreads);
    }
    
    public int getNumThreads() {
        return numThreads;
    }

    // Operations before starting to iterate
    protected abstract void init();

    protected void iterate() {

        long outerStartTime = System.nanoTime();
        int it = 0;
        int E = getNumberOfModels();
        List<Double> solveRuntimes = new ArrayList<>();
        List<Double> synchRuntimes = new ArrayList<>();
        double totalRuntime = 0.0;
        Map<Integer, Map<Integer, SolverResult>> iterateResults = new HashMap<>();

        init();

        while (!converged(it) && it < options.iter_max) {
            it++; // Starts at Iteration 1
            pre(it);
            long solveStartTime = System.nanoTime();
            if (Objects.equals(options.method, "parallel")) {
                // Parallel execution following SSA solver pattern
                Map<Integer, SolverResult> parallelResults = new ConcurrentHashMap<>();
                
                for (int e = 0; e < E; e++) {
                    final int ensembleIndex = e;
                    final int iteration = it;
                    
                    threadPool.submit(() -> {
                        try {
                            SolverResult result = analyze(iteration, ensembleIndex);
                            parallelResults.put(ensembleIndex, result);
                        } catch (Exception ex) {
                            // Log error but don't break the entire parallel execution
                            line_warning("EnsembleSolver.runAnalyzer", "Error analyzing ensemble model %d at iteration %d: %s", ensembleIndex, iteration, ex.getMessage());
                        }
                    });
                }
                
                // Wait for all tasks to complete
                threadPool.shutdown();
                try {
                    // Wait up to 10 minutes for completion, following SSA pattern
                    threadPool.awaitTermination(600, TimeUnit.SECONDS);
                } catch (InterruptedException ex) {
                    throw new RuntimeException("Parallel execution interrupted", ex);
                }
                
                // Recreate thread pool for next iteration
                this.threadPool = Executors.newFixedThreadPool(numThreads);
                
                // Store results
                iterateResults.put(it, new HashMap<>(parallelResults));
            } else {
                for (int e = 0; e < E; e++) {
                    SolverResult results_it_e = analyze(it, e);
                    if (iterateResults.containsKey(it)) {
                        iterateResults.get(it).put(e, results_it_e);
                    } else {
                        HashMap<Integer, SolverResult> modelIdxToResultMap = new HashMap<>();
                        modelIdxToResultMap.put(e, results_it_e);
                        iterateResults.put(it, modelIdxToResultMap);
                    }
                }
            }
            this.results = new HashMap<>();
            // Deepcopy iterateResults onto this.results.
            // This code is to emulate matlab behaviour
            for (int i = 1; i <= it; i++) {
                Map<Integer, SolverResult> tempMap = new HashMap<>();
                for (int e = 0; e < E; e++) {
                    tempMap.put(e, iterateResults.get(i).get(e).deepCopy());
                }
                this.results.put(i, tempMap);
            }
            double solveTime = (System.nanoTime() - solveStartTime) / 1000000000.0;
            if (options.verbose != VerboseLevel.SILENT) {
                solveRuntimes.add(solveTime);
                totalRuntime = (System.nanoTime() - outerStartTime) / 1000000000.0;
                System.out.printf("\nIter %2d. ", it);
            }

            long synchStartTime = System.nanoTime();
            post(it);
            double synchTime = (System.nanoTime() - synchStartTime) / 1000000000.0;
            synchRuntimes.add(synchTime);
            if (options.verbose != VerboseLevel.SILENT) {
                System.out.printf(" Analyze time: %.3fs. Update time: %.3fs. Runtime: %.3fs. ", solveTime, synchTime, totalRuntime);
            }
        }

        finish();

        // Print summary if verbose
        if (options.verbose != VerboseLevel.SILENT && !solveRuntimes.isEmpty()) {
            double avgSolve = solveRuntimes.stream().mapToDouble(Double::doubleValue).average().orElse(0.0);
            double avgSynch = synchRuntimes.stream().mapToDouble(Double::doubleValue).average().orElse(0.0);
            System.out.printf("\nSummary: Analyze avg time: %.3fs. Update avg time: %.3fs. Total runtime: %.3fs. ",
                       avgSolve, avgSynch, totalRuntime);
        }

        // Cleanup thread pool
        if (threadPool != null && !threadPool.isShutdown()) {
            threadPool.shutdown();
            try {
                threadPool.awaitTermination(30, TimeUnit.SECONDS);
            } catch (InterruptedException ex) {
                threadPool.shutdownNow();
            }
        }
    }

    // Operations after an iteration
    protected abstract void post(int it);

    // Operations before an iteration
    protected abstract void pre(int it);

    public void printEnsembleAvgTables() {
        int E = getNumberOfModels();
        for (int e = 0; e < E; e++) {
            solvers[e].getAvgTable();
        }
    }

    // printEnsembleAvgTables -> printEnsembleAvgTs alias
    public void printEnsembleAvgTs() {
        printEnsembleAvgTables();
    }

    // NOTE: the following LINE methods have not been migrated to JLINE
    // a) list() - duplicative purpose with getNumberOfModels()
    // b) getSolver() - solvers is public instead of using getter
    // c) setSolver() x 3 - solvers is public so no need for setters
    
    // ========== Kotlin-style Alias Methods ==========
    // Aliases for get* methods following Kotlin naming conventions
    
    public int numberOfModels() { return getNumberOfModels(); }
    protected AvgTable ensembleAvg() { return getEnsembleAvg(); }
    public int numThreads() { return getNumThreads(); }
}
