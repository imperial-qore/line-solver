/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.bench;

import jline.GlobalConstants;

import java.io.*;
import java.util.*;

/**
 * Handles benchmark regression testing by storing and comparing against baseline scores
 */
public class BenchmarkRegression {
    
    
    /**
     * Container for regression baseline data
     */
    public static class RegressionBaseline implements Serializable {
        private static final long serialVersionUID = 1L;
        public Map<String, Double> scores = new HashMap<>();
        public List<String> solverLabels = new ArrayList<>();
        public String timestamp;
        public String version;
        
        public RegressionBaseline() {
            this.timestamp = new Date().toString();
            this.version = "1.0";
        }
    }
    
    /**
     * Container for benchmark results
     */
    public static class BenchmarkResult {
        public String benchmarkName;
        public Map<String, Double> scores = new HashMap<>();
        public List<String> solverLabels = new ArrayList<>();
        
        public BenchmarkResult(String name) {
            this.benchmarkName = name;
        }
    }
    
    /**
     * Container for regression comparison results
     */
    public static class RegressionComparison {
        public boolean passed;
        public double maxError;
        public double tolerance;
        public Map<String, Double> scoreDifferences = new HashMap<>();
        public List<String> failedSolvers = new ArrayList<>();
        
        public RegressionComparison(double tolerance) {
            this.tolerance = tolerance;
        }
    }
    
    /**
     * Load regression baseline from file
     */
    public static RegressionBaseline loadBaseline(String filename) throws IOException {
        File file = new File(filename);
        if (!file.exists()) {
            return null;
        }
        
        try (ObjectInputStream ois = new ObjectInputStream(new FileInputStream(file))) {
            return (RegressionBaseline) ois.readObject();
        } catch (ClassNotFoundException e) {
            throw new IOException("Error deserializing baseline file: " + e.getMessage(), e);
        }
    }
    
    /**
     * Save regression baseline to file
     */
    public static void saveBaseline(String filename, RegressionBaseline baseline) throws IOException {
        try (ObjectOutputStream oos = new ObjectOutputStream(new FileOutputStream(filename))) {
            oos.writeObject(baseline);
        }
    }
    
    /**
     * Create baseline from current benchmark results
     */
    public static RegressionBaseline createBaseline(List<BenchmarkResult> results) {
        RegressionBaseline baseline = new RegressionBaseline();
        
        if (!results.isEmpty()) {
            // Use solver labels from first result (should be consistent across all)
            baseline.solverLabels = new ArrayList<>(results.get(0).solverLabels);
            
            // Aggregate scores from all benchmark results
            Map<String, Double> aggregatedScores = new HashMap<>();
            
            // Initialize aggregated scores
            for (String solver : baseline.solverLabels) {
                aggregatedScores.put(solver + "_ERRQ", 0.0);
                aggregatedScores.put(solver + "_ERRU", 0.0);
                aggregatedScores.put(solver + "_ERRR", 0.0);
                aggregatedScores.put(solver + "_ERRT", 0.0);
            }
            
            // Sum errors across all benchmarks
            for (BenchmarkResult result : results) {
                for (Map.Entry<String, Double> entry : result.scores.entrySet()) {
                    String key = entry.getKey();
                    Double value = entry.getValue();
                    
                    if (!Double.isInfinite(value) && !Double.isNaN(value)) {
                        aggregatedScores.put(key, aggregatedScores.get(key) + value);
                    } else {
                        // Keep infinite/NaN values as they indicate solver failures
                        aggregatedScores.put(key, value);
                    }
                }
            }
            
            baseline.scores = aggregatedScores;
        }
        
        return baseline;
    }
    
    /**
     * Compare current results against baseline using maxerronsum logic
     */
    public static RegressionComparison compareAgainstBaseline(
            List<BenchmarkResult> currentResults, 
            RegressionBaseline baseline) {
        
        double tolerance = GlobalConstants.CoarseTol;
        RegressionComparison comparison = new RegressionComparison(tolerance);
        
        // Create baseline from current results for comparison
        RegressionBaseline currentBaseline = createBaseline(currentResults);
        
        double maxError = 0.0;
        boolean passed = true;
        
        // Compare each score
        for (Map.Entry<String, Double> entry : baseline.scores.entrySet()) {
            String key = entry.getKey();
            Double baselineValue = entry.getValue();
            Double currentValue = currentBaseline.scores.get(key);
            
            if (currentValue != null) {
                double error;
                
                if (Double.isInfinite(baselineValue) && Double.isInfinite(currentValue)) {
                    // Both infinite - consider as match
                    error = 0.0;
                } else if (Double.isInfinite(baselineValue) || Double.isInfinite(currentValue)) {
                    // One infinite, one finite - major difference
                    error = Double.MAX_VALUE;
                } else if (Double.isNaN(baselineValue) && Double.isNaN(currentValue)) {
                    // Both NaN - consider as match  
                    error = 0.0;
                } else if (Double.isNaN(baselineValue) || Double.isNaN(currentValue)) {
                    // One NaN, one valid - major difference
                    error = Double.MAX_VALUE;
                } else {
                    // Normal comparison using relative error
                    if (Math.abs(baselineValue) > 1e-12) {
                        error = Math.abs(currentValue - baselineValue) / Math.abs(baselineValue);
                    } else {
                        error = Math.abs(currentValue - baselineValue);
                    }
                }
                
                comparison.scoreDifferences.put(key, error);
                maxError = Math.max(maxError, error);
                
                if (error > tolerance) {
                    passed = false;
                    comparison.failedSolvers.add(key);
                }
            } else {
                // Missing score - failure
                passed = false;
                comparison.failedSolvers.add(key);
                comparison.scoreDifferences.put(key, Double.MAX_VALUE);
            }
        }
        
        comparison.maxError = maxError;
        comparison.passed = passed;
        
        return comparison;
    }
    
    /**
     * Run regression test for a benchmark suite
     */
    public static boolean runRegressionTest(String suiteName, List<BenchmarkResult> results) {
        String baselineFilename = suiteName + "_regression.dat";
        
        try {
            // Try to load existing baseline
            RegressionBaseline baseline = loadBaseline(baselineFilename);
            
            if (baseline == null) {
                // No baseline exists - create new one
                System.out.println("No baseline file found for " + suiteName + 
                                 ". Saving current results as new baseline.");
                baseline = createBaseline(results);
                saveBaseline(baselineFilename, baseline);
                return true;
            }
            
            // Compare against baseline
            RegressionComparison comparison = compareAgainstBaseline(results, baseline);
            
            if (comparison.passed) {
                System.out.println("✓ Regression test PASSED for " + suiteName + 
                                 " (max error: " + String.format("%.2e", comparison.maxError) + 
                                 ", tolerance: " + String.format("%.2e", comparison.tolerance) + ")");
                return true;
            } else {
                System.out.println("✗ Regression test FAILED for " + suiteName);
                System.out.println("  Max error: " + String.format("%.2e", comparison.maxError));
                System.out.println("  Tolerance: " + String.format("%.2e", comparison.tolerance));
                System.out.println("  Failed solvers/metrics: " + comparison.failedSolvers.size());
                
                // Show details of failures
                for (String failedKey : comparison.failedSolvers) {
                    double error = comparison.scoreDifferences.get(failedKey);
                    System.out.println("    " + failedKey + ": error = " + String.format("%.2e", error));
                }
                
                return false;
            }
            
        } catch (IOException e) {
            System.err.println("Error in regression test for " + suiteName + ": " + e.getMessage());
            return false;
        }
    }
    
    /**
     * Convert benchmark results from Map format to BenchmarkResult
     */
    public static BenchmarkResult convertToBenchmarkResult(String benchmarkName, 
                                                          Map<String, Object> results) {
        BenchmarkResult result = new BenchmarkResult(benchmarkName);
        
        @SuppressWarnings("unchecked")
        Map<String, Double> errQ = (Map<String, Double>) results.get("ERRQ");
        @SuppressWarnings("unchecked")
        Map<String, Double> errU = (Map<String, Double>) results.get("ERRU");
        @SuppressWarnings("unchecked")
        Map<String, Double> errR = (Map<String, Double>) results.get("ERRR");
        @SuppressWarnings("unchecked")
        Map<String, Double> errT = (Map<String, Double>) results.get("ERRT");
        
        // Get solver names from one of the error maps
        if (errQ != null) {
            result.solverLabels = new ArrayList<>(errQ.keySet());
            Collections.sort(result.solverLabels);
            
            // Add all error scores
            for (String solver : result.solverLabels) {
                result.scores.put(solver + "_ERRQ", errQ.getOrDefault(solver, Double.MAX_VALUE));
                result.scores.put(solver + "_ERRU", errU != null ? errU.getOrDefault(solver, Double.MAX_VALUE) : Double.MAX_VALUE);
                result.scores.put(solver + "_ERRR", errR != null ? errR.getOrDefault(solver, Double.MAX_VALUE) : Double.MAX_VALUE);
                result.scores.put(solver + "_ERRT", errT != null ? errT.getOrDefault(solver, Double.MAX_VALUE) : Double.MAX_VALUE);
            }
        }
        
        return result;
    }
}