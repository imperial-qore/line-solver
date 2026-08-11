/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.bench;

import java.util.List;
import java.util.Map;

/**
 * Formats benchmark results in MATLAB-style tables
 */
public class BenchmarkFormatter {
    
    /**
     * Format benchmark results as a table similar to MATLAB output
     */
    public static void printBenchmarkSummary(String suiteName, List<BenchmarkRegression.BenchmarkResult> results) {
        if (results.isEmpty()) {
            System.out.println("No benchmark results to display.");
            return;
        }
        
        System.out.println("\n=== " + suiteName + " Benchmark Summary ===");
        
        // Get solver names from first result
        List<String> solvers = results.get(0).solverLabels;
        
        if (solvers.isEmpty()) {
            System.out.println("No solvers found in results.");
            return;
        }
        
        // Calculate column widths
        int nameWidth = Math.max(15, suiteName.length());
        int solverWidth = 12;
        
        // Print header
        System.out.printf("%-" + nameWidth + "s", "Benchmark");
        for (String solver : solvers) {
            System.out.printf(" %11s", solver.substring(0, Math.min(solver.length(), 11)));
        }
        System.out.println();
        
        // Print separator
        for (int i = 1; i < nameWidth; i++) {
            System.out.print("-");
        }
        for (String solver : solvers) {
            System.out.print(" ");
            for (int i = 1; i < 11; i++) {
                System.out.print("-");
            }
        }
        System.out.println();
        
        // Aggregate errors across all benchmarks
        for (String errorType : new String[]{"ERRQ", "ERRU", "ERRR", "ERRT"}) {
            System.out.printf("%-" + nameWidth + "s", errorType);
            
            for (String solver : solvers) {
                double totalError = 0.0;
                int validResults = 0;
                
                for (BenchmarkRegression.BenchmarkResult result : results) {
                    String key = solver + "_" + errorType;
                    Double error = result.scores.get(key);
                    if (error != null && !Double.isInfinite(error) && !Double.isNaN(error)) {
                        totalError += error;
                        validResults++;
                    }
                }
                
                if (validResults > 0) {
                    double avgError = totalError / validResults;
                    if (avgError < 1e-10) {
                        System.out.printf(" %11.2e", avgError);
                    } else if (avgError < 0.01) {
                        System.out.printf(" %11.6f", avgError);
                    } else {
                        System.out.printf(" %11.4f", avgError);
                    }
                } else {
                    System.out.printf(" %11s", "FAIL");
                }
            }
            System.out.println();
        }
        
        System.out.println();
    }
    
    /**
     * Print individual benchmark iteration results
     */
    public static void printIterationResult(String benchmarkName, int iteration, Map<String, Object> results) {
        @SuppressWarnings("unchecked")
        Map<String, Double> errQ = (Map<String, Double>) results.get("ERRQ");
        @SuppressWarnings("unchecked")
        Map<String, Double> errU = (Map<String, Double>) results.get("ERRU");
        @SuppressWarnings("unchecked")
        Map<String, Double> errR = (Map<String, Double>) results.get("ERRR");
        @SuppressWarnings("unchecked")
        Map<String, Double> errT = (Map<String, Double>) results.get("ERRT");
        
        if (errQ == null || errQ.isEmpty()) {
            return;
        }
        
        System.out.println("  Iteration " + iteration + ":");
        
        for (String solver : errQ.keySet()) {
            double errorQ = errQ.getOrDefault(solver, Double.MAX_VALUE);
            double errorU = errU != null ? errU.getOrDefault(solver, Double.MAX_VALUE) : Double.MAX_VALUE;
            double errorR = errR != null ? errR.getOrDefault(solver, Double.MAX_VALUE) : Double.MAX_VALUE;
            double errorT = errT != null ? errT.getOrDefault(solver, Double.MAX_VALUE) : Double.MAX_VALUE;
            
            System.out.printf("    %-12s Q: %8.2e  U: %8.2e  R: %8.2e  T: %8.2e%n", 
                            solver, errorQ, errorU, errorR, errorT);
        }
    }
    
    /**
     * Print regression test summary
     */
    public static void printRegressionSummary(String suiteName, boolean passed, double maxError, double tolerance) {
        System.out.println("\n=== Regression Test Results ===");
        System.out.println("Suite: " + suiteName);
        System.out.println("Status: " + (passed ? "PASS" : "FAIL"));
        System.out.printf("Max Error: %.2e%n", maxError);
        System.out.printf("Tolerance: %.2e%n", tolerance);
        
        if (passed) {
            System.out.println("✓ All benchmarks within tolerance");
        } else {
            System.out.println("✗ Some benchmarks exceeded tolerance");
        }
        System.out.println();
    }
}