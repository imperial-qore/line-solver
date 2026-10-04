/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.bench;

import jline.util.matrix.Matrix;
import jline.lang.Network;
import jline.lang.nodes.Node;
import jline.lang.nodes.Delay;
import jline.util.Utils;

import java.util.ArrayList;
import java.util.List;

/**
 * Utility functions for benchmarking
 */
public class BenchmarkUtils {
    
    /**
     * Calculate maximum error on sum for queue length metrics
     * Returns max absolute error relative to sum of each column
     * Similar to MATLAB: max(abs(exact-approx))/sum(exact)
     */
    public static double maxErrorOnSum(Matrix approx, Matrix exact) {
        if (approx == null || exact == null) {
            return Double.MAX_VALUE;
        }
        
        if (approx.getNumRows() != exact.getNumRows() || 
            approx.getNumCols() != exact.getNumCols()) {
            return Double.MAX_VALUE;
        }
        
        double maxError = 0.0;
        
        // Process each column (class)
        for (int j = 0; j < exact.getNumCols(); j++) {
            double colMax = 0.0;
            double colSum = 0.0;
            
            // Find max absolute difference and sum for this column
            for (int i = 0; i < exact.getNumRows(); i++) {
                double diff = Math.abs(exact.get(i, j) - approx.get(i, j));
                colMax = Math.max(colMax, diff);
                colSum += exact.get(i, j);
            }
            
            // Compute relative error for this column
            if (colSum > 0) {
                double colError = colMax / colSum;
                maxError = Math.max(maxError, colError);
            }
        }
        
        return maxError;
    }
    
    /**
     * Calculate mean error on sum for queue length metrics
     * Returns mean absolute error relative to sum of each column
     * Similar to MATLAB: mean(abs(exact-approx))/sum(exact)
     */
    public static double errOnSum(Matrix approx, Matrix exact) {
        if (approx == null || exact == null) {
            return Double.MAX_VALUE;
        }
        
        if (approx.getNumRows() != exact.getNumRows() || 
            approx.getNumCols() != exact.getNumCols()) {
            return Double.MAX_VALUE;
        }
        
        double maxError = 0.0;
        
        // Process each column (class)
        for (int j = 0; j < exact.getNumCols(); j++) {
            double colMean = 0.0;
            double colSum = 0.0;
            
            // Calculate mean absolute difference and sum for this column
            for (int i = 0; i < exact.getNumRows(); i++) {
                double diff = Math.abs(exact.get(i, j) - approx.get(i, j));
                colMean += diff;
                colSum += exact.get(i, j);
            }
            
            // Compute mean for this column
            colMean = colMean / exact.getNumRows();
            
            // Compute relative error for this column
            if (colSum > 0) {
                double colError = colMean / colSum;
                maxError = Math.max(maxError, colError);
            }
        }
        
        return maxError;
    }
    
    /**
     * Get indices of non-delay stations
     */
    public static List<Integer> getNonDelayStationIndices(Network model) {
        List<Integer> indices = new ArrayList<>();
        List<Node> stations = model.getNodes();
        
        for (int i = 0; i < stations.size(); i++) {
            if (!(stations.get(i) instanceof Delay)) {
                indices.add(i);
            }
        }
        return indices;
    }
    
    /**
     * Generate random gallery matrix (similar to MATLAB's randgallery)
     * This creates a random matrix with values in [0,1]
     */
    public static Matrix randGallery(int rows, int cols, int seed) {
        java.util.Random rand = new java.util.Random(seed);
        Matrix result = new Matrix(rows, cols);
        
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                // Generate random values similar to MATLAB's randgallery
                result.set(i, j, 0.1 + rand.nextDouble() * 0.9);
            }
        }
        return result;
    }
    
    /**
     * Calculate Mean Absolute Percentage Error (MAPE)
     * Wrapper for Utils.mape for benchmark compatibility
     */
    public static double mape(Matrix approx, Matrix exact) {
        return Utils.mape(approx, exact);
    }
    
    /**
     * Calculate utilization error for benchmarks
     * Excludes delay stations (row 0) from comparison
     */
    public static double utilizationError(Matrix approx, Matrix exact, Network model) {
        if (approx == null || exact == null) {
            return Double.MAX_VALUE;
        }
        
        // Get non-delay station indices
        List<Integer> nonDelayIndices = getNonDelayStationIndices(model);
        if (nonDelayIndices.isEmpty()) {
            return 0.0;
        }
        
        // For utilization, we typically look at specific rows (non-delay stations)
        double maxError = 0.0;
        for (int idx : nonDelayIndices) {
            if (idx < approx.getNumRows() && idx < exact.getNumRows()) {
                for (int j = 0; j < approx.getNumCols() && j < exact.getNumCols(); j++) {
                    double err = Math.abs(approx.get(idx, j) - exact.get(idx, j));
                    maxError = Math.max(maxError, err);
                }
            }
        }
        return maxError;
    }
}