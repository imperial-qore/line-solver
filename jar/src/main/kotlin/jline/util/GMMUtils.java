/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.util;

import jline.lang.processes.GMM;
import jline.util.matrix.Matrix;

import static org.apache.commons.math3.util.FastMath.*;

/**
 * Utility functions for Gaussian Mixture Model operations.
 * 
 * This class provides static methods that correspond to the MATLAB gmm_* functions.
 */
public class GMMUtils {
    
    /**
     * Calculates the minimum value with high probability (0.001 percentile).
     */
    public static double gmmMin(GMM gmm) {
        return gmmPercentile(gmm, 0.001);
    }
    
    /**
     * Calculates the maximum value with high probability (99.999 percentile).
     */
    public static double gmmMax(GMM gmm) {
        return gmmPercentile(gmm, 99.999);
    }
    
    /**
     * Calculates the median (50th percentile).
     */
    public static double gmmMedian(GMM gmm) {
        return gmmPercentile(gmm, 50.0);
    }
    
    /**
     * Calculates the specified percentile of the GMM distribution.
     * 
     * @param gmm the GMM distribution
     * @param p the percentile (0-100)
     * @return the value at the specified percentile
     */
    public static double gmmPercentile(GMM gmm, double p) {
        if (p < 0 || p > 100) {
            throw new IllegalArgumentException("Percentile must be between 0 and 100");
        }
        
        double target = p / 100.0;
        
        // Use bisection method to find the percentile
        // Start with bounds based on mean and standard deviation
        double mean = gmm.getMean();
        double std = sqrt(gmm.getVar());
        
        double lower = mean - 10 * std;
        double upper = mean + 10 * std;
        double tol = 1e-6;
        
        // Expand bounds if necessary
        while (gmm.evalCDF(lower) > target) {
            lower = lower - 10 * std;
        }
        while (gmm.evalCDF(upper) < target) {
            upper = upper + 10 * std;
        }
        
        // Bisection method
        while (upper - lower > tol) {
            double mid = (lower + upper) / 2;
            double cdf = gmm.evalCDF(mid);
            
            if (cdf < target) {
                lower = mid;
            } else {
                upper = mid;
            }
        }
        
        return (lower + upper) / 2;
    }
    
    /**
     * Calculates the mean of the minimum of k independent samples from the GMM.
     * This is an approximation using extreme value theory.
     */
    public static double gmmMeanMin(GMM gmm, int k) {
        if (k <= 0) {
            throw new IllegalArgumentException("k must be positive");
        }
        
        // For k=1, return the mean
        if (k == 1) {
            return gmm.getMean();
        }
        
        // Use approximation based on order statistics
        // For large k, the minimum approaches the lower bound
        // This is a simplified approximation
        double mean = gmm.getMean();
        double std = sqrt(gmm.getVar());
        
        // Approximate using Gumbel distribution for minima
        double gamma = 0.5772156649; // Euler-Mascheroni constant
        double scale = std / sqrt(2 * log(k));
        double location = mean - std * sqrt(2 * log(k)) + scale * gamma;
        
        return location;
    }
    
    /**
     * Calculates the mean of the maximum of k independent samples from the GMM.
     * This is an approximation using extreme value theory.
     */
    public static double gmmMeanMax(GMM gmm, int k) {
        if (k <= 0) {
            throw new IllegalArgumentException("k must be positive");
        }
        
        // For k=1, return the mean
        if (k == 1) {
            return gmm.getMean();
        }
        
        // Use approximation based on order statistics
        // For large k, the maximum approaches the upper bound
        // This is a simplified approximation
        double mean = gmm.getMean();
        double std = sqrt(gmm.getVar());
        
        // Approximate using Gumbel distribution for maxima
        double gamma = 0.5772156649; // Euler-Mascheroni constant
        double scale = std / sqrt(2 * log(k));
        double location = mean + std * sqrt(2 * log(k)) - scale * gamma;
        
        return location;
    }
    
    /**
     * Calculates the variance of the minimum of k independent samples from the GMM.
     */
    public static double gmmVarMin(GMM gmm, int k) {
        if (k <= 0) {
            throw new IllegalArgumentException("k must be positive");
        }
        
        // For k=1, return the variance
        if (k == 1) {
            return gmm.getVar();
        }
        
        // Approximate variance of minimum using extreme value theory
        double std = sqrt(gmm.getVar());
        double scale = std / sqrt(2 * log(k));
        
        // Variance of Gumbel distribution
        return (PI * PI / 6) * scale * scale;
    }
    
    /**
     * Calculates the variance of the maximum of k independent samples from the GMM.
     */
    public static double gmmVarMax(GMM gmm, int k) {
        if (k <= 0) {
            throw new IllegalArgumentException("k must be positive");
        }
        
        // For k=1, return the variance
        if (k == 1) {
            return gmm.getVar();
        }
        
        // Approximate variance of maximum using extreme value theory
        double std = sqrt(gmm.getVar());
        double scale = std / sqrt(2 * log(k));
        
        // Variance of Gumbel distribution
        return (PI * PI / 6) * scale * scale;
    }
    
    /**
     * Generates random samples from a simplified GMM representation.
     * 
     * @param sgmm Matrix where each row is [weight, mean, variance]
     * @param n Number of samples to generate
     * @return Array of random samples
     */
    public static double[] sgmmRand(Matrix sgmm, int n) {
        GMM gmm = GMM.fromMatrix(sgmm);
        return gmm.sample(n);
    }
    
    /**
     * Calculates the mean of a simplified GMM directly from matrix representation.
     * 
     * @param sgmm Matrix where each row is [weight, mean, variance]
     * @return The weighted mean
     */
    public static double sgmmMean(Matrix sgmm) {
        double mean = 0;
        for (int i = 0; i < sgmm.getNumRows(); i++) {
            mean += sgmm.get(i, 0) * sgmm.get(i, 1);  // weight * mean
        }
        return mean;
    }
    
    /**
     * Calculates the first moment (E[X]) of a simplified GMM.
     * 
     * @param sgmm Matrix where each row is [weight, mean, variance]
     * @return The first moment (mean)
     */
    public static double sgmmE1(Matrix sgmm) {
        return sgmmMean(sgmm);
    }
    
    /**
     * Calculates the second moment (E[X²]) of a simplified GMM.
     * 
     * @param sgmm Matrix where each row is [weight, mean, variance]
     * @return The second moment
     */
    public static double sgmmE2(Matrix sgmm) {
        double e2 = 0;
        for (int i = 0; i < sgmm.getNumRows(); i++) {
            double weight = sgmm.get(i, 0);
            double mean = sgmm.get(i, 1);
            double variance = sgmm.get(i, 2);
            e2 += weight * (variance + mean * mean);
        }
        return e2;
    }
    
    /**
     * Calculates the variance of a simplified GMM.
     * 
     * @param sgmm Matrix where each row is [weight, mean, variance]
     * @return The variance
     */
    public static double sgmmVar(Matrix sgmm) {
        double e1 = sgmmE1(sgmm);
        double e2 = sgmmE2(sgmm);
        return e2 - e1 * e1;
    }
    
    /**
     * Calculates the standard deviation of a simplified GMM.
     * 
     * @param sgmm Matrix where each row is [weight, mean, variance]
     * @return The standard deviation
     */
    public static double sgmmStd(Matrix sgmm) {
        return sqrt(sgmmVar(sgmm));
    }
    
    /**
     * Convolves two simplified GMMs (sum of independent random variables).
     * 
     * @param sgmm1 First GMM matrix
     * @param sgmm2 Second GMM matrix
     * @return Convolved GMM matrix
     */
    public static Matrix sgmmConvolve(Matrix sgmm1, Matrix sgmm2) {
        int n1 = sgmm1.getNumRows();
        int n2 = sgmm2.getNumRows();
        Matrix result = new Matrix(n1 * n2, 3);
        
        int idx = 0;
        for (int i = 0; i < n1; i++) {
            for (int j = 0; j < n2; j++) {
                result.set(idx, 0, sgmm1.get(i, 0) * sgmm2.get(j, 0));  // weight
                result.set(idx, 1, sgmm1.get(i, 1) + sgmm2.get(j, 1));  // mean
                result.set(idx, 2, sgmm1.get(i, 2) + sgmm2.get(j, 2));  // variance
                idx++;
            }
        }
        return result;
    }
    
    /**
     * Creates a mixture of two simplified GMMs.
     * 
     * @param sgmm1 First GMM matrix
     * @param sgmm2 Second GMM matrix  
     * @param p1 Weight for first GMM
     * @param p2 Weight for second GMM
     * @return Mixed GMM matrix
     */
    public static Matrix sgmmMixture(Matrix sgmm1, Matrix sgmm2, double p1, double p2) {
        double sum = p1 + p2;
        p1 = p1 / sum;
        p2 = p2 / sum;
        
        int n1 = sgmm1.getNumRows();
        int n2 = sgmm2.getNumRows();
        Matrix result = new Matrix(n1 + n2, 3);
        
        // Copy first GMM with scaled weights
        for (int i = 0; i < n1; i++) {
            result.set(i, 0, p1 * sgmm1.get(i, 0));
            result.set(i, 1, sgmm1.get(i, 1));
            result.set(i, 2, sgmm1.get(i, 2));
        }
        
        // Copy second GMM with scaled weights
        for (int i = 0; i < n2; i++) {
            result.set(n1 + i, 0, p2 * sgmm2.get(i, 0));
            result.set(n1 + i, 1, sgmm2.get(i, 1));
            result.set(n1 + i, 2, sgmm2.get(i, 2));
        }
        
        return result;
    }
}