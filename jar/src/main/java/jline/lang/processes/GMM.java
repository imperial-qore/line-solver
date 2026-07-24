/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import static jline.GlobalConstants.Inf;
import static jline.GlobalConstants.NegInf;

import jline.GlobalConstants;
import jline.util.Pair;
import jline.util.RandomManager;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import static jline.util.Maths.max;
import static jline.util.Maths.min;
import static org.apache.commons.math3.util.FastMath.*;
import static org.apache.commons.math3.special.Erf.erf;

/**
 * A Gaussian Mixture Model (GMM) distribution.
 * 
 * This class represents a mixture of normal distributions, where each component
 * has a weight, mean, and standard deviation.
 */
public class GMM extends ContinuousDistribution implements Serializable {
    
    private double[] weights;
    private double[] means;
    private double[] sigmas;
    private int nComponents;
    
    /**
     * Creates a GMM with specified components.
     * 
     * @param weights array of component weights (must sum to 1)
     * @param means array of component means
     * @param sigmas array of component standard deviations
     */
    public GMM(double[] weights, double[] means, double[] sigmas) {
        super("GMM", 3, new Pair<Double, Double>(NegInf, Inf));
        
        if (weights.length != means.length || weights.length != sigmas.length) {
            throw new IllegalArgumentException("All arrays must have the same length");
        }
        
        // Normalize weights
        double sum = 0;
        for (double w : weights) {
            sum += w;
        }
        
        this.nComponents = weights.length;
        this.weights = new double[nComponents];
        this.means = new double[nComponents];
        this.sigmas = new double[nComponents];
        
        for (int i = 0; i < nComponents; i++) {
            this.weights[i] = weights[i] / sum;
            this.means[i] = means[i];
            this.sigmas[i] = sigmas[i];
        }
        
        this.setParam(1, "weights", this.weights);
        this.setParam(2, "means", this.means);
        this.setParam(3, "sigmas", this.sigmas);
        
        this.mean = calculateMean();
        this.immediate = false;
    }
    
    /**
     * Creates a GMM from a simplified matrix representation.
     * Each row contains [weight, mean, sigma].
     */
    public static GMM fromMatrix(Matrix sgmm) {
        int n = sgmm.getNumRows();
        double[] weights = new double[n];
        double[] means = new double[n];
        double[] sigmas = new double[n];
        
        for (int i = 0; i < n; i++) {
            weights[i] = sgmm.get(i, 0);
            means[i] = sgmm.get(i, 1);
            sigmas[i] = sqrt(sgmm.get(i, 2)); // Convert variance to std dev
        }
        
        return new GMM(weights, means, sigmas);
    }
    
    /**
     * Converts this GMM to a simplified matrix representation.
     * Each row contains [weight, mean, variance].
     */
    public Matrix toMatrix() {
        Matrix sgmm = new Matrix(nComponents, 3);
        
        for (int i = 0; i < nComponents; i++) {
            sgmm.set(i, 0, weights[i]);
            sgmm.set(i, 1, means[i]);
            sgmm.set(i, 2, sigmas[i] * sigmas[i]); // Convert std dev to variance
        }
        
        return sgmm;
    }
    
    private double calculateMean() {
        double mean = 0;
        for (int i = 0; i < nComponents; i++) {
            mean += weights[i] * means[i];
        }
        return mean;
    }
    
    private double calculateSecondMoment() {
        double e2 = 0;
        for (int i = 0; i < nComponents; i++) {
            double var_i = sigmas[i] * sigmas[i];
            e2 += weights[i] * (var_i + means[i] * means[i]);
        }
        return e2;
    }
    
    @Override
    public double getMean() {
        return calculateMean();
    }
    
    @Override
    public double getVar() {
        double e1 = calculateMean();
        double e2 = calculateSecondMoment();
        return e2 - e1 * e1;
    }
    
    @Override
    public double getSCV() {
        double mean = getMean();
        if (abs(mean) < GlobalConstants.FineTol) {
            return Inf;
        }
        return getVar() / (mean * mean);
    }
    
    @Override
    public double getSkewness() {
        double mean = getMean();
        double var = getVar();
        double std = sqrt(var);
        
        if (std < GlobalConstants.FineTol) {
            return 0;
        }
        
        // Calculate third central moment
        double e3 = 0;
        for (int i = 0; i < nComponents; i++) {
            double mu_i = means[i];
            double sigma_i = sigmas[i];
            double diff = mu_i - mean;
            
            // E[(X_i - mean)^3] = (mu_i - mean)^3 + 3(mu_i - mean)sigma_i^2
            e3 += weights[i] * (pow(diff, 3) + 3 * diff * sigma_i * sigma_i);
        }
        
        return e3 / pow(std, 3);
    }
    
    @Override
    public double evalCDF(double x) {
        double cdf = 0;
        for (int i = 0; i < nComponents; i++) {
            cdf += weights[i] * normalCDF(x, means[i], sigmas[i]);
        }
        return cdf;
    }
    
    public double evalPDF(double x) {
        double pdf = 0;
        for (int i = 0; i < nComponents; i++) {
            pdf += weights[i] * normalPDF(x, means[i], sigmas[i]);
        }
        return pdf;
    }
    
    @Override
    public double evalLST(double s) {
        double lst = 0;
        for (int i = 0; i < nComponents; i++) {
            lst += weights[i] * exp(means[i] * s + 0.5 * sigmas[i] * sigmas[i] * s * s);
        }
        return lst;
    }
    
    @Override
    public double[] sample(int n) {
        return this.sample(n, RandomManager.getThreadRandomAsRandom());
    }
    
    @Override
    public double[] sample(int n, Random random) {
        double[] samples = new double[n];
        
        for (int i = 0; i < n; i++) {
            // Select component based on weights
            double u = random.nextDouble();
            double cumSum = 0;
            int component = 0;
            
            for (int j = 0; j < nComponents; j++) {
                cumSum += weights[j];
                if (u <= cumSum) {
                    component = j;
                    break;
                }
            }
            
            // Sample from selected component
            samples[i] = means[component] + sigmas[component] * random.nextGaussian();
        }
        
        return samples;
    }
    
    /**
     * Creates a new GMM by convolving this GMM with another GMM.
     * Convolution corresponds to the sum of two independent random variables.
     */
    public GMM convolve(GMM other) {
        int n1 = this.nComponents;
        int n2 = other.nComponents;
        int nNew = n1 * n2;
        
        double[] newWeights = new double[nNew];
        double[] newMeans = new double[nNew];
        double[] newSigmas = new double[nNew];
        
        int idx = 0;
        for (int i = 0; i < n1; i++) {
            for (int j = 0; j < n2; j++) {
                newWeights[idx] = this.weights[i] * other.weights[j];
                newMeans[idx] = this.means[i] + other.means[j];
                // Variance of sum is sum of variances for independent variables
                double var1 = this.sigmas[i] * this.sigmas[i];
                double var2 = other.sigmas[j] * other.sigmas[j];
                newSigmas[idx] = sqrt(var1 + var2);
                idx++;
            }
        }
        
        return new GMM(newWeights, newMeans, newSigmas);
    }
    
    /**
     * Creates a mixture of two GMMs with specified weights.
     */
    public static GMM mixture(GMM gmm1, GMM gmm2, double p1, double p2) {
        double sum = p1 + p2;
        p1 = p1 / sum;
        p2 = p2 / sum;
        
        int n1 = gmm1.nComponents;
        int n2 = gmm2.nComponents;
        int nNew = n1 + n2;
        
        double[] newWeights = new double[nNew];
        double[] newMeans = new double[nNew];
        double[] newSigmas = new double[nNew];
        
        // Copy components from first GMM
        for (int i = 0; i < n1; i++) {
            newWeights[i] = p1 * gmm1.weights[i];
            newMeans[i] = gmm1.means[i];
            newSigmas[i] = gmm1.sigmas[i];
        }
        
        // Copy components from second GMM
        for (int i = 0; i < n2; i++) {
            newWeights[n1 + i] = p2 * gmm2.weights[i];
            newMeans[n1 + i] = gmm2.means[i];
            newSigmas[n1 + i] = gmm2.sigmas[i];
        }
        
        return new GMM(newWeights, newMeans, newSigmas);
    }
    
    /**
     * Merges two components into one using moment matching.
     */
    public GMM merge(int i, int j) {
        if (i == j || i >= nComponents || j >= nComponents) {
            return this;
        }
        
        // Ensure i < j
        if (i > j) {
            int temp = i;
            i = j;
            j = temp;
        }
        
        // New component weight
        double pNew = weights[i] + weights[j];
        
        // New component mean (weighted average)
        double muNew = (weights[i] * means[i] + weights[j] * means[j]) / pNew;
        
        // New component variance (moment matching)
        double var1 = sigmas[i] * sigmas[i];
        double var2 = sigmas[j] * sigmas[j];
        double varNew = (weights[i] * (var1 + means[i] * means[i]) + 
                         weights[j] * (var2 + means[j] * means[j])) / pNew - muNew * muNew;
        double sigmaNew = sqrt(max(GlobalConstants.FineTol, varNew));
        
        // Create new arrays excluding component j and updating component i
        double[] newWeights = new double[nComponents - 1];
        double[] newMeans = new double[nComponents - 1];
        double[] newSigmas = new double[nComponents - 1];
        
        int idx = 0;
        for (int k = 0; k < nComponents; k++) {
            if (k == j) continue;
            
            if (k == i) {
                newWeights[idx] = pNew;
                newMeans[idx] = muNew;
                newSigmas[idx] = sigmaNew;
            } else {
                newWeights[idx] = weights[k];
                newMeans[idx] = means[k];
                newSigmas[idx] = sigmas[k];
            }
            idx++;
        }
        
        return new GMM(newWeights, newMeans, newSigmas);
    }
    
    /**
     * Automatically merges components to reduce the GMM to k components.
     */
    public GMM autoMerge(int k) {
        if (k >= nComponents) {
            return this;
        }
        
        GMM current = this;
        
        while (current.nComponents > k) {
            // Find two components with smallest combined weight
            int minI = 0, minJ = 1;
            double minWeight = current.weights[0] + current.weights[1];
            
            for (int i = 0; i < current.nComponents - 1; i++) {
                for (int j = i + 1; j < current.nComponents; j++) {
                    double combinedWeight = current.weights[i] + current.weights[j];
                    if (combinedWeight < minWeight) {
                        minWeight = combinedWeight;
                        minI = i;
                        minJ = j;
                    }
                }
            }
            
            current = current.merge(minI, minJ);
        }
        
        return current;
    }
    
    public int getNumberOfComponents() {
        return nComponents;
    }
    
    public double[] getWeights() {
        return weights.clone();
    }
    
    public double[] getMeans() {
        return means.clone();
    }
    
    public double[] getSigmas() {
        return sigmas.clone();
    }
    
    private static double normalPDF(double x, double mu, double sigma) {
        double diff = x - mu;
        return (1.0 / (sigma * sqrt(2 * PI))) * exp(-0.5 * (diff * diff) / (sigma * sigma));
    }
    
    private static double normalCDF(double x, double mu, double sigma) {
        return 0.5 * (1 + erf((x - mu) / (sigma * sqrt(2))));
    }
    
    @Override
    public MatrixCell getProcess() {
        // For GMM, return the matrix representation
        MatrixCell representation = new MatrixCell();
        representation.set(0, toMatrix());
        return representation;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder("jline.GMM(");
        sb.append(nComponents).append(" components: ");
        for (int i = 0; i < nComponents; i++) {
            if (i > 0) sb.append(", ");
            sb.append(String.format("[w=%.3f,μ=%.3f,σ=%.3f]", weights[i], means[i], sigmas[i]));
        }
        sb.append(")");
        return sb.toString();
    }
}