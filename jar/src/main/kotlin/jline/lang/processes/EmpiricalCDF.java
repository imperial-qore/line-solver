/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import static jline.GlobalConstants.Inf;
import static jline.GlobalConstants.NegInf;

import jline.util.Pair;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

/**
 * Empirical CDF for a distribution
 */
public class EmpiricalCDF extends Distribution {
    protected Matrix data;

    /**
     * Creates an EmpiricalCDF with the given data
     * @param xdata the data matrix with CDF and values
     */
    public EmpiricalCDF(Matrix xdata) {
        super("EmpiricalCdf", 2, new Pair<Double, Double>(NegInf, Inf));
        this.data = xdata;
    }

    /**
     * Creates an EmpiricalCDF with separate CDF and value data
     * @param cdfdata the CDF values
     * @param xdata the corresponding x values
     */
    public EmpiricalCDF(Matrix cdfdata, Matrix xdata) {
        super("EmpiricalCdf", 2, new Pair<Double, Double>(NegInf, Inf));
        // Combine cdfdata and xdata into unique rows
        this.data = combineAndUnique(cdfdata, xdata);
    }

    /**
     * Combine CDF and x data into unique rows
     * @param cdfdata the CDF values
     * @param xdata the x values
     * @return combined unique matrix
     */
    private Matrix combineAndUnique(Matrix cdfdata, Matrix xdata) {
        List<double[]> combined = new ArrayList<double[]>();
        
        for (int i = 0; i < Math.min(cdfdata.getNumRows(), xdata.getNumRows()); i++) {
            combined.add(new double[]{cdfdata.get(i, 0), xdata.get(i, 0)});
        }
        
        // Sort by CDF values and remove duplicates
        combined.sort((a, b) -> Double.compare(a[0], b[0]));
        
        List<double[]> unique = new ArrayList<double[]>();
        double lastCdf = NegInf;
        for (double[] row : combined) {
            if (row[0] != lastCdf) {
                unique.add(row);
                lastCdf = row[0];
            }
        }
        
        Matrix result = new Matrix(unique.size(), 2);
        for (int i = 0; i < unique.size(); i++) {
            result.set(i, 0, unique.get(i)[0]);
            result.set(i, 1, unique.get(i)[1]);
        }
        
        return result;
    }

    @Override
    public double[] sample(int n, Random random) {
        double[] Y = new double[n];
        for (int i = 0; i < n; i++) {
            Y[i] = random.nextDouble();
        }
        
        // Find indices where F(x) < 1 to remove duplicates
        List<Integer> xset = new ArrayList<Integer>();
        for (int i = 0; i < data.getNumRows(); i++) {
            if (data.get(i, 0) < 1.0) {
                xset.add(i);
            }
        }
        
        double[] result = new double[n];
        for (int i = 0; i < n; i++) {
            result[i] = interpolate(Y[i], xset);
        }
        
        return result;
    }

    /**
     * Linear interpolation for sampling
     * @param y the random value to interpolate
     * @param xset valid indices for interpolation
     * @return interpolated value
     */
    private double interpolate(double y, List<Integer> xset) {
        if (xset.isEmpty()) {
            return 0.0;
        }
        
        // Simple linear interpolation
        for (int i = 0; i < xset.size() - 1; i++) {
            int idx1 = xset.get(i);
            int idx2 = xset.get(i + 1);
            double f1 = data.get(idx1, 0);
            double f2 = data.get(idx2, 0);
            
            if (y >= f1 && y <= f2) {
                double x1 = data.get(idx1, 1);
                double x2 = data.get(idx2, 1);
                return x1 + (x2 - x1) * (y - f1) / (f2 - f1);
            }
        }
        
        // Extrapolation for edge cases
        if (y <= data.get(xset.get(0), 0)) {
            return data.get(xset.get(0), 1);
        } else {
            int lastIdx = xset.get(xset.size() - 1);
            return data.get(lastIdx, 1);
        }
    }

    @Override
    public double evalCDF(double t) {
        // Find indices where F(x) < 1 to remove duplicates
        List<Integer> xset = new ArrayList<Integer>();
        for (int i = 0; i < data.getNumRows(); i++) {
            if (data.get(i, 0) < 1.0) {
                xset.add(i);
            }
        }
        
        return interpolateCDF(t, xset);
    }

    /**
     * Interpolate CDF value at t
     * @param t the point to evaluate
     * @param xset valid indices for interpolation
     * @return CDF value at t
     */
    private double interpolateCDF(double t, List<Integer> xset) {
        if (xset.isEmpty()) {
            return 0.0;
        }
        
        // Linear interpolation for CDF
        for (int i = 0; i < xset.size() - 1; i++) {
            int idx1 = xset.get(i);
            int idx2 = xset.get(i + 1);
            double x1 = data.get(idx1, 1);
            double x2 = data.get(idx2, 1);
            
            if (t >= x1 && t <= x2) {
                double f1 = data.get(idx1, 0);
                double f2 = data.get(idx2, 0);
                return f1 + (f2 - f1) * (t - x1) / (x2 - x1);
            }
        }
        
        // Extrapolation for edge cases
        if (t <= data.get(xset.get(0), 1)) {
            return data.get(xset.get(0), 0);
        } else {
            int lastIdx = xset.get(xset.size() - 1);
            return data.get(lastIdx, 0);
        }
    }

    @Override
    public double getMean() {
        double[] moments = getMoments();
        return moments[0];
    }

    @Override
    public double getSCV() {
        double[] moments = getMoments();
        return moments[3];
    }

    @Override
    public double getSkewness() {
        double[] moments = getMoments();
        return moments[4];
    }

    /**
     * Calculate the first three moments, SCV, and skewness
     * @return array containing [m1, m2, m3, SCV, skewness]
     */
    public double[] getMoments() {
        double m1 = 0.0; // first moment
        double m2 = 0.0; // second moment  
        double m3 = 0.0; // third moment
        
        for (int i = 0; i < data.getNumRows() - 1; i++) {
            double x = (data.get(i + 1, 1) - data.get(i, 1)) / 2.0 + data.get(i, 1);
            double bin = data.get(i + 1, 0) - data.get(i, 0);
            
            m1 += x * bin;
            m2 += x * x * bin;
            m3 += x * x * x * bin;
        }
        
        double scv = (m2 / (m1 * m1)) - 1.0;
        double variance = m2 - m1 * m1;
        double skew = (m3 - 3 * m1 * variance - m1 * m1 * m1) / Math.pow(variance, 1.5);
        
        return new double[]{m1, m2, m3, scv, skew};
    }

    /**
     * Evaluate the Laplace-Stieltjes Transform at s
     * @param s the point to evaluate
     * @return LST value
     */
    public double evalLST(double s) {
        double L = 0.0;
        
        for (int i = 0; i < data.getNumRows() - 1; i++) {
            double x = (data.get(i + 1, 1) - data.get(i, 1)) / 2.0 + data.get(i, 1);
            double bin = Math.exp(-s * x) * (data.get(i + 1, 0) - data.get(i, 0));
            L += bin;
        }
        
        return L;
    }

    /**
     * Get the empirical data
     * @return the data matrix
     */
    public Matrix getData() {
        return data;
    }
}