/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.util;

import jline.util.matrix.Matrix;
import org.apache.commons.math3.util.FastMath;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

/**
 * Miscellaneous utilities
 */
public class Utils {

    /**
     * Finds the key corresponding to a target string value in a map.
     * 
     * @param map the map to search
     * @param target the string value to find
     * @return the key associated with the target value, or 0 if not found
     */
    public static int findString(Map<Integer, String> map, String target) {
        int res = 0;
        for (Map.Entry e : map.entrySet()) {
            if (e.getValue().equals(target)) res = (int) e.getKey();
        }
        return res;
    }

    /**
     * Checks if a value represents infinity.
     * Considers Double infinity and integer limit values as infinite.
     * 
     * @param val the value to check
     * @return true if the value represents infinity, false otherwise
     */
    public static boolean isInf(double val) {
        return Double.isInfinite(val) || val == Integer.MAX_VALUE || val == Integer.MIN_VALUE;
    }

    /**
     * Return mean absolute percentage error of approx with respect to exact
     *
     * @param approx The approximate values matrix
     * @param exact  The exact values matrix
     * @return The mean absolute percentage error as a percentage
     */
    public static double mape(Matrix approx, Matrix exact) {
        return mapeWithNanMean(approx, exact).mape;
    }

    /**
     * Polymorphic version that returns both MAPE and nanMean.
     * Calculates mean absolute percentage error and handles NaN/zero values appropriately.
     *
     * @param approx The approximate values matrix
     * @param exact  The exact values matrix
     * @return MapeResult containing MAPE, nanMean, and count of valid comparisons
     */
    public static MapeResult mapeWithNanMean(Matrix approx, Matrix exact) {
        int numRows = approx.getNumRows();
        double totalAbsolutePercentageError = 0;
        double totalApproxValues = 0;
        int numExactGreaterThanZero = 0;
        int numValidApprox = 0;

        for (int row = 0; row < numRows; row++) {
            double exactVal = exact.get(row, 0);
            double approxVal = approx.get(row, 0);

            // Count valid approximate values for nanMean calculation
            if (!Double.isNaN(approxVal) && !Double.isInfinite(approxVal)) {
                totalApproxValues += approxVal;
                numValidApprox++;
            }

            // Calculate MAPE only for non-zero exact values
            if (exactVal > 0) {
                if (!Double.isNaN(approxVal) && !Double.isInfinite(approxVal)) {
                    totalAbsolutePercentageError += FastMath.abs(1 - (approxVal / exactVal));
                    numExactGreaterThanZero++;
                }
            }
        }

        double mape = numExactGreaterThanZero > 0 ?
                totalAbsolutePercentageError / numExactGreaterThanZero : Double.NaN;
        double nanMean = numValidApprox > 0 ?
                totalApproxValues / numValidApprox : Double.NaN;

        return new MapeResult(mape, nanMean, numExactGreaterThanZero);
    }

    public static <T extends Object> List<T> unique(List<T> list) {
        return new ArrayList<T>(new HashSet<>(list));
    }

    /**
     * Convert double array to string representation
     * @param array The double array to convert
     * @return String representation in format [1.0,2.0,3.0]
     */
    public static String doubleArrayToString(double[] array) {
        if (array == null || array.length == 0) return "[]";
        
        StringBuilder sb = new StringBuilder("[");
        for (int i = 0; i < array.length; i++) {
            if (i > 0) sb.append(",");
            sb.append(String.format("%.15e", array[i]));
        }
        sb.append("]");
        return sb.toString();
    }

    /**
     * Parse integer array from string representation
     * @param str String representation in format [1,2,3] or 1,2,3
     * @return Parsed integer array
     */
    public static int[] parseIntArray(String str) {
        if (str.startsWith("[") && str.endsWith("]")) {
            str = str.substring(1, str.length() - 1);
        }
        String[] parts = str.split(",");
        int[] result = new int[parts.length];
        for (int i = 0; i < parts.length; i++) {
            result[i] = Integer.parseInt(parts[i].trim());
        }
        return result;
    }

    /**
     * Container class for MAPE calculation results
     */
    public static class MapeResult {
        public final double mape;
        public final double nanMean;
        public final int validCount;

        public MapeResult(double mape, double nanMean, int validCount) {
            this.mape = mape;
            this.nanMean = nanMean;
            this.validCount = validCount;
        }
    }
}
