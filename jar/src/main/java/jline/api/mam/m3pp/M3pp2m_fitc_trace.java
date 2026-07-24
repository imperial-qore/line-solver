/**
 * @file M3PP trace-based counting process fitting
 *
 * Fits a M3PP(2,m) from trace data using counting process characteristics.
 * Supports multiple fitting methods: exact_delta, approx_delta, approx_cov, approx_ag.
 *
 * @since LINE 3.0
 */
package jline.api.mam.m3pp;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.util.matrix.Matrix;

public final class M3pp2m_fitc_trace {
    private M3pp2m_fitc_trace() {}

    /** Helper to check if verbose output is enabled. */
    private static boolean isVerbose() {
        return GlobalConstants.getVerbose() != VerboseLevel.SILENT;
    }

    /**
     * Fits a M3PP(2,m) from trace data using counting process characteristics.
     */
    public static Matrix[] m3pp2m_fitc_trace(double[] T, int[] A, String method, Double t1, Double tinf) {
        // Distinct sorted labels
        List<Integer> labels = new ArrayList<Integer>();
        for (int v : A) {
            if (!labels.contains(v)) labels.add(v);
        }
        Collections.sort(labels);
        int m = labels.size();

        // Validate method for 2-class restriction
        if ("approx_cov".equals(method) && m > 2) {
            throw new IllegalArgumentException("Approximate covariance fitting only supported for two classes.");
        }

        double[] TC = cumulativeSum(T);

        // Default time scales
        double tAvg = average(T);
        double t1Used = (t1 != null) ? t1 : (10 * tAvg);
        double tinfUsed = (tinf != null) ? tinf : Math.max(10 * t1Used, (TC[TC.length - 1] - TC[0]) / 100);
        double t2 = t1Used + tAvg;
        double t3 = tinfUsed; // Controls approximation accuracy

        if (isVerbose()) System.out.println("Computing counting process at resolution " + t1Used);
        double[][] mNt1 = mtrace_iat2counts(T, A, t1Used);
        double[][] mNt2 = mNt1; // Same resolution
        if (isVerbose()) System.out.println("Computing counting process at resolution " + tinfUsed);
        double[][] mNtinf = mtrace_iat2counts(T, A, tinfUsed);
        double[][] mNt3 = mNt1;

        double[] Nt1 = sumRows(mNt1);
        double[] Nt2 = sumRows(mNt2);
        double[] Ntinf = sumRows(mNtinf);

        // Total rate
        double a = 1.0 / tAvg;

        // Per-class rates
        double[] ai = new double[m];
        for (int i = 0; i < m; i++) {
            int classLabel = labels.get(i);
            int count = 0;
            for (int v : A) {
                if (v == classLabel) count++;
            }
            ai[i] = a * count / (double) A.length;
        }

        if (isVerbose()) System.out.println("Rate: " + a);

        // Joint-process characteristics
        double bt1 = variance(Nt1) / (a * t1Used);
        double bt2 = bt1;
        double binf = variance(Ntinf) / (a * tinfUsed);

        // Third centered moment
        double meanNt2 = average(Nt2);
        double meanNt2sq = 0.0;
        double meanNt2cub = 0.0;
        for (double v : Nt2) {
            meanNt2sq += v * v;
            meanNt2cub += v * v * v;
        }
        meanNt2sq /= Nt2.length;
        meanNt2cub /= Nt2.length;
        double m3t2 = meanNt2cub - 3 * meanNt2sq * meanNt2 + 2 * meanNt2 * meanNt2 * meanNt2;

        String methodLower = method.toLowerCase();
        if ("exact_delta".equals(methodLower)) {
            double[] dvt3 = computeTraceDelta(mNt3, m, labels);
            return M3pp2m_fitc.m3pp2m_fitc(a, bt1, bt2, binf, m3t2, t1Used, t2, ai, dvt3, t3);
        } else if ("approx_delta".equals(methodLower)) {
            double[] dvt3 = computeTraceDelta(mNt3, m, labels);
            return M3pp2m_fitc_approx.m3pp2m_fitc_approx(a, bt1, bt2, binf, m3t2, t1Used, t2, ai, dvt3, t3);
        } else if ("approx_cov".equals(methodLower)) {
            double totalVar = variance(sumRows(mNt3));
            double[] vi = new double[m];
            for (int i = 0; i < m; i++) {
                vi[i] = variance(getColumn(mNt3, i));
            }
            double viSum = 0.0;
            for (double v : vi) viSum += v;
            double s = 0.5 * (totalVar - viSum);
            return M3pp22_fitc_approx_cov.m3pp22_fitc_approx_cov(a, bt1, bt2, binf, m3t2, t1Used, t2, ai, s, t3);
        } else if ("approx_ag".equals(methodLower)) {
            double[] vt3 = new double[m];
            for (int i = 0; i < m; i++) {
                vt3[i] = variance(getColumn(mNt3, i));
            }
            double[] st3 = computeTraceCovariance(mNt3, m);
            double[] gt3 = new double[m];
            for (int i = 0; i < m; i++) {
                gt3[i] = vt3[i] + st3[i];
            }
            return M3pp2m_fitc_approx_ag.m3pp2m_fitc_approx_ag(a, bt1, bt2, binf, m3t2, t1Used, t2, ai, gt3, t3);
        } else {
            throw new IllegalArgumentException("Invalid method '" + method + "'");
        }
    }

    public static Matrix[] m3pp2m_fitc_trace(double[] T, int[] A, String method) {
        return m3pp2m_fitc_trace(T, A, method, null, null);
    }

    public static Matrix[] m3pp2m_fitc_trace(double[] T, int[] A) {
        return m3pp2m_fitc_trace(T, A, "approx_delta", null, null);
    }

    /**
     * Fits a M3PP(2,m) from trace data using Matrix inputs.
     */
    public static Matrix[] m3pp2m_fitc_trace(Matrix T, Matrix A, String method, Double t1, Double tinf) {
        double[] tArray = T.toArray1D();
        double[] aDouble = A.toArray1D();
        int[] aArray = new int[A.getNumRows() * A.getNumCols()];
        for (int i = 0; i < aArray.length; i++) {
            aArray[i] = (int) aDouble[i];
        }
        return m3pp2m_fitc_trace(tArray, aArray, method, t1, tinf);
    }

    public static Matrix[] m3pp2m_fitc_trace(Matrix T, Matrix A, String method) {
        return m3pp2m_fitc_trace(T, A, method, null, null);
    }

    public static Matrix[] m3pp2m_fitc_trace(Matrix T, Matrix A) {
        return m3pp2m_fitc_trace(T, A, "approx_delta", null, null);
    }

    /**
     * Computes per-class variance difference for delta fitting from trace.
     */
    private static double[] computeTraceDelta(double[][] mNt, int m, List<Integer> labels) {
        double[] dvt = new double[m];
        for (int i = 0; i < m; i++) {
            double[] classCol = getColumn(mNt, i);
            double[] otherSum = new double[mNt.length];
            for (int row = 0; row < mNt.length; row++) {
                double sum = 0.0;
                for (int j = 0; j < m; j++) {
                    if (j != i) sum += mNt[row][j];
                }
                otherSum[row] = sum;
            }
            dvt[i] = variance(classCol) - variance(otherSum);
        }
        return dvt;
    }

    /**
     * Computes per-class covariance with complementary classes from trace.
     */
    private static double[] computeTraceCovariance(double[][] mNt, int m) {
        double[] st = new double[m];
        for (int i = 0; i < m; i++) {
            double[] classCol = getColumn(mNt, i);
            double[] otherSum = new double[mNt.length];
            for (int row = 0; row < mNt.length; row++) {
                double sum = 0.0;
                for (int j = 0; j < m; j++) {
                    if (j != i) sum += mNt[row][j];
                }
                otherSum[row] = sum;
            }
            st[i] = covariance(classCol, otherSum);
        }
        return st;
    }

    private static double[] cumulativeSum(double[] arr) {
        double[] result = new double[arr.length];
        double sum = 0.0;
        for (int i = 0; i < arr.length; i++) {
            sum += arr[i];
            result[i] = sum;
        }
        return result;
    }

    private static double[] sumRows(double[][] matrix) {
        double[] result = new double[matrix.length];
        for (int row = 0; row < matrix.length; row++) {
            double sum = 0.0;
            for (double v : matrix[row]) sum += v;
            result[row] = sum;
        }
        return result;
    }

    private static double[] getColumn(double[][] matrix, int col) {
        double[] result = new double[matrix.length];
        for (int row = 0; row < matrix.length; row++) {
            result[row] = matrix[row][col];
        }
        return result;
    }

    private static double average(double[] arr) {
        if (arr.length == 0) return 0.0;
        double sum = 0.0;
        for (double v : arr) sum += v;
        return sum / arr.length;
    }

    private static double variance(double[] arr) {
        if (arr.length == 0) return 0.0;
        double mean = average(arr);
        double sum = 0.0;
        for (double v : arr) sum += (v - mean) * (v - mean);
        return sum / arr.length;
    }

    private static double covariance(double[] x, double[] y) {
        if (x.length != y.length || x.length == 0) return 0.0;
        double meanX = average(x);
        double meanY = average(y);
        double cov = 0.0;
        for (int i = 0; i < x.length; i++) {
            cov += (x[i] - meanX) * (y[i] - meanY);
        }
        return cov / x.length;
    }

    /**
     * Converts inter-arrival times to counts at given resolution.
     */
    public static double[][] mtrace_iat2counts(double[] T, int[] A, double scale) {
        // Distinct sorted labels
        List<Integer> labels = new ArrayList<Integer>();
        for (int v : A) {
            if (!labels.contains(v)) labels.add(v);
        }
        Collections.sort(labels);
        int m = labels.size();

        // Cumulative arrival times
        double[] cumT = cumulativeSum(T);
        double totalTime = cumT[cumT.length - 1];
        int numBins = Math.max(1, (int) (totalTime / scale));

        double[][] counts = new double[numBins][m];

        // Count arrivals in each bin for each class
        for (int i = 0; i < T.length; i++) {
            double arrivalTime = cumT[i];
            int binIndex = Math.min((int) (arrivalTime / scale), numBins - 1);
            int classIndex = labels.indexOf(A[i]);
            if (classIndex >= 0 && binIndex >= 0) {
                counts[binIndex][classIndex]++;
            }
        }

        return counts;
    }
}
