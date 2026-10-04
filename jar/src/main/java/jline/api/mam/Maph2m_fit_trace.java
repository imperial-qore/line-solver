/**
 * @file Multi-class Absorbing Phase-type distribution trace-based fitting
 *
 * Fits MAPH(2,m) from empirical trace data for multiclass service time modeling.
 * Essential for data-driven phase-type distribution modeling from real measurements.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import java.util.ArrayList;
import java.util.List;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Maph2m_fit_trace {
    private Maph2m_fit_trace() {}

    /**
     * Fits a multi-class MAPH(2,m) model to trace data with class labels.
     */
    public static MatrixCell maph2m_fit_trace(double[] interArrivalTimes, int[] classLabels) {
        if (interArrivalTimes.length != classLabels.length) {
            throw new IllegalArgumentException("Inter-arrival times and class labels must have the same length");
        }

        // Determine number of classes
        int maxClass = 0;
        for (int c : classLabels) {
            if (c > maxClass) maxClass = c;
        }
        int m = maxClass + 1;

        // Compute overall moments
        double M1 = average(interArrivalTimes);
        double M2 = averageSquared(interArrivalTimes);
        double M3 = averageCubed(interArrivalTimes);

        // Compute class probabilities
        int[] classCounts = new int[m];
        for (int label : classLabels) {
            if (label >= 0 && label < m) {
                classCounts[label]++;
            }
        }
        double[] classProbs = new double[m];
        for (int i = 0; i < m; i++) {
            classProbs[i] = (double) classCounts[i] / classLabels.length;
        }

        // Compute class-specific backward moments (conditional means)
        double[] classBackwardMoments = new double[m];
        for (int i = 0; i < m; i++) {
            List<Double> classTimesIndexed = new ArrayList<Double>();
            for (int idx = 0; idx < classLabels.length; idx++) {
                if (classLabels[idx] == i) {
                    classTimesIndexed.add(interArrivalTimes[idx]);
                }
            }
            classBackwardMoments[i] = !classTimesIndexed.isEmpty() ? averageList(classTimesIndexed) : M1;
        }

        // Convert arrays to matrices
        Matrix probMatrix = new Matrix(1, classProbs.length);
        for (int i = 0; i < classProbs.length; i++) {
            probMatrix.set(0, i, classProbs[i]);
        }

        Matrix backMatrix = new Matrix(classBackwardMoments.length, 1);
        for (int i = 0; i < classBackwardMoments.length; i++) {
            backMatrix.set(i, 0, classBackwardMoments[i]);
        }

        return Maph2m_fit.maph2m_fit(M1, M2, M3, probMatrix, backMatrix);
    }

    /**
     * Fits a multi-class MAPH(2,m) model to trace data with arrival time stamps and class labels.
     */
    public static MatrixCell maph2m_fit_trace_timestamps(double[] arrivalTimes, int[] classLabels) {
        if (arrivalTimes.length != classLabels.length) {
            throw new IllegalArgumentException("Arrival times and class labels must have the same length");
        }
        if (arrivalTimes.length <= 1) {
            throw new IllegalArgumentException("Need at least 2 arrivals to compute inter-arrival times");
        }

        // Convert timestamps to inter-arrival times
        double[] interArrivalTimes = new double[arrivalTimes.length - 1];
        for (int i = 1; i < arrivalTimes.length; i++) {
            interArrivalTimes[i - 1] = arrivalTimes[i] - arrivalTimes[i - 1];
        }

        // Use the first n-1 class labels (since we have n-1 inter-arrival times)
        int[] adjustedClassLabels = new int[classLabels.length - 1];
        System.arraycopy(classLabels, 1, adjustedClassLabels, 0, classLabels.length - 1);

        return maph2m_fit_trace(interArrivalTimes, adjustedClassLabels);
    }

    private static double average(double[] arr) {
        double s = 0.0;
        for (double v : arr) s += v;
        return s / arr.length;
    }

    private static double averageSquared(double[] arr) {
        double s = 0.0;
        for (double v : arr) s += v * v;
        return s / arr.length;
    }

    private static double averageCubed(double[] arr) {
        double s = 0.0;
        for (double v : arr) s += v * v * v;
        return s / arr.length;
    }

    private static double averageList(List<Double> list) {
        double s = 0.0;
        for (double v : list) s += v;
        return s / list.size();
    }
}
