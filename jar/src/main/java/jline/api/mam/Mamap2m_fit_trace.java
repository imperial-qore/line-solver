/**
 * @file Markovian Arrival MAP with Marked arrivals trace-based fitting
 *
 * Fits MAMAP(2,m) processes from empirical trace data with inter-arrival times and class labels.
 * Essential for data-driven modeling of real multiclass arrival processes.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import java.util.HashMap;
import java.util.Map;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Mamap2m_fit_trace {
    private Mamap2m_fit_trace() {}

    /**
     * Fits a MAMAP(2,m) to trace data.
     *
     * @param interArrivalTimes Array of inter-arrival times
     * @param classLabels Array of class labels for each arrival
     * @return Fitted MAMAP(2,m) model
     */
    public static MatrixCell mamap2m_fit_trace(double[] interArrivalTimes, int[] classLabels) {
        if (interArrivalTimes.length != classLabels.length) {
            throw new IllegalArgumentException("Inter-arrival times and class labels must have same length");
        }

        // Extract basic moments
        double M1 = average(interArrivalTimes);
        double M2 = averageSquared(interArrivalTimes);
        double M3 = averageCubed(interArrivalTimes);

        // Extract class information
        int maxLabel = 0;
        for (int label : classLabels) {
            if (label > maxLabel) maxLabel = label;
        }
        int numClasses = maxLabel + 1;
        int[] classCounts = new int[numClasses];
        for (int label : classLabels) {
            if (label >= 0 && label < numClasses) {
                classCounts[label]++;
            }
        }

        double[] classProbs = new double[numClasses];
        for (int i = 0; i < numClasses; i++) {
            classProbs[i] = (double) classCounts[i] / classLabels.length;
        }
        double[] backwardMoments = computeClassBackwardMoments(interArrivalTimes, classLabels, numClasses);

        // Create proper matrices for the function call
        MatrixCell dummyMap = new MatrixCell(2);
        dummyMap.set(0, new Matrix(2, 2));
        dummyMap.set(1, new Matrix(2, 2));

        double[] F = new double[]{M2}; // Forward moments
        double[] B = backwardMoments;  // Backward moments

        // Returns a 3-tuple whose first element is the MatrixCell.
        return Mamap2m_fit_fb_multiclass.mamap2m_fit_fb_multiclass(dummyMap, classProbs, F, B).mmap;
    }

    /**
     * Compute backward moments for each class
     */
    private static double[] computeClassBackwardMoments(double[] interArrivalTimes, int[] classLabels, int numClasses) {
        double[] backwardMoments = new double[numClasses];
        Map<Integer, Double> lastArrivalTime = new HashMap<Integer, Double>();
        double currentTime = 0.0;

        for (int i = 0; i < interArrivalTimes.length; i++) {
            currentTime += interArrivalTimes[i];
            int currentClass = classLabels[i];

            if (currentClass >= 0 && currentClass < numClasses) {
                if (lastArrivalTime.containsKey(currentClass)) {
                    double backwardTime = currentTime - lastArrivalTime.get(currentClass);
                    backwardMoments[currentClass] = backwardTime;
                } else {
                    backwardMoments[currentClass] = currentTime;
                }
                lastArrivalTime.put(currentClass, currentTime);
            }
        }

        // Fill in missing values with global mean
        double globalMean = average(interArrivalTimes);
        for (int i = 0; i < backwardMoments.length; i++) {
            if (backwardMoments[i] == 0.0) {
                backwardMoments[i] = globalMean;
            }
        }

        return backwardMoments;
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
}
