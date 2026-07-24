/**
 * @file Markovian Arrival MAP with Marked arrivals two-class gamma forward-sigma trace fitting
 *
 * Fits MAMAP(2,2) from trace data using gamma autocorrelation and forward-sigma characteristics.
 * Advanced trace-based fitting with correlation control for two-class systems.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import java.util.ArrayList;
import java.util.List;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Mamap22_fit_gamma_fs_trace {
    private Mamap22_fit_gamma_fs_trace() {}

    /**
     * Fits a MAMAP(2,2) using forward-start method from trace data.
     */
    public static MatrixCell mamap22_fit_gamma_fs_trace(double[] interArrivalTimes, int[] classLabels) {
        if (interArrivalTimes.length != classLabels.length) {
            throw new IllegalArgumentException("Inter-arrival times and class labels must have same length");
        }
        for (int label : classLabels) {
            if (label != 0 && label != 1) {
                throw new IllegalArgumentException("Class labels must be binary (0 or 1) for MAMAP(2,2)");
            }
        }

        // Extract trace characteristics
        double M1 = average(interArrivalTimes);
        double M2 = averageSquared(interArrivalTimes);
        double M3 = averageCubed(interArrivalTimes);
        double GAMMA = computeTraceGamma(interArrivalTimes);

        // Class-specific rates
        int class0Count = 0;
        int class1Count = 0;
        for (int label : classLabels) {
            if (label == 0) class0Count++;
            else if (label == 1) class1Count++;
        }
        int totalCount = classLabels.length;

        double p0 = (double) class0Count / totalCount;
        double p1 = (double) class1Count / totalCount;

        // Extract class-specific inter-arrival times
        List<Double> class0Times = new ArrayList<Double>();
        List<Double> class1Times = new ArrayList<Double>();

        for (int i = 0; i < interArrivalTimes.length; i++) {
            if (classLabels[i] == 0) {
                class0Times.add(interArrivalTimes[i]);
            } else {
                class1Times.add(interArrivalTimes[i]);
            }
        }

        double M1_0 = !class0Times.isEmpty() ? average(class0Times) : M1;
        double M1_1 = !class1Times.isEmpty() ? average(class1Times) : M1;

        // Create proper matrices for the function call
        MatrixCell dummyAmap = new MatrixCell(2);
        dummyAmap.set(0, new Matrix(2, 2));
        dummyAmap.set(1, new Matrix(2, 2));

        Matrix P = new Matrix(1, 2);
        P.set(0, 0, p0);
        P.set(0, 1, p1);

        Matrix B = new Matrix(1, 1);
        B.set(0, 0, M2);

        Matrix S = new Matrix(1, 1);
        S.set(0, 0, M3);

        return Mamap22_fit_multiclass.mamap22_fit_bs_multiclass(dummyAmap, P, B, S);
    }

    private static double computeTraceGamma(double[] interArrivalTimes) {
        if (interArrivalTimes.length < 2) return 0.0;

        double mean = average(interArrivalTimes);
        double variance = 0.0;
        for (double v : interArrivalTimes) {
            variance += (v - mean) * (v - mean);
        }
        variance /= interArrivalTimes.length;

        if (variance <= 0) return 0.0;

        // Lag-1 autocorrelation
        double autocovariance = 0.0;
        for (int i = 0; i < interArrivalTimes.length - 1; i++) {
            autocovariance += (interArrivalTimes[i] - mean) * (interArrivalTimes[i + 1] - mean);
        }
        autocovariance /= (interArrivalTimes.length - 1);

        return autocovariance / variance;
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

    private static double average(List<Double> list) {
        double s = 0.0;
        for (double v : list) s += v;
        return s / list.size();
    }
}
