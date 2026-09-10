/**
 * @file Multi-class Absorbing Phase-type distribution multiclass fitting
 *
 * Fits MAPH(2,m) models to multiclass characteristics with class-specific parameters.
 * Specialized fitting for complex multiclass phase-type service time distributions.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Maph2m_fit_multiclass {
    private Maph2m_fit_multiclass() {}

    /**
     * Fits a multi-class MAPH(2,m) model to given multi-class characteristics.
     *
     * @param M1 First moment of inter-arrival times
     * @param M2 Second moment of inter-arrival times
     * @param M3 Third moment of inter-arrival times
     * @param classProbs Array of class probabilities
     * @param classRates Array of class-specific rates (may be null)
     * @param backwardMoments Array of backward moments for each class (may be null)
     * @return Fitted MAPH(2,m) model
     */
    public static MatrixCell maph2m_fit_multiclass(
            double M1,
            double M2,
            double M3,
            double[] classProbs,
            double[] classRates,
            double[] backwardMoments) {

        int m = classProbs.length;

        // Ensure class probabilities sum to 1
        double[] normalizedProbs = classProbs.clone();
        double probSum = 0.0;
        for (int i = 0; i < normalizedProbs.length; i++) {
            probSum += normalizedProbs[i];
        }
        if (Math.abs(probSum - 1.0) > 1e-6) {
            for (int i = 0; i < normalizedProbs.length; i++) {
                normalizedProbs[i] /= probSum;
            }
        }

        // If class rates not provided, distribute evenly based on probabilities
        double[] rates;
        if (classRates != null) {
            rates = classRates;
        } else {
            rates = new double[normalizedProbs.length];
            for (int i = 0; i < normalizedProbs.length; i++) {
                rates[i] = normalizedProbs[i] / M1;
            }
        }

        // If backward moments not provided, use uniform distribution
        double[] backMoments;
        if (backwardMoments != null) {
            backMoments = backwardMoments;
        } else {
            backMoments = new double[m];
            for (int i = 0; i < m; i++) {
                backMoments[i] = M1;
            }
        }

        // Convert arrays to matrices
        Matrix probMatrix = new Matrix(1, normalizedProbs.length);
        for (int i = 0; i < normalizedProbs.length; i++) {
            probMatrix.set(0, i, normalizedProbs[i]);
        }

        Matrix backMatrix = new Matrix(backMoments.length, 1);
        for (int i = 0; i < backMoments.length; i++) {
            backMatrix.set(i, 0, backMoments[i]);
        }

        return Maph2m_fit.maph2m_fit(M1, M2, M3, probMatrix, backMatrix);
    }

    public static MatrixCell maph2m_fit_multiclass(double M1, double M2, double M3, double[] classProbs) {
        return maph2m_fit_multiclass(M1, M2, M3, classProbs, null, null);
    }

    public static MatrixCell maph2m_fit_multiclass(double M1, double M2, double M3,
                                                   double[] classProbs, double[] classRates) {
        return maph2m_fit_multiclass(M1, M2, M3, classProbs, classRates, null);
    }
}
