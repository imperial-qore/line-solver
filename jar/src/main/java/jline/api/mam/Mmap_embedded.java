/**
 * @file Marked Markovian Arrival Process embedded chain analysis
 *
 * Computes embedded discrete-time Markov chain for MMAP processes.
 * Essential for analyzing state transition probabilities at marked arrival epochs.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Mmap_embedded {
    private Mmap_embedded() {}

    /**
     * Computes the embedded chain of an MMAP.
     * The embedded chain represents the transition probabilities between states
     * when arrivals occur.
     *
     * @param mmap The MMAP process
     * @return Embedded chain transition matrix
     */
    public static Matrix mmap_embedded(MatrixCell mmap) {
        Matrix D0 = mmap.get(0);
        int order = D0.getNumRows();

        // Compute total arrival rates for each state
        double[] totalArrivalRates = new double[order];
        for (int i = 1; i < mmap.size(); i++) {
            Matrix Di = mmap.get(i);
            for (int row = 0; row < order; row++) {
                for (int col = 0; col < order; col++) {
                    totalArrivalRates[row] += Di.get(row, col);
                }
            }
        }

        // Create embedded transition matrix
        Matrix P = new Matrix(order, order);

        for (int i = 0; i < order; i++) {
            if (totalArrivalRates[i] > 1e-12) {
                // For each state with positive arrival rate
                for (int j = 0; j < order; j++) {
                    double transitionRate = 0.0;

                    // Sum arrival rates from state i to state j across all classes
                    for (int c = 1; c < mmap.size(); c++) {
                        Matrix Dc = mmap.get(c);
                        transitionRate += Dc.get(i, j);
                    }

                    // Normalize by total arrival rate to get probability
                    P.set(i, j, transitionRate / totalArrivalRates[i]);
                }
            } else {
                // If no arrivals from this state, use uniform distribution
                for (int j = 0; j < order; j++) {
                    P.set(i, j, 1.0 / order);
                }
            }
        }

        return P;
    }
}
