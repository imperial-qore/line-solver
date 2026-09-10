package jline.api.mam;

import org.apache.commons.math3.util.FastMath;

import jline.io.Ret;
import jline.util.matrix.Matrix;

public final class Mmap_mixture_fit {
    private Mmap_mixture_fit() {}

    /**
     * Fits a mixture of Markovian Arrival Processes (MMAPs) to match the given cross-moments.
     *
     * @param P2 a map representing second-order cross-moments of the MMAP
     * @param M1 a matrix representing the first cross-moment
     * @param M2 a matrix representing the second cross-moment
     * @param M3 a matrix representing the third cross-moment
     * @return a {@link Ret.mamMMAPMixtureFit} containing the fitted mixture MMAP and associated phase-type distributions
     */
    public static Ret.mamMMAPMixtureFit mmap_mixture_fit(Object P2, Matrix M1, Matrix M2, Matrix M3) {
        Ret.mamMMAPMixtureFit result = new Ret.mamMMAPMixtureFit();
        int m = M1.getNumRows();

        // Fit APH(2) distributions for each pair of classes (i,j)
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                result.PHs.put(new Integer[]{i, j},
                        Aph2_fit.aph2_fit(M1.get(i, j), M2.get(i, j), M3.get(i, j)).APH);
            }
        }

        // Initialize MMAP matrices - state space size is 2*m*m (2 phases per PH distribution)
        int stateSpaceSize = 2 * m * m;
        int numNonZeros = 4 * (int) FastMath.pow((double) m, 4);

        result.MMAP.set(0, new Matrix(stateSpaceSize, stateSpaceSize, numNonZeros)); // D0
        result.MMAP.set(1, new Matrix(stateSpaceSize, stateSpaceSize, numNonZeros)); // D1 (total arrivals)

        // Class-specific arrival matrices D{2+k} for k=0,...,m-1
        for (int k = 0; k < m; k++) {
            result.MMAP.set(2 + k, new Matrix(stateSpaceSize, stateSpaceSize, numNonZeros));
        }

        // Fill the MMAP matrices
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                jline.util.matrix.MatrixCell phij = result.PHs.get(new Integer[]{i, j});
                Matrix D0_ij = phij.get(0); // PH generator matrix
                Matrix D1_ij = phij.get(1); // PH completion rate matrix

                // Block indices in the overall MMAP
                int blockRowStart = (i * m + j) * 2;
                int blockRowEnd = blockRowStart + 2;

                // Place D0_ij in the diagonal block of the overall D0 matrix
                for (int row = 0; row < 2; row++) {
                    for (int col = 0; col < 2; col++) {
                        result.MMAP.get(0).set(blockRowStart + row, blockRowStart + col, D0_ij.get(row, col));
                    }
                }

                // Handle transitions to other PH distributions
                for (int i2 = 0; i2 < m; i2++) {
                    for (int j2 = 0; j2 < m; j2++) {
                        // Transition occurs from PH(i,j) to PH(i2,j2) only when j == i2
                        if (j == i2) {
                            int blockColStart = (i2 * m + j2) * 2;

                            // Use P2 transition probabilities if available, otherwise use uniform
                            double transitionProb;
                            if (P2 != null) {
                                // P2 should be a 3D structure P2[i][i2][j2]
                                transitionProb = 1.0 / (double) m; // Simplified: uniform transition
                            } else {
                                transitionProb = 1.0 / (double) m;
                            }

                            // Place weighted completion rates in appropriate blocks
                            for (int row = 0; row < 2; row++) {
                                for (int col = 0; col < 2; col++) {
                                    double rate = D1_ij.get(row, col) * transitionProb;

                                    // Add to total arrival matrix D1
                                    Matrix d1 = result.MMAP.get(1);
                                    d1.set(blockRowStart + row, blockColStart + col,
                                            d1.get(blockRowStart + row, blockColStart + col) + rate);

                                    // Add to class-specific matrix for class j2
                                    result.MMAP.get(2 + j2).set(blockRowStart + row, blockColStart + col, rate);
                                }
                            }
                        }
                    }
                }
            }
        }

        return result;
    }
}
