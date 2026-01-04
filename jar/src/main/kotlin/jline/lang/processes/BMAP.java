/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.io.Serializable;

import static jline.api.mam.Map_pieKt.map_pie;

/**
 * A Batch Markovian Arrival Process (BMAP)
 *
 * BMAP is a point process where arrivals occur in batches.
 * Uses the standard BMAP representation:
 * - D0: infinitesimal generator for transitions without arrivals
 * - D1: rate matrix for transitions generating 1 arrival
 * - D2: rate matrix for transitions generating 2 arrivals
 * - ...
 * - Dk: rate matrix for transitions generating k arrivals
 *
 * BMAP extends MarkedMAP where each "mark" k represents a batch size k.
 */
public class BMAP extends MarkedMAP implements Serializable {

    /**
     * Default constructor creating an empty BMAP
     */
    public BMAP() {
        super();
    }

    /**
     * Construct a BMAP from a MatrixCell containing {D0, D1, D2, ..., Dk}
     *
     * @param bmap MatrixCell where:
     *             - bmap.get(0) = D0 (transitions without arrivals)
     *             - bmap.get(k) = Dk (transitions generating k arrivals) for k >= 1
     */
    public BMAP(MatrixCell bmap) {
        super(bmap);
        validateGenerator();
    }

    /**
     * Construct a BMAP from D0 and variable number of Dk matrices
     *
     * @param D0 the infinitesimal generator for transitions without arrivals
     * @param Dk variable number of matrices where Dk[i-1] represents transitions generating i arrivals
     */
    public BMAP(Matrix D0, Matrix... Dk) {
        super(createMatrixCell(D0, Dk));
        validateGenerator();
    }

    /**
     * Returns the name of this distribution type.
     *
     * @return "BMAP" to identify this as a Batch Markovian Arrival Process
     */
    @Override
    public String getName() {
        return "BMAP";
    }

    /**
     * Helper method to create MatrixCell from D0 and Dk matrices
     */
    private static MatrixCell createMatrixCell(Matrix D0, Matrix... Dk) {
        MatrixCell bmap = new MatrixCell();
        bmap.set(0, D0);

        // Compute total arrival matrix D1 (sum of all Dk)
        Matrix D1_total = new Matrix(D0.getNumRows(), D0.getNumCols());
        for (Matrix d : Dk) {
            D1_total.addEq(1.0, d);
        }
        bmap.set(1, D1_total);

        // Set individual Dk matrices (starting from index 2)
        for (int k = 0; k < Dk.length; k++) {
            bmap.set(2 + k, Dk[k]);
        }

        return bmap;
    }

    /**
     * Factory method to create BMAP from a base MAP and batch size distribution
     *
     * Given a base MAP (D0_base, D1_base) for inter-batch arrivals and a batch size
     * distribution, constructs the BMAP by scaling: Dk = D1_base * pmf[k-1]
     *
     * @param D0         the base MAP's D0 matrix (transitions without batch arrivals)
     * @param D1         the base MAP's D1 matrix (inter-batch arrival transitions)
     * @param batchSizes array of batch sizes (e.g., [1, 2, 4, 8])
     * @param pmf        probability mass function for batch sizes (must sum to 1)
     * @return a BMAP constructed from the base MAP and batch distribution
     */
    public static BMAP fromMAPWithBatchPMF(Matrix D0, Matrix D1, int[] batchSizes, double[] pmf) {
        if (batchSizes.length != pmf.length) {
            throw new IllegalArgumentException("Batch sizes and PMF must have the same length");
        }

        // Normalize PMF
        double sum = 0.0;
        for (double p : pmf) {
            sum += p;
        }
        double[] normPMF = new double[pmf.length];
        for (int i = 0; i < pmf.length; i++) {
            normPMF[i] = pmf[i] / sum;
        }

        // Find maximum batch size to determine number of D matrices needed
        int maxBatch = 0;
        for (int size : batchSizes) {
            if (size > maxBatch) {
                maxBatch = size;
            }
        }

        // Create D matrices: Dk = D1 * pmf(batchSize==k)
        Matrix[] Dk = new Matrix[maxBatch];
        for (int k = 0; k < maxBatch; k++) {
            Dk[k] = new Matrix(D0.getNumRows(), D0.getNumCols()); // Initialize to zero
        }

        // Set Dk based on batch sizes and PMF
        for (int i = 0; i < batchSizes.length; i++) {
            int batchSize = batchSizes[i];
            if (batchSize < 1 || batchSize > maxBatch) {
                throw new IllegalArgumentException("Invalid batch size: " + batchSize);
            }
            // Dk[k-1] corresponds to batch size k
            Matrix scaled = Matrix.scaleMult(D1, normPMF[i]);
            Dk[batchSize - 1].addEq(1.0, scaled);
        }

        return new BMAP(D0, Dk);
    }

    /**
     * Validates that D0 + sum(Dk) forms a proper infinitesimal generator
     * (all row sums should equal zero)
     *
     * @throws IllegalArgumentException if the row sums are not close to zero
     */
    private void validateGenerator() {
        Matrix D0 = D(0);
        int nPhases = D0.getNumRows();
        MatrixCell proc = getProcess();

        // Compute D0 + D1_total (D1_total at index 1 is already the sum of all Dk)
        Matrix generator = D0.copy();
        generator.addEq(1.0, proc.get(1)); // D1_total = sum of all Dk

        // Check row sums (should be close to zero for valid generator)
        double maxRowSum = 0.0;
        for (int i = 0; i < nPhases; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < nPhases; j++) {
                rowSum += generator.get(i, j);
            }
            maxRowSum = Math.max(maxRowSum, Math.abs(rowSum));
        }

        // Allow small numerical errors (up to 1e-6)
        if (maxRowSum > 1e-6) {
            throw new IllegalArgumentException(
                    "Invalid BMAP: D0 + sum(Dk) must form a valid infinitesimal generator " +
                    "(row sums must equal zero). Maximum row sum deviation: " + maxRowSum);
        }
    }

    /**
     * Get the maximum batch size supported by this BMAP
     *
     * @return the maximum batch size k where Dk is defined
     */
    public int getMaxBatchSize() {
        // Number of D matrices is process.size()
        // D0 is at index 0, D1_total at index 1, D1 (batch size 1) at index 2, etc.
        // So max batch size = process.size() - 2
        return getProcess().size() - 2;
    }

    /**
     * Get the mean batch size (expected number of arrivals per batch)
     *
     * Computed as weighted average of batch sizes by their rates:
     * E[batch size] = sum(k * rate_k) / sum(rate_k)
     *
     * @return the mean batch size
     */
    public double getMeanBatchSize() {
        try {
            Matrix D0 = D(0);
            MatrixCell proc = getProcess();

            // Get stationary distribution using D0 and D1_total
            Matrix pi_mat = map_pie(getProcess());
            double[] pi = pi_mat.toArray1D();

            // Compute weighted average of batch sizes
            double totalRate = 0.0;
            double weightedSum = 0.0;

            // Iterate over batch sizes (k = 1, 2, ...)
            // proc.get(0) = D0, proc.get(1) = D1_total
            // proc.get(2+k) = D_{k+1} (batch size k+1)
            for (int idx = 2; idx < proc.size(); idx++) {
                int k = idx - 1; // batch size (D1 at index 2 -> batch size 1)
                Matrix Dk = proc.get(idx);

                // Compute rate for batch size k: rate_k = pi * Dk * 1
                double rate_k = 0.0;
                for (int i = 0; i < pi.length; i++) {
                    for (int j = 0; j < Dk.getNumCols(); j++) {
                        rate_k += pi[i] * Dk.get(i, j);
                    }
                }

                totalRate += rate_k;
                weightedSum += k * rate_k;
            }

            return (totalRate > 0) ? (weightedSum / totalRate) : 0.0;
        } catch (Exception e) {
            return 0.0;
        }
    }

    /**
     * Get the arrival rates for each batch size
     *
     * @return array where result[k-1] is the rate of batch size k arrivals
     */
    public double[] getBatchRates() {
        try {
            Matrix D0 = D(0);
            MatrixCell proc = getProcess();
            int maxBatch = getMaxBatchSize();

            // Get stationary distribution using D0 and D1_total
            Matrix pi_mat = map_pie(getProcess());
            double[] pi = pi_mat.toArray1D();

            // Compute rate for each batch size
            double[] rates = new double[maxBatch];
            for (int idx = 2; idx < proc.size(); idx++) {
                int k = idx - 1; // batch size
                Matrix Dk = proc.get(idx);

                // rate_k = pi * Dk * 1
                double rate_k = 0.0;
                for (int i = 0; i < pi.length; i++) {
                    for (int j = 0; j < Dk.getNumCols(); j++) {
                        rate_k += pi[i] * Dk.get(i, j);
                    }
                }

                rates[k - 1] = rate_k;
            }

            return rates;
        } catch (Exception e) {
            return new double[getMaxBatchSize()];
        }
    }

    /**
     * Get the D_k matrix for batch size k
     *
     * @param batchSize the batch size (1, 2, 3, ...)
     * @return D_k matrix for the specified batch size, or null if out of range
     */
    public Matrix getBatchMatrix(int batchSize) {
        if (batchSize < 1 || batchSize > getMaxBatchSize()) {
            return null;
        }
        // D_k for batch size k is at index k+1 in process
        return getProcess().get(batchSize + 1);
    }

    @Override
    public String toString() {
        return "BMAP(phases=" + nPhases + ", maxBatchSize=" + getMaxBatchSize() + ")";
    }
}
