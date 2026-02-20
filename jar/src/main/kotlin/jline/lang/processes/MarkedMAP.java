/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.processes;

import jline.util.RandomManager;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.io.Serializable;
import java.util.Random;

import static jline.api.mam.Mmap_timereverseKt.mmap_timereverse;

/**
 * A Marked Markovian Arrival Process
 */
public class MarkedMAP extends Marked implements Serializable {

    public MarkedMAP() {
        super("MarkedMAP", 0);
    }

    public MarkedMAP(MatrixCell mmap) {
        super("MMAP", mmap.size());
        setProcess(mmap);
        Matrix D0 = mmap.get(0);
        Matrix D1 = mmap.get(1);
        int nPhases = D0.getNumCols();
        this.setParam(1, "D0", D0);
        this.setParam(2, "D1", D1);
        for (int i = 2; i < mmap.size(); i++) {
            this.setParam(i + 1, "D1" + i, mmap.get(i));
        }

        this.nPhases = nPhases;
    }

    /**
     * Normalizes the MarkedMAP so that D0+sum(D_i) rows form a proper infinitesimal generator.
     * Each row sum should equal zero for a valid generator.
     */
    public void normalize() {
        Matrix D0 = D(0).copy();
        MatrixCell proc = getProcess();

        // Ensure non-negative off-diagonal elements in D0 and all elements in other matrices
        for (int i = 0; i < nPhases; i++) {
            for (int j = 0; j < nPhases; j++) {
                if (i != j && D0.get(i, j) < 0) {
                    D0.set(i, j, 0.0);
                }
                // Ensure non-negative elements in all D_k matrices (k > 0)
                for (int k = 1; k < proc.size(); k++) {
                    Matrix Dk = proc.get(k);
                    if (Dk.get(i, j) < 0) {
                        Dk.set(i, j, 0.0);
                    }
                }
            }
        }

        // Adjust diagonal elements of D0 to make row sums zero
        for (int i = 0; i < nPhases; i++) {
            double rowSum = 0.0;
            // Sum off-diagonal elements of D0
            for (int j = 0; j < nPhases; j++) {
                if (i != j) {
                    rowSum += D0.get(i, j);
                }
            }
            // Sum all elements in D_k matrices (k > 0) for row i
            for (int k = 1; k < proc.size(); k++) {
                Matrix Dk = proc.get(k);
                for (int j = 0; j < nPhases; j++) {
                    rowSum += Dk.get(i, j);
                }
            }
            // Set diagonal element to make row sum zero
            D0.set(i, i, -rowSum);
        }

        // Update D0 in the process
        proc.set(0, D0);
        setProcess(proc);
    }

    @Override
    public double[] sample(int n) {
        return this.sample(n, RandomManager.getThreadRandomAsRandom());
    }

    /**
     * Generates samples from the MarkedMAP using the specified random generator.
     * Returns inter-arrival times with corresponding marks.
     *
     * @param n      the number of samples to generate
     * @param random the random number generator to use
     * @return array of inter-arrival times (marks are not returned in this basic implementation)
     */
    @Override
    public double[] sample(int n, Random random) {
        try {
            // For now, use the underlying MAP sampling by considering only D0 and sum of all D_i
            Matrix D0 = D(0);
            MatrixCell proc = getProcess();

            // Create aggregated D1 by summing all arrival matrices
            Matrix D1_agg = new Matrix(D0.getNumRows(), D0.getNumCols());
            for (int k = 1; k < proc.size(); k++) {
                D1_agg.addEq(1.0, proc.get(k));
            }

            // Use MAP sampling with aggregated arrival matrix
            return jline.api.mam.Map_sampleKt.map_sample(D0, D1_agg, n, random);
        } catch (Exception e) {
            // Fallback: generate exponential samples with the mean rate
            double[] samples = new double[n];
            double mean = getMean();
            for (int i = 0; i < n; i++) {
                samples[i] = -mean * Math.log(random.nextDouble());
            }
            return samples;
        }
    }

    public MarkedMAP toTimeReversed() {
        MatrixCell reversed = mmap_timereverse(this.process);
        return new MarkedMAP(reversed);
    }


}
